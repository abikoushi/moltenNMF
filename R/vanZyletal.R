library(Matrix)
library(readr)
library(tidyr)
library(dplyr)
library(moltenNMF)

Hdf = readRDS("H1list.rds")

ncol = n_distinct(unlist(Hdf))
SS = sapply(Hdf, length)

nrow = sum(sapply(Hdf, length))

Hmat = sparseMatrix(i=integer(0),j=integer(0),dims = c(nrow,ncol))

lev = unique(unlist(Hdf))
Hdf = unlist(Hdf)
CS = c(0L, cumsum(SS))
for(s in 1:5){
    js = match(Hdf[[s]], lev)
    Hmat[CS[s]+1:SS[s],js] = TRUE
  print(s)
}
saveRDS(Hmat, file="Hmat.rds")

df <- read_csv("dfX.csv")

nrow
dim(df)


#lev = sort(unique(Hdf$Human))
#factor(Hdf$Human)
mutate(Hdf, pos_humang = as.integer(factor(Human)))

df <- read_csv("dfX.csv")
group_by(df, species) %>% 
  summarise(n_distinct(Human))

f <- ~ species_gene + cell + sample + species
X <- moltenNMF::sparse_onehot(f, data = df)

Y = readRDS("listY.rds")
Y = do.call("c", Y)
nnzero(Y)
X0Count = readRDS("listX0Count.rds")

X0prob <- vector("list", 5)
for(i in 1:5){
  X0prob[[i]] <- bind_rows(lapply(X0Count, function(x)data.frame(x[[i]])))  
}

N0 = max(sapply(X0prob,function(X)sum(X$Freq)))
X0prob <- unlist(sapply(X0prob, function(X)X$Freq/N0))


####
# n_max = 1000
system.time({
  out <- moltenNMF:::mNMF_vb_sp(y = Y@x, X = X, L = 15L,
                                N0 = N0,probX0 = X0prob,
                                iter = 10)
})
# user  system elapsed 
# 78.874   1.830  81.087 
plot(out$ELBO[-1])
