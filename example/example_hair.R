#devtools::install_github("abikoushi/moltenPPCA", build_vignettes = TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(moltenPPCA)
library(Rcpp)

# sourceCpp("./projects/moltenPPCA/src/estimator_pois.cpp")
# source("./projects/moltenPPCA/R/Poisson.R")
# source("./projects/moltenPPCA/R/product.R")
# source("./projects/moltenPPCA/R/modelMatrix.R")

Ttitanicdf <- as.data.frame(Titanic)
ggplot(hairdf, aes(x=Hair,y=Freq))+
  geom_col()+
  facet_grid(Eye~Sex)+
  theme_classic(16)

set.seed(794)
out <- mNMF_vb(Freq~Hair+Eye+Sex, data=hairdf, L=2, iter=500)
qplot(1:500,out$ELBO, geom = "line")+
  theme_minimal(16)

f <- ~Hair+Eye+Sex
vname <- unlist(strsplit(as.character(f)[2],"[+]"))
vname <- factor(vname,levels = vname)
V <- out$shape/out$rate
VL <- vector("list", length(out$indices)-2)
for(i in 2:(length(out$indices)-1)){
  iii <- (out$indices[i]+1):out$indices[i+1]
  VL[[i]] <- mutate(gather(mutate(data.frame(V[iii,,drop=FALSE]),variable=rownames(V[iii,,drop=FALSE])),component,value,-variable),
                    facet_dummy=vname[i-1])
}
Vdf <- bind_rows(VL)

ggplot(Vdf, aes(x=component,y=variable,fill=value))+
  geom_tile(colour="gray20")+
  scale_fill_continuous(low="white", high="royalblue")+
  facet_grid(facet_dummy~., scales = "free_y", space = "free")+
  theme_classic(16)

ggplot(Vdf, aes(y=variable,x=value,fill=component))+
  geom_col(width = 1, colour="gray20")+
  facet_grid(facet_dummy~., scales = "free", space = "free")+
  scale_fill_manual(values = c("royalblue","orange2"))+
  theme_classic(16)
#ggsave("./Desktop/test.png")