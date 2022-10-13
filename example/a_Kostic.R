library(moltenNMF)
library(readxl)
library(tidyverse)
library(parallel)
library(Matrix)
library(tidytext)
library(ROCR)
#library(patchwork)
#library(Rcpp)
#sourceCpp("~/Dropbox/moltenNMF/moltenNMF.cpp")
tab0 <- read_tsv("~/data/Kostic/diabimmune_t1d_16s_otu_table.txt",skip=1)
abundance <- tab0 %>% 
  dplyr::select(starts_with("G")) %>% 
  mutate_all(as.numeric) %>% t()
meta <- read_excel("/Users/koabe/data/Kostic/diabimmune_t1d_wgs_metadata.xlsx")
dat <- read_tsv("/Users/koabe/data/Kostic/diabimmune_t1d_metaphlan_table.txt")
tab1 <- t(dat[,-1])
Taxonomylist <- strsplit(dat$Taxonomy,"\\|")
len <- sapply(Taxonomylist, length)
genuslevel <- tab1[,len==6]

colnames(genuslevel) <- sapply(Taxonomylist[len==6],function(x)x[6])
# table(floor(3*datall$Age_at_Collection/365))
datall <-as.data.frame(genuslevel) %>% 
  rownames_to_column("Gid_16S") %>% 
  left_join(meta) %>% 
  mutate(age = Age_at_Collection,
         Case_Controlcontrol = factor(Case_Control, levels = c("control","case"))) %>% 
  select(starts_with("g__"),Subject_ID,age,Case_Control,Total_Reads) %>% 
  gather(genus,abundance,-Subject_ID,-age,-Case_Control,-Total_Reads) %>% 
  mutate(abundance = abundance*Total_Reads/100) %>% 
  mutate(abundance = as.integer(round(abundance))) %>% 
  mutate(Age2 = paste(Subject_ID, age, sep = "_"))

f0 <- abundance ~ genus+Age2-1
f1 <- abundance ~ genus+Age2+Case_Control-1
sq <- seq(2,80,by=6)
# system.time({
#   out0 <- mclapply(sq,function(l){
#     mNMF_vb(f1, data = datall, L=l, iter=1000, a=0.5, b=1)},
#     mc.cores = detectCores())
#   out1 <- mclapply(sq,function(l){
#     mNMF_vb(f1, data = datall, L=l, iter=1000, a=0.5, b=1)},
#     mc.cores = detectCores())
# })
# 
# f20 <- abundance ~ genus+factor(age)+Subject_ID-1
# f21 <- abundance ~ genus+factor(age)+Subject_ID+Case_Control-1
# 
# sq <- seq(2,80,by=6)
# system.time({
#   out20 <- mclapply(sq,function(l){
#     mNMF_vb(f20, data = datall, L=l, iter=1000, a=0.5, b=1)},
#     mc.cores = detectCores())
#   out21 <- mclapply(sq,function(l){
#     mNMF_vb(f21, data = datall, L=l, iter=1000, a=0.5, b=1)},
#     mc.cores = detectCores())
# })
# 
# f3 <- abundance ~ genus+Subject_ID+paste(age, Case_Control, sep = "_")-1
# system.time({
#   out3 <- mclapply(sq,function(l){
#     mNMF_vb(f20, data = datall, L=l, iter=1000, a=0.5, b=1)},
#     mc.cores = detectCores())
# })
#save(out0, out1, out20, out21, out3, file = "moltenNMF_Kostic.Rdata")
load("./example/moltenNMF_Kostic.Rdata")
ELBO0 <- sapply(out0, function(x)x$ELBO[length(x$ELBO)])
ELBO1 <- sapply(out1, function(x)x$ELBO[length(x$ELBO)])
ELBO20 <- sapply(out20, function(x)x$ELBO[length(x$ELBO)])
ELBO21 <- sapply(out21, function(x)x$ELBO[length(x$ELBO)])
ELBO3 <- sapply(out3, function(x)x$ELBO[length(x$ELBO)])

dfelbo <- data.frame(L=sq,
           model0=ELBO0,
           model1=ELBO1,
           model20=ELBO20,
           model21=ELBO21,
           model3=ELBO3) %>% 
  pivot_longer(2:6)

ggplot(dfelbo, aes(x=L, y=-log(-value), colour=name, group=name))+
  geom_line()+
  theme_minimal(14)
ggsave("kostic_elbo_comp.pdf")

outfix <- out1[[which.max(ELBO1)]]
X1 <- sparse_model_matrix_b(f1, data = datall)
data <- model.frame(f1,  datall)
t <- terms(f1, data=data)
X <- Matrix:::model.spmatrix(trms = t, mf = data,
                             verbose = TRUE)
head(data)
X1 <- sparse.model.matrix(f1, data = datall)

prob <- xprob_mNMF(varname="Case_Control",
                   vargroup=attr(t,"term.labels")[attr(X1,"assign")],
                   V=outfix$shape/outfix$rate,
                   X=X1,Y=datall$abundance)
head(prob)
X1 <- as.matrix(X1)
colnames(X1)
df <- data.frame(p=drop(prob), x=X1[,colnames(X1)=="Case_Controlcontrol"])
head(df)
ggplot(df,aes(x=p,y=x))+
  geom_jitter(alpha=0.01, width = 0)+
  theme_minimal(16)

ggplot(df,aes(x=p, colour=x, fill=x))+
  geom_area(stat = "bin", aes(y=after_stat(density)),
            colour="black",bins=50, alpha=0.1)+
  theme_minimal(16)

pred <- prediction(df$p, df$x)
perf <- performance(pred,"tpr","fpr")

df <- data.frame(p=perf@alpha.values[[1]],
           x=perf@x.values[[1]],
           y=perf@y.values[[1]])

ggplot(df,aes(x=x,y=y))+
  geom_step()+
  geom_abline(intecept=0,slope=1,linetype=2)+
  labs(x=perf@x.name, y=perf@y.name)+
  theme_bw(14)

ggsave("./example/ROCcurve.pdf")

dsystem.time({
  ytilde <- rpredictor_mNMF(X, 100, outfix$shape, outfix$rate)
})

yq <- apply(ytilde, 1, quantile, prob=c(0.05, 0.95))
#yr <- apply(ytilde, 1, range)
head(t(yq))

fitdf <- data.frame(fit = apply(ytilde, 1, mean),
                    obs = datall$abundance,
                    lower = yq[1,], upper = yq[2,])

with(fitdf, mean(lower <= obs & obs <= upper))

ggplot(fitdf,aes(x=obs, y=fit, ymax=upper, ymin=lower))+
  geom_abline(intercept=0, slope=1, linetype=2)+
  geom_point(alpha=0.1)+
  geom_linerange(alpha=0.4)+
  theme_classic(16)
#ggsave("fit.jpg")

Vdf <- data.frame(outfix$shape/outfix$rate,
                  variable=rownames(outfix$shape),
                  facet_dummy=outfix$vargroup) %>% 
  pivot_longer(!(variable|facet_dummy), names_to = "component", 
               names_transform = list(component = readr::parse_number))
  #mutate(component=factor(component,ordered = TRUE))

controldf <- dplyr::filter(Vdf,"Case_Control" == facet_dummy) %>% 
  mutate(variable = gsub("Case_Control","",variable))

ggplot(controldf, aes(x=component, y=log(value)))+
  geom_point()+
  geom_linerange(aes(ymax=log(value), ymin=0))+
  geom_hline(yintercept = 0)+
  theme_classic()

pickcontdown <- dplyr::filter(controldf,log(value) >= 1) %>% 
  dplyr::select(component, value)

pickcontdown <- dplyr::filter(controldf,log(value) <= -8) %>% 
  dplyr::select(component, value)

pickcontdown <- dplyr::filter(controldf,log(value) <= -8) %>% 
  dplyr::select(component, value)

Agedf <- dplyr::filter(Vdf,"Age2" == facet_dummy) %>% 
  mutate(variable = gsub("Age2","",variable)) %>% 
  separate(variable, c("SubjectID","Age"), sep = "_") %>% 
  group_by(Age, SubjectID) %>% 
  mutate(value=value/sum(value)) %>% 
  mutate(Age = factor(as.integer(Age)),orderd=TRUE)

ggplot(Agedf,aes(x=Age, y = component, fill=log(value)))+
  geom_tile()+
  facet_wrap(SubjectID ~., scales = "free")+
  scale_fill_viridis_c()+
  theme_classic(14)+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(size=7))

#ggsave("Kostic_age.png")
#ggsave(filename = "~/Desktop/Kostic_age.png")

head(Vdf)
genusdf <- dplyr::filter(Vdf, "genus"==facet_dummy) %>% 
  mutate(genus=gsub("genusg__","",variable))


genusdf2 <-group_by(genusdf, component) %>% 
  dplyr::filter(component %in% pickcontdown$component) %>% 
  top_n(3,value)

log(pickcontdown$value)

ggplot(genusdf2,aes(y=reorder_within(genus,value,component), x=value))+
  geom_point()+
  geom_linerange(aes(xmin=0, xmax=value))+
  facet_wrap(component~., scales = "free")+
  scale_y_reordered()+
  labs(y="", x="")+
  theme_classic()+
  theme(strip.background = element_blank(),
        axis.text.y = element_text(colour = "black", size=9),
        axis.text.x = element_text(colour = "black", size=6))

#ggsave("Kostic_genus.png")
