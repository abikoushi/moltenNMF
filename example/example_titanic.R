#devtools::install_github("abikoushi/moltenPPCA", build_vignettes = TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)
library(moltenPPCA)
library(patchwork)

#library(Rcpp)
# sourceCpp("./projects/moltenPPCA/src/estimator_pois.cpp")
# source("./projects/moltenPPCA/R/Poisson.R")
# source("./projects/moltenPPCA/R/product.R")
# source("./projects/moltenPPCA/R/modelMatrix.R")


Titanicdf <- as.data.frame(Titanic) %>% 
  mutate(Class=factor(Class,levels=c("3rd","2nd","1st","Crew")))
ggplot(Titanicdf, aes(x=Class, y=Freq, fill=Survived))+
  geom_col(position = "fill")+
  facet_grid(Age~Sex)+
  theme_classic(16)

#ggsave("./Desktop/test.png")
set.seed(794)
out <- mNMF_vb(Freq~Survived+Class+Sex+Age-1, data=Titanicdf, L=2, iter=500)

qplot(1:length(out$ELBO),out$ELBO, geom = "line")+
  labs(x="iter",y="ELBO")+
  theme_minimal(16)

V <- out$shape/out$rate
yhat <- rowSums(prod_mPPCA(~Survived+Class+Sex+Age-1, data=Titanicdf,V))
qplot(yhat, Titanicdf$Freq, alpha=I(0.5),size=I(3))+
  labs(x="fit",y="obs")+
  geom_abline(intercept=0,slope=1,linetype=2)+
  theme_minimal(16)

Vdf <- data.frame(V,variable=rownames(V),
                  facet_dummy=out$vargroup) %>% 
  pivot_longer(!(variable|facet_dummy), names_to = "component", 
               names_transform = list(component = readr::parse_number)) %>% 
  mutate(component=factor(component))

p1 <- ggplot(Vdf, aes(y=variable, x=value, fill=component))+
  geom_col(colour="gray20")+
  facet_grid(facet_dummy~.,scales="free_y",space = "free")+
  theme_classic(16)

p2 <- ggplot(Vdf, aes(y=variable, x=value, fill=component))+
  geom_col(colour="gray20", position = "fill")+
  facet_grid(facet_dummy~.,scales="free_y",space = "free")+
  scale_x_continuous(labels=scales::percent)+
  theme_classic(16)

print(p1)
print(p2)
