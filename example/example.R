library(tidyr)
library(dplyr)
library(Matrix)
library(Rcpp)
library(ggplot2)
library(moltenPPCA)

iris_g <- mutate(iris, id=factor(row_number())) %>% 
  gather(key,value,-id,-Species)
  # group_by(key) %>%
  # mutate(value = scale(value, center = TRUE)) %>%
  # ungroup()

ggplot(iris_g,aes(x=value))+
  geom_freqpoly(bins=20)+
  facet_wrap(~key)+
  scale_color_viridis_d()+
  theme_classic(16)

set.seed(9999)
out <- mPPCA_vb(value~id+key-1,data=iris_g,L=2, iter=200)
qplot(1:200,out$logloss, geom = "line")+
  theme_minimal(16)

iii <- (out$indices[1]+1):out$indices[1+1]
df <- data.frame(out$mean[iii,], species=iris$Species)

ggplot(df, aes(X1,X2,colour=species))+
  geom_point(alpha=0.7, size=2)+
  theme_minimal(16)

#####

out <- mPPCA_Gibbs(y~id+key,data=iris_g,L=2, iter=2000)
plot(out$loglik, type="l")
iii <- (out$indices[2]+1):out$indices[2+1]
df <- data.frame(out$mean[iii,], species=iris$Species[-1])

ggplot(df, aes(X1,X2,colour=species))+
  geom_point(alpha=0.7, size=2)+
  scale_color_viridis_d()+
  theme_classic(16)


burnin <- 1:1000
mubar <- apply(out$mu[,,-burnin],1:2,mean)

iii <- (out$indices[2]+1):out$indices[2+1]
df <- data.frame(mubar[iii,], species=iris$Species[-1])

ggplot(df, aes(X1,X2,colour=species))+
  geom_point(alpha=0.7, size=2)+
  scale_color_viridis_d()+
  theme_classic(16)

# df2 <- data.frame(mubar[jjj,], key=colnames(X)[jjj]) %>% 
#   gather(l,value,-key) %>% 
#   mutate(key=sub("key","",key))
# 
# ggplot(df2, aes(key,value))+
#   geom_col()+
#   facet_wrap(~l)+
#   theme_classic(16)+
#   theme(axis.text.x = element_text(angle=90))
