library(dplyr)
library(tidyr)
library(ggplot2)
library(moltenNMF)

set.seed(1192)
n <- 10L
m <- 10L
k <- 2L
W <- matrix(rgamma(k*n,2), n, k)
W <- W/rowSums(W)
H <- matrix(rgamma(k*m,2), k, m)
H <- H/rowSums(H)
lam <- rgamma(n, 2, scale = 100)
Ns <- rpois(n,lam)
X <- matrix(0L, n, m)
for(i in 1:10){
  X[i,] <- as.vector(rmultinom(1,Ns[i],W[i,]%*%H))
}

dfX <- pivot_longer(data.frame(id=1:n,X), -id)
out <- mNMF_vb(value~name+factor(id)-1, data=dfX, L=k, iter=2500, a=0.5, b=0.1)
plot(out$ELBO, type = "l")
V <- out$shape/out$rate
yhat <- rowSums(prod_mNMF(value~name+factor(id)-1, data=dfX, V))

qplot(dfX$value, yhat, alpha=I(0.5))+
  labs(x="fit", y="obs")+
  geom_abline(intercept = 0, slope=1, linetype=2)+
  theme_minimal(16)

Vdf <- data.frame(V,variable=rownames(V),
                  facet_dummy=out$vargroup) %>% 
  pivot_longer(!(variable|facet_dummy), names_to = "component", 
               names_transform = list(component = readr::parse_number)) %>% 
  mutate(component=factor(component))

ggplot(Vdf, aes(y=variable, x=value, fill=component))+
  geom_col(width = 1)+
  facet_grid(facet_dummy~.,scales="free_y",space = "free")+
  theme_classic(16)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggplot(Vdf, aes(y=variable, x=value, fill=component))+
  geom_col(position = "fill", width = 1)+
  facet_grid(facet_dummy~.,scales="free_y",space = "free")+
  scale_x_continuous(labels=scales::percent)+
  theme_classic(16)+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


H2 <- t(V["name"==out$vargroup, ])
ord <- order(as.integer(gsub("nameX","",colnames(H2))))
H2 <- H2[,ord]
H2 <- H2/rowSums(H2)

qplot(as.vector(H[2:1,]), as.vector(H2))+
  labs(x="true", y="est")+
  geom_abline(intercept = 0, slope=1, linetype=2)+
  theme_minimal(16)
ggsave("comparisonplot_mult.png")

# plot(W[-1,2:1], V[out$vargroup=="factor(id)",])
# abline(0,1,lty=2)
