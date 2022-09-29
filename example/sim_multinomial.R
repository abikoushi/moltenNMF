library(dplyr)
library(tidyr)
library(ggplot2)
library(moltenNMF)

n <- 100L
m <- 100L
k <- 3L
perm_k <- gtools::permutations(k,k)

set.seed(555)
W <- matrix(rgamma(k*n,2), n, k)
W <- W/rowSums(W)
H <- matrix(rgamma(k*m,2), k, m)
H <- H/rowSums(H)
Ns <- 10000L
X <- matrix(0L, n, m)
for(i in 1:n){
  X[i,] <- as.vector(rmultinom(1, Ns, W[i,]%*%H))
}

dfX <- pivot_longer(data.frame(id=1:n, X), -id)
out <- mNMF_vb(value~name+factor(id)-1, data=dfX, L=k, iter=10000, a=0.5, b=1)
plot(out$ELBO, type = "l")
V <- out$shape/out$rate
yhat <- product_m(value~name+factor(id)-1, data=dfX, V)

H2 <- t(V["name"==out$vargroup, ])
ord <- order(as.integer(gsub("nameX","",colnames(H2))))
H2 <- H2[,ord]
H2 <- H2/rowSums(H2)
corv <- apply(perm_k, 1, function(x)cor(as.vector(log(H2[x,])), as.vector(log(H))))

sv <- perm_k[which.max(corv),]

qplot(dfX$value, yhat, alpha=I(0.1))+
  labs(x="fit", y="obs")+
  geom_abline(intercept = 0, slope=1, linetype=2)+
  theme_minimal(16)

Vdf <- data.frame(V,variable=rownames(V),
                  facet_dummy=out$vargroup) %>% 
  pivot_longer(!(variable|facet_dummy), names_to = "component", 
               names_transform = list(component = readr::parse_number)) %>% 
  mutate(component=factor(component))

# ggplot(Vdf, aes(y=variable, x=value, fill=component))+
#   geom_col(width = 1)+
#   facet_grid(facet_dummy~.,scales="free_y",space = "free")+
#   theme_classic(16)+
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# 
# ggplot(Vdf, aes(y=variable, x=value, fill=component))+
#   geom_col(position = "fill", width = 1)+
#   facet_grid(facet_dummy~.,scales="free_y",space = "free")+
#   scale_x_continuous(labels=scales::percent)+
#   theme_classic(16)+
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Hdf <- data.frame(true=as.vector(H[sv,]),
           est = as.vector(H2),
           component=factor(rep(1:k,ncol(H))))

ggplot(Hdf,aes(x=true,y=est, shape=component, colour=component))+
  geom_abline(intercept = 0, slope=1, linetype=2, colour="lightgrey", size=1)+
  geom_point()+
  labs(x="true value", y="estimates")+
  theme_minimal(16)
ggsave("comparisonplot_mult_h.pdf")

W2 <- V[out$vargroup=="factor(id)",]

dfW <- data.frame(component=rep(sv,each=nrow(W[-1,])),
           true=as.vector(W[-1,]),
           est=as.vector(W2[,sv]))

ggplot(dfW, aes(x=true, y=est, group=component))+
  labs(x="true value", y="estimates")+
  stat_smooth(method = "lm", formula = y~x,
              se=FALSE, linetype=2, size=1, colour="lightgrey")+
  geom_point(shape=1)+
  facet_wrap(~component, scales = "free_y", labeller = label_both)+
  theme_classic(16)
ggsave("comparisonplot_mult_w.pdf")
