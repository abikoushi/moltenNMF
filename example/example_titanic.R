#devtools::install_github("abikoushi/moltenNMF")
library(dplyr)
library(tidyr)
library(ggplot2)
library(moltenNMF)
library(Matrix)

Titanicdf <- mutate(as.data.frame(Titanic), 
                    Class=factor(Class,levels=c("3rd","2nd","1st","Crew")))

head(Titanicdf)

ggplot(Titanicdf, aes(x=Class, y=Freq, fill=Survived))+
  geom_col()+
  facet_grid(Age~Sex)+
  theme_classic(16)

#ggsave("Titanic.png")

ggplot(Titanicdf, aes(x=Class, y=Freq, fill=Survived))+
  geom_col(position = "fill")+
  facet_grid(Age~Sex)+
  theme_classic(16)

#ggsave("Titanic2.png")

set.seed(1999)
f <- Freq ~ Survived + Class+Sex+Age-1
out <- mNMF_vb(f, data=Titanicdf, L=2, iter=1000, a=0.5, b=0.1)

qplot(1:length(out$ELBO),out$ELBO, geom = "line")+
  labs(x="iter",y="ELBO")+
  theme_minimal(16)

V <- out$shape/out$rate
yhat <- product_m(f, data=Titanicdf, V)
X <- sparse_model_matrix_b(f,Titanicdf)
str(X)
dfx <- as.data.frame(which(X, arr.ind = TRUE))
lab <- rep(c("Survived","Class","Sex","Age"),diff(attr(X,"indices")))
lab <- factor(lab, levels = c("Survived","Class","Sex","Age"))
dfx <- mutate(dfx, name=lab[col])
ggplot(dfx,aes(y=row,x=col))+
  geom_tile(fill="cornflowerblue", colour="white", size=1)+
  facet_grid(~name, space = "free", scales = "free")+
  scale_y_reverse()+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 16),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
ggsave("DesMat.png")
#ytilde <- rpredictor_mrNMF(X, 2000, out$shape, out$rate, out$precision)
ytilde <- rpredictor_mNMF(X, 2000, out$shape, out$rate)
#zhat <- (Titanicdf$Freq+out$precision)/(yhat+out$precision)
ybar <- apply(ytilde, 1, mean)
yq <- apply(ytilde, 1, quantile, prob=c(0.025, 0.975))
fitdf <- data.frame(rownumber = 1:length(Titanicdf$Freq),
                    fit = ybar,
                    obs = Titanicdf$Freq,
                    lower = yq[1,], upper = yq[2,])

ggplot(fitdf,aes(x=rownumber, y=obs-fit, ymin=obs-upper, ymax = obs-lower))+
  geom_hline(yintercept=0, linetype=2)+
  geom_linerange()+
  geom_point()+
  theme_classic(16)



ggplot(fitdf,aes(x=obs-fit))+
  geom_histogram(bins = 20, alpha=0.1, colour="grey20")+
  theme_classic(16)

# ggplot(fitdf,aes(x=obs-fit))+
#   stat_ecdf()+
#   theme_classic(16)

Vdf <- data.frame(V,variable=rownames(V),
                  facet_dummy=out$vargroup) %>% 
  pivot_longer(!(variable|facet_dummy), names_to = "component", 
               names_transform = list(component = readr::parse_number)) %>% 
  mutate(component=factor(component))

ggplot(Vdf, aes(y=variable, x=value, fill=component))+
  geom_col(colour="gray20")+
  facet_grid(facet_dummy~.,scales="free_y",space = "free")+
  theme_classic(16)+
  theme(strip.text.y = element_text(angle=360),
        strip.background = element_blank(),
        legend.position = "bottom")
ggsave("TitanicV1.png")

ggplot(Vdf, aes(y=variable, x=value, fill=component))+
  geom_col(colour="gray20")+
  facet_wrap(facet_dummy~.,scales="free")+
  theme_classic(16)+
  theme(strip.text.y = element_text(angle=360),
        strip.background = element_blank(),
        legend.position = "bottom")
ggsave("TitanicV2.png")

ggplot(Vdf, aes(y=variable, x=log(value), fill=component))+
  geom_col(colour="gray20")+
  geom_vline(xintercept = 0)+
  facet_grid(facet_dummy~.,scales="free_y",space = "free")+
  theme_classic(16)+
  theme(strip.text.y = element_text(angle=360),
        strip.background = element_blank(),
        legend.position = "bottom")
ggsave("TitanicV3.png")

ggplot(Vdf, aes(y=variable, x=value, fill=component))+
  geom_col(colour="gray20", position = "fill")+
  facet_grid(facet_dummy~.,scales="free_y",space = "free")+
  scale_x_continuous(labels=scales::percent)+
  theme_classic(16)+
  theme(strip.text.y = element_text(angle=360),
        strip.background = element_blank())
ggsave("TitanicV4.png")
