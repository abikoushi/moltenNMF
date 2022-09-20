#devtools::install_github("abikoushi/moltenNMF", build_vignettes = TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)
library(moltenNMF)

Titanicdf <- as.data.frame(Titanic) %>% 
  mutate(Class=factor(Class,levels=c("3rd","2nd","1st","Crew")))

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

set.seed(1987)
out <- mrNMF_vb(Freq~Class+Sex+Age+Survived-1, data=Titanicdf, L=2, iter=1000, a=0.5, b=0.1)

qplot(1:length(out$ELBO),out$ELBO, geom = "line")+
  labs(x="iter",y="ELBO")+
  theme_minimal(16)

V <- out$shape/out$rate
yhat <- rowSums(prod_mNMF(~Class+Sex+Age+Survived-1, data=Titanicdf, V))
zhat <- (Titanicdf$Freq+out$precision)/(yhat+out$precision)

qplot(yhat*zhat, Titanicdf$Freq, alpha=I(0.5), size=I(3))+
  labs(x="fit",y="obs")+
  geom_abline(intercept=0,slope=1,linetype=2)+
  theme_minimal(16)

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

ggplot(Vdf, aes(y=variable, x=log(value), fill=component))+
  geom_col(colour="gray20")+
  geom_vline(xintercept = 0)+
  facet_grid(facet_dummy~.,scales="free_y",space = "free")+
  theme_classic(16)+
  theme(strip.text.y = element_text(angle=360),
        strip.background = element_blank(),
        legend.position = "bottom")
ggsave("TitanicV2.png")

ggplot(Vdf, aes(y=variable, x=value, fill=component))+
  #geom_vline(xintercept = 0.5, linetype=2)+
  geom_col(colour="gray20", position = "fill")+
  facet_grid(facet_dummy~.,scales="free_y",space = "free")+
  scale_x_continuous(labels=scales::percent)+
  theme_classic(16)+
  theme(strip.text.y = element_text(angle=360),
        strip.background = element_blank())
ggsave("TitanicV3.png")
