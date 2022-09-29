#devtools::install_github('EESI/themetagenomics')
library(themetagenomics)
library(tidyverse)
library(parallel)
library(tidytext)
library(moltenNMF)

dat_meta <- DAVID$META[DAVID$META$Site=="UBERON:feces",]
dat_abund <- DAVID$ABUND[DAVID$META$Site=="UBERON:feces",]
range(dat_meta$Day,na.rm = TRUE)
species <- ifelse(is.na(DAVID$TAX[,7]),DAVID$TAX[,6],DAVID$TAX[,7])
for(i in 5:1){species <- ifelse(is.na(species),DAVID$TAX[,i],species)}
species <- ifelse(is.na(species),"unknown",species)
genus <- ifelse(is.na(DAVID$TAX[,6]),DAVID$TAX[,5],DAVID$TAX[,6])
for(i in 4:1){genus <- ifelse(is.na(genus),DAVID$TAX[,i],genus)}
genus <- ifelse(is.na(genus),"unknown",genus)

dat_abund_g <- data.frame(t(dat_abund)) %>% 
  mutate(species=species) %>% 
  group_by(species) %>% 
  summarise_all(sum) %>% 
  pivot_longer(-species, "ID", values_to = "abundance") %>% 
  left_join(dat_meta, by="ID") %>% 
  group_by(ID) %>% 
  mutate(total_read=sum(abundance)) %>% 
  ungroup() %>% 
  dplyr::filter(!is.na(Day)) %>% 
  group_by(species) %>% 
  dplyr::filter(sum(abundance)!=0) %>% 
  group_by(Day, Donor) %>% 
  dplyr::filter(sum(abundance)!=0) %>% 
  ungroup() %>% 
  mutate(Day_Donor=paste(Day, Donor, sep="_"))

head(dat_abund_g)

set.seed(111)
f1 <- abundance ~ species+Day_Donor-1
# system.time({
#   out <- mclapply(2:9,function(l){mrNMF_vb(abundance~species+Day:Donor-1, data = dat_abund_g, L=l, iter=500, a=0.5, b=1)},mc.cores = detectCores())
# })
#   user  system elapsed 
#529.364 215.448 521.039 
# ELBOs <- sapply(out, function(x)x$ELBO[length(x$ELBO)])
# plot(2:9, ELBOs,type="l")
# 
# outfix <- out[[which.max(ELBOs)]]
# L <- which.max(ELBOs)+1

set.seed(111)
f1 <- abundance ~ species+Day_Donor-1
outfix <- mNMF_vb(f1, data = dat_abund_g, L=9, iter=1000, a=0.5, b=0.1)
plot(outfix$ELBO, type = "l")

X <- sparse_model_matrix_b(f1,data = dat_abund_g)
#ytilde <- rpredictor_mrNMF(X, 2000, out$shape, out$rate, out$precision)
ytilde <- rpredictor_mNMF(X, 2000, outfix$shape, outfix$rate)

ybar <- rowMeans(ytilde)
yq <- apply(ytilde, 1, quantile, prob = c(0.05, 0.995))
head(t(yq))

dffit <- data.frame(obs=dat_abund_g$abundance,
           fit=ybar,
           ymin=yq[1,],
           ymax=yq[2,])
head(dffit)

ggplot(dffit, aes(x=obs, y=fit))+
  geom_point(alpha=0.1)+
  #geom_ribbon(aes(ymin=ymin, ymax=ymax), colour="royalblue", alpha=0.1)+
  geom_abline(intercept=0, slope=1, linetype=2) +
  theme_minimal(16)

ggsave("goodnessoffit_david.png")

Vdf <- data.frame(outfix$shape/outfix$rate,variable=rownames(outfix$shape),
                  facet_dummy=outfix$vargroup) %>% 
  pivot_longer(!(variable|facet_dummy), names_to = "component", 
               names_transform = list(component = readr::parse_number)) %>% 
  mutate(component=factor(component,ordered = TRUE))
head(Vdf)
Daydf <- dplyr::filter(Vdf,"Day_Donor"== facet_dummy) %>% 
  separate(variable, c("Day","Donor"),sep=":") %>% 
  mutate(Day = as.integer(gsub("^Day_Donor|_2202$","",Day)))
tail(Daydf)

dplyr::filter(Vdf,"Donor"== facet_dummy)

ggplot(Daydf, aes(x=Day, y=value,  colour=component))+
  geom_line()+
  facet_grid(Donor~.)+
  theme_classic(16)+
  theme(strip.background = element_blank())
ggsave("david_day_line.pdf")

ggplot(Daydf, aes(x=Day, y=value, fill=component))+
  stat_summary_bin(geom="col", fun = mean, binwidth = 10, position = "fill")+
  facet_grid(Donor~.)+
  theme_classic(16)+
  theme(strip.background = element_blank())
ggsave("david_day_bar.pdf")

genusdf <- dplyr::filter(Vdf, "species"==facet_dummy) %>% 
  dplyr::mutate(variable=sub("species","",variable)) %>% 
  rename(species=variable) %>% 
  group_by(component) %>% 
  top_n(5,value) %>% 
  ungroup()

ggplot(genusdf,aes(y=reorder_within(species,value,component), x=value))+
  geom_point()+
  geom_linerange(aes(xmin=0, xmax=value))+
  facet_wrap(component~., scales = "free")+
  scale_x_continuous(n.breaks = 3)+
  scale_y_reordered()+
  labs(y="")+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black", size=12))

ggsave("david_genus.pdf")
