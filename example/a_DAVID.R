#devtools::install_github('EESI/themetagenomics')
library(themetagenomics)
library(tidyverse)
library(parallel)
library(tidytext)
library(moltenNMF)

dat_meta <- DAVID$META[DAVID$META$Site=="UBERON:feces",]
dat_abund <- DAVID$ABUND[DAVID$META$Site=="UBERON:feces",]
range(dat_meta$Day,na.rm = TRUE)
species <- ifelse(is.na(DAVID$TAX[,7]), DAVID$TAX[,6], DAVID$TAX[,7])
for(i in 5:1){species <- ifelse(is.na(species), DAVID$TAX[,i], species)}
species <- ifelse(is.na(species), "unknown", species)
genus <- ifelse(is.na(DAVID$TAX[,6]), DAVID$TAX[,5], DAVID$TAX[,6])
for(i in 4:1){genus <- ifelse(is.na(genus), DAVID$TAX[,i], genus)}
genus <- ifelse(is.na(genus),"unknown",genus)

datall <- data.frame(t(dat_abund)) %>% 
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
  ungroup()

head(datall$Donor)
datB <- dplyr::filter(datall, Donor=="2202:DonorB")

set.seed(1111)
f1 <- abundance ~ species+factor(Day)-1
sq <- seq(2,30,by=4)
system.time({
  out <- mclapply(sq,function(l){
    mNMF_vb(f1, data = datB, L=l, iter=1000, a=0.5, b=1)},
    mc.cores = detectCores())
})
# user   system  elapsed 
#600.706  28.097 390.630  
#saveRDS(out, "David_B_2_30.rds")
ELBOs <- sapply(out, function(x)x$ELBO[length(x$ELBO)])
plot(sq, ELBOs, type="b")
#
sq
outfix <- out[[which.max(ELBOs)]]
X <- sparse_model_matrix_b(f1,data = datB)
ytilde <- rpredictor_mNMF(X, 200, outfix$shape, outfix$rate)
ybar <- rowMeans(ytilde)
yq <- apply(ytilde, 1, quantile, prob = c(0.05, 0.995))
# head(t(yq))

dffit <- data.frame(obs=datB$abundance,
           fit=ybar,
           ymin=yq[1,],
           ymax=yq[2,])
head(dffit)

ggplot(dffit, aes(x=obs, y=fit))+
  geom_point(alpha=0.1)+
  geom_abline(intercept=0, slope=1, linetype=2) +
  theme_minimal(16)

ggsave("goodnessoffit_david_B.png")

Vdf <- data.frame(outfix$shape/outfix$rate,variable=rownames(outfix$shape),
                  facet_dummy=outfix$vargroup) %>% 
  pivot_longer(!(variable|facet_dummy), names_to = "component", 
               names_transform = list(component = readr::parse_number)) %>% 
  mutate(component=factor(component,ordered = TRUE))
head(Vdf$facet_dummy)
Daydf <- dplyr::filter(Vdf,"Day_Donor"== facet_dummy) %>% 
  separate(variable, c("Day","Donor"),sep=":") %>% 
  mutate(Day = as.integer(gsub("^Day_Donor|_2202$","",Day)))
tail(Daydf)
Daydf <- dplyr::filter(Vdf,"factor(Day)"== facet_dummy) %>% 
  mutate(Day = as.integer(gsub("^factor\\(Day\\)", "", variable)))

ggplot(Daydf, aes(x=Day, y=value))+
  geom_line()+
  facet_wrap(component~.)+
  theme_classic()+
  theme(strip.background = element_blank())
ggsave("david_day_line_B.pdf")

ggplot(Daydf, aes(x=Day, y=value))+
  stat_summary_bin(geom="col", fun = mean, binwidth = 10)+
  facet_wrap(component~.,scales = "free_y")+
  theme_classic()+
  theme(strip.background = element_blank())

# ggplot(Daydf, aes(x=Day, y=value, fill=component))+
#   stat_summary_bin(geom="col", fun = mean, binwidth = 10, position = "fill")+
#   #facet_grid(Donor~.)+
#   theme_classic(16)+
#   theme(strip.background = element_blank())
ggsave("david_day_bar_B.pdf")

genusdf <- dplyr::filter(Vdf, "species"==facet_dummy) %>% 
  dplyr::mutate(variable=sub("species","",variable)) %>% 
  rename(species=variable) %>% 
  group_by(component) %>% 
  mutate(value=value/sum(value)) %>% 
  top_n(7, value) %>% 
  ungroup()

genusdf2 <- dplyr::filter(genusdf,component %in% c(10,12,27))
genusdf2 <- dplyr::filter(genusdf,component %in% c(1,13,19))
genusdf2 <- dplyr::filter(genusdf,component %in% c(24,25,28))

ggplot(genusdf2, aes(y=reorder_within(species,value,component), x=value))+
  geom_point()+
  geom_linerange(aes(xmin=0, xmax=value))+
  facet_wrap(component~., scales = "free_y")+
  scale_x_continuous(n.breaks = 5)+
  scale_y_reordered()+
  labs(y="", x="relative value")+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black", size=11),
        strip.background = element_blank(),
        strip.text = element_text(size=11))

#ggsave("david_genus_B.pdf")
