library(readxl)
library(tidyverse)
library(parallel)
library(Matrix)
library(patchwork)
library(tidytext)
library(patchwork)
#library(Rcpp)
#sourceCpp("~/Dropbox/moltenNMF/moltenNMF.cpp")
tab0 <- read_tsv("~/data/Kostic/diabimmune_t1d_16s_otu_table.txt",skip=1)
abundance <- tab0 %>% 
  dplyr::select(starts_with("G")) %>% 
  mutate_all(as.numeric) %>% t
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
  mutate(Age2 = paste(age,Case_Control,sep = "_"))

head(datall)

f1 <- abundance ~ genus+Subject_ID+Age2-1
sq <- seq(10,90,by=10)
set.seed(1)
system.time({
  out <- mclapply(sq,function(l){
    mrNMF_vb(f1, data = datall, L=l, iter=500, a=0.5, b=1)},
                  mc.cores = detectCores())
})
# 437s
#saveRDS(out,file = "moltenNMF_Kostic_10_90.rds")
ELBOs <- sapply(out, function(x)x$ELBO[length(x$ELBO)])
# dim(Y)
# setEPS()
# postscript("~/Desktop/ELBOplot_Kostik.eps")
plot(sq, ELBOs, type="b", 
     xlab = "number of topics", ylab = "ELBO")
# dev.off()
wch <- which.max(ELBOs)
outfix <- out[[2]]
#L <- which.max(ELBOs)+1
plot(outfix$ELBO, type="l")

X <- sparse_model_matrix_b(f1, datall)
system.time({
  ytilde <- rpredictor_mrNMF(X, 1000, outfix$shape, outfix$rate, outfix$precision)
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
               names_transform = list(component = readr::parse_number)) %>%
  mutate(component=factor(component,ordered = TRUE))

Agedf <- dplyr::filter(Vdf,"Age2" == facet_dummy) %>% 
  mutate(variable = gsub("Age2","",variable)) %>% 
  separate(variable, c("Age", "Group"), sep = "_") %>% 
  mutate(Age = as.integer(Age))

ggplot(Agedf,aes(x=Age, y=value, fill=component, group=component))+
  stat_summary_bin(geom="col", fun = mean, binwidth = 25, position = "fill")+
  facet_grid(Group~.)+
  theme_classic(14)

ggplot(Agedf,aes(x=Age, y = value, colour = Group))+
  geom_point(alpha=0.2)+
  stat_smooth(geom = "line")+
  scale_colour_manual(values=c("orange","royalblue"))+
  facet_wrap(component ~., scales = "free_y")+
  theme_classic(14)+
  theme(strip.background = element_blank())

ggsave("Kostic_age.png")
#ggsave(filename = "~/Desktop/Kostic_age.png")

head(Vdf)
genusdf <- dplyr::filter(Vdf, "genus"==facet_dummy) %>% 
  mutate(genus=gsub("genusg__","",variable))


genusdf2 <-group_by(genusdf, component) %>% 
  top_n(3,value)

genusdf2

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

ggsave("Kostic_genus.png")
