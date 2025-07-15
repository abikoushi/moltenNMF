library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(moltenNMF)
out = readRDS("out_L20.rds")
#plot(out$ELBO[-1], type="l")
V = out$shape/out$rate
head(V)
Vcell = grepV(V, "cell")

ord = order(apply(Vcell, 2, var), decreasing = TRUE)

df = data.frame(log(Vcell[,ord])) %>% 
  rownames_to_column(var = "rowname") %>% 
  tidyr::separate(col = rowname, into = c("variable","sample", "barcode"), sep = "_")

df = mutate(df, species=if_else(grepl("^H",sample),"Human", sample)) %>% 
  mutate(species=if_else(grepl("^Pig",sample),"Pig", species)) %>% 
  mutate(species=if_else(grepl("^Mouse",sample),"Mouse", species)) %>% 
  mutate(species=if_else(grepl("^MM",sample),"MR", species)) %>% 
  mutate(species=if_else(grepl("^Maca",sample),"MC", species))

ggplot(df, aes(x=X1, y=X2, colour=species))+
  geom_point(alpha = 0.1) +
  guides(colour=guide_legend(override.aes = list(alpha=1, size=2))) +
  theme_bw()

ggsave(filename = "test_species.png")
