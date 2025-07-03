library(moltenNMF)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(umap)
library(tsne)
library(rliger)

extract_row_liger <- function(liger_obj){
  W = getMatrix(liger_obj, slot = "W")
  V = getMatrix(liger_obj, slot = "V")
  as.matrix(W + V[[1]])
}

extract_col_liger <- function(liger_obj){
  H = getMatrix(liger_obj, slot = "H")
  H = t(H[[1]])
  as.matrix(H)
}

#load
load("resMOCA_liger.Rdata")
path <- scan("datapath.txt", what = character())
cells = read_csv(path[3])
cellname = colnames(readRDS(path[1]))
head(cellname)

Vhat = extract_col_liger(liger_obj)
dim(Vhat)

resU = umap(Vhat)
saveRDS(resU, file = "resU_liger1.rds")

Vcell = mutate(data.frame(resU$layout), sample = cellname) %>% 
  left_join(cells)

head(Vcell$Main_cell_type)

# ggplot(Vcell, aes(x=Main_trajectory_umap_1, y=Main_trajectory_umap_2, colour=Main_cell_type)) +
#   geom_point(size=0.01, alpha=0.1)+
#   guides(colour=guide_legend(override.aes = list(size=2, alpha=1)))+
#   theme_bw(16)
# 
# ggsave(filename = "umap_original.png")

ggplot(Vcell, aes(x=X1, y=X2, colour=Main_trajectory)) +
  geom_point(size=0.01, alpha=0.1)+
  guides(colour=guide_legend(override.aes = list(size=2, alpha=1)))+
  theme_bw(16)
ggsave(filename = "umap_liger1.png")
