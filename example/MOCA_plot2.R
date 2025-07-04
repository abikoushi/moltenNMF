library(moltenNMF)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(umap)
library(tsne)
library(rliger)
library(ggrepel)
library(bench)

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
load("resMOCA_liger_2.Rdata")
path <- scan("datapath.txt", what = character())
cells = read_csv(path[3])

cellname = colnames(readRDS(path[1]))
head(cellname)

Vhat = extract_col_liger(liger_obj)
dim(Vhat)

resU = umap(Vhat)
saveRDS(resU, file = "resU_liger2.rds")

resU = readRDS(file = "resU_log.rds")
resU_liger = readRDS(file = "resU_liger2.rds")

Vcell = mutate(data.frame(resU$layout), sample = cellname) %>% 
  left_join(cells)

head(Vcell)
Vcell_sum = group_by(Vcell, Main_cell_type) %>% 
  summarise(X1 = mean(X1), X2 = mean(X2))
Vcell_sum

ggplot(Vcell, aes(x=X1, y=X2, colour=Main_cell_type)) +
  geom_point(size=0.01, alpha=0.1, show.legend = FALSE)+
  geom_label_repel(data = Vcell_sum, aes(label=Main_cell_type), show.legend = FALSE, max.overlaps = 16)+
  guides(colour=guide_legend(override.aes = list(size=2, alpha=1)))+
  theme_bw(16)+labs(x = "UMAP1", y = "UMAP2", title = "proposed")

ggsave(filename = "umap_prop.png", width = 9, height = 9)

Vcell = mutate(data.frame(resU_liger$layout), sample = cellname) %>% 
  left_join(cells)
Vcell_sum = group_by(Vcell, Main_cell_type) %>% 
  summarise(X1 = mean(X1), X2 = mean(X2))
Vcell_sum

p = ggplot(Vcell, aes(x=X1, y=X2, colour=Main_cell_type)) +
  geom_point(size=0.01, alpha=0.1)+
  #geom_label_repel(data = Vcell_sum, aes(label=Main_cell_type), show.legend = FALSE, max.overlaps = 200)+
  guides(colour=guide_legend(override.aes = list(size=2, alpha=1)))+
  theme_bw(16)+labs(x = "UMAP1", y = "UMAP2", title = "liger")

ggsave(p, filename = "umap_liger.png", width = 9, height = 9)


Vcell_sum = group_by(Vcell, Main_cell_type) %>% 
  summarise(Main_trajectory_umap_1 = mean(Main_trajectory_umap_1),
            Main_trajectory_umap_2 = mean(Main_trajectory_umap_2))
Vcell_sum

p = ggplot(Vcell, aes(x=Main_trajectory_umap_1, y=Main_trajectory_umap_2, colour=Main_cell_type)) +
  geom_point(size=0.01, alpha=0.1, show.legend = FALSE)+
  geom_label_repel(data = Vcell_sum, aes(label=Main_cell_type), show.legend = FALSE, max.overlaps = 16)+
  guides(colour=guide_legend(override.aes = list(size=2, alpha=1)))+
  theme_bw(16)+labs(x = "UMAP1", y = "UMAP2", title = "origin")

ggsave(p, filename = "umap_origin.png", width = 9, height = 9)

######
load("resMOCA.Rdata")
load("resMOCA_liger_2.Rdata")

Vhat = extract_col_liger(liger_obj)
Vhat = sweep(Vhat, 1, rowSums(Vhat), "/")
Vhat = Vhat[,order(apply(Vhat,2,var), decreasing = TRUE)]
colnames(Vhat) <- NULL 
Vcell =  mutate(data.frame(Vhat),sample = cellname) %>% 
  left_join(cells) %>% 
  pivot_longer(X1:X30, names_to = "factor") %>% 
  mutate(factor = as.integer(gsub("X", "",factor)))

p = ggplot(Vcell, aes(x=factor, y=value, group=sample))+
  geom_line(aes(colour=development_stage), linewidth=0.01)+
  guides(colour=guide_legend(override.aes = list(linewidth=1)))+
  scale_color_viridis_c()+
  labs(x="component", title = "liger") + 
  scale_x_continuous(n.breaks = 6)+
  theme_bw(base_size = 16)
print(p)
ggsave(plot = p, filename = "MOCA_Vcell_stage_liger.png", width = 9, height = 9)


Vhat = moltenNMF:::meanV_array(m_obj)
Vhat = moltenNMF:::rearrange_cols(Vhat, axis = 2, normalize = TRUE, decreasing = TRUE)
Vcell =  mutate(data.frame(Vhat[[2]]),sample = cellname) %>% 
  left_join(cells) %>% 
  pivot_longer(X1:X30, names_to = "factor") %>% 
  mutate(factor = as.integer(gsub("X", "", factor)))

p = ggplot(Vcell, aes(x=factor, y=value, group=sample))+
  geom_line(aes(colour=development_stage), linewidth=0.01)+
  guides(colour=guide_legend(override.aes = list(linewidth=1)))+
  scale_color_viridis_c() + 
  labs(x="component", title = "proposed") + 
  scale_x_continuous(n.breaks = 6)+
  theme_bw(base_size = 16)
#print(p)
ggsave(plot = p, filename = "MOCA_Vcell_stage.png", width = 9, height = 9)

