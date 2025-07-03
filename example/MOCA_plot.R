library(moltenNMF)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(umap)
library(tsne)

# poisloss <- function(obs, fit){
#   -mean(obs*log(fit)-fit-lgamma(obs+1)) 
# }

# rmse <- function(obs, fit){
#   sqrt(mean((as.matrix(obs) - fit)^2))
# }

fit_mnmf <- function(out){
  Vhat_t <- moltenNMF:::meanV_array(out)
  Vhat_t[[1]]%*%t(Vhat_t[[2]])
}

fit_liger <- function(liger_obj){
  W = getMatrix(liger_obj, slot = "W")
  V = getMatrix(liger_obj, slot = "V")
  H = getMatrix(liger_obj, slot = "H")
  as.matrix((W + V[[1]])%*%H[[1]])
}

extract_row_liger <- function(liger_obj){
  W = getMatrix(liger_obj, slot = "W")
  V = getMatrix(liger_obj, slot = "V")
  as.matrix(W + V[[1]])
}
####
#browseURL("https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads")

load("resMOCA.Rdata")
path <- scan("datapath.txt", what = character())
cells = read_csv(path[3])

cellname = colnames(readRDS(path[1]))
head(cellname)

head(cells$sample)
#n_distinct(cells$development_stage)
#moltenNMF:::matplot2(V = t(Vhat[[1]]), lwd = 0.005, colour = factor(cells$development_stage))

Vhat = moltenNMF:::meanV_array(m_obj)
Vhat = moltenNMF:::rearrange_cols(Vhat, axis = 2, normalize = FALSE, decreasing = TRUE)

resU = umap(log(Vhat[[2]]))
saveRDS(resU, file = "resU_log.rds")
plot(resU$layout, pch=".")

# resT = tsne(Vhat[[2]])

#colnames(cells)
# [1] "sample"                             "all_exon_count"                     "all_intron_count"                  
# [4] "all_read_count"                     "intergenic_rate"                    "embryo_id"                         
# [7] "embryo_sex"                         "nuclei_extraction_date"             "development_stage"                 
# [10] "Total_mRNAs"                        "num_genes_expressed"                "Size_Factor"                       
# [13] "Main_Cluster"                       "Main_cluster_tsne_1"                "Main_cluster_tsne_2"               
# [16] "Sub_cluster"                        "Sub_cluster_tsne_1"                 "Sub_cluster_tsne_2"                
# [19] "doublet_score"                      "detected_doublet"                   "doublet_cluster"                   
# [22] "sub_cluster_id"                     "Main_cell_type"                     "Main_trajectory"                   
# [25] "Main_trajectory_umap_1"             "Main_trajectory_umap_2"             "Main_trajectory_umap_3"            
# [28] "Main_trajectory_refined_by_cluster" "Main_trajectory_refined_umap_1"     "Main_trajectory_refined_umap_2"    
# [31] "Main_trajectory_refined_umap_3"     "Sub_trajectory_name"                "Sub_trajectory_umap_1"             
# [34] "Sub_trajectory_umap_2"              "Sub_trajectory_louvain_component"   "Sub_trajectory_Pseudotime"         



Vcell = mutate(data.frame(resU$layout), sample = cellname) %>% 
  left_join(cells)

head(Vcell$Main_cell_type)

ggplot(Vcell, aes(x=Main_trajectory_umap_1, y=Main_trajectory_umap_2, colour=Main_cell_type)) +
  geom_point(size=0.01, alpha=0.1)+
  guides(colour=guide_legend(override.aes = list(size=2, alpha=1)))+
  theme_bw(16)

ggsave(filename = "umap_original.png")

ggplot(Vcell, aes(x=X1, y=X2, colour=Main_trajectory)) +
  geom_point(size=0.01, alpha=0.1)+
  guides(colour=guide_legend(override.aes = list(size=2, alpha=1)))+
  theme_bw(16)
ggsave(filename = "umap_log.png")

Vhat = moltenNMF:::meanV_array(m_obj)
Vhat = moltenNMF:::rearrange_cols(Vhat, axis = 2, normalize = TRUE, decreasing = TRUE)
Vcell =  mutate(data.frame(Vhat[[2]]),sample = cellname) %>% 
  left_join(cells) %>% 
  pivot_longer(X1:X30, names_to = "factor") %>% 
  mutate(factor = as.integer(gsub("X", "",factor)))


head(Vcell)
p = ggplot(Vcell, aes(x=factor, y=value, group=sample))+
  geom_line(aes(colour=development_stage), linewidth=0.01)+
  guides(colour=guide_legend(override.aes = list(linewidth=1)))+
  scale_color_viridis_c()+
  theme_bw(base_size = 16)
print(p)
#ggsave(plot = p, filename = "MOCA_Vcell_stage_log.png")

####
path <- scan("datapath.txt", what = character())
genes = read_csv(path[4])
Vhat = moltenNMF:::meanV_array(m_obj)
Vhat = moltenNMF:::rearrange_cols(Vhat, axis = 2, normalize = TRUE, decreasing = TRUE)


Vgene = mutate(data.frame(Vhat[[1]]), gene_short_name = genes$gene_short_name) %>% 
  pivot_longer(X1:X30, names_to = "factor") %>% 
  mutate(factor = as.integer(gsub("X", "",factor)))

top10 = slice_max(Vgene, by = factor, n=10, order_by = value)
p = ggplot(top10, aes(x=reorder(gene_short_name,value), y=value))+
  geom_point()+
  geom_linerange(aes(ymin=0, ymax=value))+
  facet_wrap(~factor, scales="free_y")+
  theme_bw(16)+
  coord_flip()

ggsave(plot = p, filename = "top10gene.pdf", width = 20, height = 15)

###
