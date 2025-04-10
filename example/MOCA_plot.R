library(moltenNMF)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# poisloss <- function(obs, fit){
#   -mean(obs*log(fit)-fit-lgamma(obs+1)) 
# }
# 
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
Vhat = moltenNMF:::rearrange_cols(Vhat, axis = 2, normalize = TRUE, decreasing = TRUE)

Vcell = mutate(data.frame(Vhat[[2]]),sample = cellname) %>% 
  left_join(cells) %>% 
  pivot_longer(X1:X30, names_to = "factor") %>% 
  mutate(factor = as.integer(gsub("X", "",factor)))

head(Vcell)
p = ggplot(Vcell, aes(x=factor, y=value, group=sample))+
  geom_line(aes(colour=Main_cell_type), linewidth=0.01)+
  guides(colour=guide_legend(override.aes = list(linewidth=1)))+
  theme_bw()
ggsave(plot = p, filename = "MOCA_Vcell.png")


