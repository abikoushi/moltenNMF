library(rliger)
library(bench)
library(xtable)
library(dplyr)
library(knitr)
library(moltenNMF)


load("resMOCA.Rdata")
##
path <- scan("datapath.txt", what = character())
tpath <- path[2]
Vhat = moltenNMF:::meanV_array(m_obj)

system.time({
  fit_p = moltenNMF:::obsfitloss_2d_mtx(readtxt = tpath, V1 = Vhat[[1]], V2 = Vhat[[2]], n_header = 2)  
})

# ユーザ   システム       経過  
# 840.76       7.30     871.10 

saveRDS(fit_p, file = "fit_p.rds")

load("resMOCA_liger.Rdata")

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

V1 = extract_col_liger(liger_obj)
V2 = extract_col_liger(liger_obj)

system.time({
  fit_l = moltenNMF:::obsfitloss_2d_mtx(readtxt = tpath, V1 = V1, V2 = V2, n_header = 2)  
})
#ユーザ   システム       経過  
#27560.62      65.27   27985.14 
saveRDS(fit_l, file = "fit_l0.rds")

fit_l2 = readRDS("fit_l.rds")

fit_l$MSE > fit_l2$MSE

kable(log10(rbind(simplify2array(fit_p),simplify2array(fit_l))), format = "latex")
###

load("resMOCA_liger.Rdata")
load("resMOCA_liger_2.Rdata")

tab1 = bind_rows(bm2, bm_l1, bm_l2) %>% 
  mutate(method=c("poposed","liger","liger(2)")) %>% 
  dplyr::select(method,mem_alloc,total_time)

tab1
print(xtable(tab1), include.rownames = FALSE)
cat(kable(tab1, format = "latex"), file="bench.txt")


