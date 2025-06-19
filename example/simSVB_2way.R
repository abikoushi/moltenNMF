library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
library(dplyr)
library(gridExtra)

# 2-way
set_data_mf <- function(L, nrow, ncol, mu=0){
  W <- matrix(rnorm(nrow*L,0,1), ncol=L)
  H <- matrix(rnorm(L*ncol,0,1), nrow=L)
  W <- sweep(W,1,rowMeans(W)-mu)
  H <- sweep(H,2,colMeans(H)-mu)
  list(lambda=exp(W)%*%exp(H), trueW=exp(W), trueH=exp(H))
}

rearrange_winner_ord <- function(V,V_s){
  cmat = cor(V, V_s)
  ord =integer(ncol(V))
  ord[1] = which.max(cmat[1,])
  for(i in 2:L){
    ord[i] = which(cmat[i,]==max(cmat[i,-ord[1:(i-1)]]))  
  }
  list(V=V_s[,ord], cor=diag(cmat[,ord]))
}

simfunc <- function(seed, V, L, lambda, settings){
  set.seed(seed)
  Y <- matrix(rpois(length(lambda),lambda),
              nrow(lambda), ncol(lambda))
  Y <- as(Y, "TsparseMatrix")

  out_d <- moltenNMF:::NMF2D_vb(Y, rank = L, 
                                prior_shape = 1, prior_rate = 1, iter=1000,
                                display_progress = FALSE)
  V_d = moltenNMF:::meanV_array(out_d)
  resV_d = rearrange_winner_ord(rbind(V_d[[1]], V_d[[2]]), V)
  

  Vcor = vector("list", nrow(settings))
  for(set in seq_len(nrow(settings))){
    lr_param = c(settings$delay[set],
                 settings$forgetting[set])
    out_s <- moltenNMF:::NMF2D_svb(Y, rank = L,
                                   n_epochs = 200,
                                   n_baches  = settings$n_batches[set],
                                   lr_param = lr_param,
                                   lr_type = "exponential",
                                   display_progress = FALSE)
    V_s = moltenNMF:::meanV_array(out_s)
    resV_s = rearrange_winner_ord(rbind(V_s[[1]],V_s[[2]]),
                                  rbind(param$trueW,t(param$trueH)))
    
    Vcor[[set]] = resV_s$cor
  }
  list(svb=Vcor, bvb=resV_d$cor)
}


settings = expand.grid(forgetting = c(0.7,0.8,0.9),
                       delay = c(1.5,5,15),
                       n_batches = c(500,1000,2000))
#dim(settings)
#[1] 27  3
L <- 5L
ncols = c(100, 500, 1000)
for(i in 1:3){
  param <- set_data_mf(L, 100, ncols[i])
  ressim = simfunc(1, rbind(param$trueW,t(param$trueH)), L, param$lambda, settings)
  ressimdf = reshape2::melt(simplify2array(ressim$svb),
                            varnames = c("component","setid")) %>% 
    left_join(mutate(settings, setid=row_number()))
  saveRDS(ressimdf, paste0("simnmf_",ncols[i],".rds"))  
}

ressimdf1 <- readRDS("simnmf_100.rds")
ressimdf2 <- readRDS("simnmf_500.rds")
ressimdf3 <- readRDS("simnmf_1000.rds")

p1 = ggplot(ressimdf1, aes(x=n_batches, y=value, 
                           group=factor(component),
                           colour=factor(component)))+
  geom_line()+
  facet_grid(forgetting~delay, labeller = label_both)+
  scale_color_viridis_d()+
  labs(title  = paste("number of columns:", ncols[1]), colour="component")+
  theme_classic(20)+
print(p1)

p2 = ggplot(ressimdf2, aes(x=n_batches, y=value,
                           group=factor(component), colour=factor(component)))+
  geom_line()+
  facet_grid(forgetting~delay, labeller = label_both)+
  scale_color_viridis_d()+
  labs(title  = paste("number of columns:", ncols[2]), colour="component")+
  theme_classic(20)
print(p2)

p3 = ggplot(ressimdf3, aes(x=n_batches, y=value,
                           group=factor(component), colour=factor(component)))+
  geom_line()+
  facet_grid(forgetting~delay, labeller = label_both)+
  scale_color_viridis_d()+
  labs(title  = paste("number of columns:", ncols[3]), colour="component")+
  theme_classic(20)
print(p3)
pp = gridExtra::grid.arrange(p1, p2, p3, nrow=3)
ggsave("simNMF.pdf", plot = pp, width = 20, height = 20)
