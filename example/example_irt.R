# library(readr)
# library(tidyr)
# library(Matrix)
# library(dplyr)
# library(ggrepel)
# library(ggplot2)
# 
# dat <- read_csv("https://raw.githubusercontent.com/MatsuuraKentaro/RStanBook/master/chap10/input/data-shogi-player.txt",
#                 col_types = cols(Loser ="i", Winner = "i"))
# kishi_name <- read_table("https://raw.githubusercontent.com/MatsuuraKentaro/RStanBook/master/chap10/input/data-shogi-player-name.txt",
#                          col_types = cols(kid ="i", kname="c", nid = "i"))
# dat <- mutate(dat, Winner=factor(Winner), Loser=factor(Loser))
# X1 <- sparse_model_matrix_b(~Loser-1, data=dat)
# X2 <- sparse_model_matrix_b(~Winner-1, data=dat)
# X <- rbind(X1,X2)
# attr(X,"indices") <- attr(X1,"indices")
# y <- rep(0:1, each=nrow(X1))
# set.seed(777)
# system.time({
#   out3 <- mBMF_mcvb(y,X,L=2,iter=500)
# })
# #1.467   0.043   1.531
# plot(out3$loglik,type = "l")
# dfmu <- data.frame(out3$mean) %>% 
#   mutate(nid=as.integer(row_number())) %>% 
#   left_join(kishi_name)
# 
# win <- tally(group_by(dat,Winner))
# lose <- tally(group_by(dat,Loser))
# dft <- data.frame(nid=win$Winner,win=win$n,lose=lose$n,n=win$n+lose$n) %>% 
#   mutate(dft,p=win/n)
# 
# hist(dft$p)
# 
# w5 <- top_n(dft, 5, p)
# l5 <- top_n(dft, 5, -p)
# 
# ggplot(dfmu,aes(x=X1,y=X2))+
#   geom_point(alpha=0.2)+
#   geom_point(data=dfmu[l5$nid,], colour="darkred")+
#   geom_point(data=dfmu[w5$nid,], colour="darkblue")+
#   geom_text_repel(data=dfmu[w5$nid,],aes(label=kname), family="Osaka", colour="darkblue")+
#   geom_text_repel(data=dfmu[l5$nid,],aes(label=kname), family="Osaka", colour="darkred")+
#   theme_classic(16)
# 
# 
