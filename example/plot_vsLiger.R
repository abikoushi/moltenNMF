library(ggplot2)
library(dplyr)
library(bench)
load("vsliger_minibatch.Rdata")
logmu = c(-5, -2, 0)
n_rows <- c(1000, 2000, 5000)
n_cols <- c(2000)
L_rank <- c(10L, 15L, 30L)
repl <- 10 
df <- expand.grid(n_rows=n_rows,n_cols=n_cols, logmu=logmu, L_rank=L_rank)
head(df)


head(res_bm)
t1 = sapply(res_bm, function(x)x$elapsed[1])
t2 = sapply(res_bm, function(x)x$elapsed[2])
plot(t1, type="l",col=1)
lines(t2, col=2)


t1 = sapply(res_bm, function(x)x$mem_used[1])
t2 = sapply(res_bm, function(x)x$mem_used[2])

plot(t1, type="l",col=1)
lines(t2, col=2)


s1 = sapply(res_bm, function(x)x$sparsity[1])
s2 = sapply(res_bm, function(x)x$sparsity[2])
s3 = sapply(res_bm, function(x)x$sparsity[3])
s4 = sapply(res_bm, function(x)x$sparsity[4])

plot(sparsity)

df = mutate(df, con=row_number())
dfcor= reshape2::melt(simplify2array(res_cor)) %>% 
  rename(con=L1) %>% 
  left_join(df) %>% 
  mutate(Var2=c("proposed","liger")[Var2])



p = ggplot(dfcor,aes(x=Var1, y=value, linetype=factor(Var2), colour=factor(Var2), group=interaction(Var2,con)))+
  facet_grid(n_rows~L_rank+logmu, scales="free_x", labeller = label_both)+
  geom_line()+
  theme_classic(16)+labs(colour="method", linetype="method")+
  theme(strip.text.y = element_text(angle=0))

print(p)
ggsave(plot = p, filename = "vsliger.pdf", width = 18, height=6)
