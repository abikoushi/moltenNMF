library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
library(dplyr)
library(gridExtra)

ressimdf1 <- readRDS("simnmf_100.rds")
ressimdf2 <- readRDS("simnmf_500.rds")
ressimdf3 <- readRDS("simnmf_1000.rds")

settings = expand.grid(forgetting = c(0.7,0.8,0.9),
                       delay = c(1.5,5,15),
                       n_batches = c(500,1000,2000),
                       rep = 1:10)

ressimdf = reshape2::melt(simplify2array(ressimdf1$svb),
                          varnames = c("component","setid")) %>% 
  left_join(mutate(settings, setid=row_number()))

ncols = c(100, 500, 1000)

dfb = data.frame(component=1:5, value=ressimdf1$bvb)

p1 = ggplot(ressimdf, aes(x=n_batches, y=value, 
                           group=factor(component),
                           colour=factor(component),
                          shape=factor(component)))+
  stat_summary(geom = "line", fun = mean)+
  geom_point()+
  geom_hline(data = dfb, aes(yintercept = value), linetype = 3, linewidth=0.5)+
  facet_grid(forgetting~delay, labeller = label_both)+
  scale_color_viridis_d()+
  labs(title  = paste("number of columns:", ncols[1]),
       colour="component", shape="component")+
  theme_classic(20)

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
