library(ggplot2)
library(parallel)
library(moltenNMF)
library(Matrix)
library(tidyr)
library(dplyr)
library(gridExtra)

ressim1 <- readRDS("simSVB_b_100.rds")
ressim2 <- readRDS("simSVB_B_500.rds")
ressim3 <- readRDS("simSVB_b_1000.rds")
# settings = expand.grid(delay = c(0.7,0.8,0.9,1.5),
#                        forgetting = c(0.7,0.8,0.9),
#                        rep = 1:10) %>% 
#   mutate(setid=row_number())
settings = expand.grid(forgetting = c(0.6,0.7,0.8,0.9),
                       delay = c(0.7,0.8,0.9,1),
                       n_batches = c(500,1000,2000),
                       rep = 1) %>% 
  mutate(setid=row_number())


ncols = c(100, 500, 1000)

ressimdf1 = reshape2::melt(simplify2array(ressim1$svb),
                           varnames = c("component","setid")) %>% 
  left_join(settings, by="setid")

ressimdf2 = reshape2::melt(simplify2array(ressim2$svb),
                           varnames = c("component","setid"))%>% 
  left_join(settings, by="setid")

ressimdf3 = reshape2::melt(simplify2array(ressim3$svb),
                           varnames = c("component","setid"))%>% 
  left_join(settings, by="setid")

dfb1 = data.frame(component=1:5, value=ressim1$bvb)
dfb2 = data.frame(component=1:5, value=ressim1$bvb)
dfb3 = data.frame(component=1:5, value=ressim1$bvb)

p1 = ggplot(ressimdf1, aes(x=component, y=value))+
  stat_summary(geom = "line",fun = mean)+
  geom_point()+
  geom_line(data = dfb1, linetype = 3, linewidth=0.5)+
  facet_grid(delay~forgetting, labeller = label_both)+
  scale_color_viridis_d()+
  labs(title  = paste("number of columns:", ncols[1]),
       colour="n_batches", shape="n_batches")+
  theme_classic(20)
print(p1)

p2 = ggplot(ressimdf2, aes(x=component, y=value))+
  stat_summary(geom = "line",fun = mean)+
  geom_point()+
  geom_line(data = dfb2, linetype = 3, linewidth=0.5)+
  facet_grid(delay~forgetting, labeller = label_both)+
  scale_color_viridis_d()+
  labs(title  = paste("number of columns:", ncols[3]),
       colour="n_batches", shape="n_batches")+
  theme_classic(20)
print(p2)

p3 = ggplot(ressimdf3, aes(x=component, y=value))+
  stat_summary(geom = "line",fun = mean)+
  geom_point()+
  geom_line(data = dfb3, linetype = 3, linewidth=0.5)+
  facet_grid(delay~forgetting, labeller = label_both)+
  scale_color_viridis_d()+
  labs(title  = paste("number of columns:", ncols[3]),
       colour="n_batches", shape="n_batches")+
  theme_classic(20)
print(p3)

#ggsave(filename = "sim2d_nc100.pdf", plot = p1, width = 9, height = 8)
