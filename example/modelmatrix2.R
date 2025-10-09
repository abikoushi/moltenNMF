library(moltenNMF)
library(Matrix)
library(gridExtra)
library(dplyr)

df = as.data.frame(Titanic)
df = as.data.frame(HairEyeColor)
head(df)
df = expand.grid(w=1:10,x=1:10,y=1:10,z=1:10)
#f = Freq ~ Class + Sex + Age + Survived
#f = Freq ~ Hair + Eye + Sex
f = ~ w + x + y + z
X = sparse_onehot(f, data = df)
X2 = moltenNMF:::sparse_onehot2(f, data = df)
colnames(X2)
grid.arrange(image(t(X)), image(t(X2)), nrow=2)
str(X)
str(X2)

sparseMatrix(i = c(1L,214748364), p= c(0,1))
