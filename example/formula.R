d <- as.data.frame(Titanic)
str(sparse_onehot(Freq ~ Class+Sex:Age,data=d))
