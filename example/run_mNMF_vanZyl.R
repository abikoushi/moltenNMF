library(Matrix)
library(dplyr)
library(moltenNMF)
Hlist = readRDS("listH.rds")
Hlist = do.call("rbind", Hlist)
Xlist <- bind_rows(readRDS("listX.rds"))
Xlist <- mutate(Xlist, species_gene=if_else(species=="Human", NA_character_, species_gene))
f <- ~ species_gene + cell + sample + species
X <- moltenNMF::sparse_onehot(f, data = Xlist)
Y = readRDS("listY.rds")
Y = do.call("c", Y)
X = append_new(X, Hlist)

(nnzero(Y))

system.time({
  out <- moltenNMF:::mNMF_svb(y = Y, X = X,
                              L = 20L,
                              n_batches=2e7L,
                              n_epochs=10,
                              lr_type="exponential",
                              lr_param=c(15,0.9))
})

#  ユーザ   システム       経過  
# 6885.12     462.28    7412.15 

plot(out$ELBO[-1], type="l")

# system.time({
#   out <- moltenNMF:::mNMF_svb_batch(y = Y, X = X,
#                               L = 20L,
#                               n_epochs=10,
#                               lr_type="exponential",
#                               lr_param=c(1.5,0.7))
# })

saveRDS(out, file="out_L20.rds")
