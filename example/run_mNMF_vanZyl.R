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

system.time({
  out <- moltenNMF:::mNMF_svb(y = Y, X = X,
                              L = 15L,
                              n_batches=2000L,
                              n_epochs=10,
                              lr_type="exponential",
                              lr_param=c(15,0.9))
})

saveRDS(out, file="out_L15.rds")

