library(Matrix)
library(readr)
library(tidyr)
library(dplyr)
path = scan("datapath.txt", what = character())
d <- dir(path[5], full.names = TRUE)
d <- d[grep("_to_Human", d)]
species <- c("Pig","Mouse","MacaqueM","MacaF", "Human")


mtxname <- paste0(path[6],"\\", species, "_count_matrix.mtx")
gname <- paste0(path[6],"\\", species, "_count_matrix_gene.txt")
cname <- paste0(path[6],"\\", species, "_count_matrix_cell.txt")
Y <- vector("list", 5)
X0count <- vector("list", 5)
H0count = vector("list", 5)
H1list = vector("list", 5)
names(H1list) <- species
names(H0count) <- species
names(Y) <- species
for(it in 1:4){
    path_to_human <- d[grep(species[it], d)]
    Mat <- readMM(mtxname[it])
    dim(Mat)
    to_human <- read_tsv(path_to_human)
    to_human <- rename(to_human, all_of(c(gene=colnames(to_human)[1])))
    to_human <- group_by(to_human, gene) %>% 
      reframe(Human = unlist(strsplit(Human, split = "; "))) %>% 
      mutate(species=species[it])
    
    genes <- read_csv(gname[it])
    genes <- left_join(genes, to_human, by = "gene")
    
    
    cells <- read_csv(cname[it])
    cells <- tidyr::separate(cells,  cell,  into=c("sample", "code"),
                             sep = "_", 
                             extra = "merge", remove = FALSE)
    
    Y[[it]] <- as(Mat, "sparseVector")
    
    dfX = expand.grid(gene_num = seq_len(nrow(Mat)), 
                      cell_num = seq_len(ncol(Mat))) %>% 
      mutate(species_gene = genes$gene[gene_num],
             cell = cells$cell[cell_num],
             sample = cells$sample[cell_num],
             species = species[it]) %>% 
      mutate(species_gene = paste(species, species_gene, sep = ":")) 

    H = strsplit(genes$Human[dfX$gene_num], "; ")
    H0count[[it]] <- as.data.frame(table(unlist(H[-Y[[it]]@i])))
    X0count[[it]] <- lapply(dplyr::slice(dfX, -Y[[it]]@i), table)
    H1list[[it]] <- H[Y[[it]]@i]
    dfX <- dplyr::slice(dfX, Y[[it]]@i)
    
    write.table(dfX, file = "dfX.csv", sep = ",",
              row.names = FALSE,
              col.names = !(it > 1L),
              append = (it > 1L))
}
###
it = 5
Mat <- readMM(mtxname[it])
genes <- read_csv(gname[it])
genes <- left_join(genes, to_human, by = "gene")

cells <- read_csv(cname[it])
cells <- tidyr::separate(cells,  cell,  into=c("sample", "code"),
                         sep = "_", 
                         extra = "merge", remove = FALSE)

Y[[it]] <- as(Mat, "sparseVector")

dfX = expand.grid(gene_num = seq_len(nrow(Mat)), 
                  cell_num = seq_len(ncol(Mat))) %>% 
  mutate(species_gene = genes$gene[gene_num],
         cell = cells$cell[cell_num],
         sample = cells$sample[cell_num],
         species = species[it]) %>% 
  mutate(species_gene = paste(species, species_gene, sep = ":")) 

H = genes$gene[dfX$gene_num]
H1list[[it]] <- H[Y[[it]]@i]
H0count[[it]] <- as.data.frame(table(H[-Y[[it]]@i]))
X0count[[it]] <- lapply(dplyr::slice(dfX, -Y[[it]]@i), table)

dfX <- dplyr::slice(dfX, Y[[it]]@i)

write.table(dfX, file = "dfX.csv", sep = ",",
            row.names = FALSE,
            col.names = !file.exists("dfX.csv"),
            append = TRUE)

###
saveRDS(Y, file = "listY.rds")
saveRDS(X0count, file = "listX0count.rds")
saveRDS(H1list, file = "H1list.rds")
saveRDS(H0count, file = "H0cont.rds")
##

# listH <- vector("list", length(species))
# names(listH) <- species
# for(it in 1:4){
#   path_to_human <- d[grep(species[it], d)]
#   to_human <- read_tsv(path_to_human)
#   to_human <- rename(to_human, all_of(c(gene=colnames(to_human)[1])))
#   to_human <- group_by(to_human, gene) %>% 
#     reframe(Human = unlist(strsplit(Human, split = "; "))) %>% 
#     mutate(species=species[it])
#   listH[[it]] <- to_human
# }
# 
# genes <- read_csv(gname[5])
# listH[[5]] <- data.frame(gene=genes$gene,
#                          Human=genes$gene,
#                          species=species[5])
# 
# 
# 
# 
# 
