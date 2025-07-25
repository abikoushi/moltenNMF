library(Matrix)
i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(A <- sparseMatrix(i, j, x = x))  ##  8 x 10 "dgCMatrix"
writeMM(obj = A, file = "examplematrix.mtx")
moltenNMF:::size_mtx("examplematrix.mtx")
moltenNMF:::rowfilter_mtx("obs.mtx", "filtered.mtx", c(1, 5))
readMM("filtered.mtx")
moltenNMF:::rowmeanvar_mtx("obs.mtx", n_header = 2)
