library(Matrix)
i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(A <- sparseMatrix(i, j, x = x))  ##  8 x 10 "dgCMatrix"
writeMM(obj = A, file = "examplematrix.mtx")
moltenNMF:::size_mtx("examplematrix.mtx")
A
moltenNMF:::rowfilter_mtx("examplematrix.mtx", "filtered.mtx", 
                          c(1L, 4L))

readMM("filtered.mtx")

