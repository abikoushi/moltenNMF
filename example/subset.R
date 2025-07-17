library(Matrix)
i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(A <- sparseMatrix(i, j, x = x))  ##  8 x 10 "dgCMatrix"
writeMM(obj = A, file = "examplematrix.mtx")
moltenNMF:::size_mtx("examplematrix.mtx")

y = rpois(1000,0.1)
y = as(y, "sparseVector")
yv=y@x
yi=y@i


moltenNMF:::filter_y(yv, yi,
                     y@x, y@i,
                     ind=0:99)

yi
yv

spY = y[1:100]
spY@x
spY@i

X <- moltenNMF::sparse_onehot(~row+col, data=expand.grid(row=1:12, col=1:11))

id = sort(sample.int(nrow(X), 10))
subX <- X[id,]

subX@p
subX
subX@i

res = c(moltenNMF:::subset_spx(X@i, X@p, uid=id-1))
subX@p
c(res[[1]])
subX@i
c(res[[2]])
