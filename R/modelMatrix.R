#' @export
sparse_model_matrix_b <- function(object,
                                  data = environment(object),
                                  xlev = NULL,
                                  transpose = FALSE,
                                  drop.unused.levels = FALSE,
                                  row.names = TRUE,
                                  sep = "_",
                                  verbose = FALSE,
                                  na.action='na.pass'){
  data <- model.frame(object, data, xlev=xlev, na.action=na.action)
  t <- if(missing(data)) terms(object) else terms(object, data=data)
  X <- Matrix:::model.spmatrix(trms = t,
                      mf = data,
                      transpose = transpose,
                      drop.unused.levels = drop.unused.levels,
                      row.names = row.names,
                      sep = sep,
                      verbose = verbose)
  m <- attr(X, "assign")
  dropK <- logical(length = length(m))
  dropK[-1] <- TRUE
  ind <- which(diff(m)!=0)
  Y <- methods::as(X, "lgCMatrix") # equivalent to X > 0
  attr(Y, "indices") <- c(0, ind, length(m))
  attr(Y, "term.labels") <- attr(t, "term.labels")
  attr(Y, "dropK") <- dropK
  return(Y)
}
