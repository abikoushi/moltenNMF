# 2D array
NMF2D_vb <- function(Y, rank,
                     iter=100,
                     weight = NULL,
                     prior_shape = 1, prior_rate = 1,
                     Vini = NULL,
                     dims=NULL,
                     display_progress=TRUE){
  if(all(class(Y)!="dgTMatrix")){
    Y = as(Y, "TsparseMatrix") 
  }
  if(is.null(dims)){
    dims <- dim(Y)
  }
  if(is.null(Vini)){
  Vini <- list(matrix(rgamma(dims[1]*rank, 1, 1), dims[1], rank),
               matrix(rgamma(dims[2]*rank, 1, 1), dims[2], rank))
  }
  if(!is.null(weight)){
    out = doVB_pois_2D_ww(Vini,
                          y = Y@x[nn_wch],
                          rowi = Y@i[nn_wch],
                          coli = Y@j[nn_wch],
                          dims = dims,
                          L = rank, iter=iter, 
                          weight = weight,
                          a = prior_shape, b = prior_rate,
                          display_progress = display_progress)
  }else{
    naY = is.na(Y)
    if(any(naY)){
      weight <- list( (Y@Dim[1] - rowSums(naY))/Y@Dim[1],
                      (Y@Dim[2] - colSums(naY))/Y@Dim[2] )
      nn_wch = which(!is.na(Y@x))
      out = doVB_pois_2D_ww(Vini,
                            y = Y@x[nn_wch],
                            rowi = Y@i[nn_wch],
                            coli = Y@j[nn_wch],
                            dims = dims,
                            L = rank, iter=iter, 
                            weight = weight,
                            a = prior_shape, b = prior_rate,
                            display_progress = display_progress)
    }else{
      out = doVB_pois_2D(Vini,
                         y = Y@x, rowi = Y@i,  coli = Y@j,
                         dims = dims,
                         L = rank, iter=iter, 
                         a = prior_shape, b = prior_rate,
                         display_progress = display_progress)
      
    }    
  }
  return(out)
}

NMF2D_svb <- function(Y, rank,
                      n_epochs, 
                      n_baches,
                      lr_param = c(1, 0.8),
                      lr_type = "exponential",
                      dims=NULL,
                      prior_shape=1, prior_rate=1,
                      Vini = NULL,
                      weight = NULL,
                      display_progress=TRUE){
  if(all(class(Y)!="dgTMatrix")){
    Y = as(Y, "TsparseMatrix")    
  }
  if(is.null(dims)){
    dims <- dim(Y)
  }
  if(is.null(Vini)){
    Vini <- list(matrix(rgamma(dims[1]*rank, 1, 1), dims[1], rank),
                 matrix(rgamma(dims[2]*rank, 1, 1), dims[2], rank))
  }
  n_baches <- min(n_baches, length(Y@x))
  if(is.null(weight)){
    out = doVB_pois_s_2D(y = Y@x, rowi = Y@i,  coli = Y@j,
                         L = rank,
                         iter = n_epochs,
                         a = prior_shape, b = prior_rate,
                         N1 = length(Y@x),
                         Nr = dims[1],
                         Nc = dims[2],
                         bsize = n_baches,
                         lr_param = lr_param,
                         lr_type = lr_type,
                         display_progress = display_progress)
  }else{
    out = doVB_pois_s_2D_ww(y = Y@x, rowi = Y@i,  coli = Y@j,
                         L = rank,
                         iter = n_epochs,
                         weight = weight,
                         a = prior_shape, b = prior_rate,
                         N1 = length(Y@x),
                         Nr = dims[1],
                         Nc = dims[2],
                         bsize = n_baches,
                         lr_param = lr_param,
                         lr_type = lr_type,
                         display_progress = display_progress)    
  }
  return(out)
}


# NMF2D_svb_t1 <- function(Y, rank,
#                       n_epochs, 
#                       n_baches,
#                       lr_param = c(1, 0.8),
#                       lr_type = "exponential",
#                       dims=NULL,
#                       subiter = 1,
#                       prior_shape=1, prior_rate=1,
#                       Vini = NULL,
#                       weight = NULL,
#                       display_progress=TRUE){
#   if(all(class(Y)!="dgTMatrix")){
#     Y = as(Y, "TsparseMatrix")    
#   }
#   if(is.null(dims)){
#     dims <- dim(Y)
#   }
#   if(is.null(Vini)){
#     Vini <- list(matrix(rgamma(dims[1]*rank, 1, 1), dims[1], rank),
#                  matrix(rgamma(dims[2]*rank, 1, 1), dims[2], rank))
#   }
#   n_baches <- min(n_baches, length(Y@x))
#   out = doVB_pois_s_2D_t1(y = Y@x, rowi = Y@i,  coli = Y@j,
#                        L = rank,
#                        iter = n_epochs,
#                        a = prior_shape, b = prior_rate,
#                        N1 = length(Y@x),
#                        Nr = dims[1],
#                        Nc = dims[2],
#                        bsize = n_baches,
#                        lr_param = lr_param,
#                        lr_type = lr_type,
#                        display_progress = display_progress)
#   return(out)
# }

NMF2D_svb_bin <- function(filepath_x,
                         filepath_y,
                         dims,
                         rank,
                         n_epochs, 
                         n_baches,
                         lr_param,
                         lr_type = "exponential",
                         subiter = 1,
                         prior_shape=1, prior_rate=1,
                         Vini = NULL,
                         display_progress=TRUE){
  n_baches <- min(n_baches, dims[3])
  out = doVB_pois_s_2D_bin(filepath_x,
                           filepath_y,
                           L = rank,
                           iter = n_epochs,
                           a = prior_shape, b = prior_rate,
                           N1 = dims[3],
                           Nr = dims[1], Nc = dims[2],
                           bsize = n_baches,
                           lr_param = lr_param, lr_type = lr_type,
                           display_progress = display_progress)
  return(out)
}
