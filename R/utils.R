#' utils
#'
#' Functions used to run the package
#'
#' @param v The matrix from which elements are extracted
#' @param A The vector indexing the elements of the matrix v
#' @param l The length of A, given to speed up slightly
#' @param mat indicator that the indexing vector should be returned instead of of the indexed matrix
fmatch <- function(v, A, l = NULL, mat = FALSE){
  #replaces syntax of type [cbind(1:nrow(dat),match(dat[[M]], levels(dat[[M]])))]
    if(!is.factor(A)) stop("Must be factor variable")
    if(is.null(l)) l <- length(A)
    if(is.null(dim(v)) == 2) v <- cbind(1-v, v)
    if(mat){return(cbind(1:l,match(A, levels(A))))
    }else return(v[cbind(1:l,match(A, levels(A)))])
}
