#' Structural means model to compute controlled direct effects given an exposue-induced confounder of mediator
#'
#' \code{med_smean}
#'
#' Returns the controlled direct effect CDE(M) of an exposure A on an outcome Y, not
#' operating through a mediator M, in the case of exposure-induced confounding of
#' the M->Y association. Follows section 5.3.5 of Vanderweele's book on causal
#' mediation analysis.
#'
#' @param dat  a dataframe containing the exposure, outcome, mediators, and confounders
#' @param A  the exposure of interest. Can take any form
#' @param Y  the outcome, currently must be continuous
#' @param C  confounders of either X -> M and/or M -> Y. Can take any form
#' @param M  the mediators of interest
#' @param L  the exposure-induced confounders of the association of M with Y
#' @param boot  the number of bootstrap samples used to build confidence intervals
#' @param nmin  number of participants all categories of exposure must have; samples will be redrawn if this criterion is not met
#' @param mlvl  the levels of M to calculate corresponding CDE's to. Default is sample average.
#' @param quants an optional vector of quantiles for the confidence interval (95 percent by default)
#' @param mids an optional mids object to serve as template for imputations
#' @return  An S3 object of class \code{cmed_smean} containing:
#' @return  k2  the coefficient of M in regression of A, M, L and C on Y
#' @return  k3  the coefficient of A*M from regression of A, M, L and C on Y
#' @return  g1  the coefficient g1 from the regression of A and C on partial residuals of Y
#' @return  ymod1  the model of Y for the last bootstrap sample
#' @return  ymod2  the model of Y residuals for the last bootstrap sample
#' @return  cde  array story controlled direct effects from bootstrapped results
#' @examples \donttest{
#' my_list <- med_smean(data, A, Y, M, c1 + c2 + c3 + c2*c3, mlvl = quantiles(M, probs = c(0.25, 0.5, 0.75)))
#' }
#' @export
med_smean <- function(dat, A, Y, M, C = NULL, L = NULL, boot = 10, nmin = 10, mlvl = NULL, quants = c(0.025, 0.5, 0.975),
                      mids = NULL, maxit = 5){
  acol <- deparse(substitute(A)) %>% match(names(dat))
  mcol <- deparse(substitute(M)) %>% match(names(dat))
  ycol <- deparse(substitute(Y)) %>% match(names(dat))

  if(is.factor(dat[[acol]])) alen <- nlevels(dat[[acol]]) - 1 else alen <- 1
  if(is.factor(dat[[mcol]])) mlen <- nlevels(dat[[mcol]]) - 1 else mlen <- 1

  if(is.null(mlvl)){
    if(is.factor(dat[[mcol]])) mlvl <- (table(dat[[mcol]]) %>% prop.table)[-1]
    else mlvl <- mean(dat[[mcol]], na.rm = TRUE)
  }
  cde <- array(NA, dim = c(boot, alen, length(mlvl)/mlen))
  pb <- txtProgressBar(style = 3)
  if(!is.null(mids)) dat2 <- dat
  for(i in 1:boot){
     if(!is.null(mids)) {
      dat <- mice(dat2, m = 1, maxit = maxit, pred = mids$predictorMatrix, method = mids$method, print = FALSE) %>% complete
    }
    if (i == boot){ bi <- 1:nrow(dat)
    }else{
      bi <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
      while(is.factor(dat[acol]) & min(table(dat[bi,acol])) < 10){
        message("Resampling failed... retrying")
        bi <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
      }
    }
    ymod <- lm(data = dat[bi,], substitute(Y ~ A + M + A*M + C + L))
    if(i == 1) intpos <- grep(":", names(ymod$coefficients))[1]
    k3 <- ymod$coefficients[intpos:(intpos+alen*mlen-1)]
    k2 <- ymod$coefficients[(alen+2) : (alen+mlen+1)]

    if(is.factor(dat[[acol]]) && is.factor(dat[[mcol]])){
      intmat <- cbind(rep(0, alen+1), rbind(0, matrix(k3, alen, mlen)))
      yp <- dat[bi,ycol] - intmat[cbind(dat[bi, acol], dat[bi, mcol])] - c(0,k2)[as.numeric(dat[bi, mcol])]
    }else if(is.factor(dat[[acol]])) {
      yp <- dat[bi, ycol] - dat[bi, mcol] * c(0,k3)[as.numeric(dat[bi, acol])] - k2*dat[bi, mcol]
    }else if(is.factor(dat[[mcol]])) {
      yp <- dat[bi, ycol] - dat[bi, acol] * c(0,k3)[as.numeric(dat[bi, mcol])] - c(0,k2)[as.numeric(dat[bi, mcol])]
    }else yp <- dat[bi, ycol] - dat[bi, mcol]*dat[bi, acol]*k3

    # We lookup yp in parent environment to avoid collision with dataframe
    ymod2 <- lm(data = dat[bi,], substitute(eval(quote(yp), parent.frame()) ~ A + C))
    g1 <- ymod2$coefficients[1:alen + 1]


    cde[i,,] <- g1 + t(matrix(k3, alen, mlen) %*% mlvl)
    setTxtProgressBar(pb, i/boot)
  }
  out <- list(raw = list(cde = cde),
              last = list(g1 = g1, k2 = k2, k3 = k3, ymod1 = ymod, ymod2 = ymod2),
              cde = apply(cde, 2, quantile, probs = c(quants), na.rm = TRUE)
  )
  class(out) <- "cmed.smean"
  return(out)
}
