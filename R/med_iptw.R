#' Marginal structural models to compute controlled direct effects given an exposue-induced confounder of mediator
#'
#' \code{med_iptw} Computes CDE(M) for given mediator levels, in a setting with an exposure-induced confounder of
#' the mediator-outcome association. Described in chapter 5.3.1 of Tyler's book
#'
#' @param dat The original dataset
#' @param A  the exposure of interest. Must be binary or categorical
#' @param Y  the outcome, currently must be continuous
#' @param C  confounders of either X -> M and/or M -> Y. Can take any form, specified as formula
#' @param M  the mediators of interest. Must be binary or categorical
#' @param L  the exposure-induced confounders of the association of M with Y. Can take any form
#' @param mlvl a matrix or table of probability-mass functions for the mediator, to calculate CDE(M). By default, mlvl is set to
#' the observed sample distributions
#' @param boot specifies the number of bootstrap samples drawn to make the confidence intervals
#' @param nmin  number of participants all categories of exposure must have; samples will be redrawn if this criterion is not met
#' @param quants an optional vector of quantiles for the confidence interval (95 percent by default)
#' @param mids an optional mids object to serve as template for imputations
#' @return An S3 object of class \code{cmed.ipw} containing:
#' @return w  the ipw used in the marginal strucutral model
#' @return cde.int  an array where cde.int[i,,] indexes a matrix corresponding to the CDE calculated for a PMF of M given in mlvl,
#' with each row a bootsrap replicate
#' @return cde.noint a matrix of cde given no M*A interaction, with each row a bootstrap replicate
#' @return te a matrix of total effects, with each row a bootstrap replicate
#' @return ymod1 the marginal strucutral model of Y given no interaction
#' @return ymod2 the marginal strucutral model of Y allowing a A*M interaction
#' @examples \donttest{
#' my_list <- med_iptw(dat = df,  A = my_exposure, Y = my_outcome, M = my_mediator, C = a_confounder + another_confounder, boot = 1000)
#' }
#' @export
med_iptw <- function(dat, A, M, Y, C = NULL, L = NULL, regtype = "gaussian", boot = 10, nmin = 10,
                     quants = c(0.025, 0.5, 0.975), mlvl = NULL, link = logit, mids = NULL, maxit = 5){

  acol <- deparse(substitute(A)) %>% match(names(dat))
  ycol <- deparse(substitute(Y)) %>% match(names(dat))
  mcol <- deparse(substitute(M)) %>% match(names(dat))

  alen <- nlevels(dat[[acol]]) - 1
  mlen <- nlevels(dat[[mcol]]) - 1

  if(class(mlvl) == "table") mlvl <- t(as.matrix(mlvl))
  if(class(mlvl)  == "matrix"){ ilvl <- nrow(mlvl)
  }else ilvl <- 1
  if(is.null(mlvl)) mlvl <- (table(dat[[mcol]]) %>% prop.table)[-1]

  pbar <- txtProgressBar(style = 3)
  cde.int <- array(NA, dim = c(boot, alen, ilvl))
  cde.noint <- te <- matrix(NA, nrow = boot, ncol = alen)
  if(!is.null(mids)) dat2 <- dat
  for(i in 1:boot){
    if(!is.null(mids)) {
      dat <- mice(dat2, m = 1, maxit = maxit, pred = mids$predictorMatrix, method = mids$method, print = FALSE) %>% complete
    }
    if (i == boot){ bi <- 1:nrow(dat)
    }else{
      bi <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
      while(is.factor(dat[acol]) & min(table(dat[bi,acol])) < nmin){
        message("Resampling failed... retrying")
        bi <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
      }
    }

    anum <- (table(dat[bi,acol]) %>% prop.table)[as.numeric(dat[bi, acol])]
    aden.mod <- multinom(data = dat[bi,], substitute(A ~ C), trace = FALSE)
    aden <- predict(aden.mod, newdata = dat[bi,], type = "probs") %>% fmatch(dat[bi, acol])
    mnum.mod <- multinom(data = dat[bi,], substitute(M ~ A), trace = FALSE)
    mnum <- predict(mnum.mod, newdata = dat[bi,], type = "probs") %>% fmatch(dat[bi, mcol])
    mden.mod <- multinom(data = dat[bi,], substitute(M ~ A + C + L), trace = FALSE)
    mden <- predict(mden.mod, newdata = dat[bi,], type = "probs") %>% fmatch(dat[bi, mcol])

    dat$w <- anum*mnum/aden/mden ## ugh if only there was a way not to copy
    dat$wa <- anum/aden

    ymod <- glm(data = dat[bi,], substitute(Y ~ A + M), weights = w)
    cde.noint[i,] <- ymod$coefficients[1:alen + 1]

    ymod2 <- glm(data = dat[bi,], substitute(Y ~ A * M), weights = w)
    g1 <- ymod2$coefficients[1:alen + 1]
    g3 <- ymod2$coefficients %>% tail(alen*mlen)

    cde.int[i,,] <- g1 + t(matrix(g3,alen,mlen) %*% mlvl)

    tmod <- glm(data = dat[bi,], substitute(Y ~ A), weights = wa)
    te[i,] <- tmod$coefficients[1:alen + 1]

    setTxtProgressBar(pbar, i/boot)
  }
  out <- list(
    w = dat$w,
    ymod = ymod, ymod2 = ymod2,
    amod = aden.mod,
    mmod = mden.mod,
    g1 = g1, g3 = g3,
    cde.int = cde.int,
    cde.noint = cde.noint,
    te = te
  )
  class(out) <- "cmed.ipw"
  return(out)
}
