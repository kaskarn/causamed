#' Random interventional analogue effects
#'
#' \code{med_rint} computes NDEr, NIEr and TEr, the random interventional
#' analogues of the natural direct effect, natural indirect effect and total effect,
#' in the presence of exposure-induced confounding of M -> Y. Gives bootstrapped confidence
#' interval. To do: extend to multiple confounders. Might or might not work at the moment.
#'
#' The procedure is described in chapter 5.4.2 of Tyler's book
#'
#' @param dat The original dataset
#' @param A  the exposure of interest. Must be binary or categorical
#' @param Y  the outcome, currently must be continuous
#' @param C  confounders of either X -> M and/or M -> Y. Can take any form, specified as formula
#' @param M  the mediators of interest. Must be binary or categorical
#' @param L  the exposure-induced confounder of the association of M with Y. Must be binary or categorical
#' @param mlvl a matrix or table of probability-mass functions for the mediator, to calculate CDE(M). By default, mlvl is set to
#' the observed sample distributions
#' @param boot specifies the number of bootstrap samples drawn to make the confidence intervals. Default is 10 for testing purposes
#' @param nmin  number of participants all categories of exposure must have; samples will be redrawn if this criterion is not met
#' @param quants an optional vector of quantiles for the confidence interval (95 percent by default)
#' @param mids an optional mids object to serve as template for imputations
#' @return An S3 object of class \code{cmed.ipw} containing:
#' @return nde mean and 95p confidence intervals for NDEr
#' @return nie mean and 95p confidence intervals for NIEr
#' @return te mean and 95p confidence intervals for the total effect
#' @return ter mean and 95p confidence intervals for the random interventional analogue to the total effect
#' @return boots a list with nde, nie, te, and ter resutls for each bootstrap sample
#' @return raw a list of the duplicated dataset and intermediary propensity scores calculated from original data (not resampled)
#' @examples \donttest{
#' my_list <- med_rint(dat = mydat,  A = my_exposure, Y = my_outcome, M = my_mediator, C = a_confounder + another_confounder, L = my_problem, boot = 1000)
#' }
#' @export
med_rint <- function(dat, A, M, Y, C = NULL, L = NULL, astar = "astar", boot = 10, quants = c(0.025, 0.5, 0.975), nmin = 20,
                     mids = NULL, maxit = 5){

  acol <- deparse(substitute(A)) %>% match(names(dat))
  ycol <- deparse(substitute(Y)) %>% match(names(dat))
  mcol <- deparse(substitute(M)) %>% match(names(dat))
  lcol <- deparse(substitute(L)) %>% match(names(dat))
  lchar <- deparse(substitute(L))
  achar <- deparse(substitute(A))

  alen <- nlevels(dat[[acol]]) - 1
  ref <- levels(dat[[acol]])[1]

  nde <- nie <- te <- ter <- niewrong <- matrix(NA, nrow = boot, ncol = alen)
  pb <- txtProgressBar(style = 3)
  if(!is.null(mids)) dat2 <- dat
  if(!is.null(mids)) dat2 <- dat
  for(i in 1:boot){
    if (i == boot){ bi <- 1:nrow(dat)
    }else{
      bi <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
      while(is.factor(dat[acol]) & min(table(dat[bi,acol])) < nmin){
        message("Resampling failed... retrying")
        bi <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
      }
    }
    if(!is.null(mids)) {
      dat <- mice(dat[bi,], m = 1, maxit = maxit, pred = mids$predictorMatrix, method = mids$method, print = FALSE) %>% complete
      for(j in names(dat)) attr(dat[[j]], "contrasts") <- NULL
    }

    #create models from observed data
    amod <- multinom(data = dat[bi,], substitute(A ~ C), trace = FALSE)
    mmod <- multinom(data = dat[bi,], substitute(M ~ A + C + L), trace = FALSE)
    lmod <- multinom(data = dat[bi,], substitute(L ~ A + C), trace = FALSE)

    #get probability of observed A and M
    pac <- predict(amod, newdata = dat[bi,], type = "probs") %>% fmatch(dat[bi,acol])
    pmlac <- predict(mmod, newdata = dat[bi,], type = "probs") %>% fmatch(dat[bi, mcol])

    #duplicate dataset
    ddat <- do.call(rbind, lapply(levels(dat[[acol]]), function (i) mutate(dat[bi,], astar = i)))
    ddat$astar <- factor(ddat$astar, labels = levels(dat[[acol]]))

    #get probability of M under randomized intervention
    pmlasc <- do.call(cbind, lapply(levels(dat[[lcol]]),
                     function (i) predict(mmod,
                                          newdata = mutate_(ddat, .dots = setNames(list(~i, ~astar), list(lchar, achar))),
                                          type = "probs") %>% fmatch(ddat[[mcol]]))
    )
    #get probability of L under randomized intervention
    plasc <- predict(lmod, newdata = mutate_(ddat, .dots = setNames(list(~astar), achar)), type = "probs")
    if(is.null(dim(plasc))) plasc <- cbind(1-plasc, plasc)

    #marginalize over L and create weights
    num <- apply(plasc*pmlasc, 1, sum)
    ddat$w_a <- 1/pac
    ddat$w_rint <- num/pac/pmlac

    #run models and store
    nde[i,] <- lm(data = ddat[ddat$astar == ref,], substitute(Y ~ A), weights = w_rint)$coefficients[-1]
    for(j in 1:alen){
      index <- ddat[[acol]] == levels(ddat[[acol]])[j+1] & (ddat$astar %in% c(ref, levels(ddat[[acol]])[j+1]))
      nie[i,j] <- lm(data = ddat[index,], as.formula(paste0(deparse(substitute(Y)), "~ astar")), weights = w_rint)$coefficients[2]
    }

    niewrong[i,] <- lm(data = ddat[ddat[[acol]] != ref,], as.formula(paste0(deparse(substitute(Y)), "~ astar")), weights = w_rint)$coefficients[-1]
    te[i,] <- lm(data = ddat[ddat$astar == ref,], substitute(Y ~ A), weights = w_a)$coefficients[-1]
    ter[i,] <- lm(data = ddat[ddat$astar == ddat[[acol]],], as.formula(paste0(deparse(substitute(Y)), "~ astar")), weights = w_rint)$coefficients[-1]

    if(!is.null(mids)) dat <- dat2
    setTxtProgressBar(pb, i/boot)
  }

  out <- list(boots = list(nie = nie, niewrong = niewrong, te = te, ter = ter, nde = nde),
              nie = apply(nie, 2, quantile, probs = quants, na.rm = TRUE),
              te = apply(te, 2, quantile, probs = quants, na.rm = TRUE),
              ter = apply(ter, 2, quantile, probs = quants, na.rm = TRUE),
              nde = apply(nde, 2, quantile, probs = quants, na.rm = TRUE),
              nie2 = apply(niewrong, 2, quantile, probs = quants, na.rm = TRUE),
              raw = list(ddat = ddat, pac = pac, pmlac = pmlac, pmlasc = pmlasc, plasc = plasc)
              )
  class(out) <- "cmed.rint"
  return(out)
}
