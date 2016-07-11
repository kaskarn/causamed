#' \code{med_multi}
#'
#' Decomposes the effect of an exposure on an outcome into a direct effect operating through a set of intermediary
#' variables, and a direct effect involving other pathways. A key assumption is the dual lack of confounding of the
#' exposure-outcome, and the mediator-outcome associations. Described in chapter 5.2.1 of Tyler's book.
#'
#' @param dat  a dataframe containing the exposure, outcome, mediators, and confounders
#' @param A  the exposure of interest. Currently must be categorical (or binary)
#' @param Y  the outcome, currently must be continuous, ordinal or binary
#' @param C  confounders of either X -> M and/or M -> Y.
#' @param M  the mediators of interest
#' @param fam  specifies GLM link function and distribution of residuals. Default is gaussian(link = identity)
#' @param boot  number of bootstrap samples used to build the 95p confidence intervals
#' @param nmin  number of participants all categories of exposure must have; samples will be redrawn if this criterion is not met
#' @param mids an optional mids object to serve as template for imputations
#' @examples \donttest{my_list <- med_multi(dat = df,
#'  X = my_exposure,
#'  M = mediator_1 + mediator_2 + ... + mediator_1*mediator_2 + mediator_1*exposure,
#'  Y = my_continuous_outcome,
#'  C = a_confounder + another_confounder * anything,
#'  fam = gaussian(link = "identity"), boot = 1000)}
#' @export
med_multi <- function(dat, A, Y, M, C = NULL, fam = gaussian(link="identity"), boot = 10, nmin = 10, mids = NULL, maxit = 5){
  # Setup
  acol <- deparse(substitute(A)) %>% match(names(dat)) #for speedy addressing
  ycol <- deparse(substitute(Y)) %>% match(names(dat)) #moar speed
  ref <- levels(dat[[acol]])[1]

  # Initialize arrays and progress bar before bootstrap
  q2 <- q3  <- array(NA, dim = c(boot, nlevels(dat[[acol]])))
  pb <- txtProgressBar(style = 3); a <- 0
  # Run bootsrapped procedure
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

    # Run exposure and outcome models
    amod <- multinom(data = dat[bi,], substitute(A ~ C), trace = FALSE)
    ymod <- glm(data = dat[bi,], substitute(Y ~ A + M + C), family = fam)

    # Get marginal and conditional probabilities of A
    pa <- (table(dat[bi,][[acol]]) %>% prop.table)
    pac <- predict(object = amod, newdata = dat[bi,], type = "probs")

    # Get weights like VdW (strangely?) explains and calculate expectations at levels of A
    w <- fmatch(((1/pac) %*% diag(pa)), dat[bi,acol])
    q2[i,] <- sapply(levels(dat[[acol]]), function (i) weighted.mean(x = dat[bi,ycol][dat[bi,acol] == i], w = w[dat[bi,acol] == i], na.rm = TRUE))

    # Calculate Y_aMa* for all a in support of A, from model of Y
    xset <- function (i) predict(object = ymod, newdata = mutate_(dat[bi,], .dots = setNames(list(~i), names(dat)[acol])))
    ypred <- sapply( levels(dat[[acol]]), function(i) xset(i) )[dat[bi,acol] == ref,]
    q3[i,] <- apply(ypred, 2, weighted.mean, w = w[dat[bi,acol] == ref], na.rm = TRUE)

    if(!is.null(mids)) dat <- dat2
    setTxtProgressBar(pb, i/boot)
  }
  if(ymod$family[2] == "identity"){ ### additive scale outcomes
    nde <- q3 - replicate(nlevels(dat[[acol]]), q2[,1])
    nie <- q2 - q3
    pm <- nie / (nde + nie)
  }else{ ### log scales
    nde <- q3 / replicate(nlevels(dat[[acol]]), q2[,1])
    nie <- q2 / q3
    pm <- NULL #hmm
  }
  out <- list(raw = list(q2 = q2, q3 = q3, nde = nde, nie = nie, pm = pm, ypred = ypred),
              nie = apply(nie, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE),
              nde = apply(nde, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE),
              pm = apply(pm, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  class(out) <- "cmed.multi"
  return(out)
}
