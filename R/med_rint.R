#' Random interventional analogue effects
#'
#' \code{med_rint} computes NDEr, NIEr and TEr, the random interventional
#' analogues of the natural direct effect, natural indirect effect and total effect,
#' in the presence of exposure-induced confounding of M -> Y. Gives bootstrapped confidence
#' interval. To do: extend to multiple confounders. Might or might not work at the moment.
#'
#' The procedure is described in chapter 5.4.2 of Tyler's book
#'
#' @param orig_dat The original dataset
#' @param A a character string containing the name of the exposure variable. A must be categorical/binary
#' @param M a character string containing the name of the mediator variable. M must be categorical/binary
#' @param Y a character string containing the name of the outcome variable. At the moment, Y must be continuous
#' @param astar optional alternative name for (pseudo) randomized treatment level a*
#' @param C list of character strings containing the name or expression of confounder variables
#' @param quants optionally specify own confidence interval (by default: 95p CI)
#' @param boot number of bootstrap samples (default 10 for testing)
#' @return A large list containing meta-data on the procedure, along with the median and 95p confidence interval
#' for the randomized interventional analogue direct and indirect effects
#' @examples \donttest{my_list <- med_rint(dat = df,
#'  X = "my_exposure",
#'  M = "my_mediator",
#'  Y = "my_binary_outcome",
#'  C = c("a_confounder", "another_confounder"),
#'  fam = binomial(logit), boot = 1000)}
#' @export
med_rint <- function(dat, A, M, Y, C = "", L, astar = "astar", boot = 10, quants = c(0.025, 0.5, 0.975), alex = FALSE){
  alen <- levels(dat[[A]])[-1] %>% length
  mlen <- levels(dat[[M]])[-1] %>% length

  # Initialize arrays of results and progress bar
  ar_nder <- ar_pm <- ar_nier <- ar_ter <- ar_te <- array(NA, dim = c(boot, alen), dimnames = list(1:boot, levels(dat[[A]])[-1]))
  pb <- txtProgressBar(style = 3)

  # Start bootstrapping
  for(i in 1:boot)
  {
    tdat <- sample_frac(dat, 1, replace = TRUE)
    while(min(table(tdat[[A]])) < 20){
      print("resampling failed, retrying...")
      tdat <- sample_frac(dat, 1, replace = TRUE)
    }

    # Make weights
    rint_items <- rint_med.mkdata(tdat, A = A, C = C, M = M, L = L)
    setTxtProgressBar(pb, (2*i-1)/boot/2)

    # Calculate effects
    if(alex){ res <- rint_med.decompose(rint_items$an_dat %>% mutate(w = w_ak), Y = Y, A = A)
    }else res <- rint_med.decompose(rint_items$an_dat, Y = Y, A = A)

    setTxtProgressBar(pb, i/boot)
    ar_ter[i,] <- res$ter
    ar_nder[i,] <- res$nder
    ar_nier[i,] <- res$nier
    ar_pm[i,] <- res$nier / res$ter
    ar_te[i,] <- res$te
  }

  te_95 <- lapply(quants, function (i) apply(ar_te, 2, quantile, probs = i, na.rm = TRUE))
  ter_95 <- lapply(quants, function (i) apply(ar_ter, 2, quantile, probs = i, na.rm = TRUE))
  nder_95 <- lapply(quants, function (i) apply(ar_nder, 2, quantile, probs = i, na.rm = TRUE))
  pm_95 <- lapply(quants, function (i) apply(ar_pm, 2, quantile, probs = i, na.rm = TRUE))
  nier_95 <- lapply(quants, function (i) apply(ar_nier, 2, quantile, probs = i, na.rm = TRUE))

  return(list(
    arrays = list(te = ar_te, ter = ar_ter, nder = ar_nder, nier = ar_nier, pm = ar_pm),
    te = ter_95,
    ter = ter_95,
    nder = nder_95,
    pm = pm_95,
    nier = nier_95)
  )
}


#' rint_med.mkdata
#'
#' Creates the analytic dataset to compute NDEr, NIEr and TEr, the random interventional
#' analogues of the natural direct effect, natural indirect effect and total effect,
#' in the presence of exposure-induced confounding of M -> Y
rint_med.mkdata <- function(orig_dat, A, M, Y, C = "", L){
  ### replicate C with parentheses to allow arbitrary syntax in confounder specification (like equations or functions)
  C <- paste0("(",C,")")

  rint_med.mkform <- function(A, C, L, M){
    # Creates formulas to be passed to propensity score models
    a_form <- as.formula(paste0(A, "~", paste0(C, collapse = "+")))
    m_form <- as.formula(paste0(M, "~", paste(A, paste0(L, collapse = "+"), paste0(C, collapse = "+"), sep = "+")))
    l_form <- as.formula(paste0(L, "~", A, "+", paste0(C, collapse = "+")))
    return(list(a_form = a_form, m_form = m_form, l_form = l_form))
  }
  rint_med.mkmods <- function(dat, f){
    # Makes propensity score models from formulas in rint_med.mkform
    a_mod <- multinom(data = dat, f$a_form, trace = FALSE)
    m_mod <- multinom(data = dat, f$m_form, trace = FALSE)
    l_mod <- multinom(data = dat, f$l_form, trace = FALSE)
    return(list(a_mod = a_mod, m_mod = m_mod, l_mod = l_mod))
  }
  rint_med.mkden <- function(dat, mods, A, L){
    # Makes denominator
    pmlac <- predict(object = mods$m_mod, newdata = dat, type = "probs") %>% fmatch(dat[[M]])
    pac <- predict(object = mods$a_mod, newdata = dat, type = "probs") %>% fmatch(dat[[A]])
    return(list(pmlac = pmlac, pac = pac, den = pac * pmlac))
  }
  rint_med.copy_data <- function(dat, A){
    ## One copy for each value in support of A
    new <- do.call(rbind, lapply(levels(dat[[A]]), function (i) mutate(dat, astar = i)))
    new$astar <- factor(new$astar, labels = levels(dat[[A]]))
    return(new)
  }

  rint_med.mknum <- function(dat, M, L, A, mmod, lmod){
    ## get p( l | a*,c) for each value in support of L
    plasc <- predict(object = lmod, newdata = mutate_(dat, .dots = setNames(list(~astar), A)), type = "probs")
    if(is.null(dim(plasc))) plasc <- cbind(1-plasc, plasc)

    ## get p( m | l,a*,c) for each value in support of L
    pmlasc <- lapply(levels(dat[[L]]), function (i)
                      predict(object = mmod, newdata = mutate_(dat, .dots = setNames(list(~astar), A)) %>%
                                                       mutate_(.dots = setNames(list(~i), L)),
                              type = "probs") %>% fmatch(dat[[M]])
                     ) %>% do.call(cbind, .)

    ## sum p( m | l,a*,c ) * p( l | a*,c ) across levels of L
    num <- apply((plasc * pmlasc), 1, sum)


    #Same pmlasc as AK
    pmlasc_ak <- predict(object = mmod,
                       newdata = mutate_(dat, .dots = setNames(list(~astar), A)),
                       type = "probs") %>% fmatch(dat[[M]])
    num_ak <- apply((plasc * pmlasc_ak), 1, sum)
    return(list(plasc = plasc, pmlasc = pmlasc, num = num, num_ak = num_ak))
  }

  rint_med.mkweight <- function(dat){
    w <- dat$num/dat$den
    return(w)
  }
  #To get same weights as AK
  rint_med.mkweight_ak <- function(dat){
    w <- dat$num_ak/dat$den
    return(w)
  }

  ## create models
  flist <- rint_med.mkform(A, C, L, M)
  mlist <- rint_med.mkmods(orig_dat, flist)

  ## calculate denominator
  den_items <- rint_med.mkden(orig_dat, mlist, A, L)
  orig_dat$ipw_conf <- 1/den_items$pac
  orig_dat$den <- den_items$den

  ## duplicate dataset
  an_dat <- rint_med.copy_data(orig_dat, A)

  ## calculate numerator
  num_items <- rint_med.mknum(an_dat, M, L, A, mmod = mlist$m_mod, lmod = mlist$l_mod)
  an_dat$num <- num_items$num

  ## calculate weight (regular)
  an_dat$w <- rint_med.mkweight(an_dat)

  # get same weights as AK
  an_dat$num_ak <- num_items$num_ak
  an_dat$w_ak <- rint_med.mkweight_ak(an_dat)

  return(list(
    models = mlist,
    an_dat = an_dat,
    pmlac = den_items$pmlac,
    pac = den_items$pac,
    plasc = num_items$plasc,
    pmlasc = num_items$pmlasc,
    den = an_dat$den,
    num = an_dat$num,
    w = an_dat$w,

    #same weights as AK
    w_ak = an_dat$w_ak,
    num_ak = num_items$num_ak
  ))
}



#' rint_med.decomposes
#'
#' part of the functio calculating CDE(M) at different values or distributions of M, and computes weights used in
#' master function
rint_med.decompose <- function(dat, Y, A, astar = "astar"){
  ref <- levels(dat[[A]])[1]
  m_te <- lm(data = dat, as.formula(paste0(Y, "~", A)), weights = ipw_conf)
  m_ter <- lm(data = dat[dat[[astar]] == dat[[A]],], as.formula(paste0(Y, "~", astar)), weights = w)
  m_nder <- lm(data = dat[dat[[astar]] == ref,], as.formula(paste0(Y, "~", A)), weights = w)
  m_nier_prollywrong <- lm(data = dat[dat[[A]] != ref,], as.formula(paste0(Y, "~", astar)), weights = w)
  nier <- lapply(levels(dat[[A]])[-1],
                 function(i) lm(data = dat[(dat[[A]] == i) & (dat[[astar]] %in% c(ref,i)),],
                                as.formula(paste0(Y, "~", astar)),
                                weights = w)$coefficients[-1]
                ) %>% unlist

  return(list(
    models = list(m_te = m_te, m_ter = m_ter, m_nder = m_nder,
                  m_nier_prollywrong = m_nier_prollywrong),
    nder = m_nder$coefficients[-1],
    nier = nier,
    nier_prollywrong = m_nier_prollywrong$coefficients[-1],
    ter = m_ter$coefficients[-1],
    te = m_te$coefficients[-1]
  ))
}
