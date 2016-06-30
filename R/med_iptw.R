#' \code{med_iptw} Computes CDE(M) for a given mediator level, in a setting with an exposure-induced confounder of
#' the mediator-outcome association. Described in chapter 5.3.1 of Tyler's book
#'
#' @param orig_dat The original dataset
#' @param A a character string containing the name of the exposure variable. A must be categorical/binary
#' @param M a character string containing the name of the mediator variable. M must be categorical/binary
#' @param Y a character string containing the name of the outcome variable. At the moment, Y must be continuous, ordinal or binary
#' @param C vector of text strings containing the name or expression of confounder variables
#' @param L vector of text strings containing the name of exposure-induced confounders of M -> Y
#' @param noint drops the exposure-mediator interaction, in which case the CDE = NDE
#' @param mlvl a matrix of probability-mass functions for the mediator, to calculate CDE(M). By default, mlvl is set to
#' (1) the identity matrix with an extra row of zeros (everyone unexposed), and (2) p = 1/k across all k categories of the
#' mediator.
#' @param regtype specifies GLM link function and distribution of residuals. Default is gaussian(link = identity)
#' @param boot specifies the number of bootstrap samples drawn to make the confidence intervals
#' @param quants an optional vector of quantiles for the confidence interval (95 percent by default)
#' @return A list containing information about the models, the controlled direct effects computed at
#' each specified mediator level, the total effect, and the difference between those two measures at each
#' mediator level
#' @examples \donttest{my_list <- med_iptw(dat = df,
#' X = "my_exposure",
#' M = "my_mediator",
#' Y = "my_binary_outcome",
#' C = c("a_confounder", "another_confounder"),
#' regtype = binomial(logit), boot = 1000)}
#' @export
med_iptw <- function(dat, A, M, Y, C = "", L = "", regtype = "gaussian", noint = FALSE, boot = 10,
                          quants = c(0.025, 0.5, 0.975), mlvl = NULL){

  alen <- levels(dat[[A]])[-1] %>% length
  mlen <- levels(dat[[M]])[-1] %>% length
  if(is.null(mlvl)) mlvl <- rbind(rep(0, mlen), diag(mlen), rep(1/mlen, mlen))
  if(is.null(dim(mlvl))){
    if(length(mlvl) == mlen) mlvl <- t(mlvl) else stop("Non conformable mlvl")
  }

  ar_te <- array(NA, dim = c(boot, alen))
  if(noint == TRUE){ ar_cde <- ar_pm <- ar_cie <- ar_te
  }else ar_cde <- ar_pm <- ar_cie <- array(dim = c(boot, alen, nrow(mlvl)))

  pb <- txtProgressBar(style = 3)
  for(i in 1:boot)
  {
    tdat <- sample_frac(dat, 1, replace = TRUE)
    while(min(table(tdat[[A]])) < 20){
      print("resampling failed, retrying...")
      tdat <- sample_frac(dat, 1, replace = TRUE)
    }

    an_dat <- med_iptw.mkdat(tdat, A = A, M = M, Y = Y, C = C, L = L)
    setTxtProgressBar(pb, (2*i-1)/boot/2)
    res <- med_iptw.decomp(an_dat, A = A, M = M, Y = Y, regtype = regtype, noint = noint, mlvl = mlvl)
    setTxtProgressBar(pb, i/boot)

    ar_te[i,] <- res$te
    if(noint == TRUE) ar_cde[i,] <- res$cde else ar_cde[i,,] <- res$cde
    if(noint == TRUE) ar_cie[i,] <- res$cde else ar_cie[i,,] <- res$cie
    if(noint == TRUE) ar_pm[i,] <- res$pm else ar_pm[i,,] <- res$pm
  }
  if(noint == TRUE) over <- 2 else over <- c(3,2)
  te_95 <- lapply(quants, function (i) apply(ar_te, 2, quantile, probs = i, na.rm = TRUE))
  cde_95 <- lapply(quants, function (i) apply(ar_cde, over, quantile, probs = i, na.rm = TRUE))
  pm_95 <- lapply(quants, function (i) apply(ar_pm, over, quantile, probs = i, na.rm = TRUE))
  cie_95 <- lapply(quants, function (i) apply(ar_cie, over, quantile, probs = i, na.rm = TRUE))

  return(list(
    arrays = list(te = ar_te, cde = ar_cde, cie = ar_cie, pm = ar_pm),
    te = te_95,
    cde = cde_95,
    pm = pm_95,
    cie = cie_95)
  )
}
med_iptw.mkdat <- function(orig_dat, A, M, Y, C = "", L = ""){
  #replicate input with parentheses
  C <- paste0("(",C,")")
  #build formulas
  med_iptw.mkform <- function(A, M, C, L){
    a_form <- as.formula(paste0(A, "~", paste0(C, collapse = "+")))
    m_form1 <- as.formula(paste0(M, "~", A))
    m_form2 <- as.formula(paste0(M, "~", paste(A, paste(C, collapse = "+"), paste(L, collapse = "+"), sep= "+")))
    l_form <- as.formula(paste0(L, "~", paste(A, paste(C, collapse = "+"), sep = "+")))
    return(list(af = a_form , mf1 = m_form1, mf2 = m_form2, lf = l_form))
  }
  #build models
  med_iptw.mkmods <- function(dat, flist){
    a_mod <- multinom(data = dat, flist$af, trace = FALSE)
    m_mod1 <- multinom(data = dat, flist$mf1, trace = FALSE)
    m_mod2 <- multinom(data = dat, flist$mf2, trace = FALSE)
    l_mod <- multinom(data = dat, flist$lf, trace = FALSE)
    return(list(am = a_mod, mm1 = m_mod1, mm2 = m_mod2, lm = l_mod))
  }
  #get propensity scores
  med_iptw.mkprop <- function(dat, A, M, mlist){
    pa <- (table(dat[[A]])/nrow(dat))[match(dat[[A]], levels(dat[[A]]))]
    pa_c <- predict(object = mlist$am, newdata = dat, type = "probs") %>% fmatch(dat[[A]])
    pm_a <- predict(object = mlist$mm1, newdata = dat, type = "probs") %>% fmatch(dat[[M]])
    pm_lac <- predict(object = mlist$mm2, newdata = dat, type = "probs") %>% fmatch(dat[[M]])

    return(list(pa = pa, pac = pa_c, pma = pm_a, pmlac = pm_lac))
  }
  #make weight
  med_iptw.mkweight <- function(dat, plist){
    dat$ipw_conf <- (plist$pa / plist$pac) %>% unname
    dat$ipw_med <- (dat$ipw_conf * plist$pma / plist$pmlac) %>% unname
    return(dat)
  }
  flist <- med_iptw.mkform(A, M, C, L)
  mlist <- med_iptw.mkmods(orig_dat, flist)
  plist <- med_iptw.mkprop(orig_dat, A, M, mlist)
  an_dat <- med_iptw.mkweight(orig_dat, plist)
  return(an_dat)
}

#' @import dplyr
#' @importFrom nnet multinom
#' @import magrittr
med_iptw.decomp <- function(an_dat, A, M, Y, regtype = "gaussian", noint = FALSE, mlvl = NULL){
  #build models
  med_iptw.mkymod <- function(dat, regtype, noint, Y, A, M){
    #Make formula
    if(noint == TRUE){ form <- as.formula(paste0(Y, "~", paste(A, M, sep = "+")))
    }else form <- as.formula(paste0(Y, "~", paste(A, M, paste0(A, "*", M), sep = "+")))

    #models
    ymod_cde <- glm(data = an_dat, form, weights = ipw_med, family = regtype)
    ymod_te <- glm(data = an_dat, as.formula(paste0(Y, "~", A)), weights = ipw_conf, family = regtype)
    return(list(cde = ymod_cde, te = ymod_te))
  }
  #get lists of results
  med_iptw.break <- function(dat, A, M, noint, ycoef){
    alen <- levels(dat[[A]]) %>% length - 1
    mlen <- levels(dat[[M]]) %>% length - 1
    if(noint == TRUE) return(ycoef[(1:alen)+1])
    #else
    pos <- cumsum(c(1,alen,mlen,alen*mlen))
    beta <- lapply(1:3, function (i) ycoef[(pos[i]+1):pos[i+1]])
    if(mlen > 1){
      beta[[3]] <- lapply((1:mlen)-1, function (i) beta[[3]][(i*alen+1):(alen*(i+1))]) %>%
                   do.call(rbind, .) %>% as.matrix
      colnames(beta[[3]]) <- levels(dat[[A]])[-1]
      rownames(beta[[3]]) <- levels(dat[[M]])[-1]
    }
    return(beta)
  }
  #get CDE(M) for all M's
  med_iptw.ints <- function(coef_l, mlvl){
    mlen <- length(coef_l[[2]])
    fullb <- (t(replicate(nrow(mlvl), coef_l[[1]])) + mlvl%*%coef_l[[3]]) %>% as.matrix
    return(fullb)
  }

  #Default mlvl if needed
  if(noint == FALSE) if(is.null(mlvl)) mlvl <- rbind(rep(0, mlen), diag(mlen), rep(1/mlen, mlen))

  #Make model of Y
  ymods <- med_iptw.mkymod(an_dat, regtype, noint, Y, A, M)
  #Break into dummy coeffs
  cde <- med_iptw.break(an_dat, A, M, noint, ymods$cde$coefficients)
  te <- ymods$te$coefficients[-1]

  #Return NDE and TE if noint
  if(noint == TRUE) return(list(mods = ymods, cde = cde, te = te, cie = (te - cde), pm = (te-cde)/te))

  #Get CDE(M) for each M if interaction
  cde_ints <- med_iptw.ints(cde, mlvl)
  pm <- (t(replicate(nrow(cde_ints), te)) - cde_ints)/t(replicate(nrow(cde_ints), te))
  cie <- (t(replicate(nrow(cde_ints), te)) - cde_ints)

  #Return CDE(M) and TE if interaction
  return(list(mods = ymods, cde = cde_ints, te = te, pm = pm, cie = cie))
}
