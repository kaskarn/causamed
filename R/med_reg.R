#' \code{med_reg} Gives controlled direct effect, natural indirect effect, natural direct effect using regression-based formulas
#' affording the best statistical power
#'
#' @param df Data frame
#' @param X exposure, can be continuous or categorical
#' @param M mediator, must be continuous for the moment
#' @param Y outcome, must be continuous for the moment
#' @param C confounders of exposure-outcome or mediator-outcome association. Can take any form
#' @param noint supress exposure-mediator interaction term
#' @param mlvl vector of values at which to compute a controlled direct effect
#' @param delta whether the delta method should be used to get confidence intervals
#' @param ref referent level of exposure, 0 by default
#' @param treat treatment level of exposure, 1 by default
#' @export
med_reg <- function(df, X, M, Y, C = "", noint = FALSE, mlvl = NULL, delta = FALSE, ref = 0, treat = 1){
# Twoway decomposition indirect and direct effects  (VDW)
# for continuous mediator and outcome

# Sets position and lengths of coefficient vectors
  if(is.factor(df[[X]])) xlen <- levels(df[[X]]) %>% length - 1 else xlen <- 1
  if(is.factor(df[[M]])) mlen <- levels(df[[M]]) %>% length - 1 else mlen <- 1
  c_mu <- lapply(C, function (i) if(is.factor(df[[i]])) table(df[[i]])[-1]/nrow(df) else mean(df[[i]], na.rm = TRUE)) %>% unlist

#  throws error if needed
  if((xlen > 1 | mlen > 1) & delta == TRUE) stop("Standard error must be estimated by bootstrapping for polytomous treatment or mediator")

  #adds parentheses to allow flexible model specification
  for(i in c("X","M","Y","C")) assign(i, paste0("(",get(i),")"))

  #model of mediator
  m1 <- lm(data = df, as.formula(paste(M, paste(X, paste(C, collapse = "+"), sep = "+"),  sep = "~")))

  #model of outcome
  if(noint) { m2 <- lm(data = df, as.formula(paste(Y, paste(X, M, paste(C, collapse = "+"), sep = "+"),  sep = "~")))
  }else m2 <-  lm(data = df, as.formula(paste(Y, paste(X, M, paste(C, collapse = "+"), paste(X,M, sep = "*"), sep = "+"),  sep = "~")))

  #breaks up results into different variable coefficients
  p1 <- cumsum(c(1,xlen,length(c_mu)))
  p2 <- cumsum(c(1,xlen,mlen,length(c_mu),mlen*xlen))

  #creates shorthands for coefficients to be used
  b0 <- m1$coefficients[1]
  for(i in 1:4) assign(paste0("t",i), c(m2$coefficients, rep(0,xlen*mlen*noint))[(p2[i]+1):p2[i+1]] %>% unname)
  for(i in 1:2) assign(paste0("b",i), m1$coefficients[(p1[i]+1):p1[i+1]] %>% unname)

  #output results, optionally with controlled direct effects
  out <- list(
    nde = (t1 + t4*b0 + t4*b1*ref + t4*b2%*%c_mu)*(treat - ref),
    nie = (t2*b1 + t4*b1*treat)*(treat - ref),
    #throws mediated moderation in there, although not two-way decomposition
    mint = (t4*b1)*(treat - ref)
  )
  if (!is.null(mlvl)) out <- list(out, cde = do.call(rbind, lapply(mlvl, function (i) ((t1 + t4*i)*(treat - ref)) )) )

  #Delta method
  if(delta == FALSE){ return(out)
  }else{
    sigma <- cbind(
      rbind(vcov(m1), matrix(0, length(m2$coefficients), length(m1$coefficients))),
      rbind(vcov(m2), matrix(0, length(m1$coefficients), length(m2$coefficients)))
    )
    #same for Gamma
    gamma <- list(
      nde = c(t4, t4*ref, t4*c_mu, 0, 1, 0, b0 + b1*ref + b2%*%c_mu, rep(0, length(c_mu))),
      nie = c(0, t2 + t4*treat, rep(0, length(c_mu)), 0, 0, b1, b1*treat, rep(0, length(c_mu)))
    )
    if (!is.null(mlvl)) gamma <- list(gamma, cde = do.call(rbind, lapply(mlvl, function (i) c(rep(0, length(c_mu)+3), 1, 0, i, rep(0, length(c_mu))))))
    se <- lapply(c("nde", "nie"), function (i) sqrt(gamma[[i]]%*%sigma%*%gamma[[i]]) * abs(treat - ref))
  }

  #end
  return(list(out = out, se = se))
}
