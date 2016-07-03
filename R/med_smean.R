med_smean <- function(dat, A, Y, M, C = NULL, L = NULL, boot = 10, nmin = 10, mlvl = NULL){
  if("yp" %in% names(dat)) stop("ensure no variable named \"yp\" in dataframe ")

  acol <- deparse(substitute(A)) %>% match(names(dat))
  mcol <- deparse(substitute(M)) %>% match(names(dat))
  ycol <- deparse(substitute(Y)) %>% match(names(dat))

  if(is.factor(dat[[acol]])) alen <- nlevels(dat[[acol]]) - 1 else alen <- 1
  if(is.factor(dat[[mcol]])) mlen <- nlevels(dat[[mcol]]) - 1 else mlen <- 1

  ymod <- lm(data = dat, substitute(Y ~ A + M + A*M + C + L))
  intpos <- grep(":", names(ymod$coefficients))[1]
  k3 <- ymod$coefficients[intpos:(intpos+alen*mlen-1)]
  k2 <- ymod$coefficients[(alen+2) : (alen+mlen+1)]

  if(is.factor(dat[[acol]]) && is.factor(dat[[mcol]])){
    intmat <- cbind(rep(0, alen+1), rbind(0, matrix(k3, alen, mlen)))
    yp <- dat[[ycol]] - intmat[cbind(dat[[acol]], dat[[mcol]])] - c(0,k2)[as.numeric(dat[[mcol]])]
  }else if(is.factor(dat[[acol]])) {
    yp <- dat[[ycol]] - dat[[mcol]] * c(0,k3)[as.numeric(dat[[acol]])] - k2*dat[[mcol]]
  }else if(is.factor(dat[[mcol]])) {
    yp <- dat[[ycol]] - dat[[acol]] * c(0,k3)[as.numeric(dat[[mcol]])] - c(0,k2)[as.numeric(dat[[mcol]])]
  }else yp <- dat[[ycol]] - dat[[mcol]]*dat[[acol]]*k3

  #Lookup yp in parent environment to avoid collision with dataframe
  ymod2 <- lm(data = oai, substitute(eval(quote(yp), parent.frame()) ~ A + C))
  g1 <- ymod2$coefficients[1:alen + 1]

  if(is.null(mlvl)){
    if(is.factor(dat[[mcol]])) mlvl <- (table(dat[[mcol]]) %>% prop.table)[-1]
    else mlvl <- mean(dat[[mcol]], na.rm = TRUE)
  }
  cde <- g1 + matrix(k3, alen, mlen) %*% mlvl
  return(list(
    g1 = g1,
    k2 = k2,
    k3 = k3,
    ymod1 = ymod,
    ymod2 = ymod2,
    cde = cde
  ))
}

test <- med_smean(oai, bmi, hspss, obcat, age + age*p02sex, mlvl = rep(1/4, 3))


