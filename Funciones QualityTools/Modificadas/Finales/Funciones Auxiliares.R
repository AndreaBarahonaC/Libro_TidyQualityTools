library(qualityTools)
library(ggplot2)
library(patchwork)

############################## funciones .r ##############################
# .cg --------------
.cg=function(x, target, tolerance=c(-1,1), ref.interval, facCg, facCgk)
{
  if (missing(x))
    stop("x must be given as a vector")
  if (missing(target))
    target = mean(x)
  if (missing(ref.interval))
    ref.interval = pnorm(3) - pnorm(-3)
  sd = sd(x)
  mean = mean(x)
  ref.ar = qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval, mean, sd)
  if (missing(facCg))
    facCg = 0.2
  if (missing(facCgk))
    facCgk = 0.1
  if (missing(tolerance))
    stop("tolerance is missing!")
  if (length(tolerance) != 2)
    stop("tolerance has wrong length")
  Cg = (facCg * (abs(diff(tolerance))))/ref.ar
  Cgk = (facCgk * (abs(diff(tolerance))) - abs(target - mean))/(ref.ar/2)
  return(list(Cg, Cgk))
}
# internals.r --------------
.confintbeta= function(xs, thethas, varmatrix, alpha) {

  th1= thethas[[1]]
  th2 =  thethas[[2]]

  prozent=ppoints(xs)

  perzentile=qbeta(prozent, th1, th2)

  h=1e-6
  dFdth1=(qbeta(prozent, th1, th2)-qbeta(prozent, th1+h, th2))/h
  dFdth2=(qbeta(prozent, th1, th2)-qbeta(prozent, th1, th2+h))/h

  Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
  zalpha=qnorm(1-alpha/2)
  halfwidth = zalpha*sqrt(Var)

  lci=perzentile-halfwidth
  uci=perzentile+halfwidth

  bounds = list(lci, uci, perzentile)


  return (bounds)
}
.confintcauchy = function(xs, thethas, varmatrix, alpha) {

  th1= thethas[[1]]
  th2 =  thethas[[2]]

  prozent=ppoints(xs)

  perzentile=qcauchy(prozent, th1, th2)

  h=1e-6
  dFdth1=(qcauchy(prozent, th1, th2)-qcauchy(prozent, th1+h, th2))/h
  dFdth2=(qcauchy(prozent, th1, th2)-qcauchy(prozent, th1, th2+h))/h

  Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
  zalpha=qnorm(1-alpha/2)
  halfwidth = zalpha*sqrt(Var)

  lci=perzentile-halfwidth
  uci=perzentile+halfwidth

  bounds = list(lci, uci, perzentile)

  return (bounds)
}
.confintexp=function(xs, thethas, varmatrix, alpha) {
  lambda=thethas[[1]]
  prozent=ppoints(xs)
  perzentile=qexp(prozent, lambda)
  logPerzentile = log(perzentile)
  zalpha=qnorm(1-alpha/2)
  halfwidth = zalpha*sqrt(varmatrix[[1, 1]]/lambda^2)
  lci = exp(logPerzentile - halfwidth);
  uci = exp(logPerzentile + halfwidth);

  bounds = list(lci, uci, perzentile)

  return (bounds)

}
.confintgamma= function(xs, thethas, varmatrix, alpha) {
  th1= thethas[[1]]
  th2 =  thethas[[2]]
  prozent=ppoints(xs)
  perzentile=qgamma(prozent, th1, th2)

  h=1e-6
  dFdth1=(qgamma(prozent, th1, th2)-qgamma(prozent, th1+h, th2))/h
  dFdth2=(qgamma(prozent, th1, th2)-qgamma(prozent, th1, th2+h))/h

  Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
  zalpha=qnorm(1-alpha/2)
  halfwidth = zalpha*sqrt(Var)

  lci=perzentile-halfwidth
  uci=perzentile+halfwidth

  bounds = list(lci, uci, perzentile)

  return (bounds)
}
.confintlnorm=function(xs, thethas, varmatrix, alpha){
  th1= thethas[[1]]
  th2 =  thethas[[2]]

  prozent=ppoints(xs)

  perzentile=qlnorm(prozent, th1, th2)

  zp=qnorm(prozent)

  varPerzentile = varmatrix[[1, 1]]+2*varmatrix[[1, 2]]*zp+varmatrix[[2, 2]]*zp*zp

  zalpha=qnorm(1-alpha/2)
  lci=log(perzentile)-zalpha*sqrt(varPerzentile)
  uci=log(perzentile)+zalpha*sqrt(varPerzentile)

  bounds = list(exp(lci), exp(uci), perzentile)

  return (bounds)

}
.confintlogis= function(xs, thethas, varmatrix, alpha) {

  th1= thethas[[1]]
  th2 =  thethas[[2]]

  prozent=ppoints(xs)


  perzentile=qlogis(prozent, th1, th2)

  h=1e-6
  dFdth1=(qlogis(prozent, th1, th2)-qlogis(prozent, th1+h, th2))/h
  dFdth2=(qlogis(prozent, th1, th2)-qlogis(prozent, th1, th2+h))/h

  Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
  zalpha=qnorm(1-alpha/2)
  halfwidth = zalpha*sqrt(Var)

  lci=perzentile-halfwidth
  uci=perzentile+halfwidth

  bounds = list(lci, uci, perzentile)

  return (bounds)
}
.confintnorm=function(xs, thethas, varmatrix, alpha){

  prozent=ppoints(xs)

  zp=qnorm(prozent)
  perzentile=qnorm(prozent, thethas[[1]], thethas[[2]])

  varPerzentile = varmatrix[[1, 1]]+2*varmatrix[[1, 2]]*zp+varmatrix[[2, 2]]*zp*zp

  zalpha=qnorm(1-alpha/2)
  lci=perzentile-zalpha*sqrt(varPerzentile)
  uci=perzentile+zalpha*sqrt(varPerzentile)

  bounds = list(lci, uci, perzentile)

  return (bounds)
}
.confintweibull= function(xs, thethas, varmatrix, alpha) {
  th1= thethas[[1]]
  th2 =  thethas[[2]]

  prozent=ppoints(xs)

  perzentile=qweibull(prozent, th1, th2)
  q=-log(1-prozent)
  logPerzentile=log(perzentile)
  logq=log(q)
  dB=1/th2
  dA=-1/(th1^2)

  Var = varmatrix[[1, 1]]*(dA*logq)^2 + 2*varmatrix[[1, 2]]*dB*dA*logq + varmatrix[[2, 2]]*dB^2
  zalpha=qnorm(1-alpha/2)
  halfwidth = zalpha*sqrt(Var)


  lci=exp(logPerzentile-halfwidth)
  uci=exp(logPerzentile+halfwidth)

  bounds = list(lci, uci, perzentile)

  # print(data.frame(prozent, uci, perzentile, lci))

  return (bounds)
}
.gamma3 = function(data) {
  n=length(data)
  data=sort(data)

  pEmp= (seq(1:n)-0.5)/n

  weight = 1 / sqrt(pEmp*(1-pEmp))

  thld = .99*min(data)
  shape=1
  scale=1

  gammaEst = function(param) {
    return( sum(weight*(pgamma(data-param[3], shape = exp(param[1]), scale = exp(param[2]))-pEmp)^2) )
  }

  paramEst = optim(c(shape, scale, thld), gammaEst, method = "Nelder-Mead")
  paramEst = paramEst$par
  return(list(shape = exp(paramEst[1]), scale = exp(paramEst[2]), threshold = paramEst[3]))
}
.lognormal3 = function(data) {

  n=length(data)
  data=sort(data)
  #compute the empirical cumulative distribution function of the data
  pEmp= (seq(1:n)-0.5)/n
  # will minimize the weighted sum of squared distances
  # so compute weights
  weight = 1 / sqrt(pEmp*(1-pEmp))

  # initial values for optimization
  thld = .99*min(data)
  mu0 = mean(log(data-thld))
  sigma0 = sd(log(data-thld))


  lnEst = function(param) {
    return( sum(weight*(plnorm(data-param[3], meanlog = param[1], sdlog = exp(param[2]))-pEmp)^2) )
  }

  logSigma0=log(sigma0)
  # optimize gammaEst using optim function
  paramEst = optim(c(mu0,logSigma0, thld), lnEst, method = "Nelder-Mead")
  param = paramEst$par

  return(list(meanlog = param[1], sdlog = exp(param[2]), threshold = param[3]))
}
.weibull3 = function(x){
  if(any(x < 0))
    stop("x must be positive")

  n = length(x)
  x = sort(x)
  p = ((1:n)-0.5)/n
  interval = c(0.75*min(x), 0.9999*min(x))

  wb3RSquared = function(th)
  {
    return(summary(lm(log(x-th) ~ log(-log(1-p))))$r.squared)
  }

  th = (optimize(wb3RSquared, interval = interval, maximum = TRUE))$maximum

  lm.1 = lm(log(x-th) ~ log(-log(1-p)))
  estimates = list(shape = 1/coef(lm.1)[[2]], scale = exp(coef(lm.1)[[1]]), threshold = th)
  return(estimates)
}
# .charToDistFunc -----------------------------
.charToDistFunc = function(distribution, type = "q") {
  fun = NULL
  if (identical("beta", distribution))
    fun = eval(parse(text = paste(type, "beta", sep = "")))
  if (identical("cauchy", distribution))
    fun = eval(parse(text = paste(type, "cauchy", sep = "")))
  if (identical("chi-squared", distribution))
    fun = eval(parse(text = paste(type, "chisq", sep = "")))
  if (identical("exponential", distribution))
    fun = eval(parse(text = paste(type, "exp", sep = "")))
  if (identical("f", distribution))
    fun = eval(parse(text = paste(type, "f", sep = "")))
  if (identical("geometric", distribution))
    fun = eval(parse(text = paste(type, "geom", sep = "")))
  if (identical("log-normal", distribution) || identical("lognormal", distribution))         ####
    fun = eval(parse(text = paste(type, "lnorm", sep = "")))
  if (identical("log-normal3", distribution) || identical("lognormal3", distribution))       ####
    fun = eval(parse(text = paste(type, "lnorm3", sep = "")))                              ####
  if (identical("logistic", distribution))
    fun = eval(parse(text = paste(type, "logis", sep = "")))
  if (identical("negative binomial", distribution))
    fun = eval(parse(text = paste(type, "nbinom", sep = "")))
  if (identical("normal", distribution))
    fun = eval(parse(text = paste(type, "norm", sep = "")))
  if (identical("poisson", distribution))
    fun = eval(parse(text = paste(type, "pois", sep = "")))
  if (identical("t", distribution))
    fun = eval(parse(text = paste(type, "t", sep = "")))
  if (identical("weibull", distribution))
    fun = eval(parse(text = paste(type, "weibull", sep = "")))
  if (identical("weibull3", distribution))                                                   ####
    fun = eval(parse(text = paste(type, "weibull3", sep = "")))                            ####
  if (identical("gamma", distribution))
    fun = eval(parse(text = paste(type, "gamma", sep = "")))
  if (identical("gamma3", distribution))
    fun = eval(parse(text = paste(type, "gamma3", sep = "")))
  return(fun)
}
# .myADtest ---------------------
.myADTest = function(x, distribution, ...) {                                                                     ####   .MYADTESTS-FUNCTION
  #require(MASS, quietly = TRUE)
  if (missing(distribution))
    distribution = "normal"
  data.name = names(x)
  if (is.data.frame(x))
    x = x[, 1]
  dots = list()
  parameter = NULL
  smaller = NULL
  pFun = NULL
  tableValue = FALSE
  A = 0
  x <- sort(x[complete.cases(x)])
  n = length(x)
  if (n < 8)
    stop("sample size must be greater than 7")
  if (n > 40)
    warning("sample size is greater than 40")
  if (is.character(distribution)) {
    pFun = .charToDistFunc(distribution, type = "p")
    distribution = tolower(distribution)
    if (is.null(pFun))
      stop(paste(deparse(substitute(distribution)), " is not supported!"))
  }
  else {
    pFun = match.fun(distribution)
  }                                                                      ####
  if (length(dots) == 0) {
    fittedDistr = MASS::fitdistr(x, distribution)
    parameter = fittedDistr$estimate
    if (distribution == "normal") {
      parameter["mean"] = mean(x)
      parameter["sd"] = sd(x)
    }
    p = do.call(pFun, c(list(x), as.list(parameter)))
  }
  else {
    p = pFun(x, ...)
  }
  h = (2 * seq(1:n) - 1) * (log(p) + log(1 - rev(p)))
  A = -n - mean(h)
  AA = (1 + 0.75/n + 2.25/n^2) * A
  if (AA < 0.2) {
    pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
  }
  else if (AA < 0.34) {
    pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
  }
  else if (AA < 0.6) {
    pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
  }
  else {
    pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
  }
  if (identical(distribution, "cauchy")) {
    pval = NA
  }
  if (identical(distribution, "beta")) {
    pval = NA
  }
  if (identical(distribution, "chi-squared")) {
    pval = NA
  }
  if (identical(distribution, "f")) {
    pval = NA
  }
  if (identical(distribution, "t")) {
    pval = NA
  }
  if (identical(distribution, "geometric")) {
    pval = NA
  }
  if (identical(distribution, "poisson")) {
    pval = NA
  }
  if (identical(distribution, "negative-binomial")) {
    pval = NA
  }
  if (identical(distribution, "weibull")) {
    AWei = A * (1 + 1/sqrt(n))
    tableValue = TRUE
    smaller = TRUE
    if (AWei < 0.474) {
      pval = 0.25
      smaller = FALSE
    }
    if (AWei >= 0.474)
      pval = 0.25
    if (AWei >= 0.637)
      pval = 0.1
    if (AWei >= 0.757)
      pval = 0.05
    if (AWei >= 0.877)
      pval = 0.025
    if (AWei >= 1.038)
      pval = 0.01
  }
  if (identical(distribution, "exponential")) {
    AExp = A * (1 + 0.6/n)
    pval = NA
    if (0.95 < AExp) {
      pval = exp(0.731 - 3.009 * AExp + 0.15 * AExp^2)
    }
    if (0.51 < AExp & AExp < 0.95) {
      pval = exp(0.9209 - 3.353 * AExp + 0.3 * AExp^2)
    }
    if (0.26 < AExp & AExp < 0.51) {
      pval = 1 - exp(-6.1327 + 20.218 * AExp - 18.663 * AExp^2)
    }
    if (AExp < 0.26) {
      pval = 1 - exp(-12.2204 + 67.459 * AExp - 110.3 * AExp^2)
    }
  }
  if (identical(distribution, "logistic")) {
    ALogist = A * (1 + 0.25/n)
    tableValue = TRUE
    smaller = TRUE
    if (ALogist < 0.426) {
      pval = 0.25
      smaller = FALSE
    }
    if (ALogist >= 0.426) {
      pval = 0.25
    }
    if (ALogist >= 0.563) {
      pval = 0.1
    }
    if (ALogist >= 0.66) {
      pval = 0.05
    }
    if (ALogist >= 0.769) {
      pval = 0.025
    }
    if (ALogist >= 0.906) {
      pval = 0.01
    }
    if (ALogist >= 1.1) {
      pval = 0.005
    }
  }
  if (identical(distribution, "gamma")) {
    tableValue = TRUE
    gammaDF = data.frame(c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, Inf), c(0.486, 0.477, 0.475,
                                                                        0.473, 0.472, 0.472, 0.471, 0.471, 0.471, 0.47, 0.47, 0.47), c(0.657, 0.643, 0.639, 0.637,
                                                                                                                                       0.635, 0.635, 0.634, 0.633, 0.633, 0.632, 0.632, 0.631), c(0.786, 0.768, 0.762, 0.759,
                                                                                                                                                                                                  0.758, 0.757, 0.755, 0.754, 0.754, 0.754, 0.753, 0.752), c(0.917, 0.894, 0.886, 0.883,
                                                                                                                                                                                                                                                             0.881, 0.88, 0.878, 0.877, 0.876, 0.876, 0.875, 0.873), c(1.092, 1.062, 1.052, 1.048,
                                                                                                                                                                                                                                                                                                                       1.045, 1.043, 1.041, 1.04, 1.039, 1.038, 1.037, 1.035), c(1.227, 1.19, 1.178, 1.173,
                                                                                                                                                                                                                                                                                                                                                                                 1.17, 1.168, 1.165, 1.164, 1.163, 1.162, 1.161, 1.159))
    names(gammaDF) = c("m", 0.75, 0.9, 0.95, 0.975, 0.99, 0.995)
    critCheck <- gammaDF[min(which(gammaDF$m >= parameter["shape"])), 2:length(gammaDF)] > A
    if (any(critCheck)) {
      firPos <- min(which(critCheck))
    }
    else {
      firPos <- length(gammaDF)
    }
    if (firPos == 1) {
      pValue <- 1 - as.numeric(names(gammaDF)[2])
      pval = pValue
      pValue <- paste(">", pValue)
      smaller = FALSE
    }
    else {
      pValue <- 1 - as.numeric(names(gammaDF)[firPos])
      pval = pValue
      pValue <- paste("<=", pValue)
      smaller = TRUE
    }
  }
  out = list()
  out$data.name = data.name
  out$statistic = as.vector(data.frame(A = A))
  out$parameter = parameter
  out$p.value = as.vector(data.frame(p = pval))
  out$smaller = smaller
  out$tableValue = tableValue
  out$conf.int = NULL
  out$estimate = NULL
  temp = NULL
  if (is.character(distribution))
    temp = as.vector(distribution)
  else temp = deparse(substitute(distribution))
  names(temp) = "distribution"
  out$null.value = temp
  out$method = paste("Anderson Darling Test for", temp, "distribution")
  out$class = "adtest"
  invisible(out)
}
# .xyLimits -------------------
.xyLimits = function(distrCollection, lowerquantile = 0.001, upperquantile = 0.999) {
  x <- NULL
  y <- NULL
  for (i in seq_along(distrCollection$distr)) {
    object <- distrCollection$distr[[i]]
    xValues <- object$x
    parameters <- object$parameters
    distr <- object$name
    qFun <- .charToDistFunc(distr, type = "q")
    dFun <- .charToDistFunc(distr, type = "d")
    lq <- do.call(qFun, c(list(lowerquantile), as.list(parameters)))
    uq <- do.call(qFun, c(list(upperquantile), as.list(parameters)))
    x <- range(x, xValues, lq, uq)
    histObj <- hist(xValues, plot = FALSE)
    xPoints <- seq(x[1], x[2], length = 200)
    yPoints <- do.call(dFun, c(list(xPoints), as.list(parameters)))
    y <- range(y, 0, histObj$density, yPoints)
  }
  invisible(list(xlim = x, ylim = y))
}
# .sdSg, lfkp, lfrm  ---------------------
.sdSg = function(x, grouping = NULL, method = c("NOWEIGHT", "MVLUE", "RMSDF"), na.rm = TRUE, DB = TRUE) {
  DB = FALSE
  if (!is.data.frame(x) && !is.vector(x) && is.numeric(x))
    stop("x needs to be either a data.frame or a vector and numeric")
  if (is.null(grouping)) {
    if (is.data.frame(x))
      return(sd(x[, 1]))
    else return(sd(x))
  }
  else grouping = as.data.frame(grouping)
  group = unique(grouping)
  sdVec = numeric(length = length(group))
  for (i in 1:nrow(group)) {
    if (is.data.frame(x))
      temp = x[group[i, 1] == grouping[, 1], 1]
    if (is.vector(x))
      temp = x[group[i, 1] == grouping[, 1]]
    sdVec[i] = sd(temp, na.rm = T)/.c4(length(temp[!is.na(temp)]))
    if (DB) {
      print(group[i, 1])
      print(temp)
      print(length(temp[!is.na(temp)]))
    }
  }
  if (DB) {
    print(paste("std.dev: ", mean(sdVec)))
    print(sdVec)
  }
  return((mean(sdVec)))
}
.lfkp = function(wholeList, filterList) {
  if (!is.list(wholeList))
    stop(paste(deparse(substitute(wholeList)), "is not a list!"))
  if (length(wholeList) == 0)
    return(wholeList)
  if (!is.list(filterList))
    stop(paste(deparse(substitute(filterList)), "is not a list!"))
  if (length(filterList) == 0)
    return(filterList)
  logVec = lapply(names(wholeList), "%in%", names(filterList))
  filteredList = wholeList[unlist(logVec)]
  return(filteredList)
}
.lfrm = function(wholeList, filterList) {
  if (!is.list(wholeList))
    stop(paste(deparse(substitute(wholeList)), "is not a list!"))
  if (length(wholeList) == 0)
    return(wholeList)
  if (!is.list(filterList))
    stop(paste(deparse(substitute(filterList)), "is not a list!"))
  if (length(filterList) == 0)
    return(wholeList)
  logVec = lapply(names(wholeList), "%in%", names(filterList))
  filteredList = wholeList[!unlist(logVec)]
  return(filteredList)
}

############################## FUNCIONES ##############################
# Funcion Distribution ----
distribution <- function(x = NULL, distrib = "weibull", start, ...) {
  distr_coll <- DistrCollection$new()
  if (is.character(distrib))
    distrib = tolower(distrib)
  allDistr = c("beta", "cauchy", "chi-squared", "exponential", "f", "gamma", "geometric", "log-normal", "logistic", "negative binomial", "normal", "poisson",
               "t", "weibull")
  if (distrib %in% allDistr){
    distrVec = distrib
  }
  else{distrVec = c("normal")}
  if (identical(distrib, "all"))
    distrVec = allDistr
  if (identical(distrib, "quality"))
    distrVec = c("normal", "log-normal", "exponential", "weibull")
  for (i in seq(along = distrVec)) {
    temp <- suppressWarnings(FitDistr(x, densfun = distrVec[i]))
    fit <- Distr$new(x = x,
                     name = distrVec[i],
                     parameters = temp$estimate,
                     sd = temp$sd,
                     loglik = temp$loglik,
                     n = length(x))
    distr_coll$add(fit)
  }
  return(distr_coll)
}

# Funcion FitDistr ----
FitDistr <- function (x, densfun, start, ...){
  myfn <- function(parm, ...) -sum(log(dens(parm, ...)))
  mylogfn <- function(parm, ...) -sum(dens(parm, ..., log = TRUE))
  mydt <- function(x, m, s, df, log) dt((x - m)/s, df, log = TRUE) -
    log(s)
  Call <- match.call(expand.dots = TRUE)
  if (missing(start))
    start <- NULL
  dots <- names(list(...))
  dots <- dots[!is.element(dots, c("upper", "lower"))]
  if (missing(x) || length(x) == 0L || mode(x) != "numeric")
    stop("'x' must be a non-empty numeric vector")
  if (any(!is.finite(x)))
    stop("'x' contains missing or infinite values")
  if (missing(densfun) || !(is.function(densfun) || is.character(densfun)))
    stop("'densfun' must be supplied as a function or name")
  control <- list()
  n <- length(x)
  if (is.character(densfun)) {
    distname <- tolower(densfun)
    densfun <- switch(distname, beta = dbeta, cauchy = dcauchy,
                      `chi-squared` = dchisq, exponential = dexp, f = df,
                      gamma = dgamma, geometric = dgeom, `log-normal` = dlnorm,
                      lognormal = dlnorm, logistic = dlogis, `negative binomial` = dnbinom,
                      normal = dnorm, poisson = dpois, t = mydt, weibull = dweibull,
                      NULL)
    if (is.null(densfun))
      stop("unsupported distribution")
    if (distname %in% c("lognormal", "log-normal")) {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "log-Normal"), domain = NA)
      if (any(x <= 0))
        stop("need positive values to fit a log-Normal")
      lx <- log(x)
      sd0 <- sqrt((n - 1)/n) * sd(lx)
      mx <- mean(lx)
      estimate <- c(mx, sd0)
      sds <- c(sd0/sqrt(n), sd0/sqrt(2 * n))
      names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      vc <- matrix(c(sds[1]^2, 0, 0, sds[2]^2), ncol = 2,
                   dimnames = list(names(sds), names(sds)))
      names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dlnorm(x, mx, sd0, log = TRUE))), class = "FitDistr"))
    }
    if (distname == "normal") {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "Normal"), domain = NA)
      sd0 <- sqrt((n - 1)/n) * sd(x)
      mx <- mean(x)
      estimate <- c(mx, sd0)
      sds <- c(sd0/sqrt(n), sd0/sqrt(2 * n))
      names(estimate) <- names(sds) <- c("mean", "sd")
      vc <- matrix(c(sds[1]^2, 0, 0, sds[2]^2), ncol = 2,
                   dimnames = list(names(sds), names(sds)))
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dnorm(x, mx,
                                                                 sd0, log = TRUE))), class = "FitDistr"))
    }
    if (distname == "poisson") {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "Poisson"), domain = NA)
      estimate <- mean(x)
      sds <- sqrt(estimate/n)
      names(estimate) <- names(sds) <- "lambda"
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("lambda",
                                                              "lambda"))
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dpois(x, estimate,
                                                                 log = TRUE))), class = "FitDistr"))
    }
    if (distname == "exponential") {
      if (any(x < 0))
        stop("Exponential values must be >= 0")
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "exponential"), domain = NA)
      estimate <- 1/mean(x)
      sds <- estimate/sqrt(n)
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("rate",
                                                              "rate"))
      names(estimate) <- names(sds) <- "rate"
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dexp(x, estimate,
                                                                log = TRUE))), class = "FitDistr"))
    }
    if (distname == "geometric") {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "geometric"), domain = NA)
      estimate <- 1/(1 + mean(x))
      sds <- estimate * sqrt((1 - estimate)/n)
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("prob",
                                                              "prob"))
      names(estimate) <- names(sds) <- "prob"
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dgeom(x, estimate,
                                                                 log = TRUE))), class = "FitDistr"))
    }
    if (distname == "weibull" && is.null(start)) {
      if (any(x <= 0))
        stop("Weibull values must be > 0")
      lx <- log(x)
      m <- mean(lx)
      v <- var(lx)
      shape <- 1.2/sqrt(v)
      scale <- exp(m + 0.572/shape)
      start <- list(shape = shape, scale = scale)
      start <- start[!is.element(names(start), dots)]
    }
    if (distname == "gamma" && is.null(start)) {
      if (any(x < 0))
        stop("gamma values must be >= 0")
      m <- mean(x)
      v <- var(x)
      start <- list(shape = m^2/v, rate = m/v)
      start <- start[!is.element(names(start), dots)]
      control <- list(parscale = c(1, start$rate))
    }
    if (distname == "negative binomial" && is.null(start)) {
      m <- mean(x)
      v <- var(x)
      size <- if (v > m)
        m^2/(v - m)
      else 100
      start <- list(size = size, mu = m)
      start <- start[!is.element(names(start), dots)]
    }
    if (is.element(distname, c("cauchy", "logistic")) &&
        is.null(start)) {
      start <- list(location = median(x), scale = IQR(x)/2)
      start <- start[!is.element(names(start), dots)]
    }
    if (distname == "t" && is.null(start)) {
      start <- list(m = median(x), s = IQR(x)/2, df = 10)
      start <- start[!is.element(names(start), dots)]
    }
  }
  if (is.null(start) || !is.list(start))
    stop("")
  nm <- names(start)
  f <- formals(densfun)
  args <- names(f)
  m <- match(nm, args)
  if (any(is.na(m)))
    stop("'start' specifies names which are not arguments to 'densfun'")
  formals(densfun) <- c(f[c(1, m)], f[-c(1, m)])
  dens <- function(parm, x, ...) densfun(x, parm, ...)
  if ((l <- length(nm)) > 1L)
    body(dens) <- parse(text = paste("densfun(x,", paste("parm[",
                                                         1L:l, "]", collapse = ", "), ", ...)"))
  Call[[1L]] <- quote(stats::optim)
  Call$densfun <- Call$start <- NULL
  Call$x <- x
  Call$par <- start
  Call$fn <- if ("log" %in% args)
    mylogfn
  else myfn
  Call$hessian <- TRUE
  if (length(control))
    Call$control <- control
  if (is.null(Call$method)) {
    if (any(c("lower", "upper") %in% names(Call)))
      Call$method <- "L-BFGS-B"
    else if (length(start) > 1L)
      Call$method <- "BFGS"
    else Call$method <- "Nelder-Mead"
  }
  res <- eval.parent(Call)
  if (res$convergence > 0L)
    stop("optimization failed")
  vc <- solve(res$hessian)
  sds <- sqrt(diag(vc))
  structure(list(estimate = res$par, sd = sds, vcov = vc,
                 loglik = -res$value, n = n), class = "FitDistr")
}
# Funcion qqPlot ---------------------
qqPlot <- function(x, y, confbounds = TRUE, alpha, main, xlab, ylab, xlim, ylim, border = "red", bounds.col = "black", bounds.lty = 1, start, grapic = TRUE, axis.y.right = FALSE, bw.theme = FALSE,...){
  parList = list()
  if (is.null(parList[["col"]])){
    parList$col = 1:2
  }
  if (is.null(parList[["pch"]])){
    parList$pch = 19
  }
  if (is.null(parList[["lwd"]])){
    parList$lwd = 0.5
  }
  if (is.null(parList[["cex"]])){
    parList$cex = 1
  }
  if(inherits(x, "DistrCollection")){
    distList <- x$distr
    grap <- qqPlot(distList[[1]]$x, grapic = FALSE, ylab = "", xlab = "", main = paste(distList[[1]]$name,"distribution"))
    for (i in 2:length(distList)){
      aux <- qqPlot(distList[[i]]$x, grapic = FALSE, ylab = "", xlab = "", main = paste(distList[[i]]$name,"distribution"))
      grap$plot <-  grap$plot + aux$plot
    }
    show(grap$plot + plot_annotation(title = "QQ Plot for a Collection Distribution"))
    invisible()
  }
  else{
    if (missing(y))
      y = "normal"
    if(missing(alpha))
      alpha = 0.05
    if (alpha <=0 || alpha >=1)
      stop(paste("alpha should be between 0 and 1!"))
    if (missing(main))
      main = paste("QQ Plot for", deparse(substitute(y)), "distribution")
    if (missing(xlab))
      xlab = paste("Quantiles for", deparse(substitute(x)))
    if (missing(ylab))
      ylab = paste("Quantiles from", deparse(substitute(y)), "distribution")
    if (is.numeric(y)) {
      cat("\ncalling (original) qqplot from namespace stats!\n")
      return(stats::qqplot(x, y, ...))
    }
    qFun = NULL
    theoretical.quantiles = NULL
    xs = sort(x)
    distribution = tolower(y)
    distWhichNeedParameters = c("weibull", "logistic", "gamma","exponential", "f",
                                "geometric", "chi-squared", "negative binomial",
                                "poisson")

    threeParameterDistr = c("weibull3", "lognormal3", "gamma3")
    threeParameter = distribution %in% threeParameterDistr
    if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)
    if(is.character(distribution)){
      qFun = .charToDistFunc(distribution, type = "q")
      if (is.null(qFun))
        stop(paste(deparse(substitute(y)), "distribution could not be found!"))
    }
    # Puntos teoricos
    theoretical.probs = ppoints(xs)
    # Quantiles
    xq = NULL
    yq = quantile(xs, prob = c(0.25, 0.75))


    dots <- list(...)

    if(TRUE){
      fitList = .lfkp(parList, formals(qFun))
      fitList$x = xs
      fitList$densfun = distribution
      if(!missing(start))
        fitList$start = start
      if(!threeParameter){
        fittedDistr = do.call(FitDistr, fitList)
        parameter = fittedDistr$estimate

        #save the distribution parameter#
        thethas = fittedDistr$estimate
        # save the cariance-covariance matrix
        varmatrix = fittedDistr$vcov

      }else{
        parameter = do.call(paste(".",distribution, "3", sep = ""), list(xs) )    ####
        threshold = parameter$threshold
      }

      parameter = .lfkp(as.list(parameter), formals(qFun))
      params = .lfkp(parList, formals(qFun))
      parameter = .lfrm(as.list(parameter), params)
      parameter = c(parameter, params)
      theoretical.quantiles = do.call(qFun, c(list(c(theoretical.probs)), parameter))

      if(!threeParameter){
        confIntCapable = c("exponential", "log-normal", "logistic", "normal", "weibull", "gamma", "beta", "cauchy")
        getConfIntFun = .charToDistFunc(distribution, type = ".confint")
        if(confbounds == TRUE){
          if(distribution %in% confIntCapable){
            confInt = getConfIntFun(xs, thethas, varmatrix, alpha)
          }
        }
      }

      xq <- do.call(qFun, c(list(c(0.25, 0.75)), parameter))
    }
    else {
      params =.lfkp(parList, formals(qFun))
      params$p = theoretical.probs
      theoretical.quantiles = do.call(qFun, params)
      params$p = c(0.25, 0.75)
      xq = do.call(qFun, params)
    }

    params =.lfkp(parList, c(formals(plot.default), par()))

    if(!threeParameter){
      params$y = theoretical.quantiles
    }else{
      params$y = theoretical.quantiles+threshold
    }
    params$x = xs
    params$xlab = xlab
    params$ylab = ylab
    params$main = main
    if (!(is.null(params$col[1]) || is.na(params$col[1])))
      params$col = params$col[1]

    params$lwd = 1

    ############ Desde Aquí comienza la gráfica ############
    p <- ggplot(data = data.frame(x=params$x, y=params$y), mapping=aes(x=x, y=y)) +
      geom_point() + labs(x = xlab, y = ylab, title = main) + theme_minimal()

    params =.lfkp(parList, c(formals(abline), par()))
    params$a = 0
    params$b = 1
    params$col = border
    p <- p + geom_abline(intercept = params$a, slope = params$b, col = params$col, lwd = params$lwd)

    if(!threeParameter){
      if(confbounds == TRUE){
        if(distribution %in% confIntCapable){
          params =.lfkp(parList, c(formals(lines), par()))
          params$x = confInt[[3]]
          params$y = confInt[[1]]
          params$col = bounds.col
          params$lty = bounds.lty
          # La curva de abajo
          p <- p + geom_line(data = data.frame(x=params$x, y=params$y), aes(x = x, y = y),
                             col = params$col, lty = params$lty, lwd = params$lwd)
          params$x = confInt[[3]]
          params$y = confInt[[2]]
          params$col = bounds.col
          params$lty = bounds.lty
          # curva de arriba
          p <- p + geom_line(data = data.frame(x=params$x, y=params$y), aes(x = x, y = y),
                             col = params$col, lty = params$lty, lwd = params$lwd)
        }
      }
    }
    if(axis.y.right){
      p <- p + scale_y_continuous(position = "right")
    }
    if(bw.theme){
      p <- p + theme_bw()
    }
    if(main == ""){
      p <- p + labs(title = NULL)
    }
    if(grapic){
      show(p)
      invisible(list(x = theoretical.quantiles, y = xs, int = params$a, slope = params$b, plot = p))
    }
    else{
      invisible(list(x = theoretical.quantiles, y = xs, int = params$a, slope = params$b, plot = p))
    }
  }
}



# Funcion ppPlot ---------------------
ppPlot <- function (x, distribution, confbounds = TRUE, alpha, probs, main, xlab, ylab, xlim, ylim, border = "red", bounds.col = "black", bounds.lty = 1, start, grapic = TRUE, axis.y.right = FALSE, bw.theme = FALSE,...)
{
  conf.level = 0.95
  conf.lines = TRUE
  if (!(is.numeric(x) | inherits(x, "DistrCollection")))
    stop(paste(deparse(substitute(x)), " needs to be numeric or an object of class distrCollection"))
  parList = list(...)
  if (is.null(parList[["col"]]))
    parList$col = c("black", "red", "gray")
  if (is.null(parList[["pch"]]))
    parList$pch = 19
  if(is.null(parList[["lwd"]]))
    parList$lwd = 1
  if (is.null(parList[["cex"]]))
    parList$cex = 1
  qFun = NULL
  xq = NULL
  yq = NULL
  x1 = NULL
  if(missing(alpha))
    alpha = 0.05
  if (alpha <=0 || alpha >=1)
    stop(paste("alpha should be between 0 and 1!"))
  if (missing(probs))
    probs = ppoints(11)
  else if (min(probs) <= 0 || max(probs) >= 1)
    stop("probs should be values within (0,1)!")
  probs = round(probs, 2)
  if (is.numeric(x)) {
    x1 <- sort(na.omit(x))
    if (missing(xlim))
      xlim = c(min(x1) - 0.1 * diff(range(x1)), max(x1) +0.1 * diff(range(x1)))
  }
  if (missing(distribution))
    distribution = "normal"
  if (missing(ylim))
    ylim = NULL
  if (missing(main))
    main = paste("Probability Plot for", deparse(substitute(distribution)),
                 "distribution")
  if (missing(xlab))
    xlab = deparse(substitute(x))
  if (missing(ylab))
    ylab = "Probability"
  if(inherits(x, "DistrCollection")){
    distList <- x$distr
    grap <- ppPlot(distList[[1]]$x, grapic = FALSE, ylab = "", xlab = "", main = paste(distList[[1]]$name,"distribution"))
    for (i in 2:length(distList)){
      aux <- ppPlot(distList[[i]]$x, grapic = FALSE, ylab = "", xlab = "", main = paste(distList[[i]]$name,"distribution"))
      grap$plot <-  grap$plot + aux$plot
    }
    show(grap$plot + plot_annotation(title = "QQ Plot for a Collection Distribution"))
    invisible()
  }
  distWhichNeedParameters = c("weibull", "gamma", "logistic","exponential","f",
                              "geometric", "chi-squared", "negative binomial",
                              "poisson")
  # new
  threeParameterDistr = c("weibull3", "lognormal3", "gamma3")
  threeParameter = distribution %in% threeParameterDistr
  if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)
  # end new
  if (is.character(distribution)) {
    qFun = .charToDistFunc(distribution, type = "q")
    pFun = .charToDistFunc(distribution, type = "p")
    dFun = .charToDistFunc(distribution, type = "d")
    if (is.null(qFun))
      stop(paste(deparse(substitute(y)), "distribution could not be found!"))
  }
  dots <- list(...)
  if (TRUE) {
    fitList = .lfkp(parList, formals(qFun))
    fitList$x = x1
    fitList$densfun = distribution
    if (!missing(start))
      fitList$start = start
    if(!threeParameter){
      fittedDistr = do.call(FitDistr, fitList)
      parameter = fittedDistr$estimate
      #save the distribution parameter#
      thethas = fittedDistr$estimate
      # save the cariance-covariance matrix
      varmatrix = fittedDistr$vcov
    }else{
      parameter = do.call(paste(".",distribution, "3", sep = ""), list(x1) )    ####
      print(parameter[3])
      threshold = parameter$threshold
    }
    parameter = .lfkp(as.list(parameter), formals(qFun))
    params = .lfkp(parList, formals(qFun))
    parameter = .lfrm(as.list(parameter), params)
    print(parameter)
    parameter = c(parameter, params)
    # new
    if(!threeParameter){
      # array containing names of the distributions, for which conf intervals can be computed
      confIntCapable = c("exponential", "log-normal", "logistic", "normal", "weibull", "gamma", "beta", "cauchy")
      getConfIntFun = .charToDistFunc(distribution, type = ".confint")
      # if possible, compute the conf intervals
      if(confbounds == TRUE){
        if(distribution %in% confIntCapable){
          confInt = getConfIntFun(x1, thethas, varmatrix, alpha)
        }
      }# end of my code
    }
    y = do.call(qFun, c(list(ppoints(x1)), as.list(parameter)))
    yc = do.call(qFun, c(list(ppoints(x1)), as.list(parameter)))
    cv = do.call(dFun, c(list(yc), as.list(parameter)))
    print(cv)
    axisAtY = do.call(qFun, c(list(probs), as.list(parameter)))
    yq = do.call(qFun, c(list(c(0.25, 0.75)), as.list(parameter)))
    xq = quantile(x1, probs = c(0.25, 0.75))
  }
  else {
    params = .lfkp(parList, formals(qFun))
    params$p = ppoints(x1)
    y = do.call(qFun, params)
    params$p = probs
    axisAtY = do.call(qFun, params)
    params$p = c(0.25, 0.75)
    yq = do.call(qFun, params)
    xq = quantile(x1, probs = c(0.25, 0.75))
  }
  params = .lfkp(parList, c(formals(plot.default), par()))
  params$x = x1
  params$y = y
  params$xlab = xlab
  params$ylab = ylab
  params$main = main
  params$xlim = xlim
  params$axes = FALSE
  params$lwd = 1
  if (!(is.null(params$col[1]) || is.na(params$col[1])))
    params$col = params$col[1]
  # PLOT
  p <- ggplot(data.frame(x = x1, y = y), aes(x = x, y = y)) +
    geom_point(size = 1, color = "black") +
    labs(x = xlab, y = ylab, title = main) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(scale = 1/max(x1), suffix = "", accuracy = 0.01))+
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold"))
  # Regression line
  params = .lfkp(parList, c(formals(abline), par()))
  if(!threeParameter){
    params$a = 0
  }else{
    params$a = -threshold
  }
  params$b = 1
  params$col = border
  p <- p + geom_abline(intercept = params$a, slope = params$b, color = border)
  if(!threeParameter){
    if(confbounds == TRUE){
      if(distribution %in% confIntCapable){
        params =.lfkp(parList, c(formals(lines), par()))
        params$x = confInt[[3]]
        params$y = confInt[[1]]
        params$col = bounds.col
        params$lty = bounds.lty
        # La curva de abajo
        p <- p + geom_line(data = data.frame(x=params$x, y=params$y), aes(x = x, y = y),
                           col = params$col, lty = params$lty, lwd = params$lwd)
        params$x = confInt[[3]]
        params$y = confInt[[2]]
        params$col = bounds.col
        params$lty = bounds.lty
        # curva de arriba
        p <- p + geom_line(data = data.frame(x=params$x, y=params$y), aes(x = x, y = y),
                           col = params$col, lty = params$lty, lwd = params$lwd)
      }
    }
  }
  if(axis.y.right){
    p <- p + scale_y_continuous(position = "right")
  }
  if(bw.theme){
    p <- p + theme_bw()
  }
  if(main == ""){
    p <- p + labs(title = NULL)
  }
  if(grapic){
    show(p)
    invisible(list(x = x, y = y, int = params$a, slope = params$b, plot = p))
  }
  else{
    invisible(list(x = x, y = y, int = params$a, slope = params$b, plot = p))
  }
}

set.seed(1234)
x <- rnorm(20, mean = 20)
distribution <- "normal"
ppPlot(x, distribution, bounds.lty = 3, bounds.col = "red")


