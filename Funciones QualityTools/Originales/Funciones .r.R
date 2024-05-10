# https://rdrr.io/cran/qualityTools/f/
# son los archivos .r

# ORIGINALES

# -------------- sta_n.r -------------- 
.nblock = function(fdo) {
  if (class(fdo) != "facDesign") 
    stop(paste(deparse(substitute(data)), "needs to be an object of class 'facDesign'"))
  return(length(unique(fdo@block[[1]])))
}
starDesign = function(k, p = 0, alpha = c("both", "rotatable", "orthogonal"), cs, cc, data) {
  DB = FALSE
  fdo = NULL
  csFrame = NULL
  ccFrame = NULL
  starFrame = NULL
  blocks = 1
  alpha = alpha[1]
  if (DB) 
    print(alpha)
  if (missing(cc)) 
    cc = 1
  if (missing(cs)) 
    cs = 1
  if (!missing(k)) {
    nameVec = LETTERS[1:k]
  }
  if (!missing(data)) {
    fdo = data
    k = ncol(cube(fdo))
    if (class(fdo) != "facDesign") {
      stop(paste(deparse(substitute(data)), "needs to be an object of class 'facDesign'"))
    }
    if (nrow(star(fdo)) > 0) 
      stop(paste("star portion of", deparse(substitute(data)), "not empty"))
    k = length(names(fdo))
    nameVec = names(names(fdo))
    cc = nrow(centerCube(fdo))
    p = ncol(cube(fdo)) - log(nrow(unique(cube(fdo))), 2)
    blocks = .nblock(fdo) + 1
  }
  if (is.numeric(alpha)) 
    a = alpha
  if (alpha == "rotatable") 
    a = .alphaRot(k, p)
  if (alpha == "orthogonal") 
    a = .alphaOrth(k, p, cc = cc, cs = cs)
  if (alpha == "both") {
    found = FALSE
    for (i in seq(along = .rsmOrth)) {
      if (DB) {
        print(k)
        print(blocks)
        print(p)
      }
      if (.rsmOrth[[i]]$k == k) 
        if (.rsmOrth[[i]]$blocks == blocks) 
          if (.rsmOrth[[i]]$p == p) {
            found = TRUE
            cc = .rsmOrth[[i]]$cc
            cs = .rsmOrth[[i]]$cs
            p = .rsmOrth[[i]]$p
            a = .alphaOrth(k, p, cc, cs)
            break
          }
    }
    if (!found) {
      return("no starDesign with approximate rotatability and orthogonality available")
    }
  }
  starFrame = .starFrame(k, alpha = a)
  names(starFrame) = nameVec
  if (DB) 
    print(starFrame)
  if (!missing(data)) 
    star(fdo) = starFrame
  if (DB) 
    print("starFrame added")
  if (cs > 0) {
    csFrame = as.data.frame(matrix(0, nrow = cs, ncol = k))
    names(csFrame) = nameVec
    if (!missing(data)) {
      centerStar(fdo) = csFrame
      if (DB) 
        print("csFrame added")
    }
  }
  if (cc > 0) {
    ccFrame = as.data.frame(matrix(0, nrow = cc, ncol = k))
    names(ccFrame) = nameVec
    if (!missing(data)) {
      centerCube(fdo) = ccFrame
      if (DB) 
        print("ccFrame added")
    }
  }
  if (!missing(data)) 
    return(fdo)
  else return(list(star = starFrame, centerStar = csFrame, centerCube = ccFrame))
} 

# -------------- Com_t.r -------------- 
setGeneric("compPlot", function(x, main, xlab, ylab, col, cex.lab, fun = NULL, ...) standardGeneric("compPlot"))
setMethod("compPlot", signature(x = "gageRR"), function(x, main, xlab, ylab, col, cex.lab, fun = NULL, ...) {
  if (missing(xlab)) 
    xlab = ""
  if (missing(ylab)) 
    ylab = ""
  if (missing(col)) 
    col = 1
  old.par <- par(no.readonly = TRUE)
  ops = length(unique(x[, 3]))
  pts = length(unique(x[, 4]))
  n = length(x[, 5])/pts/ops
  comp = unique(x[, 4])
  val = as.data.frame(x)
  dat = val
  if (identical(fun, NULL) == FALSE) 
    means = matrix(data = NA, nrow = ops, ncol = pts)
  if (identical(fun, NULL)) 
    means = matrix(data = NA, nrow = ops, ncol = n * pts)
  if (missing(main)) 
    main = paste("Comparison Plot for", names(val)[3])
  if (missing(cex.lab)) 
    cex.lab = 10/ops
  par(mfrow = c(ops, ops))
  par(mar = c(0, 0, 0, 0))
  par(oma = c(8, 8, 2, 0))
  for (i in 1:ops) {
    dat = subset(val, val[, 3] == unique(x[, 3])[i])
    if (identical(fun, NULL) == FALSE) {
      for (j in 1:pts) means[i, j] = fun(subset(dat, dat[, 4] == unique(dat[, 4])[j])[, 5])
    }
    else means[i, ] = dat[, 5]
    for (k in 1:ops) {
      if (k == i) {
        plot(1:10, col = "transparent", axes = FALSE, xlab = "", ylab = "")
        text(5.5, 5, unique(x[, 3])[i], cex = cex.lab, font = 2, xpd = TRUE)
      }
      else {
        if (k < i) {
          plot(means[k, ], means[i, ], xlab = xlab, ylab = ylab, col = col, axes = FALSE, ...)
          box()
          if (k == 1) 
            axis(2)
          if (i == ops) 
            axis(1)
        }
        else plot(1:10, col = "transparent", axes = FALSE, xlab = "", ylab = "")
      }
    }
  }
  title(main, outer = TRUE)
  par(old.par)
  return()
}) 


# -------------- ave_t.r -------------- 
.splitDev = function(x) {
  if (x > 6) 
    dev = TRUE
  else dev = FALSE
  if (x == 1) 
    mfrow = c(1, 1)
  if (x == 2) 
    mfrow = c(1, 2)
  if (x == 3) 
    mfrow = c(2, 2)
  if (x == 4) 
    mfrow = c(2, 2)
  if (x == 5) 
    mfrow = c(2, 3)
  if (x == 6) 
    mfrow = c(2, 3) 
  if (x >= 7) 
    mfrow = c(3, 3)     
  return(list(dev, mfrow))
}
setGeneric("averagePlot", function(x, main, xlab, ylab, col, ask = TRUE, single = FALSE, ...) standardGeneric("averagePlot"))
setMethod("averagePlot", signature(x = "gageRR"), function(x, main, xlab, ylab, col, ask = TRUE, single = FALSE, ...) {
  old.par <- par(no.readonly = TRUE)
  ops = length(unique(x[, 3]))
  pts = length(unique(x[, 4]))
  n = length(x[, 5])/pts/ops
  ref = numeric(pts * n)
  val = as.data.frame(x)
  dat = val
  if (missing(xlab)) 
    xlab = "reference"
  if (missing(ylab)) 
    ylab = "values"
  if (missing(col)) 
    col = 1
  if (missing(main)) 
    main = as.vector(paste(names(val)[3], unique(x[, 3])))
  scr = TRUE
  if (identical(par("mfrow"), as.integer(c(1, 1))) == TRUE) 
    scr = FALSE
  if (single == TRUE) 
    par(mfcol = c(1, 1))
  k = 0
  m = 1
  for (j in 1:pts) {
    dat = subset(val, val[, 4] == unique(x[, 4])[j])
    ref[(k + 1):(k + n)] = mean(dat[, 5])
    k = k + n
  }
  if (single == FALSE && scr == FALSE) 
    par(mfrow = .splitDev(ops)[[2]])
  for (i in 1:ops) {
    dat = subset(val, val[, 3] == unique(x[, 3])[i])
    dat = dat[order(dat[, 4]), ]
    plot(ref, y = dat[, 5], main = main[i], xlab = xlab, ylab = ylab, col = col, ...)
    if ((i/6 - i%/%6) == 0 && .splitDev(ops)[[1]] == TRUE && scr == FALSE && single == FALSE) {
      if (ask != TRUE) 
        dev.new()
      else par(ask = TRUE)
      par(mfrow = .splitDev(ops - m * 6)[[2]])
      m = m + 1
    }
    if (scr == TRUE && (i/(prod(par("mfcol"))) - i%/%(prod(par("mfcol")))) == 0 && single == FALSE) {
      if (ask != TRUE) {
        dev.new()
        par(old.par)
      }
      else par(ask = TRUE)
    }
    if (single == TRUE && i != ops) {
      if (ask == TRUE) 
        par(ask = TRUE)
      else dev.new()
    }
  }
  return()
  par(old.par)
}) 


# -------------- myA_r.r -------------- 
.myADTest = function(x, distribution, ...) {                                                                     ####   .MYADTESTS-FUNCTION
  #require(MASS, quietly = TRUE)
  if (missing(distribution)) 
    distribution = "normal"
  data.name = names(x)
  if (is.data.frame(x)) 
    x = x[, 1]
  dots = list(...)
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
  }
  #    if (identical(distribution, "log-normal")) {                               ####
  #        x = log(x)                                                             ####
  #        distribution = "normal"                                                ####
  #    }                                                                          ####
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
  class(out) = "adtest"
  return(out)
}
print.adtest = function(x, digits = 4, quote = TRUE, prefix = "", ...) {
  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data: ", x$data.name, "\n")
  out <- character()
  if (!is.null(x$statistic)) 
    out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic, 4))))
  if (!is.null(x$parameter)) 
    out <- c(out, paste(names(x$parameter), "=", format(round(x$parameter, 3))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = digits)
    if (x$tableValue) {
      if (x$smaller) 
        out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("<=", fp)))
      else out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "=") fp else paste(">", fp)))
    }
    else {
      out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
    }
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("alternative hypothesis: ")
  if (!is.null(x$null.value)) {
    if (length(x$null.value) == 1) {
      cat("true", names(x$null.value), "is not equal to", x$null.value, "\n")
    }
    else {
      cat(x$alternative, "\nnull values:\n")
      print(x$null.value, ...)
    }
  }
  if (!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")), "percent confidence interval:\n", format(c(x$conf.int[1L], x$conf.int[2L])), "\n")
  }
  if (!is.null(x$estimate)) {
    cat("sample estimates:\n")
    print(x$estimate, ...)
  }
  cat("\n")
  invisible(x)
}

#  -------------- pcr_r -------------- 
.c4 = function(n) {
  if (n > 1 && n < 342) 
    sqrt(2/(n - 1)) * (factorial(n/2 - 1)/factorial((n - 1)/2 - 1))
  else stop("n needs to be bigger than 1 and smaller than 342")
}
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
.sdSg(1:10, grouping = c(1, 1, 1, 1, 1, 5, 5, 5, 5, 5))
pcr = function(x, distribution = "normal", lsl, usl, target, boxcox = FALSE, lambda = c(-5,                      ####   PCR-FUNCTION
                                                                                        5), main, xlim, ylim, grouping = NULL, std.dev = NULL, conf.level = 0.9973002, start, lineWidth = 1, 
               lineCol = "red", lineType = "solid", specCol = "red3", specWidth = 1, cex.text = 2, cex.val = 1.5, 
               cex.col = "darkgray", plot = TRUE, bounds.lty = 3, bounds.col = "red", ...) {                                ####
  DB = FALSE                                                                   
  data.name = deparse(substitute(x))[1]
  #require(MASS, quietly = TRUE)
  if(plot==TRUE)                                                              ####
  {                                                                           ####
    par.orig <- par(c("mar", "oma", "mfrow"))                                  ####
    on.exit(par(par.orig))                                                     ####
  }                                                                           ####
  parList = list(...)
  if (is.null(parList[["col"]])) 
    parList$col = "lightblue"
  if (is.null(parList[["border"]])) 
    parList$border = 1
  if (is.null(parList[["lwd"]])) 
    parList$lwd = 1
  if (is.null(parList[["cex.axis"]])) 
    parList$cex.axis = 1.5
  if (missing(lsl)) 
    lsl = NULL
  if (missing(usl)) 
    usl = NULL
  if (missing(target)) 
    target = NULL
  if (missing(lambda)) 
    lambda = c(-5, 5)
  if (!is.numeric(lambda)) 
    stop("lambda needs to be numeric")
  if (any(x<0) && any(distribution==c("beta","chi-squared","exponential","f","geometric","lognormal","log-normal","negative binomial","poisson","weibull","gamma")))
    stop("choosen distribution needs all values in x to be > 0!")
  if (any(x>1) && distribution=="beta")
    stop("choosen distribution needs all values in x to be between 0 and 1!")
  paramsList = vector(mode = "list", length = 0)
  estimates = vector(mode = "list", length = 0)
  varName = deparse(substitute(x))
  dFun = NULL
  pFun = NULL
  qFun = NULL
  cp = NULL
  cpu = NULL
  cpl = NULL
  cpk = NULL
  ppt = NULL
  ppl = NULL
  ppu = NULL
  xVec = numeric(0)
  yVec = numeric(0)
  if (is.vector(x)) 
    x = as.data.frame(x)
  #    if (identical(distribution, "log-normal")) {                               ####
  #        x = log(x)                                                             ####
  #        if(is.null(usl)==FALSE)                                                ####
  #         usl=log(usl)                                                          ####
  #        if(is.null(lsl)==FALSE)                                                ####
  #        lsl=log(lsl)                                                           ####
  #        distribution = "normal"                                                ####
  #        data.name = paste("log(", data.name, ")", sep = "")                    ####
  #    }                                                                          ####
  any3distr=FALSE;not3distr=FALSE                                            ####
  if(distribution=="weibull3" || distribution=="lognormal3" || distribution=="gamma3")####
    any3distr=TRUE                                                             ####
  if (distribution!="weibull3" && distribution!="lognormal3" && distribution!="gamma3")####
    not3distr=TRUE                                                                  ####
  if (boxcox) {
    distribution = "normal"
    if (length(lambda) >= 2) {
      temp = boxcox(x[, 1] ~ 1, lambda = seq(min(lambda), max(lambda), 1/10), plotit = FALSE)
      i = order(temp$y, decreasing = TRUE)[1]
      lambda = temp$x[i]
    }
    x = as.data.frame(x[, 1]^lambda)
  }
  numObs = nrow(x)
  if (!is.null(grouping)) 
    if (is.vector(grouping)) 
      grouping = as.data.frame(grouping)
  center = colMeans(x)
  if (!is.null(x) & !is.null(grouping)) {
    if (nrow(x) != nrow(grouping)) 
      stop(paste("length of ", deparse(substitute(grouping)), " differs from length of ", varName))
  }
  if (missing(main)) 
    if (boxcox) 
      main = paste("Process Capability using box cox transformation for", varName)
  else main = paste("Process Capability using", as.character(distribution), "distribution for", 
                    varName)
  if (is.null(std.dev)) {
    if (is.null(grouping)) 
      std.dev = .sdSg(x)
    else std.dev = .sdSg(x, grouping)
  }
  if (conf.level < 0 | conf.level > 1) 
    stop("conf.level must be a value between 0 and 1")
  confHigh = conf.level + (1 - conf.level)/2
  confLow = 1 - conf.level - (1 - conf.level)/2
  if (DB) {
    print(paste("confHigh:", confHigh))
    print(paste("confLow:", confLow))
  }
  distWhichNeedParameters = c("weibull", "logistic", "gamma", "exponential", "f", "geometric", 
                              "chi-squared", "negative binomial", "poisson")
  if (is.character(distribution)) {
    dis=distribution                                                           ####
    if (identical(distribution,"weibull3"))                                     ####
      dis="weibull3"                                                             ####
    if (identical(distribution,"gamma3"))                                       ####
      dis="gamma3"                                                               ####
    if (identical(distribution,"lognormal3"))                                   ####
      dis="lognormal3"                                                           ####
    qFun = .charToDistFunc(dis, type = "q")                                 ####
    pFun = .charToDistFunc(dis, type = "p")                                 ####
    dFun = .charToDistFunc(dis, type = "d")                                 ####
    if (is.null(qFun) & is.null(pFun) & is.null(dFun)) 
      stop(paste(deparse(substitute(y)), "distribution could not be found!"))
  }
  if (TRUE) {                                                                 #### distribution!="weibull3" && distribution!="lognormal3" && distribution!="gamma3"
    if (DB) 
      print("TODO: Pass the estimated parameters correctly")
    fitList = vector(mode = "list", length = 0)
    fitList$x = x[, 1]
    fitList$densfun = dis                                                   ####
    if (!missing(start)) 
      fitList$start = start
    if (not3distr)                                                          ####
    {                                                                       ####
      fittedDistr = do.call(MASS::fitdistr, fitList)
      estimates = as.list(fittedDistr$estimate)
      paramsList = estimates
    }                                                                       ####
    if (distribution=="weibull3")                                           ####
    {                                                                       ####
      paramsList= .weibull3(x[,1])                                           ####
      estimates = paramsList                                                 ####
    }                                                                       ####
    if (distribution=="lognormal3")                                         ####
    {                                                                       ####
      paramsList= .lognormal3(x[,1])                                         ####
      estimates = paramsList                                                 ####
    }                                                                       ####
    if (distribution=="gamma3")                                             ####
    {                                                                       ####
      paramsList= .gamma3(x[,1])                                             ####
      estimates = paramsList                                                 ####
    }                                                                       ####
    if (DB) 
      print(paste("parameter: ", paramsList))
  }
  paramsList = c(paramsList, .lfkp(parList, formals(qFun)))
  if (distribution == "normal") {
    paramsList$mean = center
    paramsList$sd = std.dev
    estimates = paramsList
  }
  if (boxcox) {
    if (!is.null(lsl)) 
      lsl = lsl^lambda
    if (!is.null(usl)) 
      usl = usl^lambda
    if (!is.null(target)) 
      target = target^lambda
  }
  if (is.null(lsl) && is.null(usl)) {
    paramsList$p = confLow
    lsl = do.call(qFun, paramsList)
    paramsList$p = confHigh
    usl = do.call(qFun, paramsList)
  }
  if (identical(lsl, usl)) 
    stop("lsl == usl")
  if (!is.null(lsl) && !is.null(target) && target < lsl) 
    stop("target is less than lower specification limit")
  if (!is.null(usl) && !is.null(target) && target > usl) 
    stop("target is greater than upper specification limit")
  if (!is.null(lsl) && !is.null(usl)) 
    if (lsl > usl) {
      temp = lsl
      lsl = usl
      usl = temp                       
    }
  paramsList$p = c(confLow, 0.5, confHigh)
  paramsListTemp = .lfkp(paramsList, formals(qFun))                           ####
  qs = do.call(qFun, paramsListTemp)                                         ####  
  paramsListTemp = .lfkp(paramsList, formals(pFun))                           ####
  if (!is.null(lsl) && !is.null(usl)) 
    cp = (usl - lsl)/(qs[3] - qs[1])
  if (!is.null(usl)) {
    cpu = (usl - qs[2])/(qs[3] - qs[2])
    paramsListTemp$q = usl                                                  ####
    ppu = 1 - do.call(pFun, paramsListTemp)                                 ####
  }
  if (!is.null(lsl)) {
    cpl = (qs[2] - lsl)/(qs[2] - qs[1])
    paramsListTemp$q = lsl                                                  ####
    ppl = do.call(pFun, paramsListTemp)                                     ####
  }
  cpk = min(cpu, cpl)
  ppt = sum(ppl, ppu)
  if (DB == TRUE) {
    print(cp)
    print(cpk)
    print(cpu)
    print(cpl)
    print(ppu)
    print(ppl)
    print(ppt)
  } 
  if(plot==TRUE)
  {
    if (missing(xlim)) {
      xlim <- range(x[, 1], usl, lsl)
      xlim <- xlim + diff(xlim) * c(-0.2, 0.2)
    }
    xVec <- seq(min(xlim), max(xlim), length = 200)
    dParamsList = .lfkp(paramsList, formals(dFun))
    dParamsList$x = xVec
    yVec = do.call(dFun , dParamsList)
    histObj <- hist(x[, 1], plot = FALSE)
    if (missing(ylim)) {
      ylim <- range(histObj$density, yVec)
      ylim <- ylim + diff(ylim) * c(0, 0.05)
    }
    par(mar = c(0, 0, 0, 0) + 0.1)
    par(oma = c(2, 4, 7, 4) + 0.1)
    layout(matrix(c(1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 4, 5, 5, 6, 7), nrow = 4, byrow = TRUE))
    do.call(hist, c(list(x[, 1], freq = FALSE, xlim = xlim, ylim = ylim, main = ""), parList))
    abline(h = 0, col = "gray")
    tempList = parList
    tempList$col = 1
    tempList$border = NULL
    do.call(box, tempList)
    lines(xVec, yVec, lwd = lineWidth, col = lineCol, lty = lineType)
    abline(v = usl, col = specCol, lwd = specWidth, lty = 5)
    abline(v = lsl, col = specCol, lwd = specWidth, lty = 5)
    abline(v = target, col = specCol, lwd = specWidth, lty = 5)
    if (!is.null(lsl)) 
      axis(side = 3, at = lsl, labels = paste("LSL =", format(lsl, digits = 3)), col = specCol)
    if (!is.null(usl)) 
      axis(side = 3, at = usl, labels = paste("USL =", format(usl, digits = 3)), col = specCol)
    if (!is.null(lsl) && !is.null(usl)) 
      axis(side = 3, at = c(lsl, usl), labels = c(paste("LSL =", format(lsl, digits = 3)), paste("USL =", 
                                                                                                 format(usl, digits = 3))), col = specCol)
    if (!is.null(target)) 
      text(target, max(ylim), "TARGET", pos = 1, col = cex.col, cex = cex.text)
    title(main = main, outer = TRUE)
    plot(0:5, 0:5, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
    box()
    text(2.3, 1, expression(c[p]), pos = 2, cex = cex.val)
    if (is.null(cp)) 
      text(2, 1, paste("=", "*"), pos = 4, cex = cex.val)
    else text(2, 1, paste("=", round(cp, 2)), pos = 4, cex = cex.val)
    text(2.3, 2, expression(c[pk]), pos = 2, cex = cex.val)
    if (is.null(cpk)) 
      text(2, 2, paste("=", "*"), pos = 4, cex = cex.val)
    else text(2, 2, paste("=", round(cpk, 2)), pos = 4, cex = cex.val)
    text(2.3, 3, expression(c[pkL]), pos = 2, cex = cex.val)
    if (is.null(cpl)) 
      text(2, 3, paste("=", "*"), pos = 4, cex = cex.val)
    else text(2, 3, paste("=", round(cpl, 2)), pos = 4, cex = cex.val)
    text(2.3, 4, expression(c[pkU]), pos = 2, cex = cex.val)
    if (is.null(cpu)) 
      text(2, 4, paste("=", "*"), pos = 4, cex = cex.val)
    else text(2, 4, paste("=", round(cpu, 2)), pos = 4, cex = cex.val)
    index = 1:(length(estimates) + 3)
    plot(0:5, c(0:4, max(index) + 1), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
    box()
    names(x) = data.name
    if (not3distr)                                                              ####
    {                                                                           ####
      names(x) = data.name
      adTestStats = .myADTest(x, distribution)
      A = numeric()
      p = numeric()
      if (class(adTestStats) == "adtest") {
        A = adTestStats$statistic
        p = adTestStats$p.value
      }
      else {
        A = NA
        p = NA
      }
      text(2.3, rev(index)[2], "A", pos = 2, cex = cex.val)
      text(2, rev(index)[2], paste("=", format(A, digits = 3)), pos = 4, cex = cex.val)
      text(2.3, rev(index)[3], "p", pos = 2, cex = cex.val)
      if (!is.null(adTestStats$smaller) && adTestStats$smaller) 
        text(2, rev(index)[3], paste("<", format(p, digits = 3)), pos = 4, cex = cex.val)
      if (!is.null(adTestStats$smaller) && !adTestStats$smaller) 
        text(2, rev(index)[3], paste(">=", format(p, digits = 3)), pos = 4, cex = cex.val)
      if (is.null(adTestStats$smaller)) 
        text(2, rev(index)[3], paste("=", format(p, digits = 3)), pos = 4, cex = cex.val)
      text(2.3, rev(index)[1], "n", pos = 2, cex = cex.val)
      text(2, rev(index)[1], paste("=", numObs), pos = 4, cex = cex.val)
      j = 1
      for (i in 3:(3 + length(estimates) - 1)) {
        try(text(2.3, rev(index)[i + 1], names(estimates)[[j]], pos = 2, cex = cex.val), silent = TRUE)
        try(text(2, rev(index)[i + 1], paste("=", format(estimates[[j]], digits = 3)), pos = 4, cex = cex.val), 
            silent = TRUE)
        j = j + 1
      }
    }                                                                           ####
    
    if (any3distr)                                                              ####
    {                                                                           ####
      text(2.8, rev(index)[1], "n", pos = 2, cex = cex.val)                      ####
      text(2.5, rev(index)[1], paste("=", numObs), pos = 4, cex = cex.val)       ####
      text(2.8, rev(index)[2], "A", pos = 2, cex = cex.val)                      ####
      text(2.5, rev(index)[2], paste("=", "*"), pos = 4, cex = cex.val)          ####
      text(2.8, rev(index)[3], "p", pos = 2, cex = cex.val)                      ####
      text(2.5, rev(index)[3], paste("=", "*"), pos = 4, cex = cex.val)          ####
      j = 1                                                                   ####
      for (i in 3:(3 + length(estimates) - 1)) {                                  ####
        try(text(2.8, rev(index)[i + 1], names(estimates)[[j]], pos = 2, cex = cex.val), silent = TRUE)####
        try(text(2.5, rev(index)[i + 1], paste("=", format(estimates[[j]], digits = 3)), pos = 4, cex = cex.val),#### 
            silent = TRUE)                                                      ####
        j = j + 1                                                               ####
      }                                                                           ####
    }                                                                           ####
    
    
    qqPlot(x[, 1], y = distribution, ylab = "", main = "", axes = FALSE, 
           bounds.lty = bounds.lty, bounds.col = bounds.col)     
    axis(1)
    axis(4)
    box()
    par(mar = c(0, 0, 3, 2))
    plot(c(-1, 1), c(0.5, 5), type = "n", axes = FALSE)
    box()
    text(0, 4.5, "Expected Fraction Nonconforming", cex = cex.val)
    text(-1.05, 3, expression(p[t]), pos = 4, cex = cex.val)
    text(-1.05, 2, expression(p[L]), pos = 4, cex = cex.val)
    text(-1.05, 1, expression(p[U]), pos = 4, cex = cex.val)
    text(-0.9, 3, paste("=", format(ppt, digits = 6)), pos = 4, cex = cex.val)
    if (is.null(ppl)) 
      text(-0.9, 2, paste("= 0"), pos = 4, cex = cex.val)
    else text(-0.9, 2, paste("=", format(ppl, digits = 6)), pos = 4, cex = cex.val)
    if (is.null(ppu)) 
      text(-0.9, 1, paste("= 0"), pos = 4, cex = cex.val)
    else text(-0.9, 1, paste("=", format(ppu, digits = 6)), pos = 4, cex = cex.val)
    text(0.05, 3, expression(ppm), pos = 4, cex = cex.val)
    text(0.05, 2, expression(ppm), pos = 4, cex = cex.val)
    text(0.05, 1, expression(ppm), pos = 4, cex = cex.val)
    text(0.35, 3, paste("=", format(ppt * 1e+06, digits = 6)), pos = 4, cex = cex.val)
    if (is.null(ppl)) 
      text(0.35, 2, paste("= 0"), pos = 4, cex = cex.val)
    else text(0.35, 2, paste("=", format(ppl * 1e+06, digits = 6)), pos = 4, cex = cex.val)
    if (is.null(ppu)) 
      text(0.35, 1, paste("= 0"), pos = 4, cex = cex.val)
    else text(0.35, 1, paste("=", format(ppu * 1e+06, digits = 6)), pos = 4, cex = cex.val)
    par(mar = c(0, 0, 3, 0))
    plot(c(-1, 1), c(0.5, 5), type = "n", axes = FALSE)
    box()
    text(0, 4.5, "Observed", cex = cex.val)
    obsL = 0
    obsU = 0
    if (!is.null(lsl)) {
      obsL = (sum(x < lsl)/length(x)) * 1e+06
      text(-1, 2, paste("ppm =", obsL), pos = 4, cex = cex.val)
    }
    else text(-1, 2, paste("ppm =", 0), pos = 4, cex = cex.val)
    if (!is.null(usl)) {
      obsU = (sum(x > usl)/length(x)) * 1e+06
      text(-1, 1, paste("ppm =", obsU), pos = 4, cex = cex.val)
    }
    else text(-1, 1, paste("ppm =", 0), pos = 4, cex = cex.val)
    text(-1, 3, paste("ppm =", obsL + obsU), pos = 4, cex = cex.val)   
    if(not3distr)                                                               ####
    {                                                                           ####  
      print(adTestStats)                                                         ####  
      invisible(list(lambda = lambda, cp = cp, cpk = cpk, cpl = cpl, cpu = cpu,   ####
                     ppt = ppt, ppl = ppl, ppu = ppu, A = A, usl = usl,           ####
                     lsl = lsl, target = target))                                 ####
    }                                                                           ####   
    else                                                                        ####
      invisible(list(lambda = lambda, cp = cp, cpk = cpk, cpl = cpl, cpu = cpu,   ####
                     ppt = ppt, ppl = ppl, ppu = ppu, usl = usl,                  ####
                     lsl = lsl, target = target))                                 ####
  }                                                                               ####
  invisible(list(lambda = lambda, cp = cp, cpk = cpk, cpl = cpl, cpu = cpu,        ####
                 ppt = ppt, ppl = ppl, ppu = ppu, usl = usl,                  ####
                 lsl = lsl, target = target))                                 ####
}
cp = pcr
.pcr = function(x, distribution = "normal", lsl, usl, target, boxcox = FALSE, lambda = c(-5,                     ####   .PCR-FUNCTION
                                                                                         5), main, xlim, ylim, grouping = NULL, std.dev = NULL, conf.level = 0.9973002, start, lineWidth = 1, 
                lineCol = "red", lineType = "solid", specCol = "red3", specWidth = 1, cex.text = 2, cex.val = 1.5, 
                cex.col = "darkgray", ...) {                                                                                                  
  data.name = deparse(substitute(x))
  #require(MASS, quietly = TRUE)
  par.orig <- par(c("mar", "oma", "mfrow"))                                   ####
  on.exit(par(par.orig))                                                      ####                                                                         
  parList = list(...)
  if (is.null(parList[["col"]])) 
    parList$col = "lightblue"
  if (is.null(parList[["border"]])) 
    parList$border = 1
  if (is.null(parList[["lwd"]])) 
    parList$lwd = 1
  if (is.null(parList[["cex.axis"]])) 
    parList$cex.axis = 1.5
  if (missing(lsl)) 
    lsl = NULL
  if (missing(usl)) 
    usl = NULL
  if (missing(target)) 
    target = NULL
  if (missing(lambda)) 
    lambda = c(-5, 5)
  if (!is.numeric(lambda)) 
    stop("lambda needs to be numeric")
  paramsList = vector(mode = "list", length = 0)
  estimates = vector(mode = "list", length = 0)
  varName = deparse(substitute(x))
  dFun = NULL
  pFun = NULL
  qFun = NULL
  cp = NULL
  cpu = NULL
  cpl = NULL
  cpk = NULL
  ppt = NULL
  ppl = NULL
  ppu = NULL
  xVec = numeric(0)
  yVec = numeric(0)
  if (is.vector(x)) 
    x = as.data.frame(x)
  #    if (identical(distribution, "log-normal")) {                               ####
  #        x = log(x)                                                             ####
  #        if(is.null(usl)==FALSE)                                                ####
  #         usl=log(usl)                                                          ####
  #        if(is.null(lsl)==FALSE)                                                ####
  #        lsl=log(lsl)                                                           ####
  #        distribution = "normal"                                                ####
  #        data.name = paste("log(", data.name, ")", sep = "")                    ####
  #    }                                                                          ####
  any3distr=FALSE;not3distr=FALSE                                            ####
  if(distribution=="weibull3" || distribution=="lognormal3" || distribution=="gamma3")####
    any3distr=TRUE                                                             ####
  if (distribution!="weibull3" && distribution!="lognormal3" && distribution!="gamma3")####
    not3distr=TRUE                                                                  ####
  if (boxcox) {
    distribution = "normal"
    if (length(lambda) >= 2) {
      temp = boxcox(x[, 1] ~ 1, lambda = seq(min(lambda), max(lambda), 1/10), plotit = FALSE)
      i = order(temp$y, decreasing = TRUE)[1]
      lambda = temp$x[i]
    }
    x = as.data.frame(x[, 1]^lambda)
  }
  numObs = nrow(x)
  if (!is.null(grouping)) 
    if (is.vector(grouping)) 
      grouping = as.data.frame(grouping)
  center = colMeans(x)
  if (!is.null(x) & !is.null(grouping)) {
    if (nrow(x) != nrow(grouping)) 
      stop(paste("length of ", deparse(substitute(grouping)), " differs from length of ", varName))
  }
  if (missing(main)) 
    if (boxcox) 
      main = paste("Process Capability using box cox transformation for", varName)
  else main = paste("Process Capability using", as.character(distribution), "distribution for", 
                    varName)
  if (is.null(std.dev)) {
    if (is.null(grouping)) 
      std.dev = .sdSg(x)
    else std.dev = .sdSg(x, grouping)
  }
  if (conf.level < 0 | conf.level > 1) 
    stop("conf.level must be a value between 0 and 1")
  confHigh = conf.level + (1 - conf.level)/2
  confLow = 1 - conf.level - (1 - conf.level)/2
  distWhichNeedParameters = c("weibull", "logistic", "gamma", "exponential", "f", "geometric", 
                              "chi-squared", "negative binomial", "poisson")
  if (is.character(distribution)) {
    dis=distribution                                                           ####
    if (identical(distribution,"weibull3"))                                     ####
      dis="weibull3"                                                             ####
    if (identical(distribution,"gamma3"))                                       ####
      dis="gamma3"                                                               ####
    if (identical(distribution,"lognormal3"))                                   ####
      dis="lognormal3"                                                           ####
    qFun = .charToDistFunc(dis, type = "q")                                 ####
    pFun = .charToDistFunc(dis, type = "p")                                 ####
    dFun = .charToDistFunc(dis, type = "d")                                 ####
    if (is.null(qFun) & is.null(pFun) & is.null(dFun)) 
      stop(paste(deparse(substitute(y)), "distribution could not be found!"))
  }
  if (TRUE) {                                                                 #### distribution!="weibull3" && distribution!="lognormal3" && distribution!="gamma3"
    fitList = vector(mode = "list", length = 0)
    fitList$x = x[, 1]
    fitList$densfun = dis                                                   ####
    if (!missing(start)) 
      fitList$start = start
    if (not3distr)                                                          ####
    {                                                                       ####
      fittedDistr = do.call(MASS::fitdistr, fitList)
      estimates = as.list(fittedDistr$estimate)
      paramsList = estimates
    }                                                                       ####
    if (distribution=="weibull3")                                           ####
    {                                                                       ####
      paramsList= .weibull3(x[,1])                                           ####
      estimates = paramsList                                                 ####
    }                                                                       ####
    if (distribution=="lognormal3")                                         ####
    {                                                                       ####
      paramsList= .lognormal3(x[,1])                                         ####
      estimates = paramsList                                                 ####
    }                                                                       ####
    if (distribution=="gamma3")                                             ####
    {                                                                       ####
      paramsList= .gamma3(x[,1])                                             ####
      estimates = paramsList                                                 ####
    }                                                                       ####
  }
  paramsList = c(paramsList, .lfkp(parList, formals(qFun)))
  if (distribution == "normal") {
    paramsList$mean = center
    paramsList$sd = std.dev
    estimates = paramsList
  }
  if (boxcox) {
    if (!is.null(lsl)) 
      lsl = lsl^lambda
    if (!is.null(usl)) 
      usl = usl^lambda
    if (!is.null(target)) 
      target = target^lambda
  }
  if (is.null(lsl) && is.null(usl)) {
    paramsList$p = confLow
    lsl = do.call(qFun, paramsList)
    paramsList$p = confHigh
    usl = do.call(qFun, paramsList)
  }
  if (identical(lsl, usl)) 
    stop("lsl == usl")
  if (!is.null(lsl) && !is.null(target) && target < lsl) 
    stop("target is less than lower specification limit")
  if (!is.null(usl) && !is.null(target) && target > usl) 
    stop("target is greater than upper specification limit")
  if (!is.null(lsl) && !is.null(usl)) 
    if (lsl > usl) {
      temp = lsl
      lsl = usl
      usl = temp                       
    }
  paramsList$p = c(confLow, 0.5, confHigh)
  paramsListTemp = .lfkp(paramsList, formals(qFun))                           ####
  qs = do.call(qFun, paramsListTemp)                                         ####  
  paramsListTemp = .lfkp(paramsList, formals(pFun))                           ####
  if (!is.null(lsl) && !is.null(usl)) 
    cp = (usl - lsl)/(qs[3] - qs[1])
  if (!is.null(usl)) {
    cpu = (usl - qs[2])/(qs[3] - qs[2])
    paramsListTemp$q = usl                                                  ####
    ppu = 1 - do.call(pFun, paramsListTemp)                                 ####
  }
  if (!is.null(lsl)) {
    cpl = (qs[2] - lsl)/(qs[2] - qs[1])
    paramsListTemp$q = lsl                                                  ####
    ppl = do.call(pFun, paramsListTemp)                                     ####
  }
  cpk = min(cpu, cpl)
  ppt = sum(ppl, ppu)
  
  if (missing(xlim)) {
    xlim <- range(x[, 1], usl, lsl)
    xlim <- xlim + diff(xlim) * c(-0.2, 0.2)
  }
  xVec <- seq(min(xlim), max(xlim), length = 200)
  dParamsList = .lfkp(paramsList, formals(dFun))
  dParamsList$x = xVec
  yVec = do.call(dFun , dParamsList)
  histObj <- hist(x[, 1], plot = FALSE)
  if (missing(ylim)) {
    ylim <- range(histObj$density, yVec)
    ylim <- ylim + diff(ylim) * c(0, 0.05)
  }
  par(mar = c(0, 0, 0, 0) + 0.1)
  par(oma = c(2, 4, 7, 4) + 0.1)
  do.call(hist, c(list(x[, 1], freq = FALSE, xlim = xlim, ylim = ylim, main = ""), parList))
  abline(h = 0, col = "gray")
  tempList = parList
  tempList$col = 1
  tempList$border = NULL
  do.call(box, tempList)
  lines(xVec, yVec, lwd = lineWidth, col = lineCol, lty = lineType)
  abline(v = usl, col = specCol, lwd = specWidth, lty = 5)
  abline(v = lsl, col = specCol, lwd = specWidth, lty = 5)
  abline(v = target, col = specCol, lwd = specWidth, lty = 5)
  if (!is.null(lsl)) 
    axis(side = 3, at = lsl, labels = paste("LSL =", format(lsl, digits = 3)), col = specCol)
  if (!is.null(usl)) 
    axis(side = 3, at = usl, labels = paste("USL =", format(usl, digits = 3)), col = specCol)
  if (!is.null(lsl) && !is.null(usl)) 
    axis(side = 3, at = c(lsl, usl), labels = c(paste("LSL =", format(lsl, digits = 3)), paste("USL =", 
                                                                                               format(usl, digits = 3))), col = specCol)
  if (!is.null(target)) 
    text(target, max(ylim), "TARGET", pos = 1, col = cex.col, cex = cex.text)
  title(main = main, outer = TRUE)
  
  return(list(lambda = lambda, cp = cp, cpk = cpk, cpl = cpl, cpu = cpu,        ####
              ppt = ppt, ppl = ppl, ppu = ppu, usl = usl,                  ####
              lsl = lsl, target = target))                                 ####
}

#  -------------- int_qq -------------- 
.confintbeta= function(thethas, varmatrix, alpha) {
  
  th1= thethas[[1]]
  th2 =  thethas[[2]]
  
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
  
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
.confintcauchy = function(thethas, varmatrix, alpha) {
  
  th1= thethas[[1]]
  th2 =  thethas[[2]]
  
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
  
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
.confintexp=function(thethas, varmatrix, alpha) {
  lambda=thethas[[1]]
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
  perzentile=qexp(prozent, lambda)
  logPerzentile = log(perzentile)
  zalpha=qnorm(1-alpha/2)
  halfwidth = zalpha*sqrt(varmatrix[[1, 1]]/lambda^2)
  lci = exp(logPerzentile - halfwidth);
  uci = exp(logPerzentile + halfwidth);
  
  bounds = list(lci, uci, perzentile)
  
  return (bounds)
  
}
.confintgamma= function(thethas, varmatrix, alpha) {
  th1= thethas[[1]]
  th2 =  thethas[[2]]
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
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
.confintlnorm=function(thethas, varmatrix, alpha){
  th1= thethas[[1]]
  th2 =  thethas[[2]]
  
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
  perzentile=qlnorm(prozent, th1, th2)
  
  zp=qnorm(prozent)
  
  varPerzentile = varmatrix[[1, 1]]+2*varmatrix[[1, 2]]*zp+varmatrix[[2, 2]]*zp*zp
  
  zalpha=qnorm(1-alpha/2)
  lci=log(perzentile)-zalpha*sqrt(varPerzentile)
  uci=log(perzentile)+zalpha*sqrt(varPerzentile)
  
  bounds = list(exp(lci), exp(uci), perzentile)
  
  return (bounds)
  
}
.confintlogis= function(thethas, varmatrix, alpha) {
  
  th1= thethas[[1]]
  th2 =  thethas[[2]]
  
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
  
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
.confintnorm=function(thethas, varmatrix, alpha){
  
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
  zp=qnorm(prozent)
  perzentile=qnorm(prozent, thethas[[1]], thethas[[2]])
  
  varPerzentile = varmatrix[[1, 1]]+2*varmatrix[[1, 2]]*zp+varmatrix[[2, 2]]*zp*zp
  
  zalpha=qnorm(1-alpha/2)
  lci=perzentile-zalpha*sqrt(varPerzentile)
  uci=perzentile+zalpha*sqrt(varPerzentile)
  
  bounds = list(lci, uci, perzentile)
  
  
  return (bounds)
}
.confintweibull= function(thethas, varmatrix, alpha) {
  th1= thethas[[1]]
  th2 =  thethas[[2]]
  
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
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
  #  if(any(x < 0))                                                               ####
  #    stop("x must be positive")                                                 ####
  
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

#  -------------- whi_t.r -------------- 
setGeneric("whiskersPlot", function(x, main, xlab, ylab, col, ylim, legend = TRUE, ...) standardGeneric("whiskersPlot"))
setMethod("whiskersPlot", signature(x = "gageRR"), function(x, main, xlab, ylab, col, ylim, legend = TRUE, ...) {
  ops = length(unique(x[, 3]))
  pts = length(unique(x[, 4]))
  n = length(x[, 5])/pts/ops
  val = as.data.frame(x)
  dat = val
  if (missing(xlab)) 
    xlab = paste(names(val)[4], "s", sep = "")
  if (missing(ylab)) 
    ylab = "Value"
  if (missing(main)) 
    main = paste("Whiskers Plot")
  if (missing(col)) 
    col = heat.colors(ops)
  if (length(col) < ops) 
    col = rep(col, ops)
  if (missing(ylim)) {
    max_y = numeric(pts)
    min_y = numeric(pts)
    for (i in 1:pts) {
      dat = subset(val, val[, 4] == unique(x[, 4])[i])
      max_y[i] = max(dat[, 5])
      min_y[i] = min(dat[, 5])
    }
    ylim = c(-1.25 * abs(min(min_y)), abs(max(max_y)))
  }
  plot(x = 1:((pts * n * ops) + pts + 1), y = rep(0, pts * n * ops + pts + 1), col = "transparent", xlim = c(0, pts * n * ops + pts + 1), ylim = ylim, xlab = xlab, 
       ylab = ylab, main = main, axes = FALSE)
  abline(v = seq(1, pts * n * ops + pts + 1, by = n * ops + 1), col = "gray")
  abline(h = 0, col = "gray", lty = 3)
  box()
  axis(2)
  axis(1, at = seq(1 + (n * ops/2), pts * n * ops + pts + 1 - (n * ops/2), length = pts), labels = unique(x[, 4]))
  j = 0
  for (i in 1:pts) {
    dat = subset(val, val[, 4] == unique(x[, 4])[i])
    dat = dat[order(dat[, 3]), ]
    j = j + 1
    for (k in 1:ops) {
      Max = max(dat[(((k * n) - n + 1):(k * n)), 5])
      Min = min(dat[(((k * n) - n + 1):(k * n)), 5])
      Mean = median(dat[(((k * n) - n + 1):(k * n)), 5])
      rect(j + 1, Min, j + n, Max, col = col[k])
      lines(x = c(j + 1, j + n), y = c(Mean, Mean), col = "black", lwd = 2)
      j = j + n
    }
  }
  if (legend == TRUE) 
    legend("bottomleft", legend = unique(x[, 3]), fill = col, horiz = TRUE, box.col = 1, bg = "white", title = paste(names(val)[3], "(s):", sep = ""), inset = 0.04)
  return()
}) 

#  -------------- err_t.r -------------- 
setGeneric("errorPlot", function(x, main, xlab, ylab, col, pch, type, ylim, legend = TRUE, ...) standardGeneric("errorPlot"))
setMethod("errorPlot", signature(x = "gageRR"), function(x, main, xlab, ylab, col, pch, type, ylim, legend = TRUE, ...) {
  ops = length(unique(x[, 3]))
  pts = length(unique(x[, 4]))
  n = length(x[, 5])/pts/ops
  ref = mean(response(x))
  val = as.data.frame(x)
  dat = val
  if (missing(xlab)) 
    xlab = paste(names(val)[4], "s", sep = "")
  if (missing(ylab)) 
    ylab = "Error"
  if (missing(main)) 
    main = paste("Error plot")
  if (missing(pch)) 
    pch = 1:ops
  if (length(pch) < ops) 
    pch = rep(pch, ops)
  if (missing(col)) 
    col = rainbow(ops)
  if (length(col) < ops) 
    col = rep(col, ops)
  if (missing(type)) 
    type = "o"
  if (missing(ylim)) {
    max_y = numeric(pts)
    min_y = numeric(pts)
    for (i in 1:pts) {
      dat = subset(val, val[, 4] == unique(x[, 4])[i])
      max_y[i] = max(dat[, 5] - mean(dat[, 5]))
      min_y[i] = min(dat[, 5] - mean(dat[, 5]))
    }
    ylim = c(-1.25 * abs(min(min_y)), abs(max(max_y)))
  }
  plot(x = 1:((pts * n * ops) + pts + 1), y = rep(0, pts * n * ops + pts + 1), col = "transparent", xlim = c(0, pts * n * ops + pts + 1), ylim = ylim, xlab = xlab, 
       ylab = ylab, main = main, axes = FALSE)
  abline(v = seq(1, pts * n * ops + pts + 1, by = n * ops + 1), col = "gray")
  abline(h = 0, col = "gray", lty = 3)
  box()
  axis(2)
  axis(1, at = seq(1 + (n * ops/2), pts * n * ops + pts + 1 - (n * ops/2), length = pts), labels = unique(x[, 4]))
  j = 0
  for (i in 1:pts) {
    dat = subset(val, val[, 4] == unique(x[, 4])[i])
    dat = dat[order(dat[, 3]), ]
    j = j + 1
    for (k in 1:ops) {
      points(x = (j + 1):(j + n), y = (dat[(((k * n) - n + 1):(k * n)), 5] - mean(dat[, 5])), pch = pch[k], col = col[k], type = type, ...)
      j = j + n
    }
  }
  if (legend == TRUE) 
    legend("bottomleft", legend = unique(x[, 3]), pch = pch, col = col, horiz = TRUE, box.col = 1, bg = "white", title = paste(names(val)[3], "(s):", sep = ""), 
           inset = 0.04)
  return()
})

#  -------------- ste_t.r -------------- 
setClass(Class = "steepAscent", representation = representation(name = "character", X = "data.frame", response = "data.frame"))
setMethod("response", "steepAscent", function(object) {
  out = object@response
  return(out)
})
setReplaceMethod("response", "steepAscent", function(object, value) {
  if (is.vector(value)) {
    temp = data.frame(value)
    names(temp) = deparse(substitute(value))
    if (nrow(object@X) == nrow(temp)) {
      object@response = temp
      return(object)
    }
    stop("number of rows differ!")
  }
  if (is.data.frame(value)) {
    if (nrow(object@X) == nrow(value)) {
      object@response = value
      return(object)
    }
    stop("number of rows differ!")
  }
  stop(paste(deparse(substitute(value)), " needs to be a vector or data.frame"))
})
setMethod("[", signature(x = "steepAscent", i = "ANY", j = "ANY"), function(x, i, j) {
  bound = ncol(x@X)
  if (j <= bound) 
    x@X[i, j]
  else x@response[i, j - bound]
})
setMethod("as.data.frame", "steepAscent", function(x, row.names = NULL, optional = FALSE, ...) {
  return(cbind(x@X, x@response))
})
as.data.frame.steepAscent = function(x, row.names = NULL, optional = FALSE, ...) {
  return(cbind(x@X, x@response))
}
setMethod("show", signature(object = "steepAscent"), function(object) {
  print(as.data.frame(object))
})
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "steepAscent"), function(x, y, ...) {
  Delta = (x@X)$Delta
  frame = cbind(Delta, response(x))
  names(frame) = c("Delta", names(response(x)))
  plot(frame, ...)
})
steepAscent = function(factors, response, size = 0.2, steps = 5, data) {
  DB = FALSE
  if (missing(data)) 
    return("missing an object of class 'facDesign'")
  else fdo = data
  if (missing(factors) | length(factors) < 1) 
    return("missing factors")
  if (!is.character(factors)) 
    return("factors needs to be a character")
  names(names(fdo))
  model = fits(data)[[response]]
  if (DB) 
    print(model)
  if (is.null(model)) {
    form = c(response, "~")
    for (i in seq(along = factors)) {
      if (i == 1) 
        form = c(form, factors[i])
      else form = c(form, "+", factors[i])
    }
    form = paste(form, collapse = "")
    model = lm(form, data = fdo)
  }
  if (DB) 
    print(model)
  b = numeric(length = length(factors))
  x = numeric(length = length(factors))
  names(x) = factors
  for (i in seq(along = factors)) {
    b[i] = coef(model)[factors[i]]
    if (i == 1) {
      x[i] = size * sign(b[i])
    }
    else {
      if (DB) {
        print(x[1])
        print(b[1])
        print(b[i])
      }
      x[i] = (x[1]/b[1]) * b[i]
    }
  }
  if (DB) 
    print(x)
  Run = 1:(steps + 1)
  Delta = 0:steps
  frameOut = data.frame(Run, Delta)
  initial = ncol(frameOut)
  for (i in seq(along = factors)) {
    frameOut[, i + initial] = x[i] * 0:steps
    names(frameOut)[i + initial] = paste(factors[i], ".coded", collapse = "", sep = "")
    if (DB) 
      print(factors[i])
  }
  initial = ncol(frameOut)
  for (i in seq(along = factors)) {
    frameOut[, i + initial] = code2real(lows(fdo)[[factors[i]]], highs(fdo)[[factors[i]]], x[i] * 0:steps)
    names(frameOut)[i + initial] = paste(factors[i], ".real", collapse = "", sep = "")
    if (DB) 
      print(factors[i])
  }
  soa = new("steepAscent")
  soa@X = frameOut
  response(soa) = rep(NA, times = nrow(frameOut))
  names(response(soa)) = deparse(substitute(response))
  cat("\n")
  cat(paste("Steepest Ascent for", deparse(substitute(data)), "\n"))
  cat("\n")
  print(format(frameOut, digits = 3))
  invisible(soa)
} 


#  -------------- nor_t.r -------------- 
normalPlot = function(fdo, threeWay = FALSE, na.last = NA, alpha = 0.05, response = NULL, sig.col = c("red1", "red2", "red3"), sig.pch = c(1,2,3), main, ylim, xlim, xlab, ylab, pch,  ###
                      col, border = "red", ...) {
  DB = FALSE
  #require(MASS, quietly = TRUE)
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  fdoName = deparse(substitute(fdo))                                          ###
  if(is.null(response)==FALSE)                                                ###
  {                                                                           ###
    temp=response(fdo)[response]                                               ###
    response(fdo)=temp                                                         ###
  }                                                                           ###
  parList = list(...)
  params = list()
  if (length(sig.col) < 3) 
    sig.col = as.vector(matrix(sig.col, nrow = 1, ncol = 3))
  XLIM=FALSE;YLIM=FALSE                                                       ###
  if (!(class(fdo) == "facDesign")) 
    stop(paste(deparse(substitute(fdo)), "is not an object of class facDesign"))
  mainmiss = FALSE                                                            ###
  if (missing(main))                                                          ###
    mainmiss = TRUE                                                         ###
  if (missing(ylim))                                                          ###
    YLIM=TRUE                                                               ###
  if (missing(xlim))                                                          ###
    XLIM=TRUE                                                               ###
  if (missing(xlab)) 
    xlab = "Coefficients"
  if (missing(ylab)) 
    ylab = "Theoretical Quantiles"
  if (missing(pch)) 
    pch = 19
  if (missing(col)) 
    col = "black"
  for (j in 1:ncol(response(fdo))) {
    parList = list(...)                                                     ###
    params = list()                                                         ###
    leg.col = vector()                                                      ###
    p.col = vector()                                                        ###
    p.pch = vector()                                                        ###
    leg.txt = vector()                                                      ###
    main = paste("Normal plot for", names(response(fdo))[j], "in", fdoName) ###
    if (j > 1) 
      dev.new()
    form = paste("response(fdo)[,", j, "]~")
    for (i in 1:ncol(cube(fdo))) {
      form = paste(form, names(cube(fdo))[i], sep = "")
      if (i < ncol(cube(fdo))) 
        form = paste(form, "*", sep = "")
    }
    if (DB == TRUE) 
      print(paste("form:", form))
    lm.1 = lm(as.formula(form), data = as.data.frame(fdo))
    lm.1s = summary(lm.1)
    effect = coef(lm.1s)[row.names(coef(lm.1s)) != "(Intercept)", "t value"]
    print(effect)
    if (all(is.na(effect))) 
      effect = 2 * coef(lm.1)[-pmatch("(Intercept)", names(coef(lm.1)))]      ###
    #            stop("effects could not be calculated")                            ###
    sig = summary(lm.1)$coefficients[, "Pr(>|t|)"][-pmatch("(Intercept)", names(coef(lm.1)))]
    df.resid = df.residual(lm.1)
    nc = nrow(centerCube(fdo))
    if (DB) {
      print(paste("effect:", effect))
      print(paste("df.resid:", df.resid))
      print(paste("nc:", nc))
      print(paste("sig:", sig))
      print(summary(lm.1))
    }
    tQ = ppoints(effect)
    index = order(effect)
    sQ = effect[index]
    sig = sig[index]
    if (df.resid > 0) {
      for (k in seq(along = sig)) {
        setted = FALSE
        if (abs(sig)[k] < 0.01) {
          if (!setted) {
            p.col[k] = sig.col[1]
            p.pch[k] = sig.pch[1]
            leg.txt = c(leg.txt, "p < 0.01")
            leg.col = c(leg.col, p.col)
            setted = TRUE
          }
        }
        if (abs(sig)[k] < 0.05) {
          if (!setted) {
            p.col[k] = sig.col[2]
            p.pch[k] = sig.pch[2]
            leg.txt = c(leg.txt, "p < 0.05")
            leg.col = c(leg.col, p.col)
            setted = TRUE
          }
        }
        if (abs(sig)[k] < 0.1) {
          if (!setted) {
            p.col[k] = sig.col[3]
            p.pch[k] = sig.pch[3]
            leg.txt = c(leg.txt, "p < 0.1")
            leg.col = c(leg.col, p.col)
            setted = TRUE
          }
        }
        if (abs(sig)[k] >= 0.1) {
          if (!setted) {
            p.col[k] = col
            p.pch[k] = pch 
            leg.txt = c(leg.txt, "p >= 0.1")
            leg.col = c(leg.col, p.col)
            setted = TRUE
          }
        }
      }
      leg.txt = unique(leg.txt)
      leg.col = unique(leg.col)
    }
    else                                                                    ###
    {p.col=col
    p.pch=pch}                                                              ###
    mid = round(length(tQ)/2)
    last = length(tQ)
    params$p = ppoints(effect)
    estimates = MASS::fitdistr(effect, "normal")
    params$mean = estimates$estimate[["mean"]]
    params$sd = estimates$estimate[["sd"]]                                                      
    y = do.call(qnorm, params)
    if (XLIM)                                                               ###
      xlim = range(sQ)
    if (YLIM)                                                               ###
      ylim = range(y)
    params = .lfkp(parList, c(formals(plot.default), par()))
    params$x = sQ
    params$y = y
    params$xlab = xlab
    params$ylab = ylab
    params$main = main
    params$xlim = xlim
    params$ylim = ylim
    params$pch = p.pch
    params$col = p.col
    
    do.call(plot, params)
    xp = c(qnorm(0.1), qnorm(0.99))
    yp = c(qnorm(0.1, mean = estimates$estimate[["mean"]], sd = estimates$estimate[["sd"]]), qnorm(0.99, mean = estimates$estimate[["mean"]], sd = estimates$estimate[["sd"]]))
    slope = (yp[2] - yp[1])/(xp[2] - xp[1])
    int = yp[1] - slope * xp[1]
    abline(a = int, b = slope, col = border)
    text(sQ[1:mid], y[1:mid], names(sQ)[1:mid], pos = 4)
    text(sQ[(mid + 1):last], y[(mid + 1):last], names(sQ)[(mid + 1):last], pos = 2)
    if (df.resid > 0) 
      legend("topleft", legend = leg.txt, col = leg.col, pch = p.pch, inset = 0.02)
  }
} 

#  -------------- vis_n.r -------------- 
wirePlot = function(x, y, z, data = NULL, xlim, ylim, zlim, main, xlab, ylab, border, sub, zlab, form = "fit", phi, theta, ticktype, col = 1, steps, 
                    factors, fun, plot) {
  DB = FALSE
  form = form
  fact = NULL
  if (missing(steps)) 
    steps = 25
  fdo = data
  fit = NULL
  lm.1 = NULL
  if (!is.function(col)) {
    if (identical(col, 1)) 
      col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if (identical(col, 2)) 
      col = colorRampPalette(c("blue", "white", "red"), space = "Lab")
    if (identical(col, 3)) 
      col = colorRampPalette(c("blue", "white", "orange"))
    if (identical(col, 4)) 
      col = colorRampPalette(c("gold", "white", "firebrick"))
  }
  if (is.null(data)) {
    cat("\n defaulting to persp function\n")
    return("persp")
  }
  if (class(data) != "facDesign") {
    cat("\n defaulting to persp function using formula\n")
    return("persp")
  }
  x.c = deparse(substitute(x))
  y.c = deparse(substitute(y))
  z.c = deparse(substitute(z))
  if (missing(plot)) 
    plot = TRUE
  if (missing(main)) 
    main = paste("Response Surface for", z.c)
  if (missing(ylab)) 
    ylab = paste(y.c, ": ", names(fdo)[[y.c]])
  if (missing(xlab)) 
    xlab = paste(x.c, ": ", names(fdo)[[x.c]])
  if (missing(zlab)) 
    zlab = paste(x.c, ": ", names(response(fdo)))
  if (missing(ticktype)) 
    ticktype = "detailed"
  if (missing(border)) 
    border = NULL
  if (missing(phi)) 
    phi = 30
  if (missing(theta)) 
    theta = -30
  if (missing(factors)) 
    factors = NULL
  if (missing(xlim)) 
    xlim = c(min(fdo[, x.c]), max(fdo[, x.c]))
  if (missing(ylim)) 
    ylim = c(min(fdo[, y.c]), max(fdo[, y.c]))
  allVars = c(names(names(fdo)), names(response(fdo)))
  isct = intersect(c(x.c, y.c, z.c), c(names(names(fdo)), names(response(fdo))))
  if (DB) {
    print(allVars)
    print(isct)
  }
  if (length(isct) < length(c(x.c, y.c, z.c))) {
    d = setdiff(isct, allVars)
    stop(paste(d, "could not be found\n"))
  }
  if (missing(fun)) 
    fun = NULL
  if (!is.function(fun) & !is.null(fun)) 
    if (!(fun %in% c("overall", "desirability"))) 
      stop("fun should be a function, \"overall\" or \"desirability\"")
  if (identical(fun, "desirability")) {
    obj = desires(fdo)[[z.c]]
    fun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
  }
  if (form %in% c("fit")) {
    lm.1 = fits(fdo)[[z.c]]
    if (DB) 
      print(lm.1)
    if (is.null(fit)) 
      form = "full"
  }
  if (form %in% c("quadratic", "full", "interaction", "linear")) {
  }
  if (identical(form, "interaction")) {
    form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
  }
  if (identical(form, "linear")) {
    form = paste(z.c, "~", x.c, "+", y.c)
  }
  if (identical(form, "quadratic")) {
    form = paste(z.c, "~I(", x.c, "^2) + I(", y.c, "^2)")
  }
  if (identical(form, "full")) {
    form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
    if (nrow(star(fdo)) > 0) 
      form = paste(form, "+ I(", x.c, "^2) + I(", y.c, "^2)")
    if (DB) 
      print(form)
  }
  if (is.null(form)) 
    stop(paste("invalid formula", form))
  if (is.null(lm.1)) 
    lm.1 = lm(form, data = fdo)
  if (missing(sub)) 
    sub = deparse(formula(lm.1))
  if (DB) 
    print(lm.1)
  dcList = vector(mode = "list", length = length(names(fdo)))
  names(dcList) = names(names(fdo))
  dcList[1:length(names(fdo))] = 0
  if (!is.null(factors)) {
    for (i in names(factors)) dcList[[i]] = factors[[i]][1]
  }
  if (DB) 
    print(dcList)
  help.predict = function(x, y, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = x
    dcList[[y.c]] = y
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
  }
  if (DB) {
    print(x.c)
    print(y.c)
    print(help.predict(1, 2, "A", "B", lm.1 = lm.1))
    print(help.predict(1, 2, x.c, y.c, lm.1 = lm.1))
  }
  xVec = seq(min(xlim), max(xlim), length = steps)
  yVec = seq(min(ylim), max(ylim), length = steps)
  mat = outer(xVec, yVec, help.predict, x.c, y.c, lm.1)
  if (is.function(fun)) 
    mat = try(apply(mat, c(1, 2), fun))
  if (identical(fun, "overall")) {
    main = "composed desirability"
    mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
    for (i in names(response(fdo))) {
      obj = desires(fdo)[[i]]
      fun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
      temp = outer(xVec, yVec, help.predict, x.c, y.c, fits(fdo)[[i]])
      temp = try(apply(temp, c(1, 2), fun))
      mat = mat * temp
    }
    mat = mat^(1/length(names(response(fdo))))
  }
  if (is.function(col)) {
    nrMat <- nrow(mat)
    ncMat <- ncol(mat)
    jet.colors <- colorRampPalette(c("blue", "green"))
    nbcol <- 100
    color <- col(nbcol)
    matFacet <- mat[-1, -1] + mat[-1, -ncMat] + mat[-nrMat, -1] + mat[-nrMat, -ncMat]
    facetcol <- cut(matFacet, nbcol)
  }
  else {
    color = col
    facetcol = 1
  }
  if (plot) {
    if (missing(zlim)) 
      zlim = range(mat)
    persp(xVec, yVec, mat, main = main, sub = sub, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, zlim = zlim, zlab = zlab, col = color[facetcol], 
          border = border, ticktype = ticktype, phi = phi, theta = theta)
    if (is.function(col)) {
      zlim = range(mat)
      leglevel = pretty(zlim, 6)
      legcol = col(length(leglevel))
      legpretty = as.character(abs(leglevel))
      temp = character(length(leglevel))
      temp[leglevel > 0] = "+"
      temp[leglevel < 0] = "-"
      temp[leglevel == 0] = " "
      legpretty = paste(temp, legpretty, sep = "")
      legend("topright", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
    }
  }
  invisible(list(x = xVec, y = yVec, z = mat))
}
.mfc = function(x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), 
                zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels), nlevels = 20, palette = cm.colors, col = palette(length(levels) - 1), plot.title,        ###
                plot.axes, key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, ...) {
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  mar = c(4.1, 4.1, 4.1, 4.1)
  par(mar = mar)
  plot(1, 1, type = "n", axes = FALSE, xlim, ylim, xlab = "", ylab = "", main = "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), col = col(length(levels)-1))
  leglevel = pretty(zlim, 6)
  legcol = col(length(leglevel))                                          ###
  legpretty = as.character(abs(leglevel))
  temp = character(length(leglevel))
  temp[leglevel > 0] = "+"
  temp[leglevel < 0] = "-"
  temp[leglevel == 0] = " "
  legpretty = paste(temp, legpretty, sep = "")
  legend("topright", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}
contourPlot = function(x, y, z, data = NULL, xlim, ylim, main, xlab, ylab, form = "fit", col = 1, steps, factors, fun) {
  DB = FALSE
  form = form
  fact = NULL
  if (missing(steps)) 
    steps = 25
  fdo = data
  fit = NULL
  lm.1 = NULL
  if (!is.function(col)) {
    if (identical(col, 1)) 
      col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if (identical(col, 2)) 
      col = colorRampPalette(c("blue", "white", "red"), space = "Lab")
    if (identical(col, 3)) 
      col = colorRampPalette(c("blue", "white", "orange"))
    if (identical(col, 4)) 
      col = colorRampPalette(c("gold", "white", "firebrick"))
    if (identical(col, 5)) 
      col = colorRampPalette(c("blue4", "lightblue1", "lightgreen", "green4"))
  }
  if (!is.function(col))
    col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  if (is.null(data)) {
    cat("\n defaulting to filled.contour function\n")
    return("persp")
  }
  if (class(data) != "facDesign") {
    cat("\n defaulting to filled.contour function using formula\n")
    return("persp")
  }
  x.c = deparse(substitute(x))
  y.c = deparse(substitute(y))
  z.c = deparse(substitute(z))
  if (missing(main)) 
    main = paste("Filled Contour for", z.c)
  if (missing(ylab)) 
    ylab = paste(y.c, ": ", names(fdo)[[y.c]])
  if (missing(xlab)) 
    xlab = paste(x.c, ": ", names(fdo)[[x.c]])
  if (missing(factors)) 
    factors = NULL
  if (missing(xlim)) 
    xlim = c(min(fdo[, x.c]), max(fdo[, x.c]))
  if (missing(ylim)) 
    ylim = c(min(fdo[, y.c]), max(fdo[, y.c]))
  allVars = c(names(names(fdo)), names(response(fdo)))
  isct = intersect(c(x.c, y.c, z.c), c(names(names(fdo)), names(response(fdo))))
  if (DB) {
    print(allVars)
    print(isct)
  }
  if (length(isct) < length(c(x.c, y.c, z.c))) {
    d = setdiff(isct, allVars)
    stop(paste(d, "could not be found\n"))
  }
  if (missing(fun)) 
    fun = NULL
  if (!is.function(fun) & !is.null(fun)) 
    if (!(fun %in% c("overall", "desirability"))) 
      stop("fun should be a function, \"overall\" or \"desirability\"")
  if (identical(fun, "desirability")) {
    obj = desires(fdo)[[z.c]]
    fun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
  }
  if (form %in% c("fit")) {
    lm.1 = fits(fdo)[[z.c]]
    if (DB) 
      print(lm.1)
    if (is.null(fit)) 
      form = "full"
  }
  if (form %in% c("quadratic", "full", "interaction", "linear")) {
  }
  if (identical(form, "interaction")) {
    form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
  }
  if (identical(form, "linear")) {
    form = paste(z.c, "~", x.c, "+", y.c)
  }
  if (identical(form, "quadratic")) {
    form = paste(z.c, "~I(", x.c, "^2) + I(", y.c, "^2)")
  }
  if (identical(form, "full")) {
    form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
    if (nrow(star(fdo)) > 0) 
      form = paste(form, "+ I(", x.c, "^2) + I(", y.c, "^2)")
    if (DB) 
      print(form)
  }
  if (is.null(form)) 
    stop(paste("invalid formula", form))
  if (is.null(lm.1)) 
    lm.1 = lm(form, data = fdo)
  if (DB) 
    print(lm.1)
  dcList = vector(mode = "list", length = length(names(fdo)))
  names(dcList) = names(names(fdo))
  dcList[1:length(names(fdo))] = 0
  if (!is.null(factors)) {
    for (i in names(factors)) dcList[[i]] = factors[[i]][1]
  }
  if (DB) 
    print(dcList)
  help.predict = function(x, y, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = x
    dcList[[y.c]] = y
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
  }
  if (DB) {
    print(x.c)
    print(y.c)
    print(help.predict(1, 2, "A", "B", lm.1 = lm.1))
    print(help.predict(1, 2, x.c, y.c, lm.1 = lm.1))
  }
  xVec = seq(min(xlim), max(xlim), length = steps)
  yVec = seq(min(ylim), max(ylim), length = steps)
  mat = outer(xVec, yVec, help.predict, x.c, y.c, lm.1)
  if (is.function(col)) {
    nrMat <- nrow(mat)
    ncMat <- ncol(mat)
    nbcol <- 1000
    color <- col(nbcol)
    matFacet <- mat[-1, -1] + mat[-1, -ncMat] + mat[-nrMat, -1] + mat[-nrMat, -ncMat]
    facetcol <- cut(matFacet, nbcol)
  }
  else {
    color = col
    facetcol = 1
  }
  if (is.function(fun)) 
    mat = try(apply(mat, c(1, 2), fun))
  if (identical(fun, "overall")) {
    main = "composed desirability"
    mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
    for (i in names(response(fdo))) {
      obj = desires(fdo)[[i]]
      fun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
      temp = outer(xVec, yVec, help.predict, x.c, y.c, fits(fdo)[[i]])
      temp = try(apply(temp, c(1, 2), fun))
      mat = mat * temp
    }
    mat = mat^(1/length(names(response(fdo))))
  }
  
  .mfc(xVec, yVec, mat, main = main, xlab = xlab, ylab = ylab, col = col)
  
  invisible(list(x = xVec, y = yVec, z = mat))
} 

#  -------------- distr_3.r -------------- 
#### WEIBULL3
dweibull3 <- function(x,shape,scale,threshold){
  if(missing(x))
    stop("x must be a vector")
  if(missing(threshold))
    threshold=0
  if(missing(shape))
    shape=1
  if(missing(scale))
    scale=1
  temp=function(x)
  {
    if(x>=threshold)
      return((shape/scale)*(((x-threshold)/scale)^(shape-1))*exp(-((x-threshold)/scale)^shape))
    else
      return(0)
  }
  return(unlist(lapply(x,temp)))
}
pweibull3 <- function(q,shape,scale,threshold){
  if(missing(q))
    stop("q must be a vector")
  if(missing(threshold))
    threshold=0
  if(missing(shape))
    shape=1
  if(missing(scale))
    scale=1
  temp=function(q)
  {
    if(q>=threshold)
      return(1-exp(-((q-threshold)/scale)^shape))
    else
      return(0)
  }
  return(unlist(lapply(q,temp)))
}
qweibull3 <- function(p,shape,scale,threshold,...){
  if(missing(p))
    stop("p must be a vector")
  if(missing(threshold))
    threshold=0
  if(missing(shape))
    shape=1
  if(missing(scale))
    scale=1
  myfun = function(x,p) pweibull3(q = x,
                                  threshold = threshold, scale = scale, shape = shape) - p
  temp=function(p)
  {
    return(uniroot(f=myfun,lower=threshold,upper=threshold+10000000,p=p,...)$root)       #solve myfun=0
  }
  return(unlist(lapply(p,temp)))
}
#### LOG-NORM3
dlnorm3 <- function(x,meanlog,sdlog,threshold){
  if(missing(x))
    stop("x must be a vector")
  if(missing(threshold))
    threshold=0
  if(missing(meanlog))
    meanlog=0
  if(missing(sdlog))
    sdlog=1
  temp=function(x)
  {
    if(x>threshold)
      return((1/(sqrt(2*pi)*sdlog*(x-threshold)))*exp(-(((log((x-threshold))-meanlog)^2)/(2*(sdlog)^2))))
    else
      return(0)
  }
  return(unlist(lapply(x,temp)))
}
plnorm3 <- function(q,meanlog,sdlog,threshold){
  if(missing(q))
    stop("q must be a vector")
  if(missing(threshold))
    threshold=0
  if(missing(meanlog))
    meanlog=0
  if(missing(sdlog))
    sdlog=1
  temp=function(q)
  {
    if(q>threshold)
      return(pnorm((log((q-threshold))-meanlog)/sdlog))
    else
      return(0)
  }
  return(unlist(lapply(q,temp)))
}
qlnorm3 <- function(p,meanlog,sdlog,threshold,...){
  if(missing(p))
    stop("p must be a vector")
  if(missing(meanlog))
    meanlog=0
  if(missing(threshold))
    threshold=0
  if(missing(sdlog))
    sdlog=1
  myfun = function(x, p) plnorm3(q = x,
                                 meanlog = meanlog, sdlog = sdlog, threshold = threshold) - p
  temp=function(p)
  {
    return(uniroot(f=myfun,lower=threshold,upper=threshold+10000000,p=p,...)$root)       #solve myfun=0
  }
  return(unlist(lapply(p,temp)))
}
#### GAMMA3
dgamma3 <- function(x,shape,scale,threshold){
  if(missing(x))
    stop("x must be a vector")
  if(missing(threshold))
    threshold=0
  if(missing(shape))
    shape=1
  if(missing(scale))
    scale=1
  temp=function(x)
  {
    if(x>=threshold)
      return( dgamma(x-threshold,shape,scale) )
    return(0)
  }
  return(unlist(lapply(x,temp)))
}
pgamma3 <- function(q,shape,scale,threshold){
  if(missing(q))
    stop("q must be a vector")
  if(missing(threshold))
    threshold=-2
  if(missing(shape))
    shape=1
  if(missing(scale))
    scale=1
  temp=function(q)
  {
    if(q>=threshold)
      return(pgamma(q-threshold,shape,scale))
    else
      return(0)
  }
  return(unlist(lapply(q,temp)))
}
qgamma3 <- function(p,shape,scale,threshold,...){
  if(missing(p))
    stop("p must be a vector")
  if(missing(threshold))
    threshold=0
  if(missing(shape))
    shape=1
  if(missing(scale))
    scale=1
  myfun = function(x,p) pgamma3(q = x,
                                threshold = threshold, scale = scale, shape = shape) - p
  temp=function(p)
  {
    return(uniroot(f=myfun,lower=threshold,upper=threshold+10000000,p=p,...)$root)       #solve myfun=0
  }
  return(unlist(lapply(p,temp)))
}

#  -------------- rsm_n.r -------------- 
.alphaOrth = function(k, p = 0, cc, cs) {
  alpha = sqrt((2^(k - p) * (2 * k + cs))/(2 * (2^(k - p) + cc)))
  return(alpha)
}
.alphaRot = function(k, p = 0) {
  alpha = (2^(k - p))^0.25
  return(alpha)
}
.centerPoints = function(nrow, k) {
  if (!is.numeric(nrow))
    stop("nrow must be numeric")
  if (!is.numeric(k))
    stop("ncol must be numeric")
  mat = as.data.frame(matrix(rep(0, nrow * k), nrow, ncol = k))
  return(mat)
}
.starFrame = function(k, alpha) {
  if (!is.numeric(k))
    stop("k must be numeric")
  if (!is.numeric(alpha))
    stop("alpha must be numeric")
  .starFrame = as.data.frame(matrix(0, nrow = k * 2, ncol = k))
  for (j in 1:k) {
    for (i in (2 * (j - 1) + 1):(2 * (j - 1) + 2)) {
      .starFrame[i, j] = ((-1)^i) * alpha
    }
  }
  return(.starFrame)
}
rsmDesign = function(k = 3, p = 0, alpha = "rotatable", blocks = 1, cc = 1, cs = 1, fp = 1,
                     sp = 1, faceCentered = FALSE) {
  DB = FALSE
  if (blocks > 2^(k - 1) + 1)
    stop("Blocking not possible")
  if (alpha == "rotatable")
    alpha = .alphaRot(k, p)
  if (alpha == "orthogonal")
    alpha = .alphaOrth(k, p, cc = cc, cs = cs)
  if (alpha == "both") {
    found = FALSE
    for (i in seq(along = .rsmOrth)) {
      if (DB) {
        print(k)
        print(blocks)
      }
      if (.rsmOrth[[i]]$k == k)
        if (.rsmOrth[[i]]$blocks == blocks)
          if (.rsmOrth[[i]]$p == p) {
            cc = .rsmOrth[[i]]$cc
            cs = .rsmOrth[[i]]$cs
            p = .rsmOrth[[i]]$p
            alpha = .alphaOrth(k, p, cc, cs)
            found = TRUE
            break
          }
    }
    if (!found) {
      return("no design available")
    }
  }
  if (DB) {
    print("Values")
    print(k)
    print(alpha)
    print(cc)
    print(cs)
    print(blocks)
  }
  fdo = facDesign(k = k, p = p, replicates = fp)                              ###      
  if (cc > 0) {
    temp = as.data.frame(matrix(0, nrow = cc, ncol = ncol(cube(fdo))))
    names(temp) = names(cube(fdo))
    centerCube(fdo) = temp
    if (DB)
      print("centerCube")
  }
  if (DB)
    print("star not added")
  temp = .starFrame(k, alpha)
  starportion = data.frame()
  for (i in 1:sp) {
    starportion = rbind(temp, starportion)
  }
  names(starportion) = names(cube(fdo))
  star(fdo) = starportion
  if (DB)
    print("star added")
  if (cs > 0) {
    temp = as.data.frame(matrix(0, nrow = cs, ncol = ncol(cube(fdo))))
    names(temp) = names(cube(fdo))
    centerStar(fdo) = temp
  }
  #    return(fdo)
  fdo = blocking(fdo, blocks)
  return(fdo)
}
.rsmOrth = vector(mode = "list", length = 7)
.rsmOrth[[1]] = list(k = 2, p = 0, col = 1, row = 2, blocks = 2, cc = 3, cs = 3)
.rsmOrth[[2]] = list(k = 2, p = 0, col = 1, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[3]] = list(k = 3, p = 0, col = 2, row = 3, blocks = 3, cc = 2, cs = 2)
.rsmOrth[[4]] = list(k = 3, p = 0, col = 2, row = 2, blocks = 2, cc = 2, cs = 2)
.rsmOrth[[5]] = list(k = 3, p = 0, col = 2, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[6]] = list(k = 4, p = 0, col = 3, row = 4, blocks = 5, cc = 2, cs = 2)
.rsmOrth[[7]] = list(k = 4, p = 0, col = 3, row = 3, blocks = 3, cc = 2, cs = 2)
.rsmOrth[[8]] = list(k = 4, p = 0, col = 3, row = 2, blocks = 2, cc = 2, cs = 2)
.rsmOrth[[9]] = list(k = 4, p = 0, col = 3, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[10]] = list(k = 5, p = 0, col = 4, row = 5, blocks = 9, cc = 2, cs = 4)
.rsmOrth[[11]] = list(k = 5, p = 0, col = 4, row = 4, blocks = 5, cc = 2, cs = 4)
.rsmOrth[[12]] = list(k = 5, p = 0, col = 4, row = 3, blocks = 3, cc = 2, cs = 4)
.rsmOrth[[13]] = list(k = 5, p = 0, col = 4, row = 2, blocks = 2, cc = 2, cs = 4)
.rsmOrth[[14]] = list(k = 5, p = 0, col = 4, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[15]] = list(k = 5, p = 1, col = 5, row = 4, blocks = 5, cc = 6, cs = 1)
.rsmOrth[[16]] = list(k = 5, p = 1, col = 5, row = 3, blocks = 3, cc = 6, cs = 1)
.rsmOrth[[17]] = list(k = 5, p = 1, col = 5, row = 2, blocks = 2, cc = 6, cs = 1)
.rsmOrth[[18]] = list(k = 5, p = 1, col = 5, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[19]] = list(k = 6, p = 0, col = 6, row = 6, blocks = 17, cc = 1, cs = 6)
.rsmOrth[[20]] = list(k = 6, p = 0, col = 6, row = 5, blocks = 9, cc = 1, cs = 6)
.rsmOrth[[21]] = list(k = 6, p = 0, col = 6, row = 4, blocks = 5, cc = 1, cs = 6)
.rsmOrth[[22]] = list(k = 6, p = 0, col = 6, row = 3, blocks = 3, cc = 1, cs = 6)
.rsmOrth[[23]] = list(k = 6, p = 0, col = 6, row = 2, blocks = 2, cc = 1, cs = 6)
.rsmOrth[[24]] = list(k = 6, p = 0, col = 6, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[25]] = list(k = 6, p = 1, col = 7, row = 5, blocks = 9, cc = 4, cs = 2)
.rsmOrth[[26]] = list(k = 6, p = 1, col = 7, row = 4, blocks = 5, cc = 4, cs = 2)
.rsmOrth[[27]] = list(k = 6, p = 1, col = 7, row = 3, blocks = 3, cc = 4, cs = 2)
.rsmOrth[[28]] = list(k = 6, p = 1, col = 7, row = 2, blocks = 2, cc = 4, cs = 2)
.rsmOrth[[29]] = list(k = 6, p = 1, col = 7, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[30]] = list(k = 7, p = 0, col = 8, row = 6, blocks = 17, cc = 1, cs = 11)
#.rsmOrth[[31]] = list(k = 7, p = 0, col = 8, row = 6, blocks = 17, cc = 1, cs = 11) ###
.rsmOrth[[31]] = list(k = 7, p = 0, col = 8, row = 5, blocks = 9, cc = 1, cs = 11)
.rsmOrth[[32]] = list(k = 7, p = 0, col = 8, row = 4, blocks = 5, cc = 1, cs = 11)
.rsmOrth[[33]] = list(k = 7, p = 0, col = 8, row = 3, blocks = 3, cc = 1, cs = 11)
.rsmOrth[[34]] = list(k = 7, p = 0, col = 8, row = 2, blocks = 2, cc = 1, cs = 11)
.rsmOrth[[35]] = list(k = 7, p = 0, col = 8, row = 1, blocks = 1, cc = 0, cs = 0)
.rsmOrth[[36]] = list(k = 7, p = 1, col = 9, row = 5, blocks = 9, cc = 1, cs = 4)
rsmChoose = function() {
  DB = FALSE
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  colFun = colorRampPalette(c("yellow", "red"), space = "rgb")
  colPalette = colFun(169)
  numRows = 6
  numCol = 9
  blockVals = c(1, 2, 3, 5, 9, 17)
  factorVals = c(2, 3, 4, 5, 5, 6, 6, 7, 7)
  rsmList = .rsmOrth
  plot.new()
  par(mfrow = c(6, 9))
  par(mar = c(0, 0, 0, 0))
  par(oma = c(4, 4, 4, 4))
  for (i in 1:6) for (j in 1:9) {
    par(mfg = c(i, j))
    plot(0, 0, xaxs = "i", yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n",
         xlab = "", ylab = "", bg = "red", fg = "green")
    box()
  }
  for (i in seq(along = rsmList)) {
    temp = rsmList[[i]]
    par(mfg = c(temp$row, temp$col))
    par(mfg = c(temp$row, temp$col))
    plot(0, 0, xaxs = "i", yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n",
         xlab = "", ylab = "", bg = "red", fg = "green")
    rect(0, 0, 1, 1, col = colPalette[2^((temp$k) - (temp$p))])
    text(0.1, 0.9, paste("N =", 2^((temp$k) - (temp$p)) + temp$cc * (temp$blocks - 1) + temp$cs +
                           (temp$k + temp$p) * 2), adj = c(0, 1), cex = 1.25)
    text(0.1, 0.75, paste("k =", temp$k), adj = c(0, 1), cex = 1.25)
    text(0.1, 0.6, paste("p =", temp$p), adj = c(0, 1), cex = 1.25)
    text(0.1, 0.45, ".centerPoints", adj = c(0, 1), cex = 1.25)
    text(0.1, 0.3, paste("Cube:", temp$cc), adj = c(0, 1), cex = 1.25)
    text(0.1, 0.15, paste("Axial:", temp$cs), adj = c(0, 1), cex = 1.25)
    box()
  }
  x = 1/18 + (0:8) * 2/18
  mtext(factorVals, at = x, side = 3, line = 0.5, outer = TRUE)
  mtext("number of factors k", at = 0.5, side = 3, line = 2.5, outer = TRUE)
  x = 1/12 + (5:0) * 2/12
  mtext(blockVals, at = x, side = 2, line = 0.5, outer = TRUE, las = 2)
  mtext("number of blocks", at = 0.5, side = 2, line = 2.5, outer = TRUE)
  cat("\nChoose a response surface design by clicking into the appropriate field")
  cat("\nWaiting for your selection:")
  cat("\n\n")
  flush.console()
  if (DB)
    cat("TODO: standardize the locator part in respect to possible figure region")
  x = numeric(0)
  y = numeric(0)
  xyList = locator(1)                                                         ###
  print(xyList)
  x = ceiling(xyList$x + 8)
  y = ceiling(5 - xyList$y)
  if (DB) {
    print(paste("x:", x))
    print(paste("y:", y))
  }
  if (length(x) < 1)
    return(rsmDesign(k = 2, p = 0, blocks = 2, alpha = "both"))
  if (length(y) < 1)
    return(rsmDesign(k = 2, p = 0, blocks = 2, alpha = "both"))
  #    if (!(x %in% factorVals) || !(y %in% blockVals))                           ###
  #        return(rsmDesign(k = 2, p = 0, blocks = 2, alpha = "both"))            ###
  blocks = blockVals[y]
  k = factorVals[x]
  if (x==5 || x==7 || x==9 )                                                  ###
    p = 1                                                                      ###
  else                                                                        ###
    p = 0                                                                      ###
  if (DB) {
    print(paste("blocks:", blocks))
    print(paste("k:", k))
    #      print(rsmList)                                                         ###
  }
  for (i in seq(along = rsmList)) {
    if (rsmList[[i]]$k == k)
      if (rsmList[[i]]$blocks == blocks)
        if (rsmList[[i]]$p == p)                                        ###
          return(rsmDesign(k = k, p = rsmList[[i]]$p, blocks = blocks, ###
                           alpha = "both", cc = rsmList[[i]]$cc,                 ###
                           cs = rsmList[[i]]$cs))                                ###
  }
  return(cat("\nno selection recognized\n"))
}

#  -------------- vis_2.r -------------- 
contourPlot3 = function(x, y, z, response, data = NULL, main, xlab, ylab, zlab, border, form = "linear", col = 1, col.text, cex.axis, axes = TRUE, 
                        steps, factors) {
  DB = FALSE
  out = list()
  mdo = data
  x.c = deparse(substitute(x))
  y.c = deparse(substitute(y))
  z.c = deparse(substitute(z))
  r.c = deparse(substitute(response))
  if (missing(col)) 
    col = 1
  if (missing(col.text)) 
    col.text = 1
  if (missing(cex.axis)) 
    cex.axis = 1
  if (missing(main)) 
    main = paste("Response Surface for", r.c)
  if (missing(ylab)) 
    ylab = y.c
  if (missing(xlab)) 
    xlab = x.c
  if (missing(zlab)) 
    zlab = z.c
  if (missing(border)) 
    border = "white"
  if (missing(factors)) 
    factors = NULL
  if (missing(steps)) 
    steps = 100
  col.axis = par("col.axis")
  if (!is.function(col)) {
    if (identical(col, 1)) 
      col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if (identical(col, 2)) 
      col = colorRampPalette(c("blue", "white", "red"), space = "Lab")
    if (identical(col, 3)) 
      col = colorRampPalette(c("blue", "white", "orange"))
    if (identical(col, 4)) 
      col = colorRampPalette(c("gold", "white", "firebrick"))
  }
  nameVec = names(names(mdo))
  linStrings = "-1"
  for (i in seq(along = nameVec)) linStrings = paste(linStrings, "+", nameVec[i])
  if (DB) 
    print(linStrings)
  combList = combn(nameVec, 2, simplify = FALSE)
  quadStrings = character(length = length(combList))
  for (i in seq(along = combList)) if (i == 1) 
    quadStrings[i] = paste(combList[[i]][1], ":", combList[[i]][2])
  else quadStrings[i] = paste("+", combList[[i]][1], ":", combList[[i]][2])
  quadStrings = paste(quadStrings, collapse = "")
  if (DB) 
    print(quadStrings)
  if (identical(form, "linear")) {
    form = paste(r.c, "~", linStrings)
    if (DB) 
      print(form)
  }
  if (identical(form, "quadratic")) {
    form = paste(r.c, "~", linStrings, "+", quadStrings)
    if (DB) 
      print(form)
  }
  lm.1 = lm(formula = form, data = mdo)
  if (DB) 
    print(lm.1)
  dcList = vector(mode = "list", length = length(names(mdo)))
  names(dcList) = names(names(mdo))
  dcList[1:length(names(mdo))] = 0
  if (!is.null(factors)) {
    for (i in names(factors)) dcList[[i]] = factors[[i]][1]
  }
  if (DB) 
    print(dcList)
  help.predict = function(a, b, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = 2 * b/sqrt(3)
    dcList[[y.c]] = 1 - (2 * b/sqrt(3)) - (a - b/sqrt(3))
    dcList[[z.c]] = a - b/sqrt(3)
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
  }
  a = seq(0, 1, length = steps)
  b = seq(0, sqrt(3)/2, length = steps)
  mat = outer(a, b, help.predict, x.c, y.c, lm.1)
  acc = nrow(mat)
  v = seq(acc, 1, length = acc)
  w = c(seq(2, acc, length = acc/2), seq(acc, 2, length = acc/2))
  mat[outer(w, v, `+`) <= acc] = NA
  .mfc(mat, main = main, col = col, axes = FALSE, key.axes = axis(4))
  if (axes == TRUE) {
    segments(0.5, 0, 0.5, 1, col = col.axis)
    coox1 = rep(0.49, 11)
    coox2 = rep(0.51, 11)
    cooy = seq(0, 1, length = 11)
    for (i in 2:10) {
      segments(coox1[i], cooy[i], coox2[i], cooy[i], col = col.axis)
      text(coox2[i] + 0.01, cooy[i], labels = (i - 1)/10, cex = cex.axis, col = col.text)
    }
    segments(0, 0, 0.75, 0.5, col = col.axis)
    coox1 = seq(0.745, -0.005, length = 11)
    coox2 = seq(0.755, 0.005, length = 11)
    cooy1 = seq(0.51, 0.01, length = 11)
    cooy2 = seq(0.49, -0.01, length = 11)
    for (i in 2:10) {
      segments(coox1[i], cooy1[i], coox2[i], cooy2[i], col = col.axis)
      text(coox2[i], cooy1[i] + 0.01, labels = (i - 1)/10, cex = cex.axis, col = col.text)
    }
    segments(0.25, 0.5, 1, 0, col = col.axis)
    coox1 = seq(0.245, 0.995, length = 11)
    coox2 = seq(0.255, 1.005, length = 11)
    cooy1 = seq(0.49, -0.01, length = 11)
    cooy2 = seq(0.51, 0.01, length = 11)
    for (i in 2:10) {
      segments(coox1[i], cooy1[i], coox2[i], cooy2[i], col = col.axis)
      text(coox1[i] + 0.01, cooy2[i] + 0.01, labels = (i - 1)/10, cex = cex.axis, col = col.text)
    }
  }
  segments(-0.005, 0, 0.495, 1, lwd = 5, col = border)
  segments(0.505, 1, 1.005, 0, lwd = 5, col = border)
  segments(1, -0.005, 0, -0.005, lwd = 5, col = border)
  mtext(ylab, 1, at = -0.025, cex = 1.5)
  mtext(xlab, 3, at = 0.5, cex = 1.5, line = 0.1)
  mtext(zlab, 1, at = 1.025, cex = 1.5)
  invisible(mat)
}
wirePlot3 = function(x, y, z, response, data = NULL, main, xlab, ylab, zlab, form = "linear", phi, theta, col = 1, steps, factors) {
  DB = FALSE
  out = list()
  mdo = data
  x.c = deparse(substitute(x))
  y.c = deparse(substitute(y))
  z.c = deparse(substitute(z))
  r.c = deparse(substitute(response))
  if (missing(col)) 
    col = 1
  if (missing(main)) 
    main = paste("Response Surface for", r.c)
  if (missing(ylab)) 
    ylab = y.c
  if (missing(xlab)) 
    xlab = x.c
  if (missing(zlab)) 
    zlab = z.c
  if (missing(phi)) 
    phi = 30
  if (missing(theta)) 
    theta = 30
  if (missing(factors)) 
    factors = NULL
  if (missing(steps)) 
    steps = 100
  if (!is.function(col)) {
    if (identical(col, 1)) 
      col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if (identical(col, 2)) 
      col = colorRampPalette(c("blue", "white", "red"), space = "Lab")
    if (identical(col, 3)) 
      col = colorRampPalette(c("blue", "white", "orange"))
    if (identical(col, 4)) 
      col = colorRampPalette(c("gold", "white", "firebrick"))
  }
  phi = phi%%360
  .phi = phi
  theta = theta%%360
  .theta = theta
  nameVec = names(names(mdo))
  linStrings = "-1"
  for (i in seq(along = nameVec)) linStrings = paste(linStrings, "+", nameVec[i])
  if (DB) 
    print(linStrings)
  combList = combn(nameVec, 2, simplify = FALSE)
  quadStrings = character(length = length(combList))
  for (i in seq(along = combList)) if (i == 1) 
    quadStrings[i] = paste(combList[[i]][1], ":", combList[[i]][2])
  else quadStrings[i] = paste("+", combList[[i]][1], ":", combList[[i]][2])
  quadStrings = paste(quadStrings, collapse = "")
  if (DB) 
    print(quadStrings)
  if (identical(form, "linear")) {
    form = paste(r.c, "~", linStrings)
    if (DB) 
      print(form)
  }
  if (identical(form, "quadratic")) {
    form = paste(r.c, "~", linStrings, "+", quadStrings)
    if (DB) 
      print(form)
  }
  lm.1 = lm(formula = form, data = mdo)
  if (DB) 
    print(lm.1)
  dcList = vector(mode = "list", length = length(names(mdo)))
  names(dcList) = names(names(mdo))
  dcList[1:length(names(mdo))] = 0
  if (!is.null(factors)) {
    for (i in names(factors)) dcList[[i]] = factors[[i]][1]
  }
  if (DB) 
    print(dcList)
  help.predict = function(a, b, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = 2 * b/sqrt(3)
    dcList[[y.c]] = 1 - (2 * b/sqrt(3)) - (a - b/sqrt(3))
    dcList[[z.c]] = a - b/sqrt(3)
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
  }
  a = seq(0, 1, length = steps)
  b = seq(0, sqrt(3)/2, length = steps)
  mat = outer(a, b, help.predict, x.c, y.c, lm.1)
  acc = nrow(mat)
  sca = sin(1/3 * pi)
  ncmat = ncol(mat)
  v = seq(acc, 1, length = acc)
  w = c(seq(2, acc, length = acc/2), seq(acc, 2, length = acc/2))
  mat[outer(w, v, `+`) <= acc] = NA
  if (is.function(col)) {
    nrMat <- nrow(mat)
    ncMat <- ncol(mat)
    nbcol <- 100
    color <- col(nbcol)
    matFacet = mat[-1, -1] + mat[-1, -ncmat] + mat[-acc, -1] + mat[-acc, -ncmat]
    facetcol <- cut(matFacet, nbcol)
  }
  else {
    color = col
    facetcol = 1
  }
  maxim = max(mat, na.rm = TRUE) * acc
  minim = min(mat, na.rm = TRUE) * acc
  per = persp(x = seq(0, acc, length = acc), y = seq(0, acc * sca, length = ncmat), mat * acc, phi = .phi, theta = .theta, scale = TRUE, col = "transparent", 
              border = FALSE, box = FALSE, main = main, xlab = xlab, ylab = ylab)
  lineList = contourLines(x = seq(0, acc, length = acc), y = seq(0, acc * sca, length = ncmat), mat)
  for (i in seq(along = lineList)) lines(trans3d(lineList[[i]]$x, lineList[[i]]$y, z = minim, pmat = per))
  if (.phi < 90) {
    lines(trans3d(x = seq(0, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
    lines(trans3d(x = seq(acc, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
    lines(trans3d(x = 0:acc, y = 0, z = maxim, pmat = per), lty = 2)
  }
  if (.theta > 323 || .theta < 37) {
    lines(trans3d(x = acc/2, y = acc * sca, z = minim:maxim, pmat = per), lty = 2)
    lines(trans3d(x = 0, y = 0, z = minim:maxim, pmat = per), lty = 2)
    lines(trans3d(x = acc, y = 0, z = minim:maxim, pmat = per), lty = 2)
  }
  if (.theta > 37 && .theta < 156) 
    lines(trans3d(x = 0, y = 0, z = minim:maxim, pmat = per), lty = 2)
  if (.theta > 156 && .theta < 323) {
    lines(trans3d(x = acc, y = 0, z = minim:maxim, pmat = per), lty = 2)
  }
  lines(trans3d(x = seq(0, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = minim, pmat = per), lty = 1, lwd = 2)
  lines(trans3d(x = seq(acc, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = minim, pmat = per), lty = 1, lwd = 2)
  lines(trans3d(x = 0:acc, y = 0, z = minim, pmat = per), lty = 1, lwd = 2)
  text(trans3d(x = acc/2 + acc/50, y = acc * sca + acc * sca/50, z = minim, pmat = per), labels = xlab, lwd = 2)
  text(trans3d(x = -acc/50, y = -acc * sca/50, z = minim, pmat = per), labels = ylab, lwd = 2)
  text(trans3d(x = acc + acc/50, 0, z = minim, pmat = per), labels = zlab, cex = 1, lwd = 2)
  par(new = TRUE)
  persp(x = seq(0, acc, length = acc), y = seq(0, acc * sca, length = ncmat), mat * acc, phi = .phi, theta = .theta, scale = TRUE, col = color[facetcol], 
        border = FALSE, box = FALSE)
  if (.phi > 0) {
    lines(trans3d(x = seq(0, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
    lines(trans3d(x = seq(acc, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
    lines(trans3d(x = 0:acc, y = 0, z = maxim, pmat = per), lty = 2)
  }
  if (.theta > 37 && .theta < 156) {
    lines(trans3d(x = acc/2, y = acc * sca, z = minim:maxim, pmat = per), lty = 2)
    lines(trans3d(x = acc, y = 0, z = minim:maxim, pmat = per), lty = 2)
  }
  if (.theta > 156 && .theta < 323) {
    lines(trans3d(x = acc/2, y = acc * sca, z = minim:maxim, pmat = per), lty = 2)
    lines(trans3d(x = 0, y = 0, z = minim:maxim, pmat = per), lty = 2)
  }
  if (TRUE) {
    zlim = range(mat, finite = TRUE, na.rm = TRUE)
    leglevel = pretty(zlim, 6)
    legcol = col(length(leglevel))
    legpretty = as.character(abs(leglevel))
    temp = character(length(leglevel))
    temp[leglevel > 0] = "+"
    temp[leglevel < 0] = "-"
    temp[leglevel == 0] = " "
    legpretty = paste(temp, legpretty, sep = "")
    if (.theta <= 180) 
      legend("topright", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
    if (.theta > 180) 
      legend("topleft", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
  }
  invisible(mat)
} 

#  -------------- cur_t.r -------------- 
.curvTest = function(fdo, response, DB = FALSE) {
  fullString = character(0)
  interString = character(0)
  quadString = character(0)
  quadString2 = character(0)
  y = character(0)
  nameVec = names(names(fdo))
  pValue = NA
  if (!(try(is.character(response), silent = TRUE) == TRUE)) {
    y = deparse(substitute(response))
  }
  else if (response %in% names(response(fdo))) 
    y = response
  if (length(nameVec) < 1) {
    cat("\n")
    invisible("curvTest: not enough names (factors)!")
  }
  if (nrow(centerCube(fdo)) <= 1) {
    cat("\n")
    invisible("curvTest: not enough centerPoints for a test for curvature")
  }
  if (nrow(star(fdo)) > 0) {
    cat("\n")
    invisible("curvTest: star portion exists; nothing to do")
  }
  if (length(nameVec) == 1) {
    cat("\n")
    invisible(nameVec[1])
  }
  for (i in seq(along = nameVec)) {
    if (i == 1) {
      fullString = nameVec[i]
      if (DB) 
        print(fullString)
    }
    if (length(nameVec) >= 2) {
      fullString = paste(fullString, "+", nameVec[i + 1])
      if (DB) 
        print(fullString)
    }
    if ((i + 1) >= length(nameVec)) {
      if (DB) 
        print("break")
      break
    }
  }
  for (k in 2:length(nameVec)) {
    if (DB) 
      print(k)
    temp = combn(nameVec, k, simplify = TRUE)
    if (DB) 
      print(temp)
    for (i in 1:ncol(temp)) {
      interString = character(0)
      for (j in 1:k) {
        if (j == 1) 
          interString = temp[j, i]
        else interString = paste(interString, ":", temp[j, i])
      }
      fullString = paste(fullString, "+", interString)
    }
  }
  for (i in seq(along = nameVec)) {
    if (i == 1) {
      quadString = paste("I(", nameVec[i], "^2)", sep = "")
      quadString2 = paste("I(", nameVec[i], "^2)", sep = "")
      if (DB) 
        print(quadString)
    }
    if (length(nameVec) >= 2) {
      quadString = paste(quadString, "+I(", nameVec[i + 1], "^2)", sep = "")
      quadString2 = c(quadString2, paste("I(", nameVec[i + 1], "^2)", sep = ""))
      if (DB) {
        print(quadString2)
        print(quadString)
      }
    }
    if ((i + 1) >= length(nameVec)) {
      if (DB) 
        print("break")
      break
    }
  }
  fullString = paste(fullString, "+", quadString)
  fullString = paste(y, "~", fullString)
  if (DB) 
    print(fullString)
  aov.1 = aov(formula = as.formula(fullString), data = fdo)
  if (DB) {
    print(summary(aov.1))
    print(pmatch(quadString2, row.names(summary(aov.1)[[1]])))
  }
  rows = pmatch(quadString2, row.names(summary(aov.1)[[1]]))
  if (DB) 
    print(rows)
  rows = row.names(summary(aov.1)[[1]])[rows[!is.na(rows)]]
  if (DB) 
    print(rows)
  cols = "Pr(>F)"
  tempFrame = data.frame(summary(aov.1)[[1]][rows, cols])
  if (nrow(tempFrame) > 0) {
    row.names(tempFrame) = rows
    names(tempFrame) = cols
    pValue = format(tempFrame[1, 1], digits = 3)
  }
  out = paste("Test for Curvature:  p =", pValue)
  cat("\n")
  cat(out)
  cat("\n")
  invisible(pValue)
}
summaryFits = function(fdo, lmFit = TRUE, curvTest = TRUE, origFit = TRUE) {
  summaryList = vector(mode = "list", length = 0)
  origFrame = as.data.frame(fdo)
  for (i in names(names(fdo))) origFrame[, i] = code2real(lows(fdo)[[i]], highs(fdo)[[i]], origFrame[, i])
  for (f in names(response(fdo))) {
    if (!is.null(fits(fdo)[[f]])) {
      cat(paste("----------- Summary for response '", f, "' -----------", sep = ""))
      cat("\n")
      print(summary(fits(fdo)[[f]]))
      cat("-----------")
      cat("\n")
      cat("\n")
      cat("Regression in non coded form:")
      cat("\n")
      lm.f = (lm(formula(fits(fdo)[[f]]), data = origFrame))
      coefs = coefficients(lm.f)
      coefsI = coefs[pmatch("(Intercept)", names(coefs))]
      coefsW = coefs[-pmatch("(Intercept)", names(coefs))]
      coefsW = coefsW[!is.na(coefsW)]
      temp = character(length(coefsW))
      temp[coefsW >= 0] = "+"
      temp[coefsW < 0] = "-"
      firstString = ""
      firstString = paste(firstString, format(coefsI, digits = 4))
      restString = paste(format(abs(coefsW), digits = 4), names(coefsW), sep = "*")
      restString = paste(temp, restString)
      restString = paste(restString, collapse = " ")
      fullString = paste(firstString, restString)
      fullString = paste(paste(f, " ="), fullString)
      cat("\n")
      cat(paste("  ", fullString))
      cat("\n")
      cat("\n")
      cat("-----------")
      cat("\n")
      .curvTest(fdo, f)
      cat("\n")
      cat("\n")
    }
  }
  invisible()
} 

#  -------------- dot_t.r -------------- 
dotPlot = function(x, group, xlim, ylim, col, xlab, ylab, pch, cex, breaks, stacked = TRUE, ...) {
  DB = FALSE
  pch.size = "O"
  grouped = TRUE
  parList = list(...)
  if (missing(xlab)) 
    xlab = deparse(substitute(x))
  x = x[!is.na(x)]
  if (missing(xlim)) 
    xlim = range(x)
  if (missing(ylim)) 
    ylim = c(0, 1)
  if (missing(ylab)) 
    ylab = ""
  if (missing(cex)) 
    cex = 1
  if (missing(group)) 
    group = rep(1, length(x))
  if (length(unique(group)) == 1) 
    grouped = FALSE
  if (missing(pch) || length(unique(group)) > length(pch)) 
    pch = 1:length(unique(group))
  if (missing(col) || length(unique(group)) > length(col)) 
    col = 1:length(unique(group))
  if (missing(breaks)) {
    plot(1, 1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, cex = cex, xlab = xlab, ylab = ylab)
    slotSizeX = strwidth(pch.size, units = "user", cex = cex)
    if (DB) 
      print(paste("slotSizeX:", slotSizeX))
    span = diff(range(x))
    temp1 = ppoints(2 * ceiling(span/slotSizeX))
    temp2 = numeric(length(temp1) + 2)
    temp2[2:(length(temp1) + 1)] = temp1
    temp2[1] = temp2[1] - 1.01 * diff(c(temp1[1], temp1[2]))
    temp2[length(temp2)] = rev(temp1)[1] + 1.01 * diff(c(temp1[1], temp1[2]))
    temp2 = temp2 * span + min(x)
    temp = min(x) + ppoints(span/slotSizeX) * span
    breaks = numeric(length(temp) + 2)
    breaks[2:(length(temp) + 1)] = temp
    breaks[1] = temp[1] - diff(c(temp[1], temp[2])) * 1.001
    breaks[length(breaks)] = rev(temp)[1] + diff(c(temp[1], temp[2])) * 1.001
    breaks = temp2
  }
  slotSizeY = strheight(pch.size, units = "user", cex = cex)
  if (DB) 
    print(paste("slotSizeY:", slotSizeY))
  span = diff(ylim)
  temp1 = ppoints(2 * ceiling(span/slotSizeY))
  temp2 = numeric(length(temp1) + 2)
  temp2[2:(length(temp1) + 1)] = temp1
  temp2[1] = temp2[1] - 1.01 * diff(c(temp1[1], temp1[2]))
  temp2[length(temp2)] = rev(temp1)[1] + 1.01 * diff(c(temp1[1], temp1[2]))
  yVec = temp2 * span + min(ylim)
  if (yVec[1] < 0) 
    yVec = yVec + abs(yVec[1])
  else yVec = yVec - yVec[1]
  if (DB) 
    print(paste("temp2:", temp2))
  if (DB) 
    print(paste("breaks:", breaks))
  histObj = hist(x, breaks = breaks, right = FALSE, plot = FALSE)
  hMids = histObj$mids
  hCounts = histObj$counts
  hMids = histObj$mids
  mat = matrix(NA, nrow = length(x), ncol = length(hMids))
  colMat = mat
  groupmat = mat
  numVec = 1:nrow(mat)
  cutOff = 1
  groupList = vector(mode = "list", length = length(unique(group)))
  for (k in unique(group)) {
    histObj = hist(x[group == k], breaks = breaks, plot = FALSE)
    hMids = histObj$mids
    hCounts = histObj$counts
    hMids = histObj$mids
    for (i in seq(along = hMids)) {
      value = pch[k]
      colValue = col[k]
      from = 0
      from = numVec[is.na(mat[, i])][1]
      to = from
      if (hCounts[i] == 0) 
        value = NA
      if (hCounts[i] >= 1) 
        to = to + hCounts[i] - 1
      if (to > cutOff) 
        cutOff = to
      if (DB) {
        print(paste("from:", from))
        print(paste("to:", to))
        print(paste("i:", i))
        print(paste("value:", value))
      }
      mat[from:to, i] = value
      colMat[from:to, i] = colValue
    }
    groupList[[k]] = groupmat
  }
  if (grouped && !stacked) {
    groupIndex = unique(group)
    par(mfrow = c(length(groupIndex), 1))
    for (i in groupIndex) dotPlot(x[group == i], xlim = xlim, breaks = breaks, cex = cex, xlab = xlab, ylab = ylab, col = col, pch = pch, ...)
  }
  else {
    mat = mat[1:cutOff, ]
    if (!is.matrix(mat)) 
      mat = matrix(mat, nrow = 1)
    if (DB) 
      print(mat)
    plot(1, 1, xlim = xlim, ylim = ylim, type = "n", cex = cex, xlab = xlab, ylab = ylab, ...)
    for (i in 1:nrow(mat)) {
      x = hMids[!is.na(mat[i, ])]
      y = rep(i * 0.3, times = length(x))
      y = rep(yVec[i], times = length(x))
      col = colMat[i, !is.na(mat[i, ])]
      pch = mat[i, !is.na(mat[i, ])]
      points(x, y, col = col, pch = pch, cex = cex)
    }
  }
  if (DB) 
    print(hMids)
  invisible(mat)
} 

#  -------------- min_s.r -------------- 
fracChoose = function() {
  DB = FALSE
  genList = list(6 * 9)
  genList = list(c("C = AB"), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c("D = ABC"), c("D = AB", "E = AC"), c("D = AB", 
                                                                                                                                                      "E = AC", "F = BC"), c("D = AB", "E = AC", "F = BC", "G = ABC"), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c("E = ABCD"), c("E = ABC", "F = BCD"), 
                 c("E = ABC", "F = BCD", "G = ACD"), c("E = BCD", "F = ACD", "G = ABC", "H = ABD"), c("E = ABC", "F = BCD", "G = ACD", "H = ABD", "J = ABCD"), c("E = ABC", 
                                                                                                                                                                 "F = BCD", "G = ACD", "H = ABD", "J = ABCD", "K = AB"), c("E = ABC", "F = BCD", "G = ACD", "H = ABD", "J = ABCD", "K = AB", "L = AC"), c(NULL), 
                 c(NULL), c(NULL), c("F = ABCDE"), c("F = ABCD", "G = ABDE"), c("F = ABC", "G = ABD", "H = BCDE"), c("F = BCDE", "G = ACDE", "H = ABDE", "J = ABCE"), 
                 c("F = ABCD", "G = ABCE", "H = ABDE", "J = ACDE", "K = BCDE"), c("F = ABC", "G = BCD", "H = CDE", "J = ACD", "K = AEF", "L = ADEF"), c(NULL), c(NULL), 
                 c(NULL), c(NULL), c("G = ABCDEF"), c("G = ABCD", "H = ABEF"), c("G = ABCD", "H = ACEF", "J = CDEF"), c("G = BCDF", "H = ACDF", "J = ABDE", "K = ABCE"), 
                 c("G = CDE", "H = ABCD", "J = ABF", "K = BDEF", "L = ADEF"), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c("H = ABCDEFG"), c("H = ACDFG", "J = BCEFG"), 
                 c("H = ABCG", "J = BCDE", "K = ACDF"), c("H = ABCG", "J = BCDE", "K = ACDF", "L = ABCDEFG"))
  resList = list(6 * 9)
  resList = list(c(3), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(4), c(3), c(3), c(3), c(NULL), c(NULL), c(NULL), 
                 c(NULL), c(NULL), c(NULL), c(5), c(4), c(4), c(4), c(3), c(3), c(3), c(NULL), c(NULL), c(NULL), c(6), c(4), c(4), c(4), c(4), c(4), c(NULL), c(NULL), 
                 c(NULL), c(NULL), c(7), c(5), c(4), c(4), c(4), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(8), c(6), c(5), c(5))
  facMat = matrix(rep(3:11, 6), ncol = 9, byrow = TRUE)
  runMat = matrix(c(rep(2^2, 9), rep(2^3, 9), rep(2^4, 9), rep(2^5, 9), rep(2^6, 9), rep(2^7, 9)), ncol = 9, byrow = TRUE)
  par(mfrow = c(6, 9))
  par(mar = c(0, 0, 0, 0))
  par(oma = c(4, 4, 4, 4))
  colList = vector(mode = "list")
  colList[3] = "red"
  colList[4] = "yellow"
  colList[5] = "green"
  colList[6] = "green"
  colList[7] = "green"
  colList[8] = "green"
  k = 3
  N = 2^2
  m = 0
  for (i in seq(along = genList)) {
    res = unlist(resList[[i]])
    plot(0, 0, xaxs = "i", yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "", bg = "red", fg = "green")
    box()
    if (!is.null(res)) 
      rect(0, 0, 1, 1, col = colList[[res]])
    yPos = 0.04
    xPos = 0.6
    item = rev(genList[[i]])
    for (j in seq(along = item)) {
      text(xPos, yPos + (j - 1) * 0.125, item[j], adj = c(0, 0), cex = 0.8)
    }
    p = log2((2^k)/N)
    if (!is.null(res)) {
      romNum = as.character(as.roman(res))
      text(0.1, 0.9, do.call("expression", list(substitute(2[romNum]^(k - p), list(k = k, p = p, romNum = romNum)))), adj = c(0, 1), cex = 1.5)
    }
    k = k + 1
    if ((i%%9) == 0) {
      N = N * 2
      k = 3
    }
  }
  cat("\nChoose a fractional factorial design by clicking into the appropriate field")
  cat("\nWaiting for your selection:")
  cat("\n\n")
  flush.console()
  mtext("number of runs N", side = 2, line = 2.5, outer = TRUE)
  mtext("number of variables k", side = 3, line = 2.5, outer = TRUE)
  for (numFac in 1:9) {
    mtext(numFac + 2, at = c(-0.05 + (1/9) * numFac), outer = TRUE, line = 0.5)
  }
  for (k in 1:6) {
    mtext(2^(k + 1), at = (7 - k)/6 - (0.5 * (1/6)), side = 2, outer = TRUE, line = 0.5)
  }
  if (DB) 
    cat("TODO: standardize the locator part in respect to possible figure region")
  xyList = NULL
  xyList = try(locator(1), silent = TRUE)
  x = 1
  y = 1
  if (!is.null(xyList)) {
    x = ceiling(xyList$x + 8)
    y = ceiling(6 - xyList$y)
  }
  mat = matrix(1:54, ncol = 9, byrow = TRUE)
  fdo = NULL
  if (!(x %in% 1:ncol(mat)) || !(y %in% 1:nrow(mat))) 
    return(fracDesign(k = 3, gen = NULL, replicates = 1))
  else index = mat[y, x]
  k = facMat[y, x]
  generator = genList[[index]]
  N = runMat[y, x]
  if (!is.null(generator)) {
    fdo = try(do.call("fracDesign", list(k = k, gen = generator)), silent = TRUE)
  }
  if (N >= 2^k & is.null(generator)) {
    replicates = N/(2^k)
    fdo = try(fracDesign(k = k, gen = NULL, replicates = replicates), silent = TRUE)
  }
  if (class(fdo) == "facDesign") 
    return(fdo)
  else return(genList[[mat[y, x]]])
} 

#  -------------- des_e.r -------------- 
setClass(Class = "desirability", representation = representation(response = "character", low = "numeric", high = "numeric", target = "ANY", scale = "numeric", 
                                                                 importance = "numeric"))
desirability = function(response, low, high, target = "max", scale = c(1, 1), importance = 1, constraints) {
  if (low >= high) 
    stop("the lower bound must be greater than the high bound!")
  if (any(scale <= 0)) 
    stop("the scale parameter must be greater than zero!")
  if (!is.numeric(target) & !identical(tolower(target), "min") & !identical(tolower(target), "max")) 
    stop("target needs to be \"min\", \"max\" or a numeric value")
  return(new("desirability", response = deparse(substitute(response)), low = low, high = high, target = target, scale = scale, importance = importance))
}
.desireFun = function(low, high, target = "max", scale = c(1, 1), importance = 1) {
  DB = FALSE
  if (importance > 10 | importance < 0.1) 
    stop("importance needs to be in [0.1, 10]")
  if (low >= high) 
    stop("the lower bound must be greater than the high bound!")
  if (any(scale <= 0)) 
    stop("the scale parameter must be greater than zero!")
  if (is.numeric(target)) {
    out = function(y) {
      if (DB) 
        print("target")
      flush.console()
      d = rep(0, length(y))
      d[y >= low & y <= target] = ((y[y >= low & y <= target] - low)/(target - low))^scale[1]
      d[y >= target & y <= high] = ((y[y >= target & y <= high] - high)/(target - high))^scale[2]
      return(d^importance)
    }
    return(out)
  }
  if (identical(tolower(target), "min")) {
    out = function(y) {
      if (DB) 
        print("min")
      d = rep(0, length(y))
      d[y > high] = 0
      d[y < low] = 1
      d[y >= low & y <= high] = ((y[y >= low & y <= high] - high)/(low - high))^scale[1]
      return(d^importance)
    }
    return(out)
  }
  if (identical(tolower(target), "max")) {
    out = function(y) {
      if (DB) 
        print("max")
      d = rep(0, length(y))
      d[y < low] = 0
      d[y > high] = 1
      d[y >= low & y <= high] = ((y[y >= low & y <= high] - low)/(high - low))^scale[1]
      return(d^importance)
    }
    return(out)
  }
}
setMethod("show", signature(object = "desirability"), function(object) {
  if (!is.numeric(object@target)) 
    cat("Target is to", paste(object@target, "imize", sep = ""), object@response, "\n")
  else cat("Target is ", object@target, " for", object@response, "\n")
  cat("lower Bound: ", object@low, "\n")
  cat("higher Bound: ", object@high, "\n")
  if (is.numeric(object@target)) 
    cat("Scale factor is: low =", object@scale[1], "and high =", object@scale[2], "\n")
  else if (identical("min", object@target) | identical("max", object@target)) 
    cat("Scale factor is: ", object@scale, "\n")
  cat("importance: ", object@importance, "\n")
  cat("\n")
})
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "desirability"), function(x, y, scale, main, xlab, ylab, type, col, numPoints = 500, ...) {
  xm1 = NULL
  xm2 = NULL
  ym = NULL
  y = NULL
  if (missing(main)) 
    main = paste("Desirability function for", x@response)
  if (missing(xlab)) 
    xlab = x@response
  if (missing(ylab)) 
    ylab = "Desirability"
  if (missing(type)) 
    type = "l"
  if (missing(scale)) 
    scale = x@scale
  if (missing(col)) 
    col = 1:length(scale)
  dFun = .desireFun(x@low, x@high, x@target, x@scale, x@importance)
  xVals = seq(x@low - 0.04 * diff(range(x@low, x@high)), x@high + 0.04 * diff(range(x@low, x@high)), length = numPoints)
  yVals = dFun(xVals)
  plot(xVals, yVals, main = main, xlab = xlab, ylab = ylab, type = type, col = col, ...)
  if (is.numeric(x@target)) {
    xm1 = mean(c(par("usr")[1], x@target))
    xm2 = mean(c(par("usr")[2], x@target))
    ym1 = yVals[max((1:numPoints)[xVals <= xm1])]
    ym2 = yVals[max((1:numPoints)[xVals <= xm2])]
    text(xm1 + 0.025 * diff(range(par("usr")[1:2])), ym1, paste("scale =", scale[1]), adj = c(0, 0))
    text(xm2 - 0.025 * diff(range(par("usr")[1:2])), ym2, paste("scale =", scale[2]), adj = c(1, 1))
  }
  else {
    xm1 = mean(par("usr")[c(1, 2)])
    ym1 = yVals[max((1:numPoints)[xVals <= xm1])]
    if (identical(x@target, "max")) 
      text(xm1 + 0.025 * diff(range(par("usr")[1:2])), ym1 - 0.025 * diff(range(par("usr")[3:4])), paste("scale =", scale[1]), adj = c(0, 0))
    else text(xm1 + 0.025 * diff(range(par("usr")[1:2])), ym1 + 0.025 * diff(range(par("usr")[3:4])), paste("scale =", scale[1]), adj = c(0, 1))
  }
  out = list(x = xVals, y = yVals)
  names(out) = c(x@response, "desirability")
  invisible(out)
})
overall = function(fdo, steps = 20, constraints, ...) {
  DB = FALSE
  importances = list()
  cs = list()
  if (!missing(constraints)) 
    cs = constraints
  l = vector(mode = "list", length = 0)
  fitList = fits(fdo)
  if (length(fitList) < 1) 
    stop(paste("no fits found in fits(", deparse(substitute(fdo)), ")"), sep = "")
  desList = desires(fdo)
  if (length(desList) < 1) 
    stop(paste("no desirabilities found in desires(", deparse(substitute(fdo)), ")"), sep = "")
  X = cube(fdo)
  newdata = NULL
  for (i in names(names(fdo))) {
    seqList = vector(mode = "list", length = 0)
    seqList[["length"]] = steps
    seqList[["from"]] = min(X[, i])
    seqList[["to"]] = max(X[, i])
    minC = NULL
    maxC = NULL
    if (!is.null(cs[[i]])) {
      if (length(cs[[i]]) < 2) 
        stop("length of ", names(cs[i]), "=", cs[i], " < 2 in constraints")
      minC = min(cs[[i]])
      if (!is.null(minC) & !is.na(minC)) 
        seqList[["from"]] = minC
      maxC = max(cs[[i]])
      if (!is.null(maxC) & !is.na(maxC)) 
        seqList[["to"]] = maxC
      if (maxC == minC) 
        stop(paste("equal values in constraints ", names(cs[i]), "=", cs[i]))
    }
    l[[i]] = do.call(seq, seqList)
  }
  if (DB) 
    print(l)
  newdata = expand.grid(l)
  names(newdata) = names(X)
  out = newdata
  yCharSet = intersect(names(desires(fdo)), names(fits(fdo)))
  dFrame = data.frame(matrix(NA, nrow = nrow(newdata), ncol = length(yCharSet) + 1))
  names(dFrame) = c(yCharSet, "overall")
  dFrame[, "overall"] = 1
  for (y in yCharSet) {
    obj = desList[[y]]
    dFun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
    lm.y = fitList[[y]]
    importances[[y]] = desires(fdo)[[y]]@importance
    yHat = predict(lm.y, newdata = newdata, ...)
    yDes = dFun(yHat)
    dFrame[, y] = yDes
    if (DB) {
      print(y)
      print(dFun)
      print(lm.y)
      print(dFrame)
    }
  }
  geomFac = 1/sum(unlist(importances))
  overall = apply(dFrame, 1, prod)^geomFac
  dFrame[, "overall"] = overall
  dFrame = cbind(out, dFrame)
  invisible(dFrame)
}
.desHelp = function(fdo, factors, ...) {
  if (length(factors) != length(names(fdo))) 
    stop("not enough factors specified in factors")
  if (any(is.na(factors))) 
    stop("factors contain NA")
  yCharSet = intersect(names(desires(fdo)), names(fits(fdo)))
  desList = desires(fdo)
  fitList = fits(fdo)
  yDes = vector(mode = "list")
  for (y in yCharSet) {
    obj = desList[[y]]
    dFun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
    lm.y = fitList[[y]]
    yHat = predict(lm.y, newdata = data.frame(factors), ...)
    yDes[[y]] = dFun(yHat)
  }
  return(yDes)
}
setClass(Class = "desOpt", representation = representation(facCoded = "list", facReal = "list", responses = "list", desirabilities = "list", overall = "numeric", 
                                                           all = "data.frame", fdo = "facDesign"))
as.data.frame.desOpt = function(x, row.names = NULL, optional = FALSE, ...) {
  return(x@all)
}
setMethod("as.data.frame", "desOpt", function(x, row.names = NULL, optional = FALSE, ...) {
  return(x@all)
})
.validizeConstraints = function(fdo, constraints) {
  X = as.data.frame(fdo)
  csOut = vector(mode = "list")
  for (i in names(names(fdo))) {
    csOut[[i]] = c(min(X[, i]), max(X[, i]))
  }
  if (missing(constraints)) 
    return(csOut)
  cs2 = constraints[names(names(fdo))]
  cs2 = cs2[!unlist(lapply(cs2, is.null))]
  cs2 = cs2[(unlist(lapply(cs2, length)) == 2)]
  csOut[names(cs2)] = cs2[names(cs2)]
  return(csOut)
}
.dHelp = function(model, dFun) {
  lm1 = model
  d1 = dFun
  out = function(newdata) {
    return(d1(predict(lm1, newdata = newdata)))
  }
  return(out)
}
optimum = function(fdo, constraints, steps = 25, type = "grid", start, ...) {
  DB = FALSE
  if (missing(fdo)) 
    stop("missing fdo!")
  X = as.data.frame(fdo)
  numFac = length(names(fdo))
  if (!(type %in% c("grid", "optim", "gosolnp"))) {
    warning(paste("type =", deparse(substitute(type)), "not found --> using type = \"grid\""))
    type = "grid"
  }
  constraints = .validizeConstraints(fdo, constraints)
  if (missing(start)) 
    start = as.numeric(lapply(constraints, mean))
  lower = numeric(length(constraints))
  upper = numeric(length(constraints))
  for (i in seq(along = constraints)) {
    lower[i] = min(constraints[[i]])
    upper[i] = max(constraints[[i]])
  }
  if (DB) {
    print(constraints)
    print(start)
  }
  desOpt = new("desOpt")
  desOpt@fdo = fdo
  facCoded = NA
  desirabilities = NA
  overall = NA
  setList = list()
  dList = list()
  importances = list()
  yCharSet = intersect(names(desires(fdo)), names(fits(fdo)))
  for (y in yCharSet) {
    obj = desires(fdo)[[y]]
    dFun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
    lm.y = fits(fdo)[[y]]
    importances[[y]] = desires(fdo)[[y]]@importance
    dList[[y]] = .dHelp(lm.y, dFun)
  }
  geomFac = 1/sum(unlist(importances))
  dAll = function(X) {
    newdata = data.frame(t(X))
    names(newdata) = LETTERS[1:ncol(newdata)]
    return(prod(unlist(lapply(dList, do.call, list(newdata = newdata))))^geomFac)
  }
  dAllRsolnp = function(X) {
    newdata = data.frame(t(X))
    names(newdata) = LETTERS[1:ncol(newdata)]
    return(-prod(unlist(lapply(dList, do.call, list(newdata = newdata))))^geomFac)
  }
  if (type == "optim") {
    #        print(lower)
    #       print(upper)
    temp = optim(par = start, dAll, method = "L-BFGS-B", lower = lower, upper = upper, control = list(fnscale = -1, maxit = 1000))
    facCoded = as.list(temp$par)
    names(facCoded) = names(names(fdo))
    desOpt@facCoded = facCoded
    overall = temp$value
    desirabilities = .desHelp(fdo, desOpt@facCoded)
  }
  if (type == "gosolnp") {
    #if (!require(Rsolnp, quietly = TRUE)) 
    #    stop("Package Rsolnp needs to be installed!")
    temp = Rsolnp::gosolnp(fun = dAllRsolnp, LB = lower, UB = upper)
    facCoded = as.list(temp$pars)
    names(facCoded) = names(names(fdo))
    desOpt@facCoded = facCoded
    overall = -rev(temp$values)[1]
    desirabilities = .desHelp(fdo, desOpt@facCoded)
  }
  if (type == "grid") {
    dVals = overall(fdo = fdo, constraints = constraints, steps = steps)
    index = order(dVals[, "overall"], decreasing = TRUE)[1]
    desOpt@all = dVals
    desOpt@facCoded = as.list(dVals[index, names(names(fdo))])
    desirabilities = as.list(dVals[index, names(response(fdo))])
    names(desirabilities) = names(response(fdo)) #fix for the case of having just one response 
    overall = dVals[index, "overall"]
  }
  for (i in names(desOpt@facCoded)) {
    desOpt@facReal[[i]] = code2real(lows(fdo)[[i]], highs(fdo)[[i]], desOpt@facCoded[[i]])
  }
  desOpt@desirabilities = desirabilities
  desOpt@overall = overall
  newData = do.call(data.frame, desOpt@facCoded)
  for (i in names(desOpt@desirabilities)) {
    desOpt@responses[[i]] = predict(fits(fdo)[[i]], newData)
  }
  return(desOpt)
}
setMethod("show", signature(object = "desOpt"), function(object) {
  cat(paste("\ncomposite (overall) desirability:", format(object@overall, digits = 3)))
  cat("\n")
  cat("\n")
  temp1 = do.call(data.frame, object@facCoded)
  temp2 = do.call(data.frame, object@facReal)
  facFrame = rbind(temp1, temp2)
  row.names(facFrame) = c("coded", "real")
  show(format(facFrame, digits = 3))
  temp1 = do.call(data.frame, object@responses)
  temp2 = do.call(data.frame, object@desirabilities)
  respDesFrame = rbind(temp1, temp2)
  row.names(respDesFrame) = c("Responses", "Desirabilities")
  cat("\n")
  show(format(respDesFrame, digits = 3))
}) 

#  -------------- par_t.r -------------- 
paretoPlot = function(fdo, threeWay = FALSE, abs = TRUE, decreasing = TRUE, na.last = NA, alpha = 0.05, response = NULL, xlim, ylim, xlab, ylab, main, single = TRUE, ...) {  ###
  DB = FALSE
  if(single==FALSE)                                                           ###
    par(mfrow=.splitDev(length(response(fdo)))[[2]])                           ###
  if(is.null(response)==FALSE)                                                ###
  {                                                                           ###
    temp=response(fdo)[response]                                               ###
    response(fdo)=temp                                                         ###
  }                                                                           ###
  ylimMissing = FALSE
  if (missing(ylim)) 
    ylimMissing = TRUE
  if (missing(xlab)) 
    xlab = ""
  location = "topright"
  if (decreasing == F) {
    location = "topleft"
  }
  xVals = numeric(0)
  sig.neg = NULL
  sig.pos = NULL
  effect.list = vector("list")
  for (j in 1:ncol(response(fdo))) {
    par(mar = c(5.1, 4.1, 4.1, 4.1))
    if (j > 1 && single==TRUE) {
      dev.new()
      par(mar = c(5.1, 4.1, 4.1, 4.1))
    }
    if (!any(is.na(response(fdo)[, j]))) {
      if (missing(ylab)) 
        ylabel = names(response(fdo))[j]                                ###
      else                                                                ###
        ylabel = ylab                                                   ###
      form = paste("response(fdo)[,", j, "]~")
      for (i in 1:ncol(cube(fdo))) {
        form = paste(form, names(cube(fdo))[i], sep = "")
        if (i < ncol(cube(fdo))) 
          form = paste(form, "*", sep = "")
      }
      if (DB == TRUE) 
        print(form)
      lm.1 = lm(as.formula(form), data = as.data.frame(fdo))
      coefs = coef(lm.1)[-pmatch("(Intercept)", names(coef(lm.1)))]
      df.resid = df.residual(lm.1)
      num.c = nrow(centerCube(fdo))
      if (df.resid == 0) {
        effect = 2 * coefs
        effect = effect[!is.na(effect)]
        effect.list[[j]] = effect
        if (missing(main)) 
          main = "Lenth Plot of effects"
        plt = TRUE
        limits = TRUE
        faclab = NULL
        m = length(effect)
        d = m/3
        s0 = 1.5 * median(abs(effect))
        rmedian = effect[abs(effect) < 2.5 * s0]
        PSE = 1.5 * median(abs(rmedian))
        ME = qt(1 - alpha/2, d) * PSE
        Gamma = (1 + (1 - alpha)^(1/m))/2
        SME = qt(Gamma, d) * PSE
        n = length(effect)
        if (ylimMissing) 
          if (abs) 
            ylim <- (range(c(0, abs(effect), 1.3 * ME))) * 1.1
        else ylim <- (range(c(effect, -1.3 * ME, 1.3 * ME))) * 1.1
        if (abs) {
          xVals = barplot(abs(effect), las = 2, main = main, xlab = xlab, ylim = ylim, ylab = ylabel, ...)
          abline(h = ME, col = "red")
          abline(h = SME, col = "red")
          try(axis(4, at = ME, labels = round(ME, 3), las = 2), silent = T)
          text(x = xVals[1], y = ME, "ME", pos = 3)
          try(axis(4, at = SME, labels = round(SME, 3), las = 2), silent = T)
          text(x = xVals[2], y = SME, "SME", pos = 3)
        }
        else {
          xVals = barplot(effect, las = 2, main = main, xlab = xlab, ylim = ylim, ylab = ylabel, ...)
          abline(h = c(-ME, ME), col = "red")
          abline(h = c(-SME, SME), col = "red")
          try(axis(4, at = c(-ME, ME), labels = round(c(-ME, ME), 3), las = 2), silent = T)
          text(x = xVals[1], y = c(-ME, ME), "ME", pos = c(1, 3))
          try(axis(4, at = c(-SME, SME), labels = round(c(-SME, SME), 3), las = 2), silent = T)
          text(x = xVals[2], y = c(-SME, SME), "SME", pos = c(1, 3))
        }
        if (length(xVals) >= 1) 
          for (i in 1:length(xVals)) {
            text(xVals[i], effect[i] + max(ylim) * sign(effect[i]) * 0.05, format(round(effect[i], 3)))
          }
        if (DB) 
          print(paste("MSE:", ME, "SME:", SME))
      }
      else {
        if (missing(main)) 
          main = "Standardized main effects and interactions"
        effect = ((summary(lm.1)$coefficients[-pmatch("(Intercept)", names(coef(lm.1))), 1])/(summary(lm.1)$coefficients[-pmatch("(Intercept)", names(coef(lm.1))), 
                                                                                                                         2]))
        if (all(is.na(effect))) 
          stop("effects could not be calculated")
        effect = effect[!is.na(effect)]
        effect.list[[j]] = effect
        if ((df.resid) > 0) {
          sig.pos = -qt(alpha/2, df.resid)
          sig.neg = +qt(alpha/2, df.resid)
        }
        if (ylimMissing) 
          if (abs) {
            tempVec = c(effect, sig.pos)
            tempVec = tempVec[!is.na(tempVec)]
            ylim = c(0, 1.3 * max(tempVec))
          }
        else {
          tempVec1 = c(0, effect, sig.neg, sig.pos)
          tempVec1 = tempVec1[!is.na(tempVec1)]
          tempVec2 = c(abs(effect), sig.pos, sig.neg)
          tempVec2 = tempVec2[!is.na(tempVec2)]
          ylim = c(1.3 * min(tempVec1), 1.3 * max(tempVec2))
        }
        if (DB) 
          print(paste("ylim:", ylim))
        effect = effect[order(abs(effect), na.last = TRUE, decreasing = decreasing)]
        effect = round(effect, 3)
        if (abs) {
          xVals = barplot(abs(effect), las = 2, main = main, xlab = xlab, ylim = ylim, ylab = ylabel, ...)
          if (length(xVals) >= 1) 
            for (i in 1:length(xVals)) {
              text(xVals[i], abs(effect[i] + max(ylim) * sign(effect[i]) * 0.05), format(effect[i]))
            }
        }
        else {
          xVals = barplot(effect, las = 2, main = main, xlab = xlab, ylim = ylim, ylab = ylabel, ...)
          if (length(xVals) >= 1) 
            for (i in 1:length(xVals)) {
              text(xVals[i], effect[i] + max(ylim) * sign(effect[i]) * 0.05, format(effect[i]))
            }
        }
        myDelta = diff(range(ylim)) * 0.02
        try(abline(h = sig.pos, col = "red"), silent = TRUE)
        try(axis(4, at = c(sig.pos, sig.neg), labels = round(c(sig.pos, sig.neg), 3), las = 2), silent = T)
        try(abline(h = sig.neg, col = "red"), silent = TRUE)
      }
      legend(location, legend = names(fdo), pch = paste(names(names(fdo)), sep = ""), bg = "white", inset = 0.02)
      abline(h = 0)
      box()
    }
    if (DB) {
      print(df.resid)
      print(num.c)
      print(effect)
    }
  }
  invisible(effect.list)
}
class(paretoPlot) <- "invisible" 

#  -------------- Blo_n.r -------------- 
.letToNum = function(charVec) {
  charVec = toupper(charVec)
  numVec = numeric(length = length(charVec))
  if (!is.character(charVec)) 
    stop(".letToNum: characterVector needs to be provided!")
  alphabet = LETTERS[1:26]
  for (i in seq(along = charVec)) {
    numVec[i] = match(charVec[i], alphabet)
  }
  return(numVec)
}
.letterIndex = function(char) {
  if (char %in% LETTERS[1:26]) {
    return((1:26)[LETTERS[1:26] == char])
  }
  stop("no valid LETTER specified!")
}
.isEven = function(x) {
  if (x%%2 > 0) 
    return(FALSE)
  return(TRUE)
}
.isOdd = function(x) {
  return(!.isEven(x))
}
.lociv = function(charVec) {
  lenVec = numeric(length = length(charVec))
  for (i in seq(along = charVec)) {
    lenVec[i] = length(strsplit(charVec[i], split = "")[[1]])
  }
  return(lenVec)
}
.confoundings = function(blockGenVec, lSet, DB = FALSE) {
  biVec = character(0)
  for (i in 2:length(blockGenVec)) {
    mat = combn(blockGenVec, i)
    temp = apply(mat, 2, strsplit, split = "")
    comb = lapply(temp, unlist)
    comb = lapply(comb, c, lSet)
    if (DB) {
      print("here")
      print(comb)
    }
    combFreq = sapply(comb, table)%%2
    combBool = !apply(combFreq, 2, as.logical)
    chars = row.names(combFreq)
    if (DB) 
      print(combBool)
    biTemp = character(0)
    for (j in 1:ncol(combBool)) {
      biTemp = c(biTemp, paste(chars[combBool[, j]], collapse = ""))
    }
    if (DB) 
      print(biTemp)
    biVec = c(biVec, biTemp)
  }
  return(c(blockGenVec, biVec))
}
.rsm = vector(mode = "list", length = 7)
.rsm[[1]] = list(k = 3, blocks = 2, gen = c("ABC"))
.rsm[[2]] = list(k = 3, blocks = 4, gen = c("AB", "AC"))
.rsm[[3]] = list(k = 4, blocks = 2, gen = c("ABCD"))
.rsm[[4]] = list(k = 4, blocks = 4, gen = c("ABC", "ACD"))
.rsm[[5]] = list(k = 4, blocks = 8, gen = c("AB", "BC", "CD"))
.rsm[[6]] = list(k = 5, blocks = 2, gen = c("ABCDE"))
.rsm[[7]] = list(k = 5, blocks = 4, gen = c("ABC", "CDE"))
.rsm[[8]] = list(k = 5, blocks = 8, gen = c("ABE", "BCE", "CDE"))
.rsm[[9]] = list(k = 5, blocks = 16, gen = c("AB", "AC", "CD", "DE"))
.rsm[[10]] = list(k = 6, blocks = 2, gen = c("ABCDEF"))
.rsm[[11]] = list(k = 6, blocks = 4, gen = c("ABCF", "CDEF"))
.rsm[[12]] = list(k = 6, blocks = 8, gen = c("ABEF", "ABCD", "ACE"))
.rsm[[13]] = list(k = 6, blocks = 16, gen = c("ABF", "ACF", "BDF", "DEF"))
.rsm[[14]] = list(k = 6, blocks = 32, gen = c("AB", "BC", "CD", "DE", "EF"))
.rsm[[15]] = list(k = 7, blocks = 2, gen = c("ABCDEFG"))
.rsm[[16]] = list(k = 7, blocks = 4, gen = c("ABCFG", "CDEFG"))
.rsm[[17]] = list(k = 7, blocks = 8, gen = c("ABC", "DEF", "AFG"))
.rsm[[18]] = list(k = 7, blocks = 16, gen = c("ABD", "EFG", "CDE", "ADG"))
.rsm[[19]] = list(k = 7, blocks = 32, gen = c("ABG", "BCG", "CDG", "DEG", "EFG"))
.rsm[[20]] = list(k = 7, blocks = 64, gen = c("AB", "BC", "CD", "DE", "EF", "FG"))
.blockInteractions = function(fdo, blocks = 2, useTable = "rsm") {
  DB = FALSE
  if (!(blocks %in% c(0, 1, 2, 4, 8, 16, 32, 64))) 
    stop("blocks needs to be a power of 2 up to 64!")
  gen = NULL
  if (blocks %in% c(0, 1)) {
    if (DB) 
      print("TODO: Return the Identity as generator")
    return(gen)
  }
  if (length(useTable) > 0) {
    if (!(nrow(unique(cube(fdo))) >= 2^.numFac(fdo))) 
      stop("no blocking of a fractional factorial Design --> block on replicates instead!")
    if (identical(useTable, "rsm")) {
      for (i in seq(along = .rsm)) {
        if (.rsm[[i]]$k == .numFac(fdo) & .rsm[[i]]$blocks == blocks) 
          return(.rsm[[i]]$gen)
      }
    }
    return(gen)
  }
  bgaci = matrix(nrow = 0, ncol = blocks - 1)
  if (!is.numeric(blocks)) 
    stop("blocks must be an integer")
  numCol = log2(blocks)
  blockGen = character(3)
  lSet = names(names(fdo))
  sSet = vector(mode = "list")
  for (i in length(lSet):2) {
    sSet = c(sSet, combn(lSet, i, simplify = FALSE))
  }
  if (blocks == 2) {
    index = order(sapply(sSet, length), decreasing = TRUE)[1]
    sSet = sapply(sSet, paste, collapse = "")
    return(sSet[index])
  }
  sSet = sapply(sSet, paste, collapse = "")
  if (DB) 
    print(sSet)
  possGen = combn(sSet, numCol, simplify = FALSE)
  for (i in seq(along = possGen)) {
    blockGenVec = unlist(possGen[[i]])
    if (DB) 
      print(blockGenVec)
    if (DB) 
      print(.confoundings(blockGenVec, lSet))
    newRow = .confoundings(blockGenVec, lSet)
    if (!any(newRow %in% c(lSet, ""))) 
      bgaci = rbind(bgaci, .confoundings(blockGenVec, lSet))
  }
  mat = unique(t(apply(bgaci, 1, sort)))
  temp = t(apply(mat, 1, .lociv))
  temp = t(apply(temp, 1, sort))
  ref = temp[1, ]
  index = 1
  for (i in 1:nrow(temp)) {
    if (any((ref - temp[i, ]) < 0)) {
      ref = temp[i, ]
      index = i
    }
  }
  for (i in 1:nrow(temp)) {
    if (!(any(ref - temp[i, ] > 0) | any(ref - temp[i, ] < 0))) {
      index = c(index, i)
    }
  }
  temp = unique((mat[index, ]))
  cat("\nSuggested Effects for Blocking:")
  cat("\n")
  cat(temp[1, 1:numCol])
  cat("\n")
  cat("\nInteractions Confounded with blocks:")
  cat("\n")
  cat(unique(temp[1, ]))
  cat("\n")
  cat("\n Alternate Effects for Blocking:")
  cat(temp[c(-1), 1:numCol])
  cat("\n")
  gen = temp[1, 1:numCol]
  return(gen)
}

#Reihenfolge muss noch auf Standard gesetzt werden
.blockGenCol = function(gen, fdo) {
  DB = FALSE
  blockVec = NULL
  .blockCol = NULL
  genList = gen
  genList = strsplit(genList, split = "")
  .fdo = fdo
  for (i in seq(along = genList)) {
    gen = genList[[i]]
    for (j in seq(along = gen)) {
      genTemp = .fdo[, gen[j]]
      if (j == 1) 
        blockVec = rep(1, length = length(genTemp))
      blockVec = blockVec * genTemp
      if (DB) 
        print(blockVec)
    }
    if (i == 1) 
      .blockCol = data.frame(B1 = blockVec)
    else .blockCol = cbind(.blockCol, blockVec)
  }
  names(.blockCol) = paste("B", 1:ncol(.blockCol), sep = "")
  return(.blockCol)
}
.blockCol = function(.blockGenCol) {
  DB = FALSE
  .blockCol = numeric(nrow(.blockGenCol))
  uniCol = unique(.blockGenCol)
  for (i in 1:nrow(uniCol)) {
    if (ncol(uniCol) == 1) 
      .blockCol[apply(t(as.data.frame(apply(.blockGenCol, 1, "==", uniCol[i, ]))), 2, all)] = i
    else .blockCol[apply(apply(.blockGenCol, 1, "==", uniCol[i, ]), 2, all)] = i
  }
  return(data.frame(Block = .blockCol))
}
randomize = function(fdo, random.seed, so = FALSE) {
  if (missing(random.seed)) 
    set.seed(93275938)
  else set.seed(random.seed)
  j = 1
  temp = runOrd(fdo)
  for (i in sort(unique(block(fdo)[, 1]))) {
    pos = !is.na(match(block(fdo)[, 1], i))
    count = sum(as.numeric(pos))
    if (so) {
      temp[pos, 1] = j:(j + (count - 1))
    }
    else {
      temp[pos, 1] = sample(j:(j + (count - 1)), count)
    }
    j = j + count
  }
  runOrd(fdo) = temp
  return(fdo)
}
blocking = function(fdo, blocks, BoR = FALSE, random.seed, useTable = "rsm", gen) {
  override = FALSE
  Block = data.frame(Block = rep(1, nrow(fdo))) #do not change
  block(fdo) = Block  #do not change
  fdo = randomize(fdo, so = TRUE)
  if (missing(random.seed)) {
    runif(1)
    random.seed = .Random.seed[sample(1:626, 1)]
  }
  if (missing(gen)) 
    gen = NULL
  if (blocks <= 1) {
    Block = data.frame(Block = rep(1, nrow(fdo)))
    block(fdo) = Block
    fdo = randomize(fdo, random.seed = random.seed)
    return(fdo)
  }
  if (nrow(star(fdo)) > 0 | nrow(centerStar(fdo)) > 0) {
    if (blocks == 2) {
      override = TRUE
      fdo = randomize(fdo, so = TRUE)
      numB1 = nrow(cube(fdo)) + nrow(centerCube(fdo))
      numB2 = nrow(fdo) - numB1
      block(fdo) = data.frame(Block = c(rep(1, numB1), rep(2, numB2)))
      #or by using standard order always
      
      blockGen(fdo) = data.frame(B1 = rep(NA, nrow(fdo)))
    }
    if (blocks %in% c(2, 3, 5, 9, 17)) 
      blocks = blocks - 1
    else stop("Blocking not possible")
  }
  else {
    if (!(blocks %in% c(1, 2, 4, 8, 16, 32, 64, 128))) 
      stop("Blocking not possible")
  }
  if (is.null(gen)) 
    gen = .blockInteractions(fdo, blocks, useTable)
  if (is.null(gen) & !override) {
    cat("\n")
    cat(paste("Blocking in", blocks, "blocks not possible!"))
    cat("\n")
    return(fdo)
  }
  if (!override) {
    .blockGenCol = .blockGenCol(gen, fdo)
    .blockCol = .blockCol(.blockGenCol)
    Block = .blockCol#[runOrd(fdo)[,1],] #TODO: fix this in block()
    BlockGenCol = .blockGenCol#[runOrd(fdo)[,1],] #TODO: fix this in block()
    #or by using standard order always
    block(fdo) = Block
    blockGen(fdo) = BlockGenCol
    #       block(fdo) = .blockCol
    #        blockGen(fdo) = .blockGenCol
  }
  numCC = nrow(centerCube(fdo))
  if (numCC > 0) {
    ccFrame = as.data.frame(matrix(0, nrow = numCC, ncol = ncol(cube(fdo))))
    names(ccFrame) = names(names(fdo))
    centerCube(fdo) = ccFrame
  }
  fdo = randomize(fdo, random.seed = random.seed)
  return(fdo)
} 

#  -------------- eff_t.r -------------- 
.m.interaction.plot = function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = TRUE, trace.label = deparse(substitute(trace.factor)), 
                               fixed = FALSE, xlab = deparse(substitute(x.factor)), ylab = ylabel, ylim = range(cells, na.rm = TRUE), lty = nc:1, col = 1, pch = c(1L:9, 0, letters), xpd = NULL, 
                               leg.bg = par("bg"), leg.bty = "n", xtick = FALSE, xaxt = par("xaxt"), axes = TRUE, ...) {
  ylabel <- paste(deparse(substitute(fun)), "of ", deparse(substitute(response)))
  type <- match.arg(type)
  cells <- tapply(response, list(x.factor, trace.factor), fun)
  nr <- nrow(cells)
  nc <- ncol(cells)
  xvals <- 1L:nr
  xvals = as.numeric(rownames(cells))
  if (is.ordered(x.factor)) {
    wn <- getOption("warn")
    options(warn = -1)
    xnm <- as.numeric(levels(x.factor))
    options(warn = wn)
    if (!any(is.na(xnm))) 
      xvals <- xnm
  }
  xlabs <- rownames(cells)
  ylabs <- colnames(cells)
  nch <- max(sapply(ylabs, nchar, type = "width"))
  if (is.null(xlabs)) 
    xlabs <- as.character(xvals)
  if (is.null(ylabs)) 
    ylabs <- as.character(1L:nc)
  xlim <- range(xvals)
  xleg <- xlim[2L] + 0.05 * diff(xlim)
  xlim <- xlim + c(-0.2/nr, if (legend) 0.2 + 0.02 * nch else 0.2/nr) * diff(xlim)
  matplot(xvals, cells, ..., type = type, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, xaxt = "n", col = col, lty = lty, pch = pch)
  if (axes && xaxt != "n") {
    axisInt <- function(x, main, sub, lwd, bg, log, asp, ...) axis(1, x, ...)
    mgp. <- par("mgp")
    if (!xtick) 
      mgp.[2L] <- 0
    axisInt(1, at = xvals, labels = xlabs, tick = xtick, mgp = mgp., xaxt = xaxt, ...)
  }
  if (legend) {
    yrng <- diff(ylim)
    yleg <- ylim[2L] - 0.1 * yrng
    if (!is.null(xpd) || {
      xpd. <- par("xpd")
      !is.na(xpd.) && !xpd. && (xpd <- TRUE)
    }) {
      op <- par(xpd = xpd)
      on.exit(par(op))
    }
    text(xleg, ylim[2L] - 0.05 * yrng, paste("  ", trace.label), adj = 0)
    if (!fixed) {
      ord <- sort.list(cells[nr, ], decreasing = TRUE)
      ylabs <- ylabs[ord]
      lty <- lty[1 + (ord - 1)%%length(lty)]
      col <- col[1 + (ord - 1)%%length(col)]
      pch <- pch[ord]
    }
    legend(xleg, yleg, legend = ylabs, col = col, pch = if (type %in% c("p", "b")) 
      pch, lty = if (type %in% c("l", "b")) 
        lty, bty = leg.bty, bg = leg.bg)
  }
  invisible(xvals)
}
setGeneric("effectPlot", def = function(object, factors, fun = mean, response = NULL, single = FALSE, points = FALSE, classic = FALSE, axes = TRUE, lty, xlab, ylab, ###
                                        main, ylim, ...) standardGeneric("effectPlot"))                             
setMethod(effectPlot, signature(object = "facDesign"), function(object, factors, fun = mean, response = NULL, single = FALSE, points = FALSE, classic = FALSE, axes = TRUE, ###
                                                                lty, xlab, ylab, main, ylim, ...) {
  oldMar = par("mar")
  oldOma = par("oma")
  oldMfrow = par("mfrow")
  oldMfcol = par("mfcol")
  on.exit(par(mar = oldMar, oma = oldOma, mfrow = oldMfrow, mfcol = oldMfcol))
  if(is.null(response)==FALSE)                                                ###
  {                                                                           ###
    temp=response(object)[response]                                            ###
    response(object)=temp                                                      ###
  }                                                                           ###
  ylabmiss = FALSE
  xlabmiss = FALSE
  mainmiss = FALSE
  ylimmiss = FALSE
  if (missing(ylim)) 
    ylimmiss = TRUE
  if (missing(lty)) 
    lty = 1
  X = cube(object)
  Y = as.data.frame(object@response[1:nrow(X), ])
  names(Y) = names(response(object))
  if (!missing(factors)) 
    k = length(factors)
  else #(missing(factors))                                                    ###
  {
    k = ncol(X)
    factors = names(X)
  }
  numCol = 1
  numRow = 1
  if (!single && missing(factors)) {                                          ###
    if (ncol(X) == 2) {
      numCol = 2
      numRow = 1
    }
    if (ncol(X) > 2) {
      numCol = 2
      numRow = 2
    }
  }
  if (!single && !missing(factors)) {                                         ###
    if (length(factors) == 2) {                                             ###
      numCol = 2                                                          ###
      numRow = 1                                                          ###
    }                                                                       ###
    if (length(factors) == 3) {                                             ###
      numCol = 3                                                          ###
      numRow = 1                                                          ###
    }                                                                       ###
    if (length(factors) == 4) {                                             ###
      numCol = 2                                                          ###
      numRow = 2                                                          ###
    }                                                                       ###
    if (length(factors) == 5) {                                             ###
      numCol = 3                                                          ###
      numRow = 2                                                          ###
    }                                                                       ###
    if (length(factors) == 6) {                                             ###
      numCol = 3                                                          ###
      numRow = 2                                                          ###
    }                                                                       ###
    if (length(factors) > 6) {                                              ###
      numRow = ceiling(sqrt(length(factors)))                             ###
      numCol = ceiling(sqrt(length(factors)))                             ###
    }                                                                       ###
  }                                                                           ###
  if (classic) {
    numCol = ncol(X)
    numRow = 1
  }
  if (!single) 
    par(mfrow = c(numRow, numCol))
  nextResponse = FALSE
  for (j in 1:ncol(Y)) {
    counter = 0
    cells = numeric(0)
    for (i in 1:length(factors)) {
      cells = c(cells, as.vector(tapply(Y[, j], list(X[, factors[i]], rep(0, nrow(X))), fun)))
      if (points) 
        cells = range(Y)
    }
    if (nextResponse & !single) {
      dev.new()
      par(mfrow = c(numRow, numCol))
    }
    for (i in 1:length(factors)) {
      if ((counter != 0 & counter%%(numCol * numRow) == 0) & !single) {
        dev.new()
        par(mfrow = c(numRow, numCol))
      }
      if (missing(main)) {
        main = paste("Effect Plot for", names(Y)[j])
        mainmiss = TRUE
      }
      if (mainmiss) 
        main = paste("Effect Plot for", names(Y)[j])
      if (missing(xlab)) {
        xlab = factors[i]
        xlabmiss = TRUE
      }
      if (xlabmiss) {
        if (identical(" ", names(object)[[i]])) 
          xlab = factors[i]
        else xlab = paste(factors[i], ": ", names(object)[[i]], sep = "")
      }
      if (missing(ylab)) {
        ylab = paste(deparse(substitute(fun)), "of ", names(Y)[j])
        ylabmiss = TRUE
      }
      if (ylabmiss) 
        ylab = paste(deparse(substitute(fun)), "of ", names(Y)[j])
      if (ylimmiss) 
        ylim = range(cells, na.rm = TRUE)
      if (classic & i == 1) {
        par(mar = c(5, 0, 0, 0) + 0.1)
        par(oma = c(-0.1, 4, 4, 1) + 0.1)
      }
      if (classic) {
        .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun, 
                            ylab = ylab, legend = FALSE, axes = FALSE, main = " ", ...)
        grid(NA, 2)
        axis(1, at = X[, factors[i]])
        if (i == 1) 
          axis(2)
        box()
        title(main, outer = TRUE)
      }
      else {
        .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun, 
                            ylab = ylab, legend = FALSE, axes = axes, main = main, ...)
        grid(NA, 2)
      }
      if (points) 
        points(X[, factors[i]], Y[, j], ...)
      counter = counter + 1
    }
    nextResponse = TRUE
  }
})
setMethod(effectPlot, signature(object = "taguchiDesign"), function(object, factors, fun = mean, response = NULL, single = FALSE, points = FALSE, classic = FALSE,  ###
                                                                    axes = TRUE, lty, xlab, ylab, main, ylim, ...) {
  oldMar = par("mar")
  oldOma = par("oma")
  oldMfrow = par("mfrow")
  oldMfcol = par("mfcol")
  on.exit(par(mar = oldMar, oma = oldOma, mfrow = oldMfrow, mfcol = oldMfcol))
  if(is.null(response)==FALSE)                                                ###
  {                                                                           ###
    temp=response(object)[response]                                            ###
    response(object)=temp                                                      ###
  }                                                                           ###
  ylabmiss = FALSE
  xlabmiss = FALSE
  mainmiss = FALSE
  ylimmiss = FALSE
  if (missing(ylim)) 
    ylimmiss = TRUE
  if (missing(lty)) 
    lty = 1
  X = object@design
  Y = response(object)
  if (!missing(factors)) 
    k = length(factors)
  else #(missing(factors))                                                    ###
  {
    k = ncol(X)
    factors = names(X)
  }
  numCol = 1
  numRow = 1
  if (!single && missing(factors)) {                                          ###
    if (ncol(X) == 2) {
      numCol = 2
      numRow = 1
    }
    if (ncol(X) > 2) {
      numCol = 2
      numRow = 2
    }
  }
  if (!single && !missing(factors)) {                                         ###
    if (length(factors) == 2) {                                             ###
      numCol = 2                                                          ###
      numRow = 1                                                          ###
    }                                                                       ###
    if (length(factors) == 3) {                                             ###
      numCol = 3                                                          ###
      numRow = 1                                                          ###
    }                                                                       ###
    if (length(factors) == 4) {                                             ###
      numCol = 2                                                          ###
      numRow = 2                                                          ###
    }                                                                       ###
    if (length(factors) == 5) {                                             ###
      numCol = 3                                                          ###
      numRow = 2                                                          ###
    }                                                                       ###
    if (length(factors) == 6) {                                             ###
      numCol = 3                                                          ###
      numRow = 2                                                          ###
    }                                                                       ###
    if (length(factors) > 6) {                                              ###
      numRow = ceiling(sqrt(length(factors)))                             ###
      numCol = ceiling(sqrt(length(factors)))                             ###
    }                                                                       ###
  }                                                                           ###        
  if (classic) {
    numCol = ncol(X)
    numRow = 1
  }
  if (!single) 
    par(mfrow = c(numRow, numCol))
  nextResponse = FALSE
  for (j in 1:ncol(Y)) {
    counter = 0
    cells = numeric(0)
    for (i in 1:length(factors)) {
      cells = c(cells, as.vector(tapply(Y[, j], list(X[, factors[i]], rep(0, nrow(X))), fun)))
      if (points) 
        cells = range(Y)
    }
    if (nextResponse & !single) {
      dev.new()
      par(mfrow = c(numRow, numCol))
    }
    for (i in 1:length(factors)) {
      if ((counter != 0 & counter%%(numCol * numRow) == 0) & !single) {
        dev.new()
        par(mfrow = c(numRow, numCol))
      }
      if (missing(main)) {
        main = paste("Effect Plot for", names(Y)[j])
        mainmiss = TRUE
      }
      if (mainmiss) 
        main = paste("Effect Plot for", names(Y)[j])
      if (missing(xlab)) {
        xlab = factors[i]
        xlabmiss = TRUE
      }
      if (xlabmiss) {
        if (identical(" ", names(object)[[i]])) 
          xlab = factors[i]
        else xlab = paste(factors[i], ": ", names(object)[[i]], sep = "")
      }
      if (missing(ylab)) {
        ylab = paste(deparse(substitute(fun)), "of ", names(Y)[j])
        ylabmiss = TRUE
      }
      if (ylabmiss) 
        ylab = paste(deparse(substitute(fun)), "of ", names(Y)[j])
      if (ylimmiss) 
        ylim = range(cells, na.rm = TRUE)
      if (classic & i == 1) {
        par(mar = c(5, 0, 0, 0) + 0.1)
        par(oma = c(-0.1, 4, 4, 1) + 0.1)
      }
      if (classic) {
        .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun, 
                            ylab = ylab, legend = FALSE, axes = FALSE, main = " ", ...)
        grid(NA, 2)
        axis(1, at = X[, factors[i]])
        if (i == 1) 
          axis(2)
        box()
        title(main, outer = TRUE)
      }
      else {
        .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun, 
                            ylab = ylab, legend = FALSE, axes = axes, main = main, ...)
        grid(NA, 2)
      }
      if (points) 
        points(X[, factors[i]], Y[, j], ...)
      counter = counter + 1
    }
    nextResponse = TRUE
  }
}) 

snPlot=function(object, type="nominal" , factors, fun = mean, response = NULL, 
                single = FALSE, points = FALSE, classic = FALSE, axes = TRUE, 
                lty, xlab, ylab, main, ylim, ...)
{ 
  Debugging=TRUE
  if(class(object)!="taguchiDesign")
    stop("object needs to be of class taguchiDesign") 
  Length=dim(as.data.frame(object))[1]
  resLength=dim(response(object))[2]
  temp=data.frame(attributes(object)$design)
  comp=unique(temp)
  SNi=numeric();SN=data.frame()
  m=numeric();y=numeric()
  if(missing(main))
  {
    for(k in 1:resLength)
      m[k]=paste("Effect Plot for S/N ratios of",names(response(object))[k])
    main=m
  }
  if(missing(ylab))
  {
    for(k in 1:resLength)
      y[k]=paste("means of S/N ratios for ",names(response(object))[k])
    ylab=y
  }
  if(identical(comp,temp))
    stop("taguchi design has no replicates! S/N can not be calculated!")
  for(k in 1:resLength)
  { 
    for(j in 1:dim(comp)[1])
    {  val=numeric()                                                            
    for(i in 1:Length)
    {
      if(identical(as.numeric(comp[j,]),as.numeric(temp[i,])))
      {
        val[i]=response(object)[i,k]
      }
      else
        val[i]=NA
    }
    n=Length/(dim(comp)[1])
    if(type=="nominal")
      SNi[j]=10*log10((mean(val,na.rm=TRUE)^2)/(sd(val,na.rm=TRUE)^2))
    if(type=="smaller")
      SNi[j]=-10*log10((1/n)*sum(val^2,na.rm=TRUE)) 
    if(type=="larger")
      SNi[j]=-10*log10((1/n)*sum(1/(val^2),na.rm=TRUE))
    if(Debugging==TRUE)
      print(SNi)
    for(i in 1:Length)
    {
      if(identical(as.numeric(comp[j,]),as.numeric(temp[i,])))
      {
        SN[i,k]=SNi[j]
      }
    }
    }
    tdo=object
    response(tdo)=SN[k]
    if(k>1)
      dev.new()
    effectPlot(object = tdo, factors=factors, fun = mean, response = response, 
               single = single, points = points, classic = classic, axes = axes, 
               lty = lty, xlab = xlab, ylab =ylab[k], main = main[k], ylim = ylim, ...)
  }
  names(SN)=paste("S/N",names(response(object)))
  invisible(SN)
}

#  -------------- sim_c.r -------------- 
.normalDensity2d <- function(mu1 = 150, mu2 = 165, sigma1 = 1, sigma2 = 1, rho = 0, x = seq(from = -1, to = +1, length = 10), y = seq(from = -1, 
                                                                                                                                      to = +1, length = 10)) {
  z = matrix(nrow = length(x), ncol = length(y))
  zx = 0
  zy = 0
  for (x1 in x) {
    zx = zx + 1
    for (x2 in y) {
      zy = zy + 1
      z[zx, zy] = 1/(2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) * exp(-1/(2 * (1 - rho^2)) * (((x1 - mu1)/sigma1)^2 - 2 * rho * (x1 - mu1)/sigma1 * (x2 - 
                                                                                                                                                        mu2)/sigma2 + ((x2 - mu2)/sigma2)^2))
    }
    zy = 0
  }
  return(z)
}
simProc = function(x1, x2, x3, noise = TRUE) {
  .norm2d <- function(x1, x2, mu1 = 160, mu2 = 165, rho = 0.7, sigma1 = 45, sigma2 = 22.5) {
    z = 1/(2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) * exp(-1/(2 * (1 - rho^2)) * (((x1 - mu1)/sigma1)^2 - 2 * rho * (x1 - mu1)/sigma1 * (x2 - mu2)/sigma2 + 
                                                                                       ((x2 - mu2)/sigma2)^2))
    return(z)
  }
  max_z = 0.0002200907
  min_z = 8.358082e-10
  yield = .norm2d(x1 = x1, x2 = x2)
  yield = yield - min_z
  yield = (yield/max_z) * 0.9
  if (noise) 
    yield = yield + rnorm(length(yield), mean = 0, sd = 0.007)
  return(yield)
}
print.invisible <- function(x, ...) {
  cat("")
}
class(simProc) = "invisible" 

# -------------- hel_s.r -------------- 
.charToDistFunc = function(distribution, type = "q") {                                                           ####   .CHARTODISTFUNC-FUNCTION
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

# -------------- par_l.r -------------- 
paretoChart = function(x, weight, showTable = TRUE, las = 0, main, col, border, xlab, ylab = "Frequency", percentVec, ...) {
  varName = deparse(substitute(x))[1]
  corp.col = "pink3"
  corp.border = "red3"
  if (!is.vector(x) & !is.data.frame(x) & !is.table(x)) 
    stop("x should be a vector, dataframe or a table")
  if (is.table(x)) {
    xtable = x
  }
  if (is.vector(x)) {
    if (!is.null(names(x))) 
      xtable = as.table(x)
    else xtable = table(x)
  }
  if (!missing(weight)) {
    if (!is.numeric(weight)) 
      stop("weight must be numeric!")
    if (is.null(names(weight))) 
      stop("weight is missing names for matching!")
    else {
      if (FALSE %in% (sort(names(weight)) == sort(names(xtable)))) 
        stop("names of weight and table do not match!")
      else {
        for (i in 1:length(xtable)) {
          xtable[i] = weight[names(weight) == names(xtable)[i]] * xtable[i]
        }
      }
    }
  }
  else {
    weight = FALSE
  }
  if (missing(showTable)) 
    showTable = TRUE
  if (missing(xlab)) 
    xlab = ""
  if (missing(main)) 
    main = paste("Pareto Chart for ", varName)
  if (missing(col)) 
    col = corp.col
  if (missing(border)) 
    border = corp.border
  if (missing(percentVec)) 
    percentVec = seq(0, 1, by = 0.25)
  call <- match.call(expand.dots = TRUE)
  if (length(xtable) > 1) {
    ylim = c(min(xtable), max(xtable) * 1.025)
    xtable = c(sort(xtable, decreasing = TRUE, na.last = TRUE))
    cumFreq = cumsum(xtable)
    sumFreq = sum(xtable)
    percentage = xtable/sum(xtable) * 100
    cumPerc = cumFreq/sumFreq * 100
    frameOut = data.frame(Frequency = as.numeric(xtable), Cum.Freq = cumFreq, Percentage = percentage, Cum.Perc = cumPerc)
    names(frameOut) = c(ylab, paste("Cum.", ylab), "Percentage", "Cum. Percentage")
    row.names(frameOut) = names(xtable)
    frameInt = as.data.frame(t(frameOut))
    names(frameInt) = rep(" ", dim(frameInt)[2])
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    lineHeight = par("csi") * par("lheight") * par("cex")
    tablespace = (2 * 4) * (lineHeight/par("fin")[2])
    plot.new()
    if (las == 0 | las == 1) {
      mymai = par("mai")
      mymai[1] = max(strheight(names(xtable), units = "inches")) * 3
      mymai[4] = strheight("Cumulative Percentage", units = "inches") * 8
      mymai[2] = max(strwidth(names(frameOut), units = "inches")) * 1.2
      par(mai = mymai, new = TRUE)
    }
    if (las == 2 | las == 3) {
      mymai = par("mai")
      mymai[1] = max(strwidth(names(xtable), units = "inches")) * 1.4
      mymai[4] = strheight("Cumulative Percentage", units = "inches") * 8
      mymai[2] = max(strwidth(names(frameOut), units = "inches")) * 1.2
      par(mai = mymai, new = TRUE)
    }
    if (showTable) {
      par(fig = c(0, 1, tablespace, 1))
      xValue = barplot(xtable, axes = FALSE, las = las, width = 1, space = 0.2, xlim = c(0.2, 1.2 * length(xtable)), main = main, ylim = c(0, sum(xtable) + 
                                                                                                                                             0.01 * (sum(xtable))), ylab = ylab, xlab = xlab, col = col, border = border)
      axis(1, at = xValue, labels = names(xtable), las = las)
      axis(2)
      axis(4, at = percentVec * (sumFreq), labels = percentVec)
      mtext(4, text = "Cumulative Percentage", line = 3)
      lines(xValue, cumFreq, col = corp.border)
      points(xValue, cumFreq, col = corp.border, pch = 15)
      par(fig = c(0, 1, 0, tablespace), new = TRUE)
      mymai[1] = 0
      mymai[3] = 0
      par(mai = mymai)
      plot(xValue, rep(1, length(xValue)), xlim = c(0.2, 1.2 * length(xtable)), ylim = c(0, 5), axes = FALSE, ylab = "", type = "n")
      axis(2, pos = 0.2, at = 1:4, labels = rev(c(ylab, paste("Cum.", ylab), "Percentage", "Cum. Percentage")), tick = FALSE, las = 1)
      numCol = dim(frameInt)[2]
      numRow = dim(frameInt)[1]
      for (i in 1:numCol) {
        for (j in 1:numRow) {
          text(xValue[i], numRow + 1 - j, round(frameInt[j, i]), adj = c(1, 0.5))
        }
      }
    }
    else {
      mymai[2] = mymai[4]
      par(mai = mymai)
      xValue = barplot(xtable, axes = FALSE, las = las, width = 1, space = 0.2, xlim = c(0.2, 1.2 * length(xtable)), main = main, ylim = c(0, sum(xtable) + 
                                                                                                                                             0.01 * (sum(xtable))), ylab = ylab, xlab = xlab, col = col, border = border)
      axis(1, at = xValue, labels = names(xtable), las = las)
      axis(2)
      axis(4, at = percentVec * (sumFreq), labels = percentVec)
      mtext(4, text = "Cumulative Percentage", line = 3)
      lines(xValue, cumFreq, col = corp.border)
      points(xValue, cumFreq, col = corp.border, pch = 15)
    }
  }
  else {
    warning("data should have at least two categories!")
  }
  frameOut = frameInt
  for (i in 3:nrow(frameInt)) {
    frameInt[i, ] = sprintf("%.1f%%", frameInt[i, ])
  }
  #    cat(paste("Pareto Analysis for", varName, "\n"))
  #    cat("---\n")
  print(format(frameInt, digits = 3))
  cat("\n")
  frameOut
} 

# -------------- int_t.r -------------- 
.letterPos = function(LETTER) {
  if (!(nchar(LETTER) == 1)) 
    stop("factor names should be single characters only")
  return((1:26)[LETTERS[1:26] == LETTER])
}
.testFun = function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = TRUE, trace.label = deparse(substitute(trace.factor)), 
                    fixed = FALSE, xlab = deparse(substitute(x.factor)), ylab = ylabel, ylim = range(cellNew, na.rm = TRUE), lty = nc:1, col = 1, pch = c(1L:9, 0, letters), 
                    xpd = NULL, leg.bg = par("bg"), leg.bty = "o", xtick = FALSE, xaxt = par("xaxt"), axes = TRUE, title = "", ...) {
  ylabel <- paste(deparse(substitute(fun)), "of ", deparse(substitute(response)))
  type <- match.arg(type)
  cellNew <- tapply(response, list(x.factor, trace.factor), fun)
  nr <- nrow(cellNew)
  nc <- ncol(cellNew)
  xvals <- 1L:nr
  if (is.ordered(x.factor)) {
    wn <- getOption("warn")
    options(warn = -1)
    xnm <- as.numeric(levels(x.factor))
    options(warn = wn)
    if (!any(is.na(xnm))) 
      xvals <- xnm
  }
  xlabs <- rownames(cellNew)
  ylabs <- colnames(cellNew)
  nch <- max(sapply(ylabs, nchar, type = "width"))
  if (is.null(xlabs)) 
    xlabs <- as.character(xvals)
  if (is.null(ylabs)) 
    ylabs <- as.character(1L:nc)
  xlim <- range(xvals)
  xleg <- xlim[2L] + 0.05 * diff(xlim)
  xlim <- xlim + c(-0.2/nr, if (legend) 0.2 + 0.02 * nch else 0.2/nr) * diff(xlim)
  matplot(xvals, cellNew, ..., type = type, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, xaxt = "n", col = col, lty = lty, pch = pch)
  if (axes && xaxt != "n") {
    axisInt <- function(x, main, sub, lwd, bg, log, asp, ...) axis(1, x, ...)
    mgp. <- par("mgp")
    if (!xtick) 
      mgp.[2L] <- 0
    axisInt(1, at = xvals, labels = xlabs, tick = xtick, mgp = mgp., xaxt = xaxt, ...)
  }
  if (legend) {
    legpretty = ylabs
    legend("topright", legend = legpretty, title = title, col = col, pch = if (type %in% c("p", "b")) 
      pch, lty = if (type %in% c("l", "b")) 
        lty, bty = leg.bty, bg = leg.bg, inset = 0.02)
  }
  return(list(xvals, xlabs))
}
.testFun2 = function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = TRUE, trace.label = deparse(substitute(trace.factor)), 
                     fixed = FALSE, xlab = deparse(substitute(x.factor)), ylab = ylabel, ylim = range(cells, na.rm = TRUE), lty = nc:1, col = 1, pch = c(1:9, 0, letters), xpd = NULL, 
                     leg.bg = par("bg"), leg.bty = "n", xtick = FALSE, xaxt = par("xaxt"), axes = TRUE, ...) {
  ylabel <- paste(deparse(substitute(fun)), "of ", deparse(substitute(response)))
  type <- match.arg(type)
  cells <- tapply(response, list(x.factor, trace.factor), fun)
  nr <- nrow(cells)
  nc <- ncol(cells)
  xvals <- 1:nr
  if (is.ordered(x.factor)) {
    wn <- getOption("warn")
    options(warn = -1)
    xnm <- as.numeric(levels(x.factor))
    options(warn = wn)
    if (!any(is.na(xnm))) 
      xvals <- xnm
  }
  xlabs <- rownames(cells)
  ylabs <- colnames(cells)
  nch <- max(sapply(ylabs, nchar, type = "width"))
  if (is.null(xlabs)) 
    xlabs <- as.character(xvals)
  if (is.null(ylabs)) 
    ylabs <- as.character(1:nc)
  xlim <- range(xvals)
  xleg <- xlim[2] + 0.05 * diff(xlim)
  xlim <- xlim + c(-0.2/nr, if (legend) 0.3 + 0.02 * nch else 0.3/nr) * diff(xlim)
  matplot(xvals, cells, ..., type = type, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, xaxt = "n", col = col, lty = lty, pch = pch)
  if (axes && xaxt != "n") {
    axisInt <- function(x, main, sub, lwd, bg, log, asp, ...) axis(1, x, ...)
    mgp. <- par("mgp")
    if (!xtick) 
      mgp.[2] <- 0
    axisInt(1, at = xvals, labels = xlabs, tick = xtick, mgp = mgp., xaxt = xaxt, ...)
  }
  if (legend) {
    yrng <- diff(ylim)
    yleg <- ylim[2] - 0.1 * yrng
    if (!is.null(xpd) || {
      xpd. <- par("xpd")
      !is.na(xpd.) && !xpd. && (xpd <- TRUE)
    }) {
      op <- par(xpd = xpd)
      on.exit(par(op))
    }
    text(xleg - 0.05, ylim[2] - 0.05 * yrng, paste("  ", trace.label), adj = 0)
    if (!fixed) {
      ord <- sort.list(cells[nr, ], decreasing = TRUE)
      ylabs <- ylabs[ord]
      lty <- lty[1 + (ord - 1)%%length(lty)]
      col <- col[1 + (ord - 1)%%length(col)]
      pch <- pch[ord]
    }
    legend(xleg, yleg, legend = ylabs, col = col, pch = if (type %in% c("p", "b")) 
      pch, lty = if (type %in% c("l", "b")) 
        lty, bty = leg.bty, bg = leg.bg)
  }
  return(list(xvals, xlabs))
}
#.interactionPlotOld = function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = TRUE, trace.label = deparse(substitute(trace.factor)), 
#    fixed = FALSE, xlab = deparse(substitute(x.factor)), ylab = ylabel, ylim = range(cells, na.rm = TRUE), lty = nc:1, col = 1, pch = c(1:9, 0, letters), xpd = NULL, 
#    leg.bg = par("bg"), leg.bty = "n", xtick = FALSE, xaxt = par("xaxt"), axes = TRUE, ...) {
#    ylabel <- paste(deparse(substitute(fun)), "of ", deparse(substitute(response)))
#    type <- match.arg(type)
#    cells <- tapply(response, list(x.factor, trace.factor), fun)
#    nr <- nrow(cells)
#    nc <- ncol(cells)
#    xvals <- 1:nr
#    if (is.ordered(x.factor)) {
#        wn <- getOption("warn")
#        options(warn = -1)
#        xnm <- as.numeric(levels(x.factor))
#        options(warn = wn)
#        if (!any(is.na(xnm))) 
#            xvals <- xnm
#    }
#    xlabs <- rownames(cells)
#    ylabs <- colnames(cells)
#    nch <- max(sapply(ylabs, nchar, type = "width"))
#    if (is.null(xlabs)) 
#        xlabs <- as.character(xvals)
#    if (is.null(ylabs)) 
#        ylabs <- as.character(1:nc)
#    xlim <- range(xvals)
#    xleg <- xlim[2] + 0.05 * diff(xlim)
#    xlim <- xlim + c(-0.2/nr, if (legend) 0.2 + 0.02 * nch else 0.2/nr) * diff(xlim)
#    matplot(xvals, cells, ..., type = type, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, xaxt = "n", col = col, lty = lty, pch = pch)
#    if (axes && xaxt != "n") {
#        axisInt <- function(x, main, sub, lwd, bg, log, asp, ...) axis(1, x, ...)
#        mgp. <- par("mgp")
#        if (!xtick) 
#            mgp.[2] <- 0
#        axisInt(1, at = xvals, labels = xlabs, tick = xtick, mgp = mgp., xaxt = xaxt, ...)
#    }
#    if (legend) {
#        yrng <- diff(ylim)
#        yleg <- ylim[2] - 0.1 * yrng
#        if (!is.null(xpd) || {
#            xpd. <- par("xpd")
#            !is.na(xpd.) && !xpd. && (xpd <- TRUE)
#        }) {
#            op <- par(xpd = xpd)
#            on.exit(par(op))
#        }
#        text(xleg, ylim[2] - 0.05 * yrng, paste("  ", trace.label), adj = 0)
#        if (!fixed) {
#            ord <- sort.list(cells[nr, ], decreasing = TRUE)
#            ylabs <- ylabs[ord]
#            lty <- lty[1 + (ord - 1)%%length(lty)]
#            col <- col[1 + (ord - 1)%%length(col)]
#            pch <- pch[ord]
#        }
#        legend(xleg, yleg, legend = ylabs, col = col, pch = if (type %in% c("p", "b")) 
#            pch, lty = if (type %in% c("l", "b")) 
#            lty, bty = leg.bty, bg = leg.bg)
#    }
#    invisible()
#}
interactionPlot = function(fdo, y = NULL, response = NULL, fun = mean, main, col = 1:2, ...) { ###
  DB = FALSE
  mainmiss = FALSE
  if (missing(main)) 
    mainmiss = TRUE
  if (missing(fdo) || class(fdo)!="facDesign")                                ###
    stop("fdo needs to be an object of class facDesign")                    ###
  parList = list(...)
  old.par <- par(no.readonly = TRUE)
  fdoName = deparse(substitute(fdo))                                          ###
  #    if (class(fdo) != "facDesign") {                                           ###
  #        if (any(is.null(fdo), is.null(y), is.null(response)))                  ###
  #            stop("Factors or response are not given!")                         ###
  #        .interactionPlotOld(fdo, y, response, fun = fun, main = main, ...)     ###
  #        return()                                                               ###
  #    }                                                                          ###
  #    else {                                                                     ###
  if(is.null(response)==FALSE)                                                ###
  {                                                                           ###
    temp=response(fdo)[response]                                               ###
    response(fdo)=temp                                                         ###
  }                                                                           ###
  diagNames = character(0)
  x = cube(fdo)
  runIndex = order(runOrd(fdo))
  x = x[runIndex[1:nrow(x)], ]
  y = response(fdo)[1:nrow(x), ]
  numFac = ncol(x)
  combMat = combn(names(x), 2)
  if (numFac == 2) {
    facName2 = combMat[1, 1]
    facName1 = combMat[2, 1]
    temp = with(cbind(y, x), .testFun(eval(parse(text = facName2)), eval(parse(text = facName1)), xlab = facName1, response = y, trace.label = facName1, 
                                      ylim = range(y), axes = F, fun, title = facName1, col = col, ...))
    tempList = parList
    tempList$col = 1
    tempList$lwd = 1
    tempList$side = c(2)
    do.call(axis, tempList)
    box()
    tempList$at = temp[[1]]
    tempList$labels = temp[[2]]
    tempList$side = c(1)
    do.call(axis, tempList)
    invisible()
  }
  for (r in 1:ncol(response(fdo))) {
    if (r > 1) 
      dev.new()
    y = response(fdo)[1:nrow(x), r]
    par(mfrow = c(numFac, numFac))
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0, 0, 8, 8))
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
    for (i in 1:ncol(combMat)) {
      facName1 = combMat[1, i]
      facName2 = combMat[2, i]
      rowNum = .letterPos(facName1)
      colNum = .letterPos(facName2)
      if (DB) {
        print(numFac)
        print(rowNum)
        print(colNum)
        print(facName1)
        print(facName2)
        print(y)
        print(fun)
        cat(paste(i, "\t", c(rowNum, colNum)))
      }
      par(mfg = c(rowNum, colNum))
      temp = with(cbind(x, y), .testFun(eval(parse(text = facName2)), eval(parse(text = facName1)), response = y, trace.label = facName1, ylim = range(y), 
                                        axes = F, fun = fun, title = facName1, col = col, ...))
      if (colNum == numFac) {
        tempList = parList
        tempList$col = 1
        tempList$lwd = 1
        tempList$side = 4
        do.call(axis, tempList)
      }
      if (rowNum == 1) {
        tempList = parList
        tempList$col = 1
        tempList$lwd = 1
        tempList$side = 3
        tempList$at = temp[[1]]
        tempList$labels = temp[[2]]
        do.call(axis, tempList)
      }
      box(which = "plot")
      par(mfg = c(rowNum, rowNum))
      plot(c(-1, 1), c(-1, 1), type = "n", axes = F, xlab = "", ylab = "", main = "")
      text(0, 0, facName1, cex = 4)
      diagNames = c(diagNames, facName1)
    }
    par(mfg = c(numFac, numFac))
    plot(c(-1, 1), c(-1, 1), type = "n", axes = F, xlab = "", ylab = "", main = "")
    text(0, 0, setdiff(names(x), diagNames), cex = 4)
    if (mainmiss) {
      main = paste("Interaction plot for", names(response(fdo))[r], "in", fdoName)         ###
      title(main, outer = T, ...)
    }
    else title(main[r], outer = T, ...)
  }
  #    }
  par(old.par)
  invisible()
}
class(interactionPlot) <- "invisible" 

# -------------- cg._g.r -------------- 
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

cgRunChart=function(x, target, tolerance, ref.interval, facCg, facCgk, n = 0.2, type, col, pch, xlim, ylim,main, conf.level = 0.95, cex.val = 1, cgOut = TRUE)
{
  if (missing(x)) 
    stop("x must be given as a vector")
  if (missing(target)) {
    target = mean(x)
    targetmissing = FALSE
  }
  else targetmissing = TRUE
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
    warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
  if (missing(tolerance)) {
    width = ref.ar/facCg
    tolerance = numeric(2)
    tolerance[1] = mean(x) - width/2
    tolerance[2] = mean(x) + width/2
  }
  quant1 = qnorm((1 - ref.interval)/2, mean, sd)
  quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean, sd)
  if (length(tolerance) != 2) 
    stop("tolerance has wrong length")
  if (missing(type)) 
    type = "b"
  if (missing(col)) 
    col = 1
  if (missing(pch)) 
    pch = 19
  if (missing(xlim)) 
    xlim = c(0, length(x))
  if (missing(ylim)) 
    ylim = c(min(x, target - n/2 * (abs(diff(tolerance))), quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))), quant1, quant2))
  if (missing(main))
    main="Run Chart"
  Cg = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[1]]
  Cgk = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[2]]
  par.old = par()$mar
  par(mar = c(5.1, 4.1, 4.1, 5.1))
  plot(x, type = type, col = col, main = main, pch = pch, xlim = c(xlim[1] - 0.05 * xlim[2], xlim[2]), ylim = ylim)
  abline(h = target)
  rect(xleft = min(xlim) - 0.09 * max(xlim), xright = -0.005 * max(xlim), ybottom = min(ylim), ytop = max(ylim), col = "white", border = "white")
  text(min(xlim) - 0.1 * max(xlim), target, pos = 4, srt = 0, "target", xpd = TRUE, bg = "white", cex = cex.val)
  if (targetmissing == TRUE) 
    abline(h = mean, lty = 2, col = "seagreen")
  abline(h = quant1, col = "seagreen", lty = 3)
  abline(h = quant2, col = "seagreen", lty = 3)
  abline(h = c(target + n/2 * (abs(diff(tolerance))), target - n/2 * (abs(diff(tolerance)))))
  lines(lowess(x), col = "red")
  temp = axis(4, at = c(target + n/2 * (abs(diff(tolerance))), target - n/2 * (abs(diff(tolerance))), mean(x)), labels = FALSE)
  text(max(xlim) + 0.075 * max(xlim), temp[3], srt = 0, adj = 0, substitute(x[tar] + a %.% b, list(a = round(n/2, 4), b = "T")), xpd = TRUE, cex = cex.val)
  text(max(xlim) + 0.075 * max(xlim), temp[2], srt = 0, adj = 0, expression(bar(x)), xpd = TRUE, cex = cex.val, col = "seagreen")
  text(max(xlim) + 0.075 * max(xlim), temp[1], srt = 0, adj = 0, substitute(x[tar] - a %.% b, list(a = round(n/2, 4), b = "T")), xpd = TRUE, cex = cex.val)
  lt1 = round(((1 - ref.interval)/2) * 100, 3)
  lt2 = round(((ref.interval + (1 - ref.interval)/2)) * 100, 3)
  text(max(xlim) + 0.075 * max(xlim), quant1, srt = 0, adj = 0, substitute(x[a * b], list(a = lt1, b = "%")), xpd = TRUE, cex = cex.val, col = "seagreen")
  text(max(xlim) + 0.075 * max(xlim), quant2, srt = 0, adj = 0, substitute(x[a * b], list(a = lt2, b = "%")), xpd = TRUE, cex = cex.val, col = "seagreen")
  if(cgOut == TRUE)
  {
    legend("topright",legend=(c(paste("Cg: ",Cg),paste("Cgk :",Cgk))),inset=c(0,0.06),bty="n")
  }
  par(mar = par.old)
  invisible(list(Cg, Cgk))
}

cgHist=function(x, target, tolerance, ref.interval, facCg, facCgk, n = 0.2, col, xlim, ylim, main, conf.level = 0.95, cex.val = 1, cgOut = TRUE)
{
  if (missing(x)) 
    stop("x must be given as a vector")
  if (missing(target)) {
    target = mean(x)
    targetmissing = FALSE
  }
  else targetmissing = TRUE
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
    warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
  if (missing(tolerance)) {
    width = ref.ar/facCg
    tolerance = numeric(2)
    tolerance[1] = mean(x) - width/2
    tolerance[2] = mean(x) + width/2
  }
  quant1 = qnorm((1 - ref.interval)/2, mean, sd)
  quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean, sd)
  if (length(tolerance) != 2) 
    stop("tolerance has wrong length")
  if (missing(col)) 
    col = "lightblue"
  if (missing(xlim)) 
    xlim = c(0, length(x))
  if (missing(ylim)) 
    ylim = c(min(x, target - n/2 * (abs(diff(tolerance))), quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))), quant1, quant2))
  if (missing(main))
    main=paste("Histogram of", deparse(substitute(x)), "- target")
  Cg = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[1]]
  Cgk = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[2]]
  
  x.c = x - target
  temp = hist(x.c, plot = FALSE)
  #   par.old = par()$mar
  #    par(mar = c(3.1, 4.1, 0.55, 2.1))
  hist(x.c, freq = FALSE, ylim = c(0, max(density(x.c)$y, temp$density) + 0.4 * max(density(x.c)$y, temp$density)), col = col, main=main)
  test = t.test(x.c, mu = 0, conf.level = conf.level)
  lines(x = c(0, 0), y = c(0, max(density(x.c)$y, temp$density) + 0.3 * max(density(x.c)$y, temp$density)), col = "red")
  lines(x = c(test$conf.int[1], test$conf.int[1]), y = c(0, max(density(x.c)$y, temp$density) + 0.3 * max(density(x.c)$y, temp$density)), col = "blue", lty = 3)
  lines(x = c(test$conf.int[2], test$conf.int[2]), y = c(0, max(density(x.c)$y, temp$density) + 0.3 * max(density(x.c)$y, temp$density)), col = "blue", lty = 3)
  
  if(cgOut == TRUE)
  {
    legend("topright", legend = c("conf.int",paste("Cg: ",Cg),paste("Cgk :",Cgk)), col = c("blue",-1,-1), lty = c(3,-1,-1), inset = c(0.01, 0.06))
  }
  else
    legend("topright", legend = c("conf.int"), col = c("blue"), lty = c(3), inset = c(0.01, 0.06))
  #    par(mar = par.old)
  legend("topleft", legend = c(expression(paste(H[0], " : Bias = 0")), paste("t-value: ", round(test$statistic, 3)), paste("p-value: ", round(test$p.value, 
                                                                                                                                              3))), inset = c(-0.01, 0.04), bty = "n")
  lines(density(x.c))
  #   par(mar = par.old)
  #   par(mar = par.old)
  invisible(list(Cg, Cgk))
}

cgToleranceView=function(x, target, tolerance, ref.interval, facCg, facCgk, n = 0.2, type, col, pch, xlim, ylim, main, conf.level = 0.95, cex.val = 1, cgOut = TRUE)
{
  if (missing(x)) 
    stop("x must be given as a vector")
  if (missing(target)) {
    target = mean(x)
    targetmissing = FALSE
  }
  else targetmissing = TRUE
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
    warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
  if (missing(tolerance)) {
    width = ref.ar/facCg
    tolerance = numeric(2)
    tolerance[1] = mean(x) - width/2
    tolerance[2] = mean(x) + width/2
  }
  quant1 = qnorm((1 - ref.interval)/2, mean, sd)
  quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean, sd)
  if (length(tolerance) != 2) 
    stop("tolerance has wrong length")
  if (missing(type)) 
    type = "b"
  if (missing(col)) 
    col = 1
  if (missing(pch)) 
    pch = 19
  if (missing(xlim)) 
    xlim = c(0, length(x))
  if (missing(ylim)) 
    ylim = c(min(x, target - n/2 * (abs(diff(tolerance))), quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))), quant1, quant2))
  if (missing(main))
    main="Tolerance View"   
  Cg = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[1]]
  Cgk = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[2]]    
  plot(x, type=type, col = col, main = main, pch = pch, xlim = xlim, ylim = c(min(x, tolerance[2], tolerance[1], target + n/2 * (tolerance[2] - 
                                                                                                                                   tolerance[1]), target - n/2 * (tolerance[2] - tolerance[1])), max(x, tolerance[2], tolerance[1], target + n/2 * (tolerance[2] - tolerance[1]), target - 
                                                                                                                                                                                                       n/2 * (tolerance[2] - tolerance[1]))))
  abline(h = c(tolerance[1], tolerance[2]), lty = 2, col = "red")
  abline(h = target)
  abline(h = c(target + n/2 * (tolerance[2] - tolerance[1]), target - n/2 * (tolerance[2] - tolerance[1])))
  if(cgOut == TRUE)
  {
    legend("topright",legend=(c(paste("Cg: ",Cg),paste("Cgk :",Cgk))),inset=c(0,0.06),bty="n")
  }
  invisible(list(Cg, Cgk))
}

cg = function(x, target, tolerance, ref.interval, facCg, facCgk, n = 0.2, type, col, pch, xlim, ylim, conf.level = 0.95, cex.val = 1.5) {
  old.par <- par(no.readonly = TRUE)                                             #save old par settings
  if (missing(x)) 
    stop("x must be given as a vector")
  if (missing(target)) {
    target = mean(x)
    targetmissing = FALSE
  }
  else targetmissing = TRUE
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
    warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
  if (missing(tolerance)) {
    width = ref.ar/facCg
    tolerance = numeric(2)
    tolerance[1] = mean(x) - width/2
    tolerance[2] = mean(x) + width/2
  }
  quant1 = qnorm((1 - ref.interval)/2, mean, sd)
  quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean, sd)
  if (length(tolerance) != 2) 
    stop("tolerance has wrong length")
  if (missing(type)) 
    type = "b"
  if (missing(col)) 
    col = 1
  if (missing(pch)) 
    pch = 19
  if (missing(xlim)) 
    xlim = c(0, length(x))
  if (missing(ylim)) 
    ylim = c(min(x, target - n/2 * (abs(diff(tolerance))), quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))), quant1, quant2))
  Cg = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[1]]
  Cgk = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[2]]
  layout(matrix(data = c(1, 1, 1, 1, 1, 1, 2, 3, 4), ncol = 3, nrow = 3))
  #  par.old = par()$mar
  #  par(mar = c(5.1, 4.1, 4.1, 10.1))
  cgRunChart(x=x, target=target, tolerance=tolerance, ref.interval=ref.interval, 
             facCg=facCg, facCgk=facCgk, n = n, type=type, col=col, pch=pch, 
             xlim=xlim, ylim=ylim, main="Run Chart", conf.level = conf.level, 
             cex.val = cex.val, cgOut = FALSE)
  #   par(mar = par.old)
  plot(0:6, 0:6, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
  text(2.5, 5, expression(bar(x)), pos = 2, cex = cex.val)
  text(2.5, 5, paste("=", round(mean, 2)), pos = 4, cex = cex.val)
  text(2.5, 4, expression(s), pos = 2, cex = cex.val)
  text(2.5, 4, paste("=", round(sd, 2)), pos = 4, cex = cex.val)
  text(2.5, 3, expression(target), pos = 2, cex = cex.val)
  text(2.5, 3, paste("=", round(target, 5)), pos = 4, cex = cex.val)
  text(2.5, 2, expression(bold(C[g])), pos = 2, cex = cex.val)
  text(2.5, 2, paste("=", round(Cg, 2)), pos = 4, cex = cex.val)
  text(2.5, 1, expression(bold(C[gk])), pos = 2, cex = cex.val)
  text(2.5, 1, paste("=", round(Cgk, 2)), pos = 4, cex = cex.val)
  rect(0, 0, 6, 6)
  par(mar = c(3.1, 4.1, 0.55, 2.1))
  cgHist(x=x, target=target, tolerance=tolerance, ref.interval=ref.interval, 
         facCg=facCg, facCgk=facCgk, n = n, col="lightblue",xlim=xlim, 
         ylim=ylim, main=paste("Histogram of", deparse(substitute(x)), "- target"),
         conf.level = conf.level, cex.val = cex.val, cgOut = FALSE)
  #    par(mar = par.old)
  par(mar = c(3.1, 4.1, 1.55, 2.1))
  cgToleranceView(x=x, target=target, tolerance=tolerance, ref.interval=ref.interval, 
                  facCg=facCg, facCgk=facCgk, n = n, type=type, col=col, pch=pch, 
                  xlim=xlim, ylim=ylim, main="Tolerance View", conf.level = conf.level, cex.val = cex.val, cgOut = FALSE)
  par(mar = old.par)
  invisible(list(Cg, Cgk))
} 

#set.seed(123)
#temp=rnorm(125,mean = 10.01 ,sd = 0.1)
#a=cg(temp, target = 10, tolerance = c(8,12),type="b")
#b=.cg(temp, target = 10, tolerance = c(8,12))
#c=cgRunChart(temp, target = 10, tolerance = c(8,12),n = 0.2, type="b", col=4, pch=3,conf.level = 0.95, cex.val = 1,main="Hello")
#d=cgHist(temp, target = 10, tolerance = c(8,12))
#e=cgToleranceView(temp, target = 10, tolerance = c(8,12))




# -------------- gag_n.r -------------- 
setClass("gageRR", representation = representation(X = "data.frame", ANOVA = "aov", RedANOVA = "aov", method = "character", Estimates = "list", Varcomp = "list", 
                                                   Sigma = "numeric", GageName = "character", GageTolerance = "numeric", DateOfStudy = "character", PersonResponsible = "character", Comments = "character", 
                                                   b = "factor", a = "factor", y = "numeric", facNames = "character", numO = "numeric", numP = "numeric", numM = "numeric"))
setMethod("show", signature(object = "gageRR"), function(object) {
  print(as.data.frame(object))
})
setMethod("[", signature(x = "gageRR", i = "ANY", j = "ANY"), function(x, i, j) {
  x@X[i, j]
})
setMethod("summary", signature(object = "gageRR"), function(object) {
  if (all(is.na(object@X$Measurement))) 
    return(print(as.data.frame(object)))
  cat("\n")
  cat(paste("Operators:\t", object@numO, "\tParts:\t", object@numP))
  cat("\n")
  cat(paste("Measurements:\t", object@numM, "\tTotal:\t", nrow(object@X)))
  cat("\n")
  cat("----------")
  cat("\n")
  return(gageRR(object, method = object@method))
})
setMethod("response", "gageRR", function(object) {
  out = object@X$Measurement
  return(out)
})
setReplaceMethod("response", "gageRR", function(object, value) {
  object@X$Measurement = value
  return(object)
})
setMethod("names", signature(x = "gageRR"), function(x) {
  return(names(as.data.frame(x)))
})
setMethod("as.data.frame", "gageRR", function(x, row.names = NULL, optional = FALSE, ...) {
  return(x@X)
})
as.data.frame.gageRR = function(x, row.names = NULL, optional = FALSE, ...) {
  return(x@X)
}
setGeneric("tolerance", function(x) standardGeneric("tolerance"))
setGeneric("tolerance<-", function(x, value) standardGeneric("tolerance<-"))
setMethod("tolerance", "gageRR", function(x) unlist(x@GageTolerance))
setReplaceMethod("tolerance", "gageRR", function(x, value) {
  if (!is.numeric(value)) 
    stop(paste(deparse(substitute(value)), "needs to be numeric"))
  x@GageTolerance = value
  return(x)
})
setGeneric("sigma", function(x) standardGeneric("sigma"))
setGeneric("sigma<-", function(x, value) standardGeneric("sigma<-"))
setMethod("sigma", "gageRR", function(x) unlist(x@Sigma))
setReplaceMethod("sigma", "gageRR", function(x, value) {
  if (!is.numeric(value)) 
    stop(paste(deparse(substitute(value)), "needs to be numeric"))
  x@Sigma = value
  return(x)
})
.aip = function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = FALSE, trace.label = deparse(substitute(trace.factor)), 
                fixed = FALSE, xlab = deparse(substitute(x.factor)), ylab = ylabel, ylim = range(cellNew, na.rm = TRUE), lty = nc:1, col = 1, pch = c(1L:9, 0, letters), 
                xpd = NULL, leg.bg = par("bg"), leg.bty = "o", xtick = FALSE, xaxt = par("xaxt"), axes = TRUE, title = "", ...) {
  ylabel <- paste(deparse(substitute(fun)), "of ", deparse(substitute(response)))
  type <- match.arg(type)
  cellNew <- tapply(response, list(x.factor, trace.factor), fun)
  nr <- nrow(cellNew)
  nc <- ncol(cellNew)
  xvals <- 1L:nr
  if (is.ordered(x.factor)) {
    wn <- getOption("warn")
    options(warn = -1)
    xnm <- as.numeric(levels(x.factor))
    options(warn = wn)
    if (!any(is.na(xnm))) 
      xvals <- xnm
  }
  xlabs <- rownames(cellNew)
  ylabs <- colnames(cellNew)
  nch <- max(sapply(ylabs, nchar, type = "width"))
  if (is.null(xlabs)) 
    xlabs <- as.character(xvals)
  if (is.null(ylabs)) 
    ylabs <- as.character(1L:nc)
  xlim <- range(xvals)
  xleg <- xlim[2L] + 0.05 * diff(xlim)
  xlim <- xlim + c(-0.2/nr, if (legend) 0.2 + 0.02 * nch else 0.2/nr) * diff(xlim)
  matplot(xvals, cellNew, ..., type = type, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, xaxt = "n", col = col, lty = lty, pch = pch)
  if (axes && xaxt != "n") {
    axisInt <- function(x, main, sub, lwd, bg, log, asp, ...) axis(1, x, ...)
    mgp. <- par("mgp")
    if (!xtick) 
      mgp.[2L] <- 0
    axisInt(1, at = xvals, labels = xlabs, tick = xtick, mgp = mgp., xaxt = xaxt, ...)
  }
  if (legend) {
    legend("topright", legend = ylabs, title = title, col = col, pch = if (type %in% c("p", "b")) 
      pch, lty = if (type %in% c("l", "b")) 
        lty, bty = leg.bty, bg = leg.bg, inset = 0.02)
  }
  legend("topright", legend = ylabs, title = title, col = col, pch = if (type %in% c("p", "b")) 
    pch, lty = if (type %in% c("l", "b")) 
      lty, bty = leg.bty, bg = leg.bg, inset = c(-0.2, 0), xpd = TRUE)
  invisible()
}
gageRRDesign = function(Operators = 3, Parts = 10, Measurements = 3, method = "crossed", sigma = 6, randomize = TRUE) {
  method = method
  opvec = factor()
  partvec = factor()
  yName = aName = bName = abName = NA
  yName = "Measurement"
  aName = "Operator"
  bName = "Part"
  abName = "Operator:Part"
  Operators = unique(Operators)
  Parts = unique(Parts)
  if (!is.numeric(sigma)) 
    stop("sigma needs to be numeric")
  if (method != "nested" && method != "crossed") 
    warning("unknown method specified --> defaulting to \"method = crossed\"")
  if (!is.numeric(Measurements)) 
    stop("Number of Measurements per Part not specified!")
  else Measurements = round(Measurements[1])
  if (!is.numeric(Operators) && !is.character(Operators)) 
    stop("Operator needs to be numeric 'Operator = 3' or character 'Operator = c(\"A\",\"B\", \"C\"")
  if (is.numeric(Operators)) 
    opvec = factor(LETTERS[1:Operators[1]])
  if (is.character(Operators)) 
    opvec = factor(Operators)
  if (length(unique(opvec)) > 26) 
    stop("To many Operators!")
  if (length(unique(opvec)) < 2) 
    stop("Not enough Operators")
  if (!is.numeric(Parts) && !is.character(Parts)) 
    stop("Parts needs to be numeric 'Parts = 3' or character 'Parts = c(\"A\",\"B\", \"C\"")
  if (is.numeric(Parts)) 
    partvec = factor(LETTERS[1:Parts[1]])
  if (is.character(Parts)) 
    partvec = factor(Parts)
  if (length(unique(partvec)) > 26) 
    stop("To many Parts!")
  if (length(unique(partvec)) < 2) 
    stop("To few Parts")
  Measurement = rep(NA, (length(opvec) * length(partvec) * Measurements))
  outFrame = data.frame()
  if (method == "crossed") {
    temp = expand.grid(opvec, partvec)
    o = rep(temp[, 1], Measurements)
    p = rep(temp[, 2], Measurements)
  }
  else {
    p = rep(sort(rep(partvec, length(opvec))), Measurements)
    o = (rep(opvec, length(Measurement)/length(opvec)))
    p = p[order(o,p)]
    o = o[order(o,p)]
  }
  if (randomize) 
    outFrame = data.frame(StandardOrder = 1:length(Measurement), RunOrder = sample(1:length(Measurement), length(Measurement)), Operator = factor(o), Part = factor(p), 
                          Measurement)
  else outFrame = data.frame(StandardOrder = 1:length(Measurement), RunOrder = 1:length(Measurement), Operator = factor(o), Part = factor(p), Measurement)
  outFrame = outFrame[order(outFrame$RunOrder), ]
  gageRRObj = new("gageRR")
  gageRRObj@facNames = c(yName, aName, bName, abName)
  names(gageRRObj@facNames) = c("yName", "aName", "bName", "abName")
  gageRRObj@Sigma = sigma
  gageRRObj@method = method
  gageRRObj@a = factor(o)
  gageRRObj@b = factor(p)
  gageRRObj@y = as.numeric(Measurement)
  gageRRObj@method = method
  gageRRObj@Sigma = sigma
  gageRRObj@X = outFrame
  return(gageRRObj)
}
gageRR = function(gdo, method = "crossed", sigma = 6, alpha = 0.25, DM = NULL, HM = NULL, tolerance = NULL, dig = 3, ...) {
  if (method %in% c("crossed", "nested")) 
    method = method
  else method = gdo@method
  yName = names(gdo)[5]
  aName = names(gdo)[3]
  bName = names(gdo)[4]
  if(method == "crossed")
    abName = paste(aName, ":", bName, sep = "")
  if(method == "nested")
    abName = paste(bName, "(", aName, ")", sep = "")
  bTobName = paste(bName, "to", bName, sep = " ")
  a = gdo@X[, aName]
  b = gdo@X[, bName]
  y = gdo@X[, yName]
  nestedFormula = as.formula(paste(yName, "~", aName, "/", bName))
  crossedFormula = as.formula(paste(yName, "~", aName, "*", bName))
  reducedFormula = as.formula(paste(yName, "~", aName, "+", bName))
  if (!is.null(tolerance)) 
    tolerance(gdo) = tolerance
  if (is.na(y) || !is.numeric(y)) 
    stop("Measurements need to be numeric")
  if (method == "nested") {
    numA <- nlevels(a[, drop = T])
    numB <- nlevels(b[, drop = T])
    numMPP <- length(y)/((numB) * numA)
    gdo@numO = numA
    gdo@numP = numB
    gdo@numM = numMPP
    fit = aov(nestedFormula, data = gdo)
    meanSq <- anova(fit)[, 3]
    gdo@ANOVA = fit
    gdo@method = "nested"
    MSa = meanSq[1]
    MSab = meanSq[2]
    MSe = meanSq[3]
    Cerror = MSe
    Cb = (MSab - MSe)/numMPP
    Ca = (MSa - MSab)/(numB * numMPP)
    if (Ca <= 0) 
      Ca = 0
    if (Cb <= 0) 
      Cb = 0
    Cab = 0
    totalRR = Ca + Cab + Cerror
    repeatability = Cerror
    reproducibility = Ca
    bTob = Cb
    totalVar = Cb + Ca + Cab + Cerror
    estimates = list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
    varcomp = list(totalRR = totalRR, repeatability = repeatability, reproducibility = reproducibility, bTob = bTob, totalVar = totalVar)
    gdo@Estimates = estimates
    gdo@Varcomp = varcomp
  }
  if (method == "crossed") {
    numA <- nlevels(a[, drop = T])
    numB <- nlevels(b[, drop = T])
    numMPP <- length(a)/(numA * numB)
    gdo@numO = numA
    gdo@numP = numB
    gdo@numM = numMPP
    fit = aov(crossedFormula, data = gdo)
    model <- anova(fit)
    gdo@ANOVA = fit
    gdo@method = "crossed"
    MSb = MSa = MSab = MSe = 0
    if (bName %in% row.names(model)) 
      MSb = model[bName, "Mean Sq"]
    else warning(paste("missing factor", bName, "in model"))
    if (aName %in% row.names(model)) 
      MSa = model[aName, "Mean Sq"]
    else warning(paste("missing factor", aName, "in model"))
    if (abName %in% row.names(model)) 
      MSab = model[abName, "Mean Sq"]
    else warning(paste("missing interaction", abName, "in model"))
    if ("Residuals" %in% row.names(model)) 
      MSe = model["Residuals", "Mean Sq"]
    else warning("missing Residuals in model")
    Cb = Ca = Cab = Cerror = 0
    Cb = (MSb - MSab)/(numA * numMPP)
    Ca = (MSa - MSab)/(numB * numMPP)
    Cab = (MSab - MSe)/(numMPP)
    Cerror = (MSe)
    gdo@RedANOVA = gdo@ANOVA
    if ((Cab < 0) || (model[abName, "Pr(>F)"] >= alpha)) {
      redFit <- aov(reducedFormula, data = gdo)
      model <- anova(redFit)
      MSb = MSa = MSab = MSe = 0
      if (bName %in% row.names(model)) 
        MSb = model[bName, "Mean Sq"]
      else warning(paste("missing factor", bName, "in model"))
      if (aName %in% row.names(model)) 
        MSa = model[aName, "Mean Sq"]
      else warning(paste("missing factor", aName, "in model"))
      if ("Residuals" %in% row.names(model)) 
        MSe = model["Residuals", "Mean Sq"]
      else warning("missing Residuals in model")
      Cb = Ca = Cab = Cerror = 0
      Cb = (MSb - MSe)/(numA * numMPP)
      Ca = (MSa - MSe)/(numB * numMPP)
      Cab = 0
      Cerror = (MSe)
      gdo@RedANOVA = redFit
    }
    gdo@method = "crossed"
    Ca = max(0, Ca)
    Cb = max(0, Cb)
    Cab = max(0, Cab)
    totalRR = Ca + Cab + Cerror
    repeatability = Cerror
    reproducibility = Ca + Cab
    bTob = max(0, Cb)
    totalVar = Cb + Ca + Cab + Cerror
    estimates = list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
    varcomp = list(totalRR = totalRR, repeatability = repeatability, reproducibility = reproducibility, a = Ca, a_b = Cab, bTob = bTob, totalVar = totalVar)
    gdo@Estimates = estimates
    gdo@Varcomp = varcomp
  }
  cat("\n")
  cat(paste("AnOVa Table - ", gdo@method, "Design\n"))
  print(summary(gdo@ANOVA))
  cat("\n")
  cat("----------\n")
  if (!identical(gdo@RedANOVA, gdo@ANOVA) && gdo@method == "crossed") {
    cat(paste("AnOVa Table Without Interaction - ", gdo@method, "Design\n"))
    print(summary(gdo@RedANOVA))
    cat("\n")
    cat("----------\n")
  }
  Source = names(gdo@Varcomp)
  Source[Source == "repeatability"] = " repeatability"
  Source[Source == "reproducibility"] = " reproducibility"
  Source[Source == "a_b"] = paste("  ", abName)
  Source[Source == "a"] = paste("  ", aName)
  Source[Source == "bTob"] = bTobName
  VarComp = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))]), 3)
  Contribution = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))])/as.numeric(gdo@Varcomp[length(gdo@Varcomp)]), 3)
  VarComp = t(data.frame(gdo@Varcomp))
  VarCompContrib = VarComp/gdo@Varcomp$totalVar
  Stdev = sqrt(VarComp)
  StudyVar = Stdev * gdo@Sigma
  StudyVarContrib = StudyVar/StudyVar["totalVar", ]
  SNR = 1
  ptRatio = NULL
  temp = NULL
  if ((length(gdo@GageTolerance) > 0) && (gdo@GageTolerance > 0)) {
    ptRatio = StudyVar/gdo@GageTolerance
    temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib, ptRatio)
    names(temp)[6] = c("P/T Ratio")
    row.names(temp) = c(Source)
  }
  else {
    temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib)
    row.names(temp) = c(Source)
  }
  cat("\n")
  cat("Gage R&R\n")
  tempout = temp
  print(format(tempout, digits = dig))
  cat("\n")
  cat("---\n")
  cat(" * Contrib equals Contribution in %\n")
  SNRTemp = sqrt(2) * (temp[bTobName, "Stdev"]/temp["totalRR", "Stdev"])
  if (SNRTemp > 1) 
    SNR = SNRTemp
  cat(paste(" **Number of Distinct Categories (truncated signal-to-noise-ratio) =", floor(SNR), "\n"))
  cat("\n")
  invisible(gdo)
}
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "gageRR"), function(x, y, main, xlab, ylab, col, lwd, fun = mean, ...) {
  horiz = FALSE
  parList = list(...)
  gdo = x
  yName = names(gdo)[5]
  aName = names(gdo)[3]
  bName = names(gdo)[4]
  abName = paste(aName, ":", bName, sep = "")
  if (missing(col)) 
    col = 2:(length(unique(gdo[, 3])) + 1)
  if (missing(lwd)) 
    lwd = 1
  par(mfrow = c(3, 2))
  temp = NULL
  Source = names(gdo@Varcomp)
  VarComp = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))]), 3)
  Contribution = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))])/as.numeric(gdo@Varcomp[length(gdo@Varcomp)]), 3)
  VarComp = t(data.frame(gdo@Varcomp))
  VarCompContrib = VarComp/gdo@Varcomp$totalVar
  Stdev = sqrt(VarComp)
  StudyVar = Stdev * gdo@Sigma
  StudyVarContrib = StudyVar/StudyVar["totalVar", ]
  if ((length(gdo@GageTolerance) > 0) && (gdo@GageTolerance > 0)) {
    ptRatio = StudyVar/gdo@GageTolerance
    temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib, ptRatio)
    contribFrame = data.frame(VarCompContrib, StudyVarContrib, ptRatio)
    names(temp)[6] = c("P/T Ratio")
    row.names(temp) = c(Source)
    SNR = sqrt(2 * (temp["bTob", "VarComp"]/temp["totalRR", "VarComp"]))
  }
  else {
    temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib)
    contribFrame = data.frame(VarCompContrib, StudyVarContrib)
  }
  bTob = paste(bName, "To", bName, sep = "")
  Source[Source == "bTob"] = bTob
  row.names(contribFrame) = Source
  if (gdo@method == "crossed") 
    contribFrame = contribFrame[-match(c("totalVar", "a", "a_b"), row.names(temp)), ]
  else contribFrame = contribFrame[-match(c("totalVar"), row.names(temp)), ]
  numBars = ncol(contribFrame)
  ymax = max(max(contribFrame))
  main1 = NA
  if (missing(main) || is.na(main[1])) 
    main1 = "Components of Variation"
  else main1 = main[1]
  xlab1 = NA
  if (missing(xlab) || is.na(xlab[1])) 
    xlab1 = "component"
  else xlab1 = xlab[1]
  ylab1 = NA
  if (missing(ylab) || is.na(ylab[1])) 
    ylab1 = ""
  else ylab1 = ylab[1]
  argList = list(...)
  redList = argList[names(argList) != "cex"]
  mybp = do.call(barplot, c(list(t(contribFrame), xlab = xlab1, ylab = ylab1, main = main1, names.arg = rep("", 4), axes = F, beside = T, ylim = c(0, 1.3 * 
                                                                                                                                                     ymax), col = col[1:numBars]), redList))
  axis(1, at = colMeans(mybp), labels = names(as.data.frame(t(contribFrame))), ...)
  axis(2, ...)
  box()
  legend("topright", names(contribFrame), col = col[1:numBars], pch = c(15, 15), horiz = horiz, inset = 0.02)
  if (gdo@method == "crossed") {
    main2 = NA
    if (missing(main) || is.na(main[2])) 
      main2 = paste(yName, "by", bName)
    else main2 = main[2]
    xlab2 = NA
    if (missing(xlab) || is.na(xlab[2])) 
      xlab2 = bName
    else xlab2 = xlab[2]
    ylab2 = NA
    if (missing(ylab) || is.na(ylab[2])) 
      ylab2 = yName
    else ylab2 = ylab[2]
    boxplot(split(gdo[, yName], gdo[, bName]), xlab = xlab2, ylab = ylab2, main = main2, ...)
    mByPa = split(gdo[, 5], as.numeric(gdo[, 4]))
    lines(sort(as.numeric(gdo[, 4])), lapply(mByPa, median)[sort(as.numeric(gdo[, 4]))], lwd = lwd)
    points(sort(as.numeric(gdo[, 4])), lapply(mByPa, median)[sort(as.numeric(gdo[, 4]))], lwd = lwd, pch = 13, cex = 2)
    main3 = NA
    if (missing(main) || is.na(main[3])) 
      main3 = paste(yName, "by", aName)
    else main3 = main[3]
    xlab3 = NA
    if (missing(xlab) || is.na(xlab[3])) 
      xlab3 = aName
    else xlab3 = xlab[3]
    ylab3 = NA
    if (missing(ylab) || is.na(ylab[3])) 
      ylab3 = yName
    else ylab3 = ylab[3]
    colVec = .mapping(gdo[, 3], sort(unique(gdo[, 3])), col[1:length(unique(gdo[, 3]))])
    boxplot(split(gdo[, yName], gdo[, aName]), col = colVec, xlab = xlab3, ylab = ylab3, main = main3, ...)
    mByOp = split(gdo[, 5], as.numeric(gdo[, 3]))
    lines(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, mean)[sort(names(mByOp))], lwd = lwd)
    points(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, median)[sort(names(mByOp))], lwd = lwd, pch = 13, cex = 2)
    agg = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = mean)
    tab = table(agg[, 2])
    sgSize = tab[1]
    aggSd = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = sd)
    tab = table(aggSd[, 2])
    sm = mean(aggSd[, 3])
    aggMean = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = mean)
    xm = mean(agg[, 3])
    UCL = xm + ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
    LCL = xm - ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
    values = c(UCL, xm, LCL)
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))
    plot(agg[, 3], type = "n", axes = FALSE, xlab = aName, ylab = expression(bar(x)), main = expression(paste(bar(x), " Chart")))
    box()
    abline(h = xm, col = 3)
    abline(h = UCL, col = 2)
    abline(h = LCL, col = 2)
    axis(2)
    axis(4, at = c(xm, UCL, LCL), labels = c("", "", ""))
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, xm, substitute(bar(x) == xm, list(xm = round(xm, 2))), adj = 0, srt = 0, xpd = TRUE)
    par(mar = old.par)
    j = 1
    for (i in 1:length(tab)) {
      lines(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
      points(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
      if (i < length(tab)) 
        abline(v = j + tab[i] - 1 + 0.5, lty = 2)
      axis(1, at = j, labels = names(tab[i]))
      j = j + tab[i]
    }
    main4 = NA
    if (missing(main) || is.na(main[4])) 
      main4 = paste("Interaction", abName)
    else main4 = main[4]
    xlab4 = NA
    if (missing(xlab) || is.na(xlab[4])) 
      xlab4 = names(gdo)[4]
    else xlab4 = xlab[4]
    ylab4 = NA
    if (missing(ylab) || is.na(ylab[4])) 
      ylab4 = paste(as.character(body(match.fun(fun)))[2], "of", names(gdo)[5])
    else ylab4 = ylab[4]
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))
    .aip(gdo[, 4], gdo[, 3], response = gdo[, 5], xlab = xlab4, ylab = ylab4, main = main4, col = col, type = "b", title = names(gdo)[3], ...)
    par(mar = old.par)
    D3 = c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.284, 0.308, 0.329, 0.348)
    D4 = c(0, 3.267, 2.574, 2.282, 2.115, 2.004, 1.924, 1.864, 1.816, 1.777, 1.744, 1.716, 1.692, 1.671, 1.652)
    helpRange = function(x) {
      return(diff(range(x)))
    }
    aggForLimits = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = helpRange)
    Rm = mean(aggForLimits[, 3])
    UCL = D4[sgSize] * Rm
    LCL = D3[sgSize] * Rm
    agg = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = helpRange)
    tab = table(agg[, 2])
    sgSize = tab[1]
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))
    plot(agg[, 3], ylim = c(0, max(max(agg[, 3]), UCL)), type = "n", xlab = aName, ylab = "R", axes = FALSE, main = "R Chart")
    axis(2)
    axis(4, at = c(Rm, UCL, LCL), labels = c("", "", ""))
    box()
    abline(h = Rm, col = 3)
    abline(h = UCL, col = 2)
    abline(h = LCL, col = 2)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, Rm, substitute(bar(R) == Rm, list(Rm = round(Rm, 2))), adj = 0, srt = 0, xpd = TRUE)
    j = 1
    for (i in 1:length(tab)) {
      lines(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
      points(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
      if (i < length(tab)) 
        abline(v = j + tab[i] - 1 + 0.5, lty = 2)
      axis(1, at = j, labels = names(tab[i]))
      j = j + tab[i]
    }
    par(mar = old.par)
  }
  if(gdo@method == "nested")
  {
    main2 = NA
    if (missing(main) || is.na(main[2])) 
      main2 = paste(yName, "By", bName, "Within", aName)
    else main2 = main[2]
    xlab2 = NA
    if (missing(xlab) || is.na(xlab[2])) 
      xlab2 = NA
    else xlab2 = xlab[2]
    ylab2 = NA
    if (missing(ylab) || is.na(ylab[2])) 
      ylab2 = yName
    else ylab2 = ylab[2]
    agg = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = function(x) {
      return(x)
    })
    plot(1:nrow(agg), main = main2, xlab = xlab2, ylab = ylab2, ylim = range(agg[, 3]), axes = FALSE)
    axis(2)
    box()
    label2 = ""
    for (i in 1:nrow(agg)) {
      points(rep(i, length(agg[i, 3])), agg[i, 3])
      axis(1, at = i, labels = agg[i, 1])
      if (agg[i, 2] != label2) {
        axis(1, at = i, labels = agg[i, 2], line = 1, tick = FALSE)
        label2 = agg[i, 2]
      }
    }
    aggm = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = mean)
    lines(aggm[, 3])
    points(aggm[, 3], pch = 13, cex = 2)
    main3 = NA
    if (missing(main) || is.na(main[3])) 
      main3 = paste(yName, "by", aName)
    else main3 = main[3]
    xlab3 = NA
    if (missing(xlab) || is.na(xlab[3])) 
      xlab3 = aName
    else xlab3 = xlab[3]
    ylab3 = NA
    if (missing(ylab) || is.na(ylab[3])) 
      ylab3 = yName
    else ylab3 = ylab[3]
    colVec = .mapping(gdo[, 3], sort(unique(gdo[, 3])), col[1:length(unique(gdo[, 3]))])
    boxplot(split(gdo[, yName], gdo[, aName]), col = colVec, xlab = xlab3, ylab = ylab3, main = main3, ...)
    mByOp = split(gdo[, 5], as.numeric(gdo[, 3]))
    lines(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, mean)[sort(names(mByOp))], lwd = lwd)
    points(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, mean)[sort(names(mByOp))], lwd = lwd, pch = 13, cex = 2)
    
    agg = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = mean)
    tab = table(agg[, 2])
    sgSize = tab[1]
    aggSd = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = sd)
    tab = table(aggSd[, 2])
    sm = mean(aggSd[, 3])
    aggMean = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = mean)
    xm = mean(agg[, 3])
    UCL = xm + ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
    LCL = xm - ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
    values = c(UCL, xm, LCL)
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))
    plot(agg[, 3], type = "n", axes = FALSE, xlab = aName, ylab = expression(bar(x)), main = expression(paste(bar(x), " Chart")))
    box()
    abline(h = xm, col = 3)
    abline(h = UCL, col = 2)
    abline(h = LCL, col = 2)
    axis(2)
    axis(4, at = c(xm, UCL, LCL), labels = c("", "", ""))
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, xm, substitute(bar(x) == xm, list(xm = round(xm, 2))), adj = 0, srt = 0, xpd = TRUE)
    par(mar = old.par)
    j = 1
    for (i in 1:length(tab)) {
      lines(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
      points(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
      if (i < length(tab)) 
        abline(v = j + tab[i] - 1 + 0.5, lty = 2)
      axis(1, at = j, labels = names(tab[i]))
      j = j + tab[i]
    }
    
    par(mar = old.par)
    D3 = c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.284, 0.308, 0.329, 0.348)
    D4 = c(0, 3.267, 2.574, 2.282, 2.115, 2.004, 1.924, 1.864, 1.816, 1.777, 1.744, 1.716, 1.692, 1.671, 1.652)
    helpRange = function(x) {
      return(diff(range(x)))
    }
    aggForLimits = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = helpRange)
    Rm = mean(aggForLimits[, 3])
    UCL = D4[sgSize] * Rm
    LCL = D3[sgSize] * Rm
    agg = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = helpRange)
    tab = table(agg[, 2])
    sgSize = tab[1]
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))
    plot(agg[, 3], ylim = c(0, max(max(agg[, 3]), UCL)), type = "n", xlab = aName, ylab = "R", axes = FALSE, main = "R Chart")
    axis(2)
    axis(4, at = c(Rm, UCL, LCL), labels = c("", "", ""))
    box()
    abline(h = Rm, col = 3)
    abline(h = UCL, col = 2)
    abline(h = LCL, col = 2)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, Rm, substitute(bar(R) == Rm, list(Rm = round(Rm, 2))), adj = 0, srt = 0, xpd = TRUE)
    j = 1
    for (i in 1:length(tab)) {
      lines(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
      points(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
      if (i < length(tab)) 
        abline(v = j + tab[i] - 1 + 0.5, lty = 2)
      axis(1, at = j, labels = names(tab[i]))
      j = j + tab[i]
      
    }
    
    
    
    #        plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
    #        plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
    plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
  }
}) 


# -------------- mul_t.r -------------- 
.xlimcalc = function(x) {
  r = range(x)
  d = diff(range(x))
  xlim = c(r[1] - 0.04 * d, r[2] + 0.04 * d)
  return(xlim)
}
.mapping = function(x, oldValues, newValues) {
  if (length(oldValues) != length(newValues)) {
    print("old and new")
    print(oldValues)
    print(newValues)
    warning(paste("unequal length of", deparse(substitute(oldValues)), "and", deparse(substitute(newValues))))
  }
  out = numeric(length(x))
  for (i in seq(along = newValues)) {
    index = (x == oldValues[i])
    out[index] = newValues[i]
  }
  return(out)
}
mvPlot = function(y, factor1, factor2, factor3, factor4, fun = mean, points = TRUE, connect = TRUE, col = c(1, 2, 3, 4), pch = c(1, 2, 3, 4), xlim, 
                  ylim, main, main.sub, horiz = FALSE, lwd.b = 1, lwd.w = 1, pch.b = 15, pch.w = 17, col.w = 2, col.b = 1, ...) {
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  if (missing(ylim)) 
    ylim = range(y)
  if (missing(xlim)) 
    xlim = NA
  if (missing(factor4)) 
    if (missing(factor3)) {
      if (!is.vector(factor1) | !is.vector(factor2)) 
        stop(paste(deparse(substitute(factor1)), "and", deparse(substitute(factor2)), "must be vectors!"))
      if (missing(main.sub)) 
        main.sub = ""
      temp = list(factor1)
      names(temp) = deparse(substitute(factor1))
      factor1 = temp
      temp = list(factor2)
      names(temp) = deparse(substitute(factor2))
      factor2 = temp
      temp = list(y)
      names(temp) = deparse(substitute(y))
      y = temp
      print(factor2)
      print(factor1)
      .mv2Plot(y = y, factor1, factor2, fun = fun, points = points, connect = connect, col = col, pch = pch, xlim = xlim, ylim = ylim, main, main.sub = main.sub, 
               horiz = horiz, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
    }
  else {
    if (!is.vector(factor1) | !is.vector(factor2) | !is.vector(factor3)) 
      stop(paste(deparse(substitute(factor1)), "and", deparse(substitute(factor2)), "and", deparse(substitute(factor3)), "must be vectors!"))
    temp = list(factor1)
    names(temp) = deparse(substitute(factor1))
    factor1 = temp
    temp = list(factor2)
    names(temp) = deparse(substitute(factor2))
    factor2 = temp
    temp = list(factor3)
    names(temp) = deparse(substitute(factor3))
    factor3 = temp
    temp = list(y)
    names(temp) = deparse(substitute(y))
    y = temp
    .mv3Plot(y, factor1, factor2, factor3, fun = fun, points = points, connect = connect, col = col, pch = pch, xlim = xlim, ylim, horiz = horiz, main, 
             main.sub, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
  }
  else {
    if (!is.vector(factor1) | !is.vector(factor2) | !is.vector(factor3) | !is.vector(factor4)) 
      stop(paste(deparse(substitute(factor1)), "and", deparse(substitute(factor2)), "and", deparse(substitute(factor3)), "and", deparse(substitute(factor4)), 
                 "must be vectors!"))
    temp = list(factor1)
    names(temp) = deparse(substitute(factor1))
    factor1 = temp
    temp = list(factor2)
    names(temp) = deparse(substitute(factor2))
    factor2 = temp
    temp = list(factor3)
    names(temp) = deparse(substitute(factor3))
    factor3 = temp
    temp = list(factor4)
    names(temp) = deparse(substitute(factor4))
    factor4 = temp
    temp = list(y)
    names(temp) = deparse(substitute(y))
    y = temp
    .mv4Plot(y, factor1, factor2, factor3, factor4, fun = fun, points = points, connect = connect, col = col, pch = pch, xlim = xlim, ylim, horiz = horiz, 
             main, main.sub, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
  }
  invisible()
}
.mv4Plot = function(y, factor1, factor2, factor3, factor4, fun, points, connect, col, pch, xlim, ylim, horiz, main, main.sub, lwd.b, lwd.w, pch.b, 
                    pch.w, col.w, col.b, DB = FALSE, ...) {
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  colNum = length(unique(factor1[[1]]))
  rowNum = length(unique(factor2[[1]]))
  par(mfrow = c(rowNum, colNum))
  plot.new()
  if (missing(ylim)) {
    if (points) {
      ylim = range(y)
    }
    else {
      if (DB) 
        print("TODO: Calculation of ylim")
      ylim = range(y)
    }
  }
  fullMat = data.frame(factor1[[1]], factor2[[1]], factor3[[1]], factor4[[1]], y[[1]])
  split1Mat = split(fullMat, factor1)
  if (DB) 
    print(str(split1Mat))
  for (i in seq(along = split1Mat)) {
    split12Mat = split(split1Mat[[i]], factor2[[1]])
    if (DB) 
      print(str(split12Mat))
    for (j in seq(along = split12Mat)) {
      temp = split12Mat[[j]]
      if (DB) {
        print(paste("split12Mat:", i, "and", j))
        print(temp)
      }
      par(mfg = c(i, j))
      constFac1 = unique(split1Mat[[i]][, 1])
      constFac2 = unique(split12Mat[[j]][, 2])
      if (missing(main.sub) || is.na(main.sub)) 
        main.sub = paste(names(factor1), "=", constFac1, "|", names(factor2), "=", constFac2)
      yPart = list(split1Mat[[j]][, 5])
      names(yPart) = names(y)
      f1Part = list(split1Mat[[j]][, 3])
      names(f1Part) = names(factor3)
      f2Part = list(split1Mat[[j]][, 4])
      names(f2Part) = names(factor4)
      .mv2Plot(y = yPart, factor1 = f1Part, factor2 = f2Part, fun = fun, points = points, connect = connect, main, main.sub, col = col, pch = pch, xlim = xlim, 
               ylim = ylim, horiz = horiz, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
      main.sub = NA
    }
  }
  invisible()
}
.mv3Plot = function(y, factor1, factor2, factor3, fun, points, connect, col, pch, xlim, ylim, main, main.sub, lwd.b, lwd.w, pch.b, pch.w, col.w, 
                    col.b, DB = FALSE, ...) {
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  colNum = length(unique(factor1[[1]]))
  par(mfcol = c(1, colNum))
  plot.new()
  if (missing(ylim)) {
    if (points) {
      ylim = range(y)
    }
    else {
      if (DB) 
        print("TODO: Calculation of ylim")
      ylim = range(y)
    }
  }
  fullMat = data.frame(factor1[[1]], factor2[[1]], factor3[[1]], y[[1]])
  split1Mat = split(fullMat, fullMat[, 1])
  if (DB) {
    print(fullMat)
    print(str(split1Mat))
  }
  for (i in seq(along = split1Mat)) {
    par(mfg = c(1, i))
    constFac1 = unique(split1Mat[[i]][, 1])
    if (missing(main.sub) || is.na(main.sub)) 
      main.sub = paste(names(factor1), "=", constFac1)
    yPart = list(split1Mat[[i]][, 4])
    names(yPart) = names(y)
    f1Part = list(split1Mat[[i]][, 2])
    names(f1Part) = names(factor2)
    f2Part = list(split1Mat[[i]][, 3])
    names(f2Part) = names(factor3)
    .mv2Plot(y = yPart, factor1 = f1Part, factor2 = f2Part, fun = fun, points = points, connect = connect, col = col, pch = pch, ylim = ylim, xlim = xlim, 
             main, main.sub = main.sub, lwd.b = lwd.b, lwd.w = lwd.w, pch.b = pch.b, pch.w = pch.w, col.w = col.w, col.b = col.b, ...)
    main.sub = NA
  }
  invisible()
}
.mv2Plot = function(y, factor1, factor2, fun, points, connect, horiz, main, main.sub, col, pch, xlim, ylim, xlab, ylab, cex = 1, lwd.b, lwd.w, pch.b, 
                    pch.w, col.w, col.b, DB = FALSE, ...) {
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  f1Name = NA
  f2Name = NA
  f1 = NA
  f2c = NA
  if (is.list(y)) {
    if (missing(ylab)) 
      ylab = names(y)
    y = y[[1]]
  }
  else {
    if (missing(ylab)) 
      ylab = deparse(substitute(y))
  }
  if (is.list(factor1)) {
    f1Name = names(factor1)
    factor1 = factor1[[1]]
    if (DB) 
      print(factor1)
  }
  else f1Name = deparse(substitute(factor1))
  if (is.list(factor2)) {
    if (missing(xlab)) 
      xlab = names(factor2)
    f2Name = names(factor2)
    factor2 = factor2[[1]]
    if (DB) 
      print(factor2)
  }
  else {
    if (missing(xlab)) 
      xlab = deparse(substitute(factor2))
  }
  if (missing(main)) 
    main = paste("Multi Vari Plot for", f1Name, "and", f2Name)
  levFac1 = sort(unique(factor1))
  levFac2 = sort(unique(factor2))
  if (missing(ylim)) 
    ylim = range(y)
  if (missing(col)) 
    col = 1:length(levFac1)
  else {
    if (length(col) < length(levFac1)) 
      col = 1:length(levFac1)
    else col = col[1:length(levFac1)]
  }
  if (missing(pch)) 
    pch = 1:length(levFac1)
  else {
    if (length(pch) < length(levFac1)) 
      pch = 1:length(levFac1)
    else pch = pch[1:length(levFac1)]
  }
  if (DB) 
    print(factor2)
  numTimes = 0
  matComplete = data.frame()
  matWithin = data.frame()
  matBetween = data.frame()
  mat = data.frame(f1 = factor1, f2 = factor2, y = y)
  if (DB) 
    print(mat)
  for (i in seq(along = levFac1)) {
    numTimes = numTimes + 1
    matPart = subset(mat, f1 == levFac1[i])
    levFac2 = unique(matPart$f2)
    levFac2Coded = numTimes:(numTimes + length(levFac2) - 1)
    numTimes = rev(levFac2Coded)[1] + 1
    matPart = cbind(matPart, f2c = .mapping(matPart$f2, levFac2, levFac2Coded))
    matComplete = rbind(matComplete, matPart)
    names(matComplete) = names(matPart)
    matBetween = rbind(matBetween, data.frame(f2c = mean(levFac2Coded), yFun = with(matPart, fun(matPart$y))))
    names(matBetween) = c("f2c", "yFun")
    if (DB) {
      print(paste("levFac2:", levFac2))
      print(paste("levFac1", levFac1[i], "\n"))
      print(paste("i:", i))
      print(matPart)
      print(paste("levFac2Coded:", levFac2Coded))
    }
    for (j in seq(along = matPart$f2c)) {
      matSubPart = subset(matPart, f2c == f2c[j])
      yFun = fun(as.numeric(matSubPart$y))
      matWithin = rbind(matWithin, data.frame(f1 = levFac1[i], f2c = matPart$f2c[j], yFun = yFun))
    }
  }
  f1pch = with(matComplete, .mapping(f1, levFac1, pch))
  f1col = with(matComplete, .mapping(f1, levFac1, col))
  if (missing(xlim) || is.na(xlim)) 
    xlim = .xlimcalc(matComplete$f2c)
  if (horiz) {
    ylim = c(min(ylim), 1.15 * max(ylim))
    with(matComplete, plot(f2c, y, main = main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = f1col, pch = f1pch, axes = FALSE, ...))
    legend("top", title = f1Name, legend = levFac1, pch = pch, col = col, horiz = TRUE)
  }
  else {
    par(mar = c(5, 4, 4, 6) + 0.1)
    with(matComplete, plot(f2c, y, main = main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = f1col, pch = f1pch, axes = FALSE, ...))
    xCoord = par("usr")[2] + 0.02 * diff(range(par("usr")[2], par("usr")[1]))
    legend(xCoord, par("usr")[4], title = f1Name, legend = levFac1, pch = pch, col = col, xpd = TRUE)
  }
  title(main.sub, line = 0.5)
  axis(2, ...)
  axis(1, at = matComplete$f2c, labels = matComplete$f2, ...)
  box()
  with(matWithin, points(f2c, yFun, col = col.w, pch = pch.w))
  if (DB) 
    print("WATCH HERE")
  for (i in seq(along = levFac1)) {
    subMat = subset(matWithin, f1 == levFac1[i])
    temp = duplicated(subMat)
    with(subMat[!temp, ], lines(f2c, yFun, col = col.w, lwd = lwd.w))
    if (DB) 
      print(subMat[!temp, ])
  }
  with(matBetween, points(f2c, yFun, col = col.b, pch = pch.b))
  with(matBetween, lines(f2c, yFun, col = col.b, lwd = lwd.b))
  if (DB) 
    print(matComplete)
  invisible()
} 

# -------------- internals.r -------------- 
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

# -------------- tag_e.r -------------- 
.NAMES = LETTERS[c(1:8, 10:26)]
setClass(Class = "taguchiDesign", representation = representation(name = "character", factors = "list", design = "data.frame", designType = "character", 
                                                                  replic = "data.frame", response = "data.frame", Type = "data.frame", block = "data.frame", runOrder = "data.frame", standardOrder = "data.frame", desireVal = "list", 
                                                                  desirability = "list", fits = "data.frame"))
setClass("taguchiFactor", representation = representation(values = "ANY", name = "character", unit = "character", type = "character"), prototype = prototype(values = NA, 
                                                                                                                                                             name = " ", unit = " ", type = "numeric"))
setGeneric("values", function(object) standardGeneric("values"))
setGeneric("values<-", function(object, value) standardGeneric("values<-"))
setMethod("values", "taguchiDesign", function(object) {
  listOut = vector(mode = "list")
  for (i in names(object@design)) {
    listOut[[i]] = .values(object@factors[[i]])
  }
  return(listOut)
})
setReplaceMethod("values", "taguchiDesign", function(object, value) {
  for (i in names(value)) {
    if (i %in% names(object@design)) 
      if (length(value[[i]]) == length(unique(object@design[, i]))) 
        .values(object@factors[[i]]) = value[[i]]
    else stop("Number of values greater or less than number of factor settings!")
  }
  return(object)
})
setGeneric(".values", function(object) standardGeneric(".values"))
setGeneric(".values<-", function(object, value) standardGeneric(".values<-"))
setMethod(".values", "taguchiFactor", function(object) object@values)
setReplaceMethod(".values", "taguchiFactor", function(object, value) {
  object@values <- value
  return(object)
})
setMethod(".unit", "taguchiFactor", function(object) object@unit)
setReplaceMethod(".unit", "taguchiFactor", function(x, value) {
  x@unit <- value
  x
})
setMethod("units", "taguchiDesign", function(x) {
  return(sapply(factors(x), .unit))
})
setMethod("units<-", "taguchiDesign", function(x, value) {
  for (i in 1:length(x@factors)) if (length(value) > 1) 
    .unit(x@factors[[i]]) = as.character(value[i])
  else .unit(x@factors[[i]]) = as.character(value[1])
  x
})
setMethod("factors", "taguchiDesign", function(x) x@factors)
setReplaceMethod("factors", "taguchiDesign", function(x, value) {
  if (length(value) != ncol(x@design)) 
    stop("\nNumber of factors doesn't match with number of columns for factorial Design\n")
  x@factors <- value
  x
})
setMethod("names", "taguchiDesign", function(x) {
  return(sapply(x@factors, names))
})
setReplaceMethod("names", "taguchiDesign", function(x, value) {
  for (i in 1:length(x@factors)) names(x@factors[[i]]) = as.character(value[i])
  x
})
setMethod("names", "taguchiFactor", function(x) {
  x@name
})
setReplaceMethod("names", "taguchiFactor", function(x, value) {
  x@name <- value
  x
})
setMethod("as.data.frame", "taguchiDesign", function(x, row.names = NULL, optional = FALSE) {
  frameOut = cbind(x@standardOrder, x@runOrder, x@replic, x@design, x@response)
  return(frameOut)
})
as.data.frame.taguchiDesign = function(x, row.names = NULL, optional = FALSE, ...) {
  frameOut = cbind(x@standardOrder, x@runOrder, x@replic, x@design, x@response)
  return(frameOut)
}
setMethod("show", signature(object = "taguchiDesign"), function(object) {
  print(format(as.data.frame(object), digits = 4))
})
setMethod("response", "taguchiDesign", function(object) {
  return(object@response)
})
setReplaceMethod("response", "taguchiDesign", function(object, value) {
  #    print(deparse(substitute(value)))                                          ###
  if (!is.numeric(value) & !is.data.frame(value)) 
    stop("vector or data.frame must be given")
  if (is.numeric(value)) {
    if (length(value) != nrow(object@design)) 
      stop("differing lengths")
    temp = data.frame(value)
    names(temp) = deparse(substitute(value))[1]
    value = temp
  }
  if (is.data.frame(value)) {
    if (nrow(value) != nrow(object@design)) 
      stop("differing number of rows")
  }
  object@response = value
  return(object)
})
setMethod(".nfp", "taguchiDesign", function(object) {
  x = factors(object)
  DB = FALSE
  if (is.list(x) && length(x[[1]]) > 0) {
    numAttr = length(attributes(x[[1]])) - 1
    .numFac = length(x)
    len = 0
    for (i in names(x)) if (length(x[[i]]@values) > len) 
      len = length(x[[i]]@values)
    numAttr = numAttr + len
    numrows = numAttr - 1
    frameOut = data.frame(matrix(NA, ncol = .numFac, nrow = numrows))
    names(frameOut) = names(x)
    rownames(frameOut) = c(paste("value", 1:len), "name", "unit", "type")
    for (i in names(x)) {
      vin = 1:length(x[[i]]@values)
      frameOut[vin, i] = x[[i]]@values
      frameOut[numrows - 2, i] = x[[i]]@name
      frameOut[numrows - 1, i] = x[[i]]@unit
      frameOut[numrows, i] = x[[i]]@type
    }
    print(frameOut)
  }
})
setMethod("show", signature(object = "taguchiDesign"), function(object) {
  print(as.data.frame(object))
})
setMethod("summary", signature(object = "taguchiDesign"), function(object) {
  cat(paste("Taguchi", toupper(object@designType), "Design"))
  cat("\n")
  cat("Information about the factors:\n\n")
  .nfp(object)
  cat("\n")
  cat("-----------\n")
  cat("\n")
  print(as.data.frame(object))
  cat("\n")
  cat("-----------\n")
  cat("\n")
})
taguchiDesign = function(design, randomize = TRUE, replicates = 1) {
  DB = FALSE
  odo = NA
  type = "single"
  for (i in seq(along = .oaList)) {
    pmatch(design, .oaList[[i]]$id)
    if (!is.na(pmatch(design, .oaList[[i]]$id))) {
      if (DB) 
        print(design)
      temp = .oaList[[i]]
      design = temp$design
      repVec = rep(1, nrow(design))
      if (replicates > 1) {
        X = temp$design
        for (i in 1:(replicates - 1)) {
          design = rbind(design, X)
          repVec = c(repVec, rep(i + 1, times = nrow(X)))
        }
      }
      Replicate = data.frame(Replicate = as.numeric(repVec))
      if (DB) 
        print(Replicate)
      odo = new("taguchiDesign")
      odo@design = design
      names(odo@design) = .NAMES[1:ncol(design)]
      odo@name = temp$id
      odo@designType = temp$type
      odo@replic = Replicate
      StandOrder = 1:nrow(odo@design)
      RunOrder = StandOrder
      if (randomize) {
        RunOrder = sample(1:nrow(odo@design), nrow(odo@design), replace = FALSE, prob = NULL)
      }
      odo@design = odo@design[order(RunOrder), ]
      odo@replic = data.frame(Replicate = odo@replic[order(RunOrder), 1])
      row.names(odo@design) = odo@design$RunOrder
      odo@runOrder = data.frame(RunOrder = data.frame(RunOrder = RunOrder)[order(RunOrder), ])
      odo@standardOrder = data.frame(StandOrder = data.frame(StandOrder = StandOrder)[order(RunOrder), ])
      odo@response = data.frame(y = rep(NA, nrow(odo@design)))
      tfList = vector("list", ncol(design))
      for (i in seq(along = tfList)) tfList[i] = new("taguchiFactor")
      names(tfList) = names(odo@design)
      factors(odo) = tfList
      valList = list(length = length(names(odo)))
      for (i in names(names(odo))) valList[[i]] = sort(unique(odo@design[, i]))
      values(odo) = valList
      return(odo)
    }
  }
  return(NA)
}
oaChoose = function(factors1, factors2, level1, level2, ia) {
  params = list(factors1 = 0, factors2 = 0, level1 = 0, level2 = 0, ia = 0)
  if (!missing(ia)) 
    params$ia = ia
  if (!missing(factors2)) 
    params$factors2 = factors2
  if (!missing(level2)) 
    params$level2 = level2
  do.call(taguchiChoose, params)
}
taguchiChoose = function(factors1 = 0, factors2 = 0, level1 = 0, level2 = 0, ia = 0) {
  if (factors1 == 0 & factors2 == 0 & level1 == 0 & level2 == 0 & ia == 0) {
    temp = vector(mode = "character", length = length(.oaList))
    for (i in 1:length(.oaList)) temp[i] = .oaList[[i]]$id
    temp = c(temp, rep(" ", (length(temp)%/%6 + 1) * 6 - length(temp)))
    mat = data.frame(matrix(temp, ncol = 6, byrow = TRUE))
    names(mat) = rep(" ", ncol(mat))
    print(mat)
    cat("\n")
    cat("Choose a design using e.g. taguchiDesign(\"L4_2\")")
    cat("\n")
  }
  else {
    DB = FALSE
    if (factors2 <= 0) 
      level2 = 0
    if (DB) 
      print(level2)
    Anzahl_Spalten = factors1 + factors2 + ia
    ss = list()
    for (i in seq(along = .oaList)) {
      li = .oaList[[i]]
      if (li$factors1 >= factors1 & li$factors2 >= factors2 & (li$levels1 == level1 | li$levels1 == level2) & (li$levels2 == level2 | li$levels2 == level1) & 
          li$anzahl_spalten >= Anzahl_Spalten) 
        ss[i] = li$id
    }
    out = as.character(ss)
    out = out[out != "NULL"]
    if (length(out) > 0) {
      cat(paste(factors1, "factors on", level1, "levels and", factors2, "factors on", level2, "levels with", ia, "desired interactions to be estimated\n"))
      cat("\n")
      cat("Possible Designs:\n")
      cat("\n")
      cat(paste(out, sep = " | "))
      cat("\n")
      cat("\n")
      cat(paste("Use taguchiDesign(\"", out[1], "\") or different to create a taguchi design object\n", sep = ""))
    }
    else {
      cat("No Design Found\n")
      cat("\n")
      out = NA
    }
    invisible(out)
  }
}
.replace2s = function(x) {
  if (!is.data.frame(x)) 
    stop(paste(deparse(substitute(x)), "needs to be a data.frame"))
  for (i in 1:ncol(x)) x[x[, i] == 2, i] = -1
  return(x)
}
.helpAliasTable = function(fdo, k, degree = 3) {
  if (degree > k) {
    degree = k
  }
  if (class(fdo) == "facDesign") 
    X = unique(cube(fdo))
  if (class(fdo) == "taguchiDesign") {
    X = unique(fdo@design)
    X = .replace2s(X)
  }
  N = nrow(X)
  columns = names(X[, 1:k])
  X1 = matrix(1, nrow = N, ncol = 1)
  nameVec = c("Identity")
  for (i in 1:degree) {
    temp = combn(columns, i)
    for (j in 1:ncol(temp)) {
      if (class(fdo) == "facDesign") 
        index = names(names(fdo)) %in% temp[, j]
      if (class(fdo) == "taguchiDesign") 
        index = names(X) %in% temp[, j]
      if (length((1:length(index))[index]) == 1) {
        X1 = cbind(X1, X[, index])
        nameVec = c(nameVec, temp[, j])
      }
      else {
        X1 = cbind(X1, apply(X[, index], 1, prod))
        nameVec = c(nameVec, paste(temp[, j], sep = "", collapse = ""))
      }
    }
    X1 = data.frame(X1)
    names(X1) = nameVec
  }
  return(X1)
}
aliasTable = function(fdo, degree, show = TRUE) {
  if (class(fdo) == "facDesign") {
    X = unique(cube(fdo))
    N = nrow(X)
    k = log2(N)
    kPlusP = ncol(X)
    if (missing(degree)) 
      degree = min(c(4, k + 1))
    X1 = .helpAliasTable(fdo, k, degree = degree - 1)
    X2 = .helpAliasTable(fdo, k = kPlusP, degree)
  }
  if (class(fdo) == "taguchiDesign") {
    if (length(table(as.numeric(as.matrix(fdo@design)))) != 2) 
      stop("calculation of an alias table for mixed designs is not supported")
    k = ncol(fdo@design)
    if (missing(degree)) 
      degree = min(c(3, k))
    X1 = unique(fdo@design)
    X1 = .replace2s(X1)
    X2 = .helpAliasTable(fdo, k, degree)
    X1 = cbind(data.frame(Identity = rep(1, times = nrow(X1))), X1)
  }
  logVec = !(names(X2) %in% names(X1))
  X2 = X2[, logVec]
  X1 = as.matrix(X1)
  X2 = as.matrix(X2)
  alias.matrix = solve(t(X1) %*% X1) %*% t(X1) %*% X2
  if (show) 
    print(round(alias.matrix, 2))
  invisible(alias.matrix)
}
setMethod("identity", signature(x = "taguchiDesign"), function(x) {
  identity = character(0)
  identityList = vector(mode = "list", length = 0)
  resolution = numeric(0)
  temp = NULL
  A = aliasTable(x)
  if (any(dim(A) == 0)) 
    return(identityList)
  temp = as.matrix(A["Identity", ])
  boolTemp = apply(temp, 2, as.logical)
  identity = row.names(temp)[boolTemp[, 1]]
  if (length(identity) > 0) {
    charList = strsplit(toupper(identity), split = "")
    identityList = lapply(charList, match, LETTERS[1:26])
    names(identityList) = identity
  }
  cat("Defining relations:\n")
  if (length(identityList) > 0) {
    for (i in 1:length(identityList)) {
      identLen = length((strsplit(names(identityList)[i], split = character(0))[[1]]))
      if (length(resolution) == 0 || identLen > resolution) 
        resolution = c(resolution, identLen)
      cat("I = ", names(identityList)[i], "\t\tColumns:", identityList[[i]], "\n")
    }
    cat("\nResolution: ", as.character(as.roman(min(resolution))), "\n")
  }
  invisible(identityList)
}) 

####
#.NAMES = LETTERS[c(1:8, 10:26)]
setClass(Class = "pbDesign", representation = representation(name = "character", factors = "list", design = "data.frame", designType = "character",
                                                             replic = "data.frame", response = "data.frame", Type = "data.frame", block = "data.frame", runOrder = "data.frame", standardOrder = "data.frame", desireVal = "list",
                                                             desirability = "list", fits = "data.frame"))

setClass("pbFactor", representation = representation(values = "ANY", name = "character", unit = "character", type = "character"), prototype = prototype(values = NA,
                                                                                                                                                        name = " ", unit = " ", type = "numeric"))

#setGeneric("values", function(object) standardGeneric("values"))
#setGeneric("values<-", function(object, value) standardGeneric("values<-"))
setMethod("values", "pbDesign", function(object) {
  listOut = vector(mode = "list")
  for (i in names(object@design)) {
    listOut[[i]] = .values(object@factors[[i]])
  }
  return(listOut)
})
setReplaceMethod("values", "pbDesign", function(object, value) {
  for (i in names(value)) {
    if (i %in% names(object@design))
      if (length(value[[i]]) == length(unique(object@design[, i])))
        .values(object@factors[[i]]) = value[[i]]
    else stop("Number of values greater or less than number of factor settings!")
  }
  return(object)
})

#setGeneric(".values", function(object) standardGeneric(".values"))
#setGeneric(".values<-", function(object, value) standardGeneric(".values<-"))
setMethod(".values", "pbFactor", function(object) object@values)
setReplaceMethod(".values", "pbFactor", function(object, value) {
  object@values <- value
  return(object)
})
setMethod(".unit", "pbFactor", function(object) object@unit)
setReplaceMethod(".unit", "pbFactor", function(x, value) {
  x@unit <- value
  x
})
setMethod("names", "pbFactor", function(x) {
  x@name
})
setReplaceMethod("names", "pbFactor", function(x, value) {
  x@name <- value
  x
})

setMethod("units", "pbDesign", function(x) {
  return(sapply(factors(x), .unit))
})
setMethod("units<-", "pbDesign", function(x, value) {
  for (i in 1:length(x@factors)) if (length(value) > 1)
    .unit(x@factors[[i]]) = as.character(value[i])
  else .unit(x@factors[[i]]) = as.character(value[1])
  x
})
setMethod("factors", "pbDesign", function(x) x@factors)
setReplaceMethod("factors", "pbDesign", function(x, value) {
  if (length(value) != ncol(x@design))
    stop("\nNumber of factors doesn't match with number of columns for factorial Design\n")
  x@factors <- value
  x
})
setMethod("names", "pbDesign", function(x) {
  return(sapply(x@factors, names))
})
setReplaceMethod("names", "pbDesign", function(x, value) {
  for (i in 1:length(x@factors)) names(x@factors[[i]]) = as.character(value[i])
  x
})
setMethod("names", "pbFactor", function(x) {
  x@name
})

setMethod("as.data.frame", "pbDesign", function(x, row.names = NULL, optional = FALSE) {
  frameOut = cbind(x@standardOrder, x@runOrder, x@replic, x@design, x@response)
  return(frameOut)
})
as.data.frame.pbDesign = function(x, row.names = NULL, optional = FALSE, ...) {
  frameOut = cbind(x@standardOrder, x@runOrder, x@replic, x@design, x@response)
  return(frameOut)
}
setMethod("show", signature(object = "pbDesign"), function(object) {
  print(format(as.data.frame(object), digits = 4))
})
setMethod("response", "pbDesign", function(object) {
  return(object@response)
})
setReplaceMethod("response", "pbDesign", function(object, value) {
  #    print(deparse(substitute(value)))                                          ###
  if (!is.numeric(value) & !is.data.frame(value))
    stop("vector or data.frame must be given")
  if (is.numeric(value)) {
    if (length(value) != nrow(object@design))
      stop("differing lengths")
    temp = data.frame(value)
    names(temp) = deparse(substitute(value))[1]
    value = temp
  }
  if (is.data.frame(value)) {
    if (nrow(value) != nrow(object@design))
      stop("differing number of rows")
  }
  object@response = value
  return(object)
})
setMethod(".nfp", "pbDesign", function(object) {
  x = factors(object)
  DB = FALSE
  if (is.list(x) && length(x[[1]]) > 0) {
    numAttr = length(attributes(x[[1]])) - 1
    .numFac = length(x)
    len = 0
    for (i in names(x)) if (length(x[[i]]@values) > len)
      len = length(x[[i]]@values)
    numAttr = numAttr + len
    numrows = numAttr - 1
    frameOut = data.frame(matrix(NA, ncol = .numFac, nrow = numrows))
    names(frameOut) = names(x)
    rownames(frameOut) = c(paste("value", 1:len), "name", "unit", "type")
    for (i in names(x)) {
      vin = 1:length(x[[i]]@values)
      frameOut[vin, i] = x[[i]]@values
      frameOut[numrows - 2, i] = x[[i]]@name
      frameOut[numrows - 1, i] = x[[i]]@unit
      frameOut[numrows, i] = x[[i]]@type
    }
    print(frameOut)
  }
})
setMethod("show", signature(object = "pbDesign"), function(object) {
  print(as.data.frame(object))
})
setMethod("summary", signature(object = "pbDesign"), function(object) {
  cat(paste("Plackett-Burman", toupper(object@designType), "Design"))
  cat("\n")
  cat("Information about the factors:\n\n")
  .nfp(object)
  cat("\n")
  cat("-----------\n")
  cat("\n")
  print(as.data.frame(object))
  cat("\n")
  cat("-----------\n")
  cat("\n")
})


.pbDesign=function(n)
{
  k=n-1
  if(k==27)
  {
    X=matrix(c(+1, -1, +1, +1, +1, +1, -1, -1, -1,
               +1, +1, -1, +1, +1, +1, -1, -1, -1,
               -1, +1, +1, +1, +1, +1, -1, -1, -1,
               -1, -1, -1, +1, -1, +1, +1, +1, +1,
               -1, -1, -1, +1, +1, -1, +1, +1, +1,
               -1, -1, -1, -1, +1, +1, +1, +1, +1,
               +1, +1, +1, -1, -1, -1, +1, -1, +1,
               +1, +1, +1, -1, -1, -1, +1, +1, -1,
               +1, +1, +1, -1, -1, -1, -1, +1, +1),nrow=9,ncol=9,byrow=TRUE)
    Y=matrix(c(-1, +1, -1, -1, -1, +1, -1, -1, +1,
               -1, -1, +1, +1, -1, -1, +1, -1, -1,
               +1, -1, -1, -1, +1, -1, -1, +1, -1,
               -1, -1, +1, -1, +1, -1, -1, -1, -1,
               +1, -1, -1, -1, -1, +1, +1, -1, -1,
               -1, +1, -1, +1, -1, -1, -1, +1, -1,
               -1, -1, +1, -1, -1, +1, -1, +1, -1,
               +1, -1, -1, +1, -1, -1, -1, -1, +1,
               -1, +1, -1, -1, +1, -1, +1, -1, -1),nrow=9,ncol=9,byrow=TRUE)
    Z=matrix(c(+1, +1, -1, +1, -1, +1, +1, -1, +1,
               -1, +1, +1, +1, +1, -1, +1, +1, -1,
               +1, -1, +1, -1, +1, +1, -1, +1, +1,
               +1, -1, +1, +1, +1, -1, +1, -1, +1,
               +1, +1, -1, -1, +1, +1, +1, +1, -1,
               -1, +1, +1, +1, -1, +1, -1, +1, +1,
               +1, -1, +1, +1, -1, +1, +1, +1, -1,
               +1, +1, -1, +1, +1, -1, -1, +1, +1,
               -1, +1, +1, -1, +1, +1, +1, -1, +1),nrow=9,ncol=9,byrow=TRUE)
    design=data.frame(rbind(cbind(X,Y,Z),cbind(X,Y,Z),cbind(X,Y,Z),rep(-1,27)))
  }
  else
  {
    if(k<3)
      stop("k needs to be grater than three!")
    if(k==3)
      firstRow=c(+1, -1, +1)[1:k]
    if(k>3 && k<=7)
      firstRow=c(+1, +1, +1, -1, +1, -1, -1)[1:k]
    if(k>=8 && k<=11)
      firstRow=c(+1, +1, -1, +1, +1, +1, -1, -1, -1, +1, -1)[1:k]
    if(k>=12 && k<=15)
      firstRow=c(+1, +1, +1, +1, -1, +1, -1, +1, +1, -1, -1, +1, -1, -1, -1)[1:k]
    if(k>=16 && k<=19)
      firstRow=c(+1, +1, -1, -1, +1, +1, +1, +1, -1, +1, -1, +1, -1, -1, -1, -1, +1, +1, -1)[1:k]
    if(k>=20 && k<=23)
      firstRow=c(+1, +1, +1, +1, +1, -1, +1, -1, +1, +1, -1, -1, +1, -1, -1, +1, -1, +1, -1, -1, -1, -1)[1:k]
    if(k>=24 && k<=26)
      firstRow=c(-1, -1, -1, -1, +1, -1, +1, -1, +1, +1, +1, -1, +1, +1, -1, -1, -1, +1, +1, +1, +1, +1, -1, -1, +1, +1, -1, +1, -1, -1, +1)[1:k]
    if(k==27)
      print("insert exception here!")
    if(k>=29 && k<=31)
      firstRow=c(-1, -1, -1, -1, +1, -1, +1, -1, +1, +1, +1, -1, +1, +1, -1, -1, -1, +1, +1, +1, +1, +1, -1, -1, +1, +1, -1, +1, -1, -1, +1)[1:k]
    if(k>=32 && k<=35)
      firstRow=c(-1, +1, -1, +1, +1, +1, -1, -1, -1, +1, +1, +1, +1, +1, -1, +1, +1, +1, -1, -1, +1, -1, -1, -1, -1, +1, -1, +1, -1, +1, +1, -1, -1, +1, -1)[1:k]
    if(k>=36 && k<=39)
      firstRow=c(+1, +1, -1, -1, +1, +1, +1, +1, -1, +1, -1, +1, -1, -1, -1, -1, +1, +1, -1, -1, +1, +1, -1, -1, +1, +1, +1, +1, -1, +1, -1, +1, -1, -1, -1, -1, +1, +1, -1)[1:k]
    if(k>=40 && k<=43)
      firstRow=c(+1, +1, -1, -1, +1, -1, +1, -1, -1, +1, +1, +1, -1, +1, +1, +1, +1, +1, -1, -1, -1, +1, -1, +1, +1, +1, -1, -1, -1, -1, -1, +1, -1, -1, -1, +1, +1, -1, +1, -1, +1, +1, -1)[1:k]
    if(k>=44 && k<=47)
      firstRow=c(+1, +1, +1, +1, +1, -1, +1, +1, +1, +1, -1, -1, +1, -1, +1, -1, +1, +1, +1, -1, -1, +1, -1, -1, +1, +1, -1, +1, +1, -1, -1, -1, +1, -1, +1, -1, +1, +1, -1, -1, -1, -1, +1, -1, -1, -1, -1)[1:k]
    if(k>=48 && k<=59)
      firstRow=c(+1, +1, -1, +1, +1, +1, -1, +1, -1, +1, -1, -1, +1, -1, -1, +1, +1, +1, -1, +1, +1, +1, +1, -1, -1, +1, +1, +1, +1, +1, -1, -1, -1, -1, -1, +1, +1, -1, -1, -1, -1, +1, -1, -1, -1, +1, +1, -1, +1, +1, -1, +1, -1, +1, -1, -1, -1, +1, -1)[1:k]
    
    design=matrix(NA,nrow=k+1,ncol=k)
    design[1,]=firstRow
    nextRow=firstRow
    
    for(i in 2:k)
    {
      nextRow=nextRow[c(k,1:k-1)]
      design[i,]=nextRow
    }
    lastRow=rep(-1,k)
    design[k+1,]=lastRow
    design=data.frame(design)
  }
  return(design)
}

pbDesign = function(n, k , randomize = TRUE, replicates = 1) {
  if(missing(n)&&missing(k))
    stop("Either n or k must be set!")
  if(missing(n)==FALSE && missing(k)==FALSE && k!=n-1 )
    stop("Wrong combination of n and k")
  if(missing(n))
    n=k+1
  if(missing(k))
    k=n-1  
  DB = FALSE
  odo = NA
  if (DB)
    print(n)
  design = .pbDesign(n)
  repVec = rep(1, nrow(design))
  if (replicates > 1) {
    X = .pbDesign(n)
    for (i in 1:(replicates - 1)) {
      design = rbind(design, X)
      repVec = c(repVec, rep(i + 1, times = nrow(X)))
    }
  }
  Replicate = data.frame(Replicate = as.numeric(repVec))
  if (DB)
    print(Replicate)
  odo = new("pbDesign")
  odo@design = design
  #  names(odo@design) = .NAMES[1:ncol(design)]
  odo@replic = Replicate
  StandOrder = 1:nrow(odo@design)
  RunOrder = StandOrder
  if (randomize) {
    RunOrder = sample(1:nrow(odo@design), nrow(odo@design), replace = FALSE, prob = NULL)
  }
  odo@design = odo@design[order(RunOrder), ]
  odo@replic = data.frame(Replicate = odo@replic[order(RunOrder), 1])
  row.names(odo@design) = odo@design$RunOrder
  odo@runOrder = data.frame(RunOrder = data.frame(RunOrder = RunOrder)[order(RunOrder), ])
  odo@standardOrder = data.frame(StandOrder = data.frame(StandOrder = StandOrder)[order(RunOrder), ])
  odo@response = data.frame(y = rep(NA, nrow(odo@design)))
  tfList = vector("list", ncol(design))
  for (i in seq(along = tfList)) tfList[i] = new("pbFactor")
  names(tfList) = names(odo@design)
  factors(odo) = tfList
  valList = list(length = length(names(odo)))
  for (i in names(names(odo))) valList[[i]] = sort(unique(odo@design[, i]))
  values(odo) = valList
  return(odo)
}




















# -------------- dis_t.r -------------- 
setClass("distr", representation(x = "vector", name = "character", parameters = "numeric", sd = "numeric", n = "numeric", loglik = "numeric"))
setClass("distrCollection", representation(distr = "list"))
setMethod("[", signature(x = "distrCollection", i = "ANY"), function(x, i, drop = missing) {
  x@distr[i]
})
setMethod("show", signature(object = "distrCollection"), function(object) {
  distrList = object@distr
  cat("\n")
  for (i in seq(along = distrList)) {
    temp = distrList[[i]]
    cat("\n")
    cat("fitted distribution is", temp@name, ":\n")
    print(temp@parameters)
    cat("\n")
  }
})
setMethod("summary", signature(object = "distrCollection"), function(object) {
  numDist = length(object@distr)
  gofMatrix = data.frame(matrix(nrow = numDist, ncol = 3))
  names(gofMatrix) = c("Distribution", "A", "p.value")
  cat("\n------ Fitted Distribution and estimated parameters ------\n")
  for (i in seq(along = object@distr)) {
    distrObj = object@distr[[i]]
    x = distrObj@x
    distribution = distrObj@name
    parameters = distrObj@parameters
    statistic = NA
    p.value = NA
    temp = .myADTest(x, distribution)
    try(statistic <- as.numeric(temp$statistic), silent = TRUE)
    try(p.value <- as.numeric(temp$p.value), silent = TRUE)
    gofMatrix[i, ] = c(distribution, as.numeric(statistic), as.numeric(p.value))
    cat("\n")
    cat("fitted distribution is", distribution, ":\n")
    print(parameters)
  }
  cat("\n")
  cat("\n------ Goodness of Fit - Anderson Darling Test ------\n")
  cat("\n")
  gofMatrixPrint = gofMatrix
  gofMatrixPrint[, 2] = signif(as.numeric(gofMatrixPrint[, 2]), 4)
  gofMatrixPrint[, 3] = signif(as.numeric(gofMatrixPrint[, 3]), 4)
  print(gofMatrixPrint)
})
distribution = function(x, distribution = "weibull", start, ...) {
  #if (!require(MASS, quietly = TRUE)) 
  #   stop("Package MASS needs to be installed!")
  if (is.character(distribution)) 
    distribution = tolower(distribution)
  allDistr = c("beta", "cauchy", "chi-squared", "exponential", "f", "gamma", "geometric", "log-normal", "logistic", "negative binomial", "normal", "poisson", 
               "t", "weibull")
  if (distribution %in% allDistr) 
    distrVec = distribution
  else distrVec = c("normal")
  if (identical(distribution, "all")) 
    distrVec = allDistr
  if (identical(distribution, "quality")) 
    distrVec = c("normal", "log-normal", "exponential", "weibull")
  distrColl = new("distrCollection")
  for (i in seq(along = distrVec)) {
    fit = new("distr")
    temp = fitdistr(x, distrVec[i], ...)
    fit@x = x
    fit@name = distrVec[i]
    fit@parameters = temp$estimate
    fit@sd = temp$sd
    fit@loglik = temp$loglik
    distrColl@distr[i] = fit
  }
  return(distrColl)
}
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "distr"), function(x, y, main, xlab, xlim, ylim, ylab, line.col, line.width, box = TRUE, ...) {
  object = x
  xVals = object@x
  parameters = object@parameters
  lq = NULL
  uq = NULL
  y = NULL
  if (missing(line.col)) 
    line.col = "red"
  if (missing(line.width)) 
    line.width = 1
  if (missing(main)) 
    main = object@name
  if (missing(xlab)) 
    xlab = "x"
  if (missing(ylab)) 
    ylab = "Density"
  distr = object@name
  qFun = .charToDistFunc(distr, type = "q")
  dFun = .charToDistFunc(distr, type = "d")
  adTestStats = .myADTest(xVals, distr)
  print(adTestStats)
  if (class(adTestStats) == "adtest") {
    A = adTestStats$statistic
    p = adTestStats$p.value
  }
  else {
    A = NA
    p = NA
  }
  histObj = hist(xVals, plot = FALSE)
  if (missing(xlim)) {
    lq = do.call(qFun, c(list(1e-04), as.list(parameters)))
    uq = do.call(qFun, c(list(0.9999), as.list(parameters)))
    xlim = range(lq, uq, xVals)
  }
  xPoints = seq(xlim[1], xlim[2], length = 200)
  yPoints = do.call(dFun, c(list(xPoints), as.list(parameters)))
  if (missing(ylim)) {
    ylim = range(0, histObj$density, yPoints)
  }
  hist(xVals, freq = FALSE, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, main = main, ...)
  lines(xPoints, yPoints, col = line.col, lwd = line.width)
  abline(h = 0)
  legend("topright", c(paste(c(names(parameters), "A", "p"), ": ", c(format(parameters, digits = 3), format(A, digits = 3), format(p, digits = 3))), sep = " "), 
         inset = 0.02)
  if (box) {
    box()
  }
})
.xylimits = function(distrCollection, lowerquantile = 0.001, upperquantile = 0.999) {
  x = NULL
  y = NULL
  for (i in seq(along = distrCollection@distr)) {
    object = distrCollection@distr[[i]]
    xValues = object@x
    parameters = object@parameters
    distr = object@name
    qFun = .charToDistFunc(distr, type = "q")
    dFun = .charToDistFunc(distr, type = "d")
    lq = do.call(qFun, c(list(lowerquantile), as.list(parameters)))
    uq = do.call(qFun, c(list(upperquantile), as.list(parameters)))
    x = range(xValues, x, lq, uq)
    histObj = hist(xValues, plot = FALSE)
    xPoints = seq(x[1], x[2], length = 200)
    yPoints = do.call(dFun, c(list(xPoints), as.list(parameters)))
    y = range(y, 0, histObj$density, yPoints)
  }
  invisible(list(xlim = x, ylim = y))
}
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "distrCollection"), function(x, y, xlab, ylab, xlim, ylim, line.col, line.width, ...) {
  y = NULL
  object = x
  distrList = object@distr
  numDist = length(object@distr)
  numColWin = ceiling(numDist/2)
  if (missing(xlim)) 
    xlim = .xylimits(object)$xlim
  if (missing(ylim)) 
    ylim = .xylimits(object)$ylim
  if (missing(line.col)) 
    line.col = "red"
  if (missing(line.width)) 
    line.width = 1
  lapply(distrList, plot, xlim = xlim, ylim = ylim, line.col = line.col, line.width = line.width, ...)
  cat(paste("Total of", numDist, "plots created"))
  cat("\n")
  cat(paste("Use par(mfrow = c(2,", numColWin, ") to see all of them!", sep = ""))
  cat("\n")
})
qqPlot <- function (x, y, confbounds = TRUE, alpha, main, xlab, ylab, xlim, ylim, border = "red", bounds.col = "black", bounds.lty = 1,
                    start, ...) {
  DB = FALSE
  parList = list(...)
  if (is.null(parList[["col"]])) 
    parList$col = 1:2
  if (is.null(parList[["pch"]])) 
    parList$pch = 19
  if (is.null(parList[["lwd"]])) 
    parList$lwd = 1
  if (is.null(parList[["cex"]])) 
    parList$cex = 1
  
  #if (!require(MASS)) 
  #    stop("Package MASS needs to be installed!")
  
  if (class(x) == "distrCollection") {
    distList = x@distr
    for (i in 1:length(distList)) {
      d = distList[[i]]
      do.call(qqPlot, c(list(x = d@x, y = d@name), parList))
    }
    invisible()
  }
  if (missing(y)) 
    y = "normal"
  if(missing(alpha))
    alpha = 0.05
  if (alpha <=0 || alpha >=1) 
    stop(paste("alpha should be between 0 and 1!"))		
  if (missing(main)) 
    main = paste("Q-Q Plot for", deparse(substitute(y)), 
                 "distribution")
  if (missing(xlab)) 
    xlab = paste("Quantiles for", deparse(substitute(x)))
  if (missing(ylab)) 
    ylab = paste("Quantiles from", deparse(substitute(y)), 
                 "distribution")
  if (is.numeric(y)) {
    cat("\ncalling (original) qqplot from namespace stats!\n")
    return(stats::qqplot(x, y, ...))
  }
  qFun = NULL
  theoretical.quantiles = NULL
  xs = sort(x)
  distribution = tolower(y)
  distWhichNeedParameters = c("weibull", "logistic", "gamma", 
                              "exponential", "f", "geometric", "chi-squared", "negative binomial", 
                              "poisson")
  
  
  # new
  threeParameterDistr = c("weibull3", "lognormal3", "gamma3")                
  threeParameter = distribution %in% threeParameterDistr
  if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)
  # end new
  
  if (is.character(distribution)) {
    qFun = .charToDistFunc(distribution, type = "q")
    if (is.null(qFun)) 
      stop(paste(deparse(substitute(y)), "distribution could not be found!"))
  }
  theoretical.probs = ppoints(xs)
  
  xq = NULL
  yq = quantile(xs, prob = c(0.25, 0.75))
  dots <- list(...)
  if (TRUE) {
    if (DB) 
      print("TODO: Pass the estimated parameters correctly")
    fitList = .lfkp(parList, formals(qFun))
    fitList$x = xs
    fitList$densfun = distribution
    if (!missing(start)) 
      fitList$start = start
    if (DB) {
      print(fitList)
      print("Ende")
    }
    # new
    if(!threeParameter){
      fittedDistr = do.call(fitdistr, fitList)
      parameter = fittedDistr$estimate
      
      #save the distribution parameter#
      thethas = fittedDistr$estimate
      # save the cariance-covariance matrix
      varmatrix = fittedDistr$vcov
      # end of my code
      
      # new code for three parameter
    } else {
      parameter = do.call(paste(".",distribution, "3", sep = ""), list(xs) )    ####
      threshold = parameter$threshold
    }
    
    parameter = .lfkp(as.list(parameter), formals(qFun))
    params = .lfkp(parList, formals(qFun))
    parameter = .lfrm(as.list(parameter), params)
    parameter = c(parameter, params)
    theoretical.quantiles = do.call(qFun, c(list(c(theoretical.probs)), 
                                            parameter))
    
    # new
    if(!threeParameter){		
      # array containing names of the distributions, for which conf intervals can be computed
      confIntCapable = c("exponential", "log-normal", "logistic", "normal", "weibull", "gamma", "beta", "cauchy")
      getConfIntFun = .charToDistFunc(distribution, type = ".confint")
      # if possible, compute the conf intervals
      if(confbounds == TRUE){
        if(distribution %in% confIntCapable){
          confInt = getConfIntFun(xs, thethas, varmatrix, alpha)
        }
      }# end of my code
    }
    
    xq <- do.call(qFun, c(list(c(0.25, 0.75)), parameter))
    if (DB) {
      print(paste("parameter: ", parameter))
      print(xq)
    }
  }
  else {
    params = .lfkp(parList, formals(qFun))
    params$p = theoretical.probs
    theoretical.quantiles = do.call(qFun, params)
    params$p = c(0.25, 0.75)
    xq = do.call(qFun, params)
  }
  
  params = .lfkp(parList, c(formals(plot.default), par()))	
  
  if(!threeParameter){
    params$y = theoretical.quantiles
  }  else {
    params$y = theoretical.quantiles+threshold
  }
  params$x = xs
  params$xlab = xlab
  params$ylab = ylab
  params$main = main
  if (!(is.null(params$col[1]) || is.na(params$col[1]))) 
    params$col = params$col[1]
  if (!missing(xlim)) 
    params$xlim = xlim
  if (!missing(ylim)) 
    params$ylim = ylim
  params$lwd = 1
  do.call(plot, params)
  pParams = params
  pParams = .lfkp(pParams, list(x = 1, y = 1, col = 1, cex = 1))
  do.call(points, pParams)
  params = .lfkp(parList, c(formals(abline), par()))
  params$a = 0
  params$b = 1
  params$col = border
  do.call(abline, params)
  
  if(!threeParameter){
    # plot the confInt if available
    if(confbounds == TRUE){
      if(distribution %in% confIntCapable){
        params = .lfkp(parList, c(formals(lines), par()))	
        params$x = confInt[[3]]
        params$y = confInt[[1]]
        params$col = bounds.col
        params$lty = bounds.lty
        do.call(lines, params)
        
        params$x = confInt[[3]]
        params$y = confInt[[2]]
        params$col = bounds.col
        params$lty = bounds.lty
        do.call(lines, params)
      }
    } #end of my function
  }
  
  invisible(list(x = theoretical.quantiles, y = xs, int = params$a, 
                 slope = params$b))
}

ppPlot <- function (x, distribution, confbounds = TRUE, alpha, probs, main, xlab, ylab, xlim, ylim, 
                    border = "red", bounds.col = "black", bounds.lty = 1, grid = TRUE, box = TRUE, stats = TRUE, start, 
                    ...) {
  DB = FALSE
  conf.level = 0.95
  conf.lines = TRUE
  #if (!require(MASS)) 
  #    stop("Package MASS needs to be installed!")
  if (!(is.numeric(x) | (class(x) == "distrCollection"))) 
    stop(paste(deparse(substitute(x)), " needs to be numeric or an object of class distrCollection"))
  parList = list(...)
  if (is.null(parList[["col"]])) 
    parList$col = c("black", "red", "gray")
  if (is.null(parList[["pch"]])) 
    parList$pch = 19
  if (is.null(parList[["lwd"]])) 
    parList$lwd = 1
  if (is.null(parList[["cex"]])) 
    parList$cex = 1
  if (DB) 
    print(parList)
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
      xlim = c(min(x1) - 0.1 * diff(range(x1)), max(x1) + 
                 0.1 * diff(range(x1)))
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
  if (class(x) == "distrCollection") {
    distList = x@distr
    for (i in 1:length(distList)) {
      d = distList[[i]]
      do.call(ppPlot, c(list(x = d@x, distribution = d@name), 
                        parList))
    }
    invisible()
  }
  distWhichNeedParameters = c("weibull", "gamma", "logistic", 
                              "exponential", "f", "geometric", "chi-squared", "negative binomial", 
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
    if (DB) 
      print("TODO: Pass the estimated parameters correctly")
    fitList = .lfkp(parList, formals(qFun))
    fitList$x = x1
    fitList$densfun = distribution
    if (!missing(start)) 
      fitList$start = start
    
    if(!threeParameter){
      fittedDistr = do.call(fitdistr, fitList)
      parameter = fittedDistr$estimate
      #save the distribution parameter#
      thethas = fittedDistr$estimate
      # save the cariance-covariance matrix
      varmatrix = fittedDistr$vcov		
    } else {
      parameter = do.call(paste(".",distribution, "3", sep = ""), list(x1) )    ####
      print(parameter[3])
      threshold = parameter$threshold
    }
    parameter = .lfkp(as.list(parameter), formals(qFun))
    params = .lfkp(parList, formals(qFun))
    parameter = .lfrm(as.list(parameter), params)
    print(parameter)
    parameter = c(parameter, params)
    if (DB) {
      print(qFun)
      print(as.list(parameter))
      print(list(probs))
    }
    
    
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
    if (DB) {
      print(paste("parameter: ", parameter))
      print(xq)
    }
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
  do.call(plot, params)
  pParams = params
  params = .lfkp(parList, list(cex.main = 1, cex.axis = 1, 
                               cex.lab = 1))
  params$side = 1
  axisAtX = do.call(axis, params)
  params$side = 2
  params$at = axisAtY
  params$labels = probs
  params$las = 2
  do.call(axis, params)
  if (grid) {
    params = .lfkp(parList, c(formals(abline), list(lwd = 1, 
                                                    col = 1)))
    params$h = axisAtY
    params$v = axisAtX
    if (!(is.null(params$col[3]) || is.na(params$col[3]))) 
      params$col = params$col[3]
    else params$col = 1
    if (!(is.null(params$lwd[2]) || is.na(params$lwd[2]))) 
      params$lwd = params$lwd[2]
    else params$lwd = 1
    do.call(abline, params)
  }
  pParams = .lfkp(pParams, list(x = 1, y = 1, col = 1, cex = 1))
  do.call(points, pParams)
  params = .lfkp(parList, c(formals(abline), par()))
  if(!threeParameter){
    params$a = 0
  }  else {
    params$a = -threshold
  }
  params$b = 1
  params$col = border
  do.call(abline, params)
  
  if(!threeParameter){
    # plot the confInt if available
    if(confbounds == TRUE){
      if(distribution %in% confIntCapable){
        params = .lfkp(parList, c(formals(lines), par()))	
        params$x = confInt[[3]]
        params$y = confInt[[2]]
        params$col = bounds.col
        params$lty = bounds.lty
        do.call(lines, params)
        
        params$x = confInt[[3]]
        params$y = confInt[[1]]
        params$col = bounds.col
        params$lty = bounds.lty
        do.call(lines, params)		
      }
    } #end of my function
  }
  if (box) 
    box()
  invisible(list(x = x, y = y, int = params$a, slope = params$b))
}

# -------------- mix_n.r -------------- 
setClass(Class = "mixDesign", representation = representation(name = "character", factors = "list", total = "numeric", lower = "numeric", design = "data.frame", 
                                                              designType = "character", pseudo = "data.frame", response = "data.frame", Type = "data.frame", block = "data.frame", runOrder = "data.frame", standardOrder = "data.frame", 
                                                              desireVal = "list", desirability = "list", fits = "data.frame"))
setMethod("factors", "mixDesign", function(x) x@factors)
setReplaceMethod("factors", "mixDesign", function(x, value) {
  if (length(value) != ncol(x@pseudo)) 
    stop("\nNumber of factors doesn't match with number of columns for factorial Design\n")
  x@factors <- value
  x
})
setMethod("names", "mixDesign", function(x) {
  return(sapply(x@factors, names))
})
setReplaceMethod("names", "mixDesign", function(x, value) {
  for (i in 1:length(x@factors)) names(x@factors[[i]]) = as.character(value[i])
  x
})
setMethod("as.data.frame", "mixDesign", function(x, row.names = NULL, optional = FALSE) {
  frameOut = cbind(x@standardOrder, x@runOrder, x@Type, x@pseudo, x@response)
  return(frameOut)
})
as.data.frame.mixDesign = function(x, row.names = NULL, optional = FALSE, ...) {
  frameOut = cbind(x@standardOrder, x@runOrder, x@Type, x@pseudo, x@response)
  return(frameOut)
}
setMethod("show", signature(object = "mixDesign"), function(object) {
  print(format(as.data.frame(object), digits = 4))
})
setMethod("response", "mixDesign", function(object) {
  return(object@response)
})
setReplaceMethod("response", "mixDesign", function(object, value) {
  print(deparse(substitute(value)))
  if (!is.numeric(value) & !is.data.frame(value)) 
    stop("vector or data.frame must be given")
  if (is.numeric(value)) {
    if (length(value) != nrow(object@pseudo)) 
      stop("differing lengths")
    temp = data.frame(value)
    names(temp) = deparse(substitute(value))[1]
    value = temp
  }
  if (is.data.frame(value)) {
    if (nrow(value) != nrow(object@pseudo)) 
      stop("differing number of rows")
  }
  object@response = value
  return(object)
})
.npp = function(mdo) {
  pseudo = mdo@pseudo
  Type = mdo@Type
  temp = as.character(as.data.frame(mdo)$Type)
  tab = table(temp)
  nums = data.frame(matrix(tab, nrow = 4, ncol = length(tab), byrow = TRUE))
  nums[1:2, ] = NA
  nums[4, ] = c(nrow(Type), rep("", length(tab) - 1))
  row.names(nums) = c("Unique", "Replicates", "Sub Total", "Total")
  names(nums) = names(tab)
  for (i in names(tab)) {
    sSet = pseudo[Type == i, ]
    usSet = unique(sSet)
    nums["Unique", i] = nrow(usSet)
    if (nrow(usSet) == 1) {
      nums["Replicates", i] = nrow(sSet)
    }
    else {
      for (j in 1:nrow(usSet)) {
        uCount1 = sum(apply(apply(sSet, 1, "==", usSet[j, ]), 2, "all") * 1)
        if (j == 1) {
          uCount2 = uCount1
          nums["Replicates", i] = uCount1
        }
        if (uCount2 != uCount1) 
          nums["Replicates", i] = -1
      }
    }
  }
  cat("Information about the Design Points:")
  cat("\n")
  cat("\n")
  print(nums)
}
setMethod("show", signature(object = "mixDesign"), function(object) {
  print(as.data.frame(object))
})
setMethod(".nfp", "mixDesign", function(object) {
  x = factors(object)
  if (is.list(x) && length(x[[1]]) > 0) {
    numAttr = length(attributes(x[[1]]))
    .numFac = length(x)
    frameOut = data.frame(matrix(ncol = .numFac, nrow = (numAttr - 1)))
    for (i in 1:(numAttr - 1)) {
      charVec = character(0)
      for (j in 1:.numFac) {
        charVec = c(charVec, names(attributes(x[[1]])[i]), "\t\t")
        frameOut[i, j] = attributes(x[[j]])[[i]]
      }
    }
    names(frameOut) = names(x)
    rownames(frameOut) = names(attributes(x[[1]]))[1:(numAttr - 1)]
  }
  else {
    stop("no list given or length of list < 1")
  }
  print(frameOut)
})
setMethod("summary", signature(object = "mixDesign"), function(object) {
  cat(paste("Simplex", toupper(object@designType), "Design"))
  cat("\n")
  cat("Information about the factors:\n\n")
  .nfp(object)
  cat("\n-----------\n")
  cat("\n")
  .npp(object)
  cat("\n-----------\n")
  cat("\n")
  cat("Information about the constraints:\n\n")
  lower = object@lower
  temp = character(0)
  for (i in seq(along = lower)) temp = c(temp, paste(LETTERS[i], ">=", lower[i]))
  cat(temp)
  cat("\n")
  cat("\n-----------\n")
  cat("\n")
  times = nrow(object@pseudo)
  pseudo = format(object@pseudo, digits = 2)
  design = format(object@design, digits = 2)
  amount = design
  if (object@total[2] != 1) 
    amount = format(object@design * object@total[2], digits = 2)
  temp = c("                             ", "PseudoComponent", "_|_", "Proportion", "_|_", "Amount")
  cat(temp)
  cat("\n")
  cat("\n")
  temp = cbind(pseudo, `_` = rep(" ", times = times), `|` = rep("|", times = times), `_` = rep(" ", times = times), design)
  temp = cbind(temp, `_` = rep(" ", times = times), `|` = rep("|", times = times), `_` = rep(" ", times = times), amount)
  temp = cbind(object@standardOrder, object@runOrder, object@Type, `|` = rep("|", times = times), temp, `|` = rep("|", times = times), object@response)
  show(temp)
  cat("\n-----------\n")
  cat("\n")
  cat(paste("Mixture Total:", object@total[1], "equals", object@total[2]))
  cat("\n")
  cat("\n")
  invisible(as.data.frame(object))
})
setMethod("units", "mixDesign", function(x) {
  return(sapply(factors(x), .unit))
})
setMethod("units<-", "mixDesign", function(x, value) {
  for (i in 1:length(x@factors)) if (length(value) > 1) 
    .unit(x@factors[[i]]) = as.character(value[i])
  else .unit(x@factors[[i]]) = as.character(value[1])
  x
})
setMethod("highs", "mixDesign", function(object) {
  listOut = vector(mode = "list")
  for (i in names(factors(object))) {
    listOut[i] = .high(object@factors[[i]])
  }
  return(listOut)
})
setReplaceMethod("highs", "mixDesign", function(object, value) {
  for (i in seq(along = object@factors)) if (length(value) > 1) 
    .high(object@factors[[i]]) = value[i]
  else .high(object@factors[[i]]) = value[1]
  return(object)
})
setMethod("lows", "mixDesign", function(object) {
  listOut = vector(mode = "list")
  for (i in names(factors(object))) {
    listOut[i] = .low(object@factors[[i]])
  }
  return(listOut)
})
setReplaceMethod("lows", "mixDesign", function(object, value) {
  for (i in seq(along = object@factors)) {
    if (length(value) > 1) 
      .low(object@factors[[i]]) = value[i]
    else .low(object@factors[[i]]) = value[1]
  }
  return(object)
})
.permut = function(x) {
  if (any(is.na(x))) 
    stop(paste(deparse(substitute(x)), "contains NA"))
  x = sort(x, decreasing = FALSE)
  n = length(x)
  num = 1:n
  frameOut = matrix(NA, nrow = 1, ncol = n)
  frameOut[1, ] = x
  while (TRUE) {
    highest = NA
    for (j in 1:(n - 1)) {
      if (x[j] < x[j + 1]) 
        highest = j
    }
    if (is.na(highest)) 
      return(frameOut)
    else {
      l = max((num)[x[highest] < x])
      temp = x[l]
      x[l] = x[highest]
      x[highest] = temp
      x[(highest + 1):n] = rev(x[(highest + 1):n])
    }
    frameOut = rbind(frameOut, x, deparse.level = 2)
  }
}
.simplexCentroid = function(p) {
  if (p <= 1 | !is.numeric(p)) 
    stop("invalid value for p")
  frameOut = NA
  for (i in 1:p) {
    initial = rep(0, times = p)
    initial[1:i] = 1/i
    mat = .permut(initial)
    if (i == 1) 
      frameOut = mat
    else frameOut = rbind(frameOut, mat, deparse.level = 2)
  }
  frameOut = data.frame(frameOut, row.names = 1:(2^p - 1))
  names(frameOut) = LETTERS[1:p]
  return(frameOut)
}
.simplexCentroid(4)
mixDesign = function(p, n = 3, type = "lattice", center = TRUE, axial = FALSE, delta, replicates = 1, lower, total = 1, randomize, seed) {
  DB = FALSE
  frameOut = NA
  out = new("mixDesign")
  if (missing(p)) 
    stop("the number of factors p must be given")
  if (p <= 1 | !is.numeric(p)) 
    stop("invalid value for p")
  if (!(type %in% c("lattice", "centroid"))) 
    stop("type needs to be \"lattice\" or \"centroid\"")
  out@designType = type
  if (missing(delta)) 
    delta = (p - 1)/(2 * p)
  if (missing(randomize)) 
    randomize = TRUE
  if (!missing(seed)) 
    set.seed(seed)
  if (!is.numeric(total)) 
    stop("total needs to be a numeric vector with <= 2 arguments")
  else {
    if (total[1] > 1 || total[1] <= 0) 
      stop("totall[1] needs to be within (0,1]")
    if (is.na(total[2])) 
      total[2] = 1
    if (total[2] <= 0) 
      stop("total[2] needs to be > 0")
  }
  out@total = total
  if (!is.numeric(replicates)) 
    stop("replicates need to be numeric")
  if (delta > (p - 1)/p) {
    delta = (p - 1)/(2 * p)
    warning(paste("delta was reset to:", delta))
  }
  if (missing(lower)) 
    lower = 0
  if (length(lower) == 1) 
    lower = rep(lower, p)
  if (length(lower) > 1) {
    initial = rep(0, p)
    if (length(lower) < p) 
      lower[, (length(lower) + 1):p] = 0
  }
  out@lower = lower
  repTemp = list(center = 1, axial = 1)
  for (i in 1:(p - 1)) repTemp[[as.character(i)]] = 1
  if (length(replicates) > 1) {
    for (i in 1:length(replicates)) repTemp[[i]] = replicates[[i]]
    replicates = repTemp
  }
  if (length(replicates) == 1) {
    for (i in 1:length(repTemp)) repTemp[[i]] = replicates
    replicates = repTemp
  }
  if (DB) 
    print(replicates)
  N = factorial(p + n - 1)/(factorial(n) * factorial(p - 1))
  if (identical(type, "lattice")) {
    j = 1
    x = numeric(p)
    j = 1
    x[1] = n
    for (i in 1:N) {
      x[j + 1] = n - (sum(x[1:j]))
      if (j < (p - 1)) 
        x[(j + 2):p] = 0
      if (i == 1) 
        frameOut = data.frame(matrix(x, ncol = p, nrow = 1))
      else frameOut = rbind(frameOut, x)
      if (DB) 
        print(x)
      logVec = rep(FALSE, p)
      logVec[1:(p - 1)] = x[1:p - 1] > 0
      if (any(logVec)) {
        j = max((1:p)[logVec])
        x[j] = x[j] - 1
      }
    }
    frameOut = (frameOut/(n))
    names(frameOut) = LETTERS[1:p]
  }
  if (identical(type, "centroid")) {
    frameOut = .simplexCentroid(p)
  }
  frameOutCopy = frameOut[1, ]
  Type = data.frame(Type = "1")
  temp = apply(ceiling(frameOut), 2, ">", 0) * 1
  temp = apply(temp, 1, sum)
  for (i in 1:nrow(frameOut)) {
    typ = as.character(temp[i])
    times = replicates[[typ]]
    if (is.null(times)) 
      times = 0
    else times = times - 1
    repFrame = frameOut[i, ]
    if (all(frameOut[i, ] == frameOut[i, 1])) 
      typFrame = data.frame(Type = "center")
    else typFrame = data.frame(Type = paste(typ, "-blend", sep = ""))
    if (times >= 1) {
      for (j in 1:times) {
        repFrame = rbind(repFrame, frameOut[i, ])
        typFrame = rbind(typFrame, data.frame(Type = paste(typ, "-blend", sep = "")))
      }
    }
    frameOutCopy = rbind(frameOutCopy, repFrame)
    Type = rbind(Type, typFrame)
    if (i == 1) {
      frameOutCopy = frameOutCopy[c(-1), ]
      Type = data.frame(Type = Type[c(-1), ])
    }
  }
  frameOut = frameOutCopy
  if (DB) {
    print(Type)
    print(frameOutCopy)
  }
  keepIndex = (1:nrow(frameOut))[!apply(Type, 1, "==", "center")]
  Type = data.frame(Type = Type[keepIndex, ])
  frameOut = frameOut[keepIndex, ]
  if (center) {
    center = data.frame(matrix(1/p, nrow = 1, ncol = p))
    if (DB) 
      print(center)
    times = replicates$center
    if (is.null(times)) 
      times = 0
    else times = times - 1
    if (n == p) 
      times = times - 1
    if ((n%%p) == 0) {
      numCenter = n%/%p
    }
    temp = center
    if (times >= 1) {
      for (i in 1:times) center = rbind(center, temp)
    }
    names(center) = names(frameOut)
    frameOut = rbind(frameOut, center)
    Type = rbind(Type, data.frame(Type = rep("center", times + 1)))
    if (DB) 
      print(frameOut)
  }
  if (axial) {
    temp = rep(NA, p)
    axial = data.frame(matrix(NA, ncol = p, nrow = p))
    for (i in 1:p) {
      temp[i] = delta + 1/p
      temp[c(-i)] = (1 - temp[i])/(p - 1)
      axial[i, ] = temp
    }
    times = replicates$axial
    if (is.null(times)) 
      times = 0
    else times = times - 1
    temp = axial
    if (times >= 1) {
      for (i in 1:times) axial = rbind(axial, temp)
    }
    names(axial) = names(frameOut)
    frameOut = rbind(frameOut, axial)
    Type = rbind(Type, data.frame(Type = rep("axial", (times + 1) * p)))
    if (DB) {
      print(frameOut)
      print(Type)
    }
  }
  StandOrder = 1:nrow(frameOut)
  RunOrder = StandOrder
  if (randomize) {
    RunOrder = sample(1:nrow(frameOut), nrow(frameOut), replace = FALSE, prob = NULL)
  }
  frameOut = frameOut[order(RunOrder), ]
  row.names(frameOut) = frameOut$RunOrder
  out@pseudo = frameOut
  out@runOrder = data.frame(RunOrder = RunOrder)
  out@standardOrder = data.frame(StandOrder = StandOrder)
  out@Type = data.frame(Type = Type[order(RunOrder), ])
  out@response = data.frame(y = rep(NA, nrow(out@pseudo)))
  design = frameOut
  design[, ] = NA
  for (i in 1:ncol(frameOut)) {
    design[, i] = frameOut[, i] * (total[1] - sum(lower)) + lower[i]
  }
  out@design = design
  listFac = vector("list", p)
  for (i in seq(along = listFac)) listFac[i] = new("doeFactor")
  names(listFac) = LETTERS[1:p]
  factors(out) = listFac
  if (out@total[2] != 1) {
    lows(out) = lower * out@total[2]
    highs(out) = 1 * out@total[2]
    units(out) = "NA"
  }
  else if (any(out@lower != 0)) {
    lows(out) = out@lower
    highs(out) = 1 * out@total[1]
    units(out) = "%"
  }
  else {
    lows(out) = 0
    highs(out) = 1 * out@total[1]
    units(out) = "%"
  }
  return(out)
} 

# -------------- ort_s.r --------------
.L4_2 = list(id = "L4_2", type = "single", niv_max = 2, runs = 4, anzahl_spalten = 3, levels1 = 2, levels2 = 0, factors1 = 3, factors2 = 0, colOrder = c(1,2, 3), r4 = c(1, 2), r2 = c(1, 2, 3), design = data.frame(`1` = c(1, 1, 2, 2), `2` = c(1, 2, 1, 2), `3` = c(1, 2, 2, 1)))
.L8_2 = list(id = "L8_2", type = "single", niv_max = 2, runs = 8, anzahl_spalten = 7, levels1 = 2, levels2 = 0, factors1 = 7, factors2 = 0, colOrder = c(1, 2, 4, 7, 3, 5, 6), r4 = c(1, 2, 4), r2 = c(1, 2, 4, 7), r1 = c(1, 2, 4, 7, 3, 5, 6), 
             design = data.frame(`1` = c(1, 1, 1, 1, 2, 2, 2, 2), `2` = c(1, 1, 2, 2, 1, 1, 2, 2), `3` = c(1, 1, 2, 2, 2, 2, 1, 1), 
                                 `4` = c(1, 2, 1, 2, 1, 2, 1, 2), `5` = c(1, 2, 1, 2, 2, 1, 2, 1), `6` = c(1, 2, 2, 1, 1, 2, 2, 1), 
                                 `7` = c(1, 2, 2, 1, 2, 1, 1, 2)),
             ia_table = data.frame(`1` = c(1, "NA", "NA", "NA", "NA", "NA", "NA"), `2` = c(3, 2, "NA", "NA", "NA", "NA", "NA"), 
                                   `3` = c(2, 1, 3, "NA", "NA", "NA", "NA"), `4` = c(5, 6, 7, 4, "NA", "NA", "NA"),
                                   `5` = c(4, 7, 6, 1, 5, "NA", "NA"), `6` = c(7, 4, 5, 2, 3, 6, "NA"), 
                                   `7` = c(6, 5, 4,3, 2, 1, 7)))
.L9_3 = list(id = "L9_3", type = "single", niv_max = 3, runs = 9, anzahl_spalten = 4, levels1 = 3, levels2 = 0, factors1 = 4, factors2 = 0, 
             colOrder = c(1, 2, 3, 4), r4 = c(1, 2), r1 = c(1, 2, 3, 4), 
             design = data.frame(`1` = c(1, 1, 1, 2, 2, 2, 3, 3, 3), `2` = c(1, 2, 3, 1, 2, 3, 1, 2, 3), 
                                 `3` = c(1, 2, 3, 2, 3, 1, 3, 1, 2), `4` = c(1, 2, 3, 3, 1, 2, 2, 3, 1)), ia_table = data.frame(`2` = c("3,4")))
.L12_2 = list(id = "L12_2", type = "single", niv_max = 2, runs = 12, anzahl_spalten = 11, levels1 = 2, levels2 = 0, factors1 = 11, factors2 = 0, 
              colOrder = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), r1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), 
              design = data.frame(`1` = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2), `2` = c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2), 
                                  `3` = c(1, 1, 2, 1, 2, 2, 2, 2, 1, 2, 1, 1), `4` = c(1, 1, 2, 2, 1, 2, 2, 1, 2, 1, 2, 1), `5` = c(1, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2), 
                                  `6` = c(1, 2, 1, 1, 2, 2, 1, 2, 2, 1, 2, 1), `7` = c(1, 2, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2), `8` = c(1, 2, 1, 2, 2, 1, 2, 1, 2, 2, 1, 1), 
                                  `9` = c(1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1, 2), `10` = c(1, 2, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2), `11` = c(1, 2, 2, 2, 1, 1, 1, 2, 1, 2, 2, 1)),
              ia_table = data.frame(`1` = c("ohne interaction table")))
.L16_2 = list(id = "L16_2", type = "single", niv_max = 2, runs = 15, anzahl_spalten = 15, levels1 = 2, levels2 = 0, factors1 = 15, factors2 = 0, 
              colOrder = c(1, 2, 4, 7, 8, 11, 13, 14, 3, 5, 6, 9, 10, 12, 15), r4 = c(1, 2, 4, 8), r3 = c(1, 2, 4, 8, 15),
              r2 = c(1, 2, 4, 7, 8, 11, 13, 14), r1 = c(1,2, 4, 7, 8, 11, 13, 14, 3, 5, 6, 9, 10, 12, 15), 
              design = data.frame(`1` = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2), `2` = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2), 
                                  `3` = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1), `4` = c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2), 
                                  `5` = c(1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1), `6` = c(1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1), 
                                  `7` = c(1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2), `8` = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2), 
                                  `9` = c(1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1), `10` = c(1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1)), 
              ia_table = data.frame(`1` = c(1, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), `2` = c(3, 2, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                                    `3` = c(2, 1, 3, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"),`4` = c(5, 6, 7, 4, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                                    `5` = c(4, 7, 6, 1, 5, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), `6` = c(7, 4, 5, 2, 3, 6, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                                    `7` = c(6, 5, 4, 3, 2, 1, 7, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), `8` = c(9, 10, 11, 12, 13, 14, 15, 8, "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                                    `9` = c(8, 11, 10, 13, 12, 15, 14, 1, 9, "NA", "NA", "NA", "NA", "NA", "NA"), `10` = c(11, 8, 9, 14, 15, 12, 13, 2, 3, 10, "NA", "NA", "NA", "NA", "NA"), 
                                    `11` = c(10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 11, "NA", "NA", "NA", "NA"), `12` = c(13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 12, "NA", "NA",  "NA"),
                                    `13` = c(12, 15, 14, 9, 8, 11, 10, 5, 4, 7, 6, 1, 13, "NA", "NA"), `14` = c(15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 14, "NA"), `15` = c(14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 15)))
.L16_4 = list(id = "L16_4", type = "single", niv_max = 4, runs = 16, anzahl_spalten = 5, levels1 = 4, levels2 = 0, factors1 = 5, factors2 = 0, 
              design = data.frame(`1` = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), `2` = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4), 
                                  `3` = c(1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1), `4` = c(1, 2, 3, 4, 3, 4, 1, 2, 4, 3, 2, 1, 2, 1, 4, 3), `5` = c(1, 2, 3, 4, 4, 3, 2, 1, 2, 1, 4, 3, 3, 4, 1, 2)), ia_table = data.frame(`2` = c("3,4,5")))
.L18_2_3 = list(id = "L18_2_3", type = "mixed", niv_max = 3, runs = 18, anzahl_spalten = 8, levels1 = 2, levels2 = 3, factors1 = 1, factors2 = 7, 
                design = data.frame(`1` = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2), `2` = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3), 
                                    `3` = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), `4` = c(1, 2, 3, 1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2), 
                                    `5` = c(1, 2, 3, 2, 3, 1, 1, 2, 3, 3, 1, 2, 3, 1, 2, 2, 3, 1), `6` = c(1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2), 
                                    `7` = c(1, 2, 3, 3, 1, 2, 2, 3, 1, 2, 3, 1, 3, 1, 2, 1, 2, 3), `8` = c(1, 2, 3, 3, 2, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 2, 3, 1)), 
                ia_table = data.frame(`1` = c("ohne interaction table")))
.L25_5 = list(id = "L25_5", type = "single", niv_max = 5, runs = 25, anzahl_spalten = 6, levels1 = 5, levels2 = 0, factors1 = 6, factors2 = 0, 
              design = data.frame(`1` = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5), `2` = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
                                  `3` = c(1, 2, 3, 4, 5, 2, 3, 4, 5, 1, 3, 4, 5, 1, 2, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4), `4` = c(1, 2, 3, 4, 5, 3, 4, 5, 1, 2, 5, 1, 2, 3, 4, 2, 3, 4, 5, 1, 4, 5, 1, 2, 3), 
                                  `5` = c(1, 2, 3, 4, 5, 4, 5, 1, 2, 3, 2, 3, 4, 5, 1, 5, 1, 2, 3, 4, 3, 4, 5, 1, 2), `6` = c(1, 2, 3, 4, 5, 5, 1, 2, 3, 4, 4, 5, 1, 2, 3, 3, 4, 5, 1, 2, 2, 3, 4, 5, 1)), 
              ia_table = data.frame(`2` = c("3,4,5,6")))
.L27_3 = list(id = "L27_3", type = "single", niv_max = 3, runs = 27, anzahl_spalten = 13, levels1 = 3, levels2 = 0, factors1 = 13, factors2 = 0, 
              colOrder = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), r4 = c(1, 2, 5), r2 = c(1, 2, 5, 9, 10, 12, 13), r1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), 
              design = data.frame(
                `1` = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3), `2` = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3),
                `3` = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 2, 2, 2, 3, 3, 3, 1, 1, 1, 3, 3, 3, 1, 1, 1, 2, 2, 2), `4` = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 1, 1, 1),
                `5` = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), `6` = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2),
                `7` = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 3, 1, 2, 3, 1, 2, 3, 1, 2, 2, 3, 1, 2, 3, 1, 2, 3, 1), `8` = c(1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2),
                `9` = c(1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2, 1, 2, 3, 3, 1, 2, 1, 2, 3, 2, 3, 1), `10` = c(1, 2, 3, 2, 3, 1, 3, 1, 2, 3, 1, 2, 1, 2, 3, 2, 3, 1, 2, 3, 1, 3, 1, 2, 1, 2, 3),
                `11` = c(1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3, 1), `12` = c(1, 2, 3, 3, 1, 2, 2, 3, 1, 2, 3, 1, 1, 2, 3, 3, 1, 2, 3, 1, 2, 2, 3, 1, 1, 2, 3),
                `13` = c(1, 2, 3, 3, 1, 2, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 2, 3, 1, 1, 2, 3, 3, 1, 2)), 
              ia_table = data.frame(
                `1` = c(1,"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"), `2` = c("3,4",2,"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"),
                `3` = c("2,4","1,4", 3,"NA","NA","NA", "NA","NA","NA","NA","NA","NA","NA"),`4` = c("2,3","1,3","1,2",4,"NA","NA","NA","NA","NA","NA","NA","NA","NA"),
                `5` = c( "6,7","8,11","9,13","10,12","5","NA","NA","NA", "NA","NA","NA","NA","NA"),`6` = c( "5,7","9,12","10,11","8,13","1,7",6,"NA","NA","NA", "NA","NA", "NA", "NA"),
                `7` = c("5,6", "10,13", "8,12", "9,11", "1,6", "1,5", 7, "NA", "NA", "NA", "NA", "NA", "NA"),`8` = c("9,10","5,11", "7,12","6,13","2,11","4,13","3,12",8,"NA","NA","NA","NA","NA"),
                `9` = c("8,10","6,12","5,13","7,11","3,13","2,12","4,11","1,10",9,"NA","NA","NA","NA"),`10` = c("8,9","7,13","6,11","5,12","4,12","3,11","2,13","1,9","1,8",10,"NA","NA","NA"),
                `11` = c("12,13","5,8","6,10", "7,9","2,8","3,10","4,9","2,5","4,7","3,6", 11,"NA","NA"),`12` = c("11,13","6,9","7,8","5,10","4,10","2,9","3,8","3,7","2,6","4,5","1,13",12,"NA"),
                `13` = c("11,12","7,10","5,9","6,8","3,9","4,8", "2,10","4,6","3,5","2,7","1,12","1,11",13)))
.L32_2 = list(id = "L32_2", type = "single", niv_max = 2, runs = 32, anzahl_spalten = 31, levels1 = 2, levels2 = 0, factors1 = 31, factors2 = 0, 
              colOrder = c(1, 2, 4, 7, 8, 11, 13, 14, 16, 19, 21, 22, 25, 26, 28, 31, 3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30), r4 = c(1, 2, 4, 8, 16), 
              r3 = c(1, 2, 4, 8, 16, 31), r2 = c(1, 2, 4, 8, 16, 31, 7, 11, 13, 14, 19, 21, 22, 25, 26, 28), r1 = c(1, 2, 4, 7, 8, 11, 13, 14, 16, 19, 21, 22, 25, 26, 28, 31, 3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30),
              design = data.frame(`1` = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), `2` = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2), 
                                  `3` = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1), `4` = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2), 
                                  `5` = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1), `6` = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1), 
                                  `7` = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2), `8` = c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2), 
                                  `9` = c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1), `10` = c(1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1), 
                                  `11` = c(1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2), `12` = c(1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1), 
                                  `13` = c(1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2), `14` = c(1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2), 
                                  `15` = c(1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1), `16` = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2), 
                                  `17` = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1), `18` = c(1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1), 
                                  `19` = c(1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2), `20` = c(1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1), 
                                  `21` = c(1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2), `22` = c(1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2), 
                                  `23` = c(1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1), `24` = c(1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1), 
                                  `25` = c(1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2), `26` = c(1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2), 
                                  `27` = c(1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1), `28` = c(1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2), 
                                  `29` = c(1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1), `30` = c(1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1), 
                                  `31` = c(1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2)))
.ia_table = data.frame(`1` = c(1, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `2` = c(3, 2, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `3` = c(2, 1, 3, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `4` = c(5, 6, 7, 4, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",  "NA"), 
                       `5` = c(4, 7, 6, 1, 5, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `6` = c(7, 4, 5, 2, 3, 6, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `7` = c(6, 5, 4, 3, 2, 1, 7, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"),
                       `8` = c(9, 10, 11, 12, 13, 14, 15, 8, "NA", "NA", "NA", "NA", "NA", "NA", "NA",  "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `9` = c(8, 11, 10, 13, 12, 15, 14, 1, 9, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `10` = c(11, 8, 9, 14, 15, 12, 13, 2, 3, 10, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `11` = c(10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 11, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `12` = c(13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 12, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `13` = c(12, 15, 14, 9, 8, 11, 10, 5, 4, 7, 6, 1, 13, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"),
                       `14` = c(15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 14, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"),
                       `15` = c(14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 15, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `16` = c(17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `17` = c(16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30, 1, 17, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"),
                       `18` = c(19, 16, 17, 22, 23, 20, 21, 26, 27, 24, 25, 30, 31, 28, 29, 2, 3, 18, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `19` = c(18, 17, 16, 23, 22, 21, 20, 27, 26, 25, 24, 31, 30, 29, 28, 3, 2, 1, 19, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `20` = c(21, 22, 23, 16, 17, 18, 19, 28, 29, 30, 31, 24, 25, 26, 27, 4, 5, 6, 7, 20, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `21` = c(20, 23, 22, 17, 16, 19, 18, 29, 28, 31, 30, 25, 24, 27, 26, 5, 4, 7, 6, 1, 21, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `22` = c(23, 20, 21, 18, 19, 16, 17, 30, 31, 28, 29, 26, 27, 24, 25, 6, 7, 4, 5, 2, 3, 22, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `23` = c(22, 21, 20, 19, 18, 17, 16, 31, 30, 29, 28, 27, 26, 25, 24, 7, 6, 5, 4, 3, 2, 1, 23, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"),
                       `24` = c(25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23, 8, 9, 10, 11, 12, 13, 14, 15, 24, "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                       `25` = c(24, 27, 26, 29, 28, 31, 30, 17, 16, 19, 18, 21, 20, 23, 22, 9, 8, 11, 10, 13, 12, 15, 14, 1, 25, "NA", "NA", "NA", "NA", "NA", "NA"),
                       `26` = c(27, 24, 25, 30, 31, 28, 29, 18, 19, 16, 17, 22, 23, 20, 21, 10, 11, 8, 9, 14, 15, 12, 13, 2, 3, 26, "NA", "NA", "NA", "NA", "NA"), 
                       `27` = c(26, 25, 24, 31, 30, 29, 28, 19, 18, 17, 16, 23, 22, 21, 20, 11, 10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 27, "NA", "NA", "NA", "NA"), 
                       `28` = c(29, 30, 31, 24, 25, 26, 27, 20, 21, 22, 23, 16, 17, 18, 19, 12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 28, "NA", "NA", "NA"), 
                       `29` = c(28, 31, 30, 25, 24, 27, 26, 21, 20, 23, 22, 17, 16, 19, 18, 13, 12, 15, 14, 9, 8, 11, 10, 5, 4, 7, 6, 1, 29, "NA", "NA"),
                       `30` = c(31, 28, 29, 26, 27, 24, 25, 22, 23, 20, 21, 18, 19, 16, 17, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 30, "NA"), 
                       `31` = c(30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 31))
.L32_2$ia_table = .ia_table
.L32_2_4 = list(id = "L32_2_4", type = "mixed", niv_max = 4, runs = 32, anzahl_spalten = 10, levels1 = 2, levels2 = 4, factors1 = 1, factors2 = 9, 
                design = data.frame(`1` = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), `2` = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), 
                                    `3` = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4), `4` = c(1, 2, 3, 4, 1, 2, 3, 4, 2, 1, 4, 3, 2, 1, 4, 3, 4, 3, 2, 1, 4, 3, 2, 1, 3, 4, 1, 2, 3, 4, 1, 2), 
                                    `5` = c(1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1, 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1), `6` = c(1, 2, 3, 4, 2, 1, 4, 3, 4, 3, 2, 1, 3, 4, 1, 2, 4, 3, 2, 1, 3, 4, 1, 2, 1, 2, 3, 4, 2, 1, 4, 3), 
                                    `7` = c(1, 2, 3, 4, 3, 4, 1, 2, 1, 2, 3, 4, 3, 4, 1, 2, 2, 1, 4, 3, 4, 3, 2, 1, 2, 1, 4, 3, 4, 3, 2, 1), `8` = c(1, 2, 3, 4, 3, 4, 1, 2, 2, 1, 4, 3, 4, 3, 2, 1, 3, 4, 1, 2, 1, 2, 3, 4, 4, 3, 2, 1, 2, 1, 4, 3), 
                                    `9` = c(1, 2, 3, 4, 4, 3, 2, 1, 3, 4, 1, 2, 2, 1, 4, 3, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1, 1, 2, 3, 4), `10` = c(1, 2, 3, 4, 4, 3, 2, 1, 4, 3, 2, 1, 1, 2, 3, 4, 3, 4, 1, 2, 2, 1, 4, 3, 2, 1, 4, 3, 3, 4, 1, 2)), 
                ia_table = data.frame(`1` = c("linearer Graph")))
.L36_2_3_a = list(id = "L36_2_3_a", type = "mixed", niv_max = 3, runs = 36, anzahl_spalten = 23, levels1 = 2, levels2 = 3, factors1 = 11, factors2 = 12, 
                  design = data.frame(`1` = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), `2` = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2), 
                                      `3` = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1), `4` = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1), 
                                      `5` = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2), `6` = c(1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1), 
                                      `7` = c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2), `8` = c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1), 
                                      `9` = c(1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2), `10` = c(1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
                                      `11` = c(1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1), `12` = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), 
                                      `13` = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2), `14` = c(1, 2, 3, 1, 2, 3, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 1, 2, 3, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 1, 2, 3), 
                                      `15` = c(1, 2, 3, 1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 2, 3, 1, 3, 1, 2, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1), `16` = c(1, 2, 3, 2, 3, 1, 1, 2, 3, 1, 2, 3, 3, 1, 2, 1, 2, 3, 3, 1, 2, 3, 1, 2, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2), 
                                      `17` = c(1, 2, 3, 2, 3, 1, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 1, 2, 3, 3, 1, 2, 1, 2, 3, 3, 1, 2, 2, 3, 1), `18` = c(1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2), 
                                      `19` = c(1, 2, 3, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 2, 3, 1, 2, 3, 1, 1, 2, 3, 1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3), `20` = c(1, 2, 3, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3,  1, 1, 2, 3, 2, 3, 1), 
                                      `21` = c(1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 2, 3, 3, 1, 2, 1, 2, 3, 3, 1, 2, 2, 3, 1, 2, 3, 1), `22` = c(1, 2, 3, 3, 1, 2, 2, 3, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 1, 2, 3, 3, 1, 2), 
                                      `23` = c(1, 2, 3, 3, 1, 2, 3, 1, 2, 2, 3, 1, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3, 1, 2, 3, 1, 3, 1, 2, 1, 2, 3, 1, 2, 3)), ia_table = data.frame(`1` = c("ohne linearen Graph")))
.L36_2_3_b = list(id = "L36_2_3_b", type = "mixed", niv_max = 3, runs = 36, anzahl_spalten = 16, levels1 = 2, levels2 = 3, factors1 = 3, factors2 = 13, 
                  design = data.frame(`1` = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2), `2` = c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2), 
                                      `3` = c(1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1), `4` = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3), 
                                      `5` = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), `6` = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2), 
                                      `7` = c(1, 2, 3, 1, 2, 3, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 1, 2, 3, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 1, 2, 3), `8` = c(1, 2, 3, 1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 2, 3, 1, 3, 1, 2, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1), 
                                      `9` = c(1, 2, 3, 2, 3, 1, 1, 2, 3, 1, 2, 3, 3, 1, 2, 1, 2, 3, 3, 1, 2, 3, 1, 2, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2), `10` = c(1, 2, 3, 2, 3, 1, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 1, 2, 3, 3, 1, 2, 1, 2, 3, 3, 1, 2, 2, 3, 1), 
                                      `11` = c(1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2), `12` = c(1, 2, 3, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 2, 3, 1, 2, 3, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3, 1, 1, 2, 3), 
                                      `13` = c(1, 2, 3, 3, 1, 2, 1, 2, 3, 2, 3, 1, 3, 1, 2, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 2, 3, 1), `14` = c(1, 2, 3, 3, 1, 2, 2, 3, 1, 1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 2, 3, 3, 1, 2, 1, 2, 3, 3, 1, 2, 2, 3, 1, 2, 3, 1), 
                                      `15` = c(1, 2, 3, 3, 1, 2, 2, 3, 1, 3, 1, 2, 1, 2, 3, 2, 3, 1, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 1, 2, 3, 3, 1, 2), `16` = c(1, 2, 3, 3, 1, 2, 3, 1, 2, 2, 3, 1, 2, 3, 1, 1, 2, 3, 3, 1, 2, 2, 3, 1, 2, 3, 1, 3, 1, 2, 1, 2, 3, 1, 2, 3)), 
                  ia_table = data.frame(`1` = c("linearer Graph")))
.L50_2_5 = list(id = "L50_2_5", type = "mixed", niv_max = 5, runs = 15, anzahl_spalten = 12, levels1 = 2, levels2 = 5, factors1 = 1, factors2 = 11, 
                design = data.frame(`1` = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), 
                                    `2` = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5), 
                                    `3` = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5), 
                                    `4` = c(1, 2, 3, 4, 5, 2, 3, 4, 5, 1, 3, 4, 5, 1, 2, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 5, 1, 3, 4, 5, 1, 2, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4), 
                                    `5` = c(1, 2, 3, 4, 5, 3, 4, 5, 1, 2, 5, 1, 2, 3, 4, 2, 3, 4, 5, 1, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 1, 2, 3, 4, 5, 3, 4, 5, 1, 2, 5, 1, 2, 3, 4, 2, 3, 4, 5, 1), 
                                    `6` = c(1, 2, 3, 4, 5, 4, 5, 1, 2, 3, 2, 3, 4, 5, 1, 5, 1, 2, 3, 4, 3, 4, 5, 1, 2, 5, 1, 2, 3, 4, 3, 4, 5, 1, 2, 1, 2, 3, 4, 5, 4, 5, 1, 2, 3, 2, 3, 4, 5, 1), 
                                    `7` = c(1, 2, 3, 4, 5, 5, 1, 2, 3, 4, 4, 5, 1, 2, 3, 3, 4, 5, 1, 2, 2, 3, 4, 5, 1, 4, 5, 1, 2, 3, 3, 4, 5, 1, 2, 2, 3, 4, 5, 1, 1, 2, 3, 4, 5, 5, 1, 2, 3, 4), 
                                    `8` = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4, 4, 5, 1, 2, 3, 3, 4, 5, 1, 2, 2, 3, 4, 5, 1, 5, 1, 2, 3, 4, 2, 3, 4, 5, 1, 3, 4, 5, 1, 2), 
                                    `9` = c(1, 2, 3, 4, 5, 2, 3, 4, 5, 1, 1, 2, 3, 4, 5, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 2, 3, 4, 5, 1, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 4, 5, 1, 2, 3), 
                                    `10` = c(1, 2, 3, 4, 5, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 1, 2, 3, 4, 5, 2, 3, 4, 5, 1, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 4, 5, 1, 2, 3, 2, 3, 4, 5, 1, 4, 5, 1, 2, 3), 
                                    `11` = c(1, 2, 3, 4, 5, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4, 4, 5, 1, 2, 3, 1, 2, 3, 4, 5, 2, 3, 4, 5, 1, 5, 1, 2, 3, 4, 2, 3, 4, 5, 1, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2), 
                                    `12` = c(1, 2, 3, 4, 5, 5, 1, 2, 3, 4, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 5, 1, 2, 3, 4, 3, 4, 5, 1, 2, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 3, 4, 5, 1, 2, 1, 2, 3, 4, 5)), ia_table = data.frame(`1` = c("linearer Graph")))
.L8_4_2 = list(id = "L8_4_2", type = "mixed", niv_max = 4, runs = 8, anzahl_spalten = 5, levels1 = 4, levels2 = 2, factors1 = 1, factors2 = 4, design = data.frame(`1` = c(1,  1, 2, 2, 3, 3, 4, 4), `2` = c(1, 2, 1, 2, 1, 2, 1, 2), `3` = c(1, 2, 1, 2, 2, 1, 2, 1), `4` = c(1, 2, 2, 1, 1, 2, 2, 1), `5` = c(1, 2, 2, 1, 2, 1, 1, 2)), 
               ia_table = data.frame(`1` = c("ohne interaction table")))
.L16_4_2_a = list(id = "L16_4_2_a", type = "mixed", niv_max = 4, runs = 16, anzahl_spalten = 13, levels1 = 4, levels2 = 2, factors1 = 1, factors2 = 12, 
                  design = data.frame(`1` = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), `2` = c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2), `3` = c(1, 1, 2, 2,1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1), `4` = c(1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1), `5` = c(1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2), 
                                      `6` = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2), `7` = c(1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1), `8` = c(1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1), `9` = c(1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2), `10` = c(1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1),
                                      `11` = c(1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2), `12` = c(1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2), `13` = c(1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1)), ia_table = data.frame(`1` = c("ohne interaction table")))
.L16_4_2_b = list(id = "L16_4_2_b", type = "mixed", niv_max = 4, runs = 16, anzahl_spalten = 11, levels1 = 4, levels2 = 2, factors1 = 2, factors2 = 9, 
                  design = data.frame(`1` = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), `2` = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4), `3` = c(1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1), `4` = c(1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1), `5` = c(1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2), 
                                      `6` = c(1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1), `7` = c(1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1), `8` = c(1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2), `9` = c(1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2), `10` = c(1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2), 
                                      `11` = c(1,  2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1)), ia_table = data.frame(`1` = c("ohne interaction table")))
.L16_4_2_c = list(id = "L16_4_2_c", type = "mixed", niv_max = 4, runs = 16, anzahl_spalten = 9, levels1 = 4, levels2 = 2, factors1 = 3, factors2 = 6, 
                  design = data.frame(`1` = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), `2` = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4), `3` = c(1, 2, 3, 4,2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1), `4` = c(1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1), `5` = c(1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2), 
                                      `6` = c(1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1), `7` = c(1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2), `8` = c(1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2), `9` = c(1, 2, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2)), ia_table = data.frame(`1` = c("ohne interaction table")))
.L16_4_2_d = list(id = "L16_4_2_d", type = "mixed", niv_max = 4, runs = 16, anzahl_spalten = 7, levels1 = 4, levels2 = 2, factors1 = 4, factors2 = 3, 
                  design = data.frame(`1` = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), `2` = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4), `3` = c(1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1), `4` = c(1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1), `5` = c(1, 2, 3, 4, 1, 2, 4, 3, 4, 3, 2, 1, 2, 1, 4, 3), 
                                      `6` = c(1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2), `7` = c(1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2)), ia_table = data.frame(`1` = c("ohne interaction table")))
.L18_6_3 = list(id = "L18_6_3", type = "mixed", niv_max = 6, runs = 18, anzahl_spalten = 7, levels1 = 3, levels2 = 6, factors1 = 6, factors2 = 1, 
                design = data.frame(`1` = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6), `2` = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), `3` = c(1, 2, 3, 1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2), `4` = c(1, 2, 3, 2, 3, 1, 1, 2, 3, 3, 1, 2, 3, 1, 2, 2, 3, 1), `5` = c(1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 3, 1, 1, 2, 3, 3, 1, 2), 
                                    `6` = c(1, 2, 3, 3, 1, 2, 2, 3, 1, 2, 3, 1, 3, 1, 2, 1, 2, 3), `7` = c(1, 2, 3, 3, 1, 2, 3, 1, 2, 1, 2, 3, 2, 3, 1, 2, 3, 1)), ia_table = data.frame(`1` = c("ohne interaction table")))
.oaList = list(.L4_2, .L8_2, .L9_3, .L12_2, .L16_2, .L16_4, .L18_2_3, .L25_5, .L27_3, .L32_2, .L32_2_4, .L36_2_3_a, .L36_2_3_b, .L50_2_5, .L8_4_2, 
               .L16_4_2_a, .L16_4_2_b, .L16_4_2_c, .L16_4_2_d, .L18_6_3) 









