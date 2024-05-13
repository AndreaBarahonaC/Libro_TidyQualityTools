library(R6)
library(MASS)
library(patchwork)

# Class Distr ----
Distr <- R6Class("Distr",
                 public = list(
                   x = NULL,
                   name = NULL,
                   parameters = NULL,
                   sd = NULL,
                   n = NULL,
                   loglik = NULL,

                   initialize = function(x, name, parameters, sd, n, loglik) {
                     self$x <- x
                     self$name <- name
                     self$parameters <- parameters
                     self$sd <- sd
                     self$n <- n
                     self$loglik <- loglik
                   },

                   plot = function(main = NULL, xlab = NULL, xlim = NULL, ylim = NULL, ylab = NULL, line.col = "red", line.width = 1, box = TRUE, ...) {
                     object <- self
                     xVals <- object$x
                     parameters <- object$parameters
                     lq <- NULL
                     uq <- NULL
                     y <- NULL

                     if (missing(line.col)) {
                       line.col <- "red"
                     }
                     if (missing(line.width)) {
                       line.width <- 1
                     }
                     if (missing(main)) {
                       main <- object$name
                     }
                     if (missing(xlab)) {
                       xlab <- "x"
                     }
                     if (missing(ylab)) {
                       ylab <- "Density"
                     }

                     distr <- object$name
                     qFun <- .charToDistFunc(distr, type = "q")
                     dFun <- .charToDistFunc(distr, type = "d")
                     adTestStats <- .myADTest(xVals, distr)

                     if (class(adTestStats) == "adtest") {
                       A <- adTestStats$statistic
                       p <- adTestStats$p.value
                     } else {
                       A <- NA
                       p <- NA
                     }

                     histObj <- hist(xVals, plot = FALSE)
                     df <- data.frame(
                       mid = histObj$mids,
                       density = histObj$density
                     )
                     width <- diff(df$mid)[1]

                     if (missing(xlim)) {
                       lq <- do.call(qFun, c(list(1e-04), as.list(parameters)))
                       uq <- do.call(qFun, c(list(0.9999), as.list(parameters)))
                       xlim <- range(lq, uq, xVals)
                     }

                     xPoints <- seq(xlim[1], xlim[2], length = 200)
                     yPoints <- do.call(dFun, c(list(xPoints), as.list(parameters)))

                     if (missing(ylim)) {
                       ylim <- range(0, histObj$density, yPoints)
                     }

                     # Histograma
                     p1 <- ggplot(df, aes(x = mid, y = density)) +
                       geom_bar(stat = "identity", width = width, fill = "lightblue", color = "black", alpha = 0.5) +
                       labs(y = "Density", x = "x", title = main) +
                       theme_minimal() + theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
                       guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
                       geom_line(data = data.frame(x = xVec, y = yVec), aes(x = x, y = y), color = "red", linewidth = 0.5) + # densidad
                       theme(legend.position = "none")

                     # Caja de Info
                     p2 <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
                       theme_bw() +
                       theme(
                         axis.text = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()
                       ) +
                       xlim(c(0.25,0.26)) + ylim(c(0.19, 0.36))
                     {
                       # n y A
                       p2 <- p2 +
                         annotate('text', x = 0.25, y = 0.25,
                                  label = paste("A==", round(as.numeric(A), digits = 3)),
                                  parse = TRUE, size = 3, hjust = 0)

                       # p
                       if (!is.null(adTestStats$smaller) && adTestStats$smaller){
                         p2 <- p2 +
                           annotate('text',x = 0.25,y = 0.20,
                           label = paste("p<", round(as.numeric(p), digits =3)),
                           parse = TRUE,size = 3,hjust = 0
                         )
                       }
                       if (!is.null(adTestStats$smaller) && !adTestStats$smaller){
                         p2 <- p2 +
                           annotate('text',x = 0.25, y = 0.20,
                           label = paste("p>=", round(as.numeric(p),digits = 3)),
                           parse = TRUE,size = 3,hjust = 0
                         )
                       }
                       if (is.null(adTestStats$smaller)){
                         p2 <- p2 +
                           annotate('text',x = 0.25,y = 0.20,
                           label = paste("p==", round(as.numeric(p), digits = 3)),
                           parse = TRUE,size = 3,hjust = 0
                         )
                       }

                       # mean y sd
                       p2 <- p2 + annotate('text', x = 0.25, y = 0.35,
                                           label = paste("mean==", round(parameters[[1]], digits = 3)),
                                           parse = TRUE, size = 3, hjust = 0) +
                         annotate('text', x = 0.25, y = 0.30,
                                  label = paste("sd==", round(parameters[[2]], digits = 3)),
                                  parse = TRUE, size = 3, hjust = 0)
                       }

                     p1 + inset_element(p2, left = 0.7, right = 1, top = 1, bottom = 0.60)
                   }
                 )
)

# Class DistrCollection ----
DistrCollection <- R6::R6Class("DistrCollection",
                               public = list(
                                 distr = NULL,
                                 initialize = function() {
                                   self$distr <- list()
                                 },
                                 add = function(distr) {
                                   self$distr <- append(self$distr, list(distr))
                                 },
                                 get = function(i) {
                                   self$distr[[i]]
                                 },
                                 show = function() {
                                   cat("\n")
                                   for (i in seq_along(self$distr)) {
                                     temp <- self$distr[[i]]
                                     cat("\n")
                                     cat("fitted distribution is", temp$name, ":\n")
                                     print(temp$parameters)
                                     cat("\n")
                                   }
                                 },
                                 summary = function() {
                                   numDist <- length(self$distr)
                                   gofMatrix <- data.frame(matrix(nrow = numDist, ncol = 3))
                                   names(gofMatrix) <- c("Distribution", "A", "p.value")
                                   cat("\n------ Fitted Distribution and estimated parameters ------\n")
                                   for (i in seq_along(self$distr)) {
                                     distrObj <- self$distr[[i]]
                                     x <- distrObj$x
                                     distribution <- distrObj$name
                                     parameters <- distrObj$parameters
                                     statistic <- NA
                                     p.value <- NA
                                     temp <- .myADTest(x, distribution)
                                     try(statistic <- as.numeric(temp$statistic), silent = TRUE)
                                     try(p.value <- as.numeric(temp$p.value), silent = TRUE)
                                     gofMatrix[i, ] <- c(distribution, as.numeric(statistic), as.numeric(p.value))
                                     cat("\n")
                                     cat("fitted distribution is", distribution, ":\n")
                                     print(parameters)
                                   }
                                   cat("\n")
                                   cat("\n------ Goodness of Fit - Anderson Darling Test ------\n")
                                   cat("\n")
                                   gofMatrixPrint <- gofMatrix
                                   gofMatrixPrint[, 2] <- signif(as.numeric(gofMatrixPrint[, 2]), 4)
                                   gofMatrixPrint[, 3] <- signif(as.numeric(gofMatrixPrint[, 3]), 4)
                                   print(gofMatrixPrint)
                                 },
                                 plot = function(xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, line.col = "red", line.width = 1, ...) {
                                   distrList <- self$distr
                                   numDist <- length(self$distr)
                                   numColWin <- ceiling(numDist/2)
                                   if (missing(xlim)) {
                                     xlim <- .xyLimits(self)$xlim
                                   }
                                   if (missing(ylim)) {
                                     ylim <- .xyLimits(self)$ylim
                                   }
                                   if (missing(line.col)) {
                                     line.col <- "red"
                                   }
                                   if (missing(line.width)) {
                                     line.width <- 1
                                   }
                                   lapply(distrList, function(d) plot(d, xlim = xlim, ylim = ylim, line.col = line.col, line.width = line.width, ...))
                                   cat(paste("Total of", numDist, "plots created"))
                                   cat("\n")
                                   cat(paste("Use par(mfrow = c(2,", numColWin, ") to see all of them!", sep = ""))
                                   cat("\n")
                                 }
                               )
)

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
                            vcov = vc, n = n, loglik = sum(dlnorm(x, mx, sd0, log = TRUE))), class = "fitdistr"))
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
                                                                 sd0, log = TRUE))), class = "fitdistr"))
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
                                                                 log = TRUE))), class = "fitdistr"))
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
                                                                log = TRUE))), class = "fitdistr"))
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
                                                                 log = TRUE))), class = "fitdistr"))
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
                 loglik = -res$value, n = n), class = "fitdistr")
}
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
# Funcion .xyLimits ----
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



# Funcion qqPlot falta cambiaaaar desde aquiiii ---------------------
qqPlot_o <- function(x, y, confbounds = TRUE, alpha, main, xlab, ylab, xlim, ylim, border = "red", bounds.col = "black", bounds.lty = 1, start, ...) {
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

  distr_coll <- DistrCollection$new()
  if (inherits(x, "DistrCollection")) {
    aux <- FitDistr(x[,1],distribution)
    new_dis <- Distr$new(x[,1],
                         name = distribution,
                         parameters = aux$estimate,
                         sd = aux$sd,
                         n = aux$n,
                         loglik = aux$loglik
    )
    distr_coll$add(new_dis)
    for (i in seq_along(distr_coll$distr)) {
      d <- distr_coll$distr[[i]]
      do.call(qqPlot, c(list(x = d$x[,1], y = d$name), parList))
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
  xs = sort(x[,1])
  distribution = tolower(y)
  distWhichNeedParameters = c("weibull", "logistic", "gamma",
                              "exponential", "f", "geometric", "chi-squared", "negative binomial",
                              "poisson")

  threeParameterDistr = c("weibull3", "lognormal3", "gamma3")
  threeParameter = distribution %in% threeParameterDistr
  if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)

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
    if(!threeParameter){
      fittedDistr = do.call(fitdistr, fitList)
      parameter = fittedDistr$estimate

      thethas = fittedDistr$estimate
      varmatrix = fittedDistr$vcov

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
    if (DB) {
      print(paste("parameter: ", parameter))
      print(xq)
    }
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
  pParams =.lfkp(pParams, list(x = 1, y = 1, col = 1, cex = 1))
  do.call(points, pParams)
  params =.lfkp(parList, c(formals(abline), par()))
  params$a = 0
  params$b = 1
  params$col = border
  do.call(abline, params)

  if(!threeParameter){
    if(confbounds == TRUE){
      if(distribution %in% confIntCapable){
        params =.lfkp(parList, c(formals(lines), par()))
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
                    ...)
{
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
