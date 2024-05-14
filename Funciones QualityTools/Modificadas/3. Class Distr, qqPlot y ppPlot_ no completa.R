library(R6)
library(patchwork)
library(ggplot2)

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

                   plot = function(main = NULL, xlab = NULL, xlim = NULL, ylim = NULL, ylab = NULL, line.col = "red", box=TRUE,line.width = 1, ...)
                   {
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
                       labs(y = ylab, x = xlab, title = main) + xlim(xlim) + ylim(ylim) +
                       theme_minimal() + theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
                       guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
                       geom_line(data = data.frame(x = xPoints, y = yPoints), aes(x = x, y = y), color = line.col, linewidth = line.width) + # densidad
                       theme(legend.position = "none")

                     # Caja de Info
                     if (box==FALSE) {
                       p1
                     }
                     else {
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
                                      parse = TRUE,size = 3,hjust = 0)
                         }
                         if (!is.null(adTestStats$smaller) && !adTestStats$smaller){
                           p2 <- p2 +
                             annotate('text',x = 0.25, y = 0.20,
                                      label = paste("p>=", round(as.numeric(p),digits = 3)),
                                      parse = TRUE,size = 3,hjust = 0)
                         }
                         if (is.null(adTestStats$smaller)){
                           p2 <- p2 +
                             annotate('text',x = 0.25,y = 0.20,
                                      label = paste("p==", round(as.numeric(p), digits = 3)),
                                      parse = TRUE,size = 3,hjust = 0)
                         }

                         # mean y sd
                         p2 <- p2 + annotate('text', x = 0.25, y = 0.35,
                                             label = paste("mean==", round(self$parameters[[1]], digits = 3)),
                                             parse = TRUE, size = 3, hjust = 0) +
                           annotate('text', x = 0.25, y = 0.30,
                                    label = paste("sd==", round(self$parameters[[2]], digits = 3)),
                                    parse = TRUE, size = 3, hjust = 0)
                         }

                       p1 + inset_element(p2, left = 0.7, right = 1, top = 1, bottom = 0.60)
                     }
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
                                 plot = function(xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, line.col = "red", line.width = 1, box = TRUE , ...) {
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
                                   p <- distrList[[1]]$plot(xlab = xlab, ylab = ylab, line.col = line.col, line.width = line.width, box = box)
                                   for (i in 2:length(distrList)) {
                                     p <- p+distrList[[i]]$plot(xlab = xlab, ylab = ylab, line.col = line.col, line.width = line.width, box = box)
                                   }
                                   p
                                 }
                               )
)


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

# Funciones auxiliares
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
# Funcion qqPlot ---------------------
qqPlot <- function(x, y, confbounds = TRUE, alpha, main, xlab, ylab, xlim, ylim, border = "red", bounds.col = "black", bounds.lty = 1, start, grapic = TRUE, axis.y.right = FALSE, bw.theme = FALSE,...){
  DB = FALSE
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
      if (DB)
        print("TODO: Pass the estimated parameters correctly")
      fitList = .lfkp(parList, formals(qFun))
      fitList$x = xs
      fitList$densfun = distribution
      if(!missing(start))
        fitList$start = start
      if(DB){
        print(fitList)
        print("Ende")
      }
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
  if(missing(main))
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
      }
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
