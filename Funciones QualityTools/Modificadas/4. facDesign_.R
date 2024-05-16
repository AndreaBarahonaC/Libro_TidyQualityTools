library(R6)
library(MASS)

####Necesito .helpAliasTable##################################################

#Funcion .replace2s#
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
  if (class(fdo)[1] == "facDesign")
    X = unique(fdo$cube)
  if (class(fdo)[1] == "taguchiDesign") {
    X = unique(fdo$design)
    X = .replace2s(X)
  }
  N = nrow(X)
  columns = names(X[, 1:k])
  X1 = matrix(1, nrow = N, ncol = 1)
  nameVec = c("Identity")
  for (i in 1:degree) {
    temp = combn(columns, i)
    for (j in 1:ncol(temp)) {
      if (class(fdo)[1] == "facDesign")
        index = names(X) %in% temp[, j]
      if (class(fdo)[1] == "taguchiDesign")
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

####Necesito aliasTable#######################################################
aliasTable <- function (fdo, degree, show = TRUE)
{
  if (class(fdo)[1] == "facDesign") {
    X = unique(fdo$cube)
    N = nrow(X)
    k = log2(N)
    kPlusP = ncol(X)
    if (missing(degree))
      degree = min(c(4, k + 1))
    X1 = .helpAliasTable(fdo, k, degree = degree - 1)
    X2 = .helpAliasTable(fdo, k = kPlusP, degree)
  }
  if (class(fdo)[1] == "taguchiDesign") {
    if (length(table(as.numeric(as.matrix(fdo$design)))) !=
        2)
      stop("calculation of an alias table for mixed designs is not supported")
    k = ncol(fdo$design)
    if (missing(degree))
      degree = min(c(3, k))
    X1 = unique(fdo$design)
    X1 = .replace2s(X1)
    X2 = .helpAliasTable(fdo, k, degree)
    X1 = cbind(data.frame(Identity = rep(1, times = nrow(X1))),
               X1)
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

###Necesito .fdoOrth y .NAMES#######################################
.fdoOrth = vector(mode = "list", length = 3)
###
.fdoOrth[[1]] = list(k = 3, gen = "C=AB", p = 1)
###
.fdoOrth[[2]] = list(k = 4, gen = "D=ABC", p = 1)
###
.fdoOrth[[3]] = list(k = 5, gen = c("D=AB","E=AC"), p = 2)
###
.fdoOrth[[4]] = list(k = 6, gen = c("D=AB","E=AC","F=BC"), p = 3)
###
.fdoOrth[[5]] = list(k = 7, gen = c("D=AB","E=AC","F=BC","G=ABC"), p = 4)
###
.fdoOrth[[6]] = list(k = 5, gen = "E=ABCD", p = 1)
###
.fdoOrth[[7]] = list(k = 6, gen = c("E=ABC","F=BCD"), p = 2)
###
.fdoOrth[[8]] = list(k = 7, gen = c("E=ABC","F=BCD","G=ACD"), p = 3)
###
.fdoOrth[[9]] = list(k = 8, gen = c("E=BCD","F=ACD","G=ABC","H=ABD"), p = 4)
###
.fdoOrth[[10]] = list(k = 9, gen = c("E=ABC","F=BCD","G=ACD","H=ABD","J=ABCD"), p = 5)
###
.fdoOrth[[11]] = list(k = 10, gen = c("E=ABC","F=BCD","G=ACD","H=ABD","J=ABCD","K=AB"), p = 6)
###
.fdoOrth[[12]] = list(k = 11, gen = c("E=ABC","F=BCD","G=ACD","H=ABD","J=ABCD","K=AB","L=AC"), p = 7)
###
.fdoOrth[[13]] = list(k = 6, gen = "F=ABCDE", p = 1)
###
.fdoOrth[[14]] = list(k = 7, gen = c("F=ABCD","G=ABDE"), p = 2)
###
.fdoOrth[[15]] = list(k = 8, gen = c("F=ABC","G=ABD","H=BCDE"), p = 3)
###
.fdoOrth[[16]] = list(k = 9, gen = c("F=BCDE","G=ACDE","H=ABDE","J=ABCE"), p = 4)
###
.fdoOrth[[17]] = list(k = 10, gen = c("F=ABCD","G=ABCE","H=ABDE","J=ACDE","K=BCDE"), p = 5)
###
.fdoOrth[[18]] = list(k = 11, gen = c("F=ABC","G=BCD","H=CDE","J=ACD","K=AEF","L=ADEF"), p = 6)
###
.fdoOrth[[19]] = list(k = 7, gen = "G=ABCDEF", p = 1)
###
.fdoOrth[[20]] = list(k = 8, gen = c("G=ABCD","H=ABEF"), p = 2)
###
.fdoOrth[[21]] = list(k = 9, gen = c("G=ABCD","H=ABEF","J=CDEF"), p = 3)
###
.fdoOrth[[22]] = list(k = 10, gen = c("G=BCDF","H=ACDF","J=ABDE","K=ABCE"), p = 4)
###
.fdoOrth[[23]] = list(k = 11, gen = c("G=CDE","H=ABCD","J=ABF","K=BDEF","L=ADEF"), p = 5)
###
.fdoOrth[[24]] = list(k = 8, gen = "H=ABCDEFG", p = 1)
###
.fdoOrth[[25]] = list(k = 9, gen = c("H=ACDFG","J=BCEFG"), p = 2)
###
.fdoOrth[[26]] = list(k = 10, gen = c("H=ABCG","J=BCDE","K=ACDF"), p = 3)
###
.fdoOrth[[27]] = list(k = 11, gen = c("H=ABCG","J=BCDE","K=ACDF","L=ABCDEFG"), p = 4)
###
.NAMES = LETTERS[c(1:8, 10:26)]



####Necesito .m.interaction.plot###########################
.m.interaction.plot <- function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = TRUE, trace.label = deparse(substitute(trace.factor)),
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

##Necesito clase facDesign.c######################################################
facDesign.c <- R6Class("facDesign", public = list(name = NULL,
                                                 factors = NULL,
                                                 cube = data.frame(),
                                                 star = data.frame(),
                                                 centerCube = data.frame(),
                                                 centerStar = data.frame(),
                                                 generator = NULL,
                                                 response = NULL,
                                                 block = NULL,
                                                 blockGen = NULL,
                                                 runOrder = NULL,
                                                 standardOrder = NULL,
                                                 desireVal = NULL,
                                                 desirability = NULL,
                                                 fits = NULL,

                                                 names = function(value){
                                                   if(missing(value)){
                                                     n <- c()
                                                     for (i in 1:length(self$factors)) {
                                                       n[i] <- self$factors[[i]]$name
                                                     }
                                                     return(n)
                                                   }
                                                   else {
                                                     for (i in 1:length(self$factors)){
                                                       self$factors[[i]]$name = as.character(value[i])
                                                     }
                                                     invisible(self)
                                                   }

                                                 },

                                                 as.data.frame = function(row.names = NULL, optional = FALSE, ...) {
                                                   if (!is.null(self$cube)) {
                                                     frameOut = self$cube
                                                   }
                                                   else return(NULL)
                                                   if (!is.null(self$centerCube))
                                                     frameOut = rbind(frameOut, self$centerCube)
                                                   if (!is.null(self$star))
                                                     frameOut = rbind(frameOut, self$star)
                                                   if (!is.null(self$centerStar))
                                                     frameOut = rbind(frameOut, self$centerStar)
                                                   if (!is.null(self$factors) && length(self$factors) == dim(frameOut)[2]) {
                                                     names(frameOut) = as.character(names(self$cube))
                                                   }
                                                   if (!is.null(self$blockGen) && nrow(self$blockGen) > 0) {
                                                     frameOut = cbind(self$blockGen, frameOut)
                                                   }
                                                   if (!is.null(self$block) && nrow(self$block) > 0) {
                                                     frameOut = cbind(self$block, frameOut)
                                                   }
                                                   if (!is.null(self$runOrder) && nrow(self$runOrder) > 0) {
                                                     frameOut = cbind(self$runOrder, frameOut)
                                                   }
                                                   if (!is.null(self$standardOrder) && nrow(self$standardOrder) > 0) {
                                                     frameOut = cbind(self$standardOrder, frameOut)
                                                   }
                                                   if (!is.null(self$response) && nrow(frameOut) == nrow(self$response))
                                                     frameOut = cbind(frameOut, self$response)
                                                   else {
                                                     temp = as.data.frame(matrix(NA, nrow = nrow(frameOut), ncol = ncol(self$response)))
                                                     names(temp) = names(self$response)
                                                     frameOut = cbind(frameOut, temp)
                                                   }
                                                   runIndex = order(self$runOrder[,1])
                                                   out = frameOut[runIndex, ]
                                                   return(out)
                                                 },

                                                 nrow = function(){
                                                   return(nrow(self$as.data.frame()))
                                                 },

                                                 get = function(i,j){
                                                   return(self$as.data.frame()[i, j])
                                                 },

                                                 lows = function(value){
                                                   if (missing(value)) {
                                                     listOut = vector(mode = "list")
                                                     for (i in seq(along = self$factors)) {
                                                       listOut[self$factors[[i]]$name] = self$factors[[i]]$.low()
                                                     }
                                                     return(listOut)
                                                   }
                                                   else {
                                                     for (i in seq(along = self$factors)) {
                                                       self$factors[[i]]$.low(value[i])
                                                     }
                                                     invisible(self)
                                                   }
                                                 },

                                                 highs = function(value){
                                                   if (missing(value)) {
                                                     listOut = vector(mode = "list")
                                                     for (i in seq(along = self$factors)) {
                                                       listOut[self$factors[[i]]$name] = self$factors[[i]]$.high()
                                                     }
                                                     return(listOut)
                                                   }
                                                   else {
                                                     for (i in seq(along = self$factors)) {
                                                       self$factors[[i]]$.high(value[i])
                                                     }
                                                     invisible(self)
                                                   }
                                                 },

                                                 .nfp = function(){
                                                   x = self$factors
                                                   atr <- c('low','high','name','unit','type')
                                                   if (is.list(x) && length(x[[1]]) > 0) {
                                                     numAttr = length(x[[1]]$attributes())
                                                     .numFac = length(x)
                                                     frameOut = data.frame(matrix(ncol = .numFac, nrow = numAttr ))
                                                     for (i in 1:numAttr ) {
                                                       charVec = character(0)
                                                       for (j in 1:.numFac) {
                                                         charVec = c(charVec, atr[i], "\t\t")
                                                         frameOut[i, j] = x[[j]]$attributes()[i]
                                                       }
                                                     }
                                                     names(frameOut) = self$names()
                                                     rownames(frameOut) = atr[1:numAttr ]
                                                   }
                                                   else {
                                                     stop("no list given or length of list < 1")
                                                   }
                                                   print(frameOut)
                                                 },

                                                 identity = function(){
                                                   identity = character(0)
                                                   identityList = vector(mode = "list", length = 0)
                                                   resolution = numeric(0)
                                                   temp = NULL
                                                   A = aliasTable(self, show = FALSE)
                                                   if (any(dim(A) == 0))
                                                     return(identityList)
                                                   temp = as.matrix(A["Identity", ])
                                                   boolTemp = apply(temp, 2, as.logical)
                                                   identity = row.names(temp)[boolTemp[, 1]]
                                                   if (length(identity) > 0) {
                                                     charList = strsplit(toupper(identity), split = "")
                                                     identityList = lapply(charList, match, .NAMES[1:25])
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
                                                 },

                                                 summary = function(){
                                                   doeFactors = self$factors
                                                   cat("Information about the factors:\n\n")
                                                   self$.nfp()
                                                   cat("-----------\n")
                                                   print(self$as.data.frame())
                                                   temp = aliasTable(self, show = FALSE)
                                                   if (ncol(temp) > 0) {
                                                     cat("\n---------\n\n")
                                                     identity(self)
                                                     cat("\n")
                                                   }
                                                   invisible(self$as.data.frame())
                                                 },

                                                 .response = function(value){
                                                   if(missing(value)){
                                                     iIntern <- order(self$runOrder[,1])
                                                     out <- data.frame(self$response[iIntern,])
                                                     names(out) <- names(self$response)
                                                     return(out)
                                                   }
                                                   else{
                                                     index = order(self$runOrder[,1])
                                                     if (!is.vector(value) && !is.data.frame(value))
                                                       stop("vector or data.frame expected!")
                                                     if (is.vector(value) && (is.numeric(value) || is.na(value))) {
                                                       if (nrow(self$response) != length(value))
                                                         stop(paste("Number of rows for Design does not equal length of vector ", nrow(object),
                                                                    " != ", length(value), " "))
                                                       self$response <- data.frame(value)
                                                       self$response[index, ] <- value
                                                       names(self$response) <- make.names(deparse(substitute(value)))
                                                       invisible(self)
                                                     }
                                                     if (is.data.frame(value)) {
                                                       self$response <- value
                                                       self$response[index, ] <- value
                                                       invisible(self)
                                                     }

                                                   }

                                                 },

                                                 effectPlot = function(factors, fun = mean, response = NULL, single = FALSE, points = FALSE, classic = FALSE, axes = TRUE, ###
                                                                       lty, xlab, ylab, main, ylim, ...) {
                                                   oldMar = par("mar")
                                                   oldOma = par("oma")
                                                   oldMfrow = par("mfrow")
                                                   oldMfcol = par("mfcol")
                                                   on.exit(par(mar = oldMar, oma = oldOma, mfrow = oldMfrow, mfcol = oldMfcol))
                                                   if(is.null(response)==FALSE)                                                ###
                                                   {                                                                           ###
                                                     temp=self$.response()[response]                                            ###
                                                     self$.response(temp)                                                      ###
                                                   }                                                                           ###
                                                   ylabmiss = FALSE
                                                   xlabmiss = FALSE
                                                   mainmiss = FALSE
                                                   ylimmiss = FALSE
                                                   if (missing(ylim))
                                                     ylimmiss = TRUE
                                                   if (missing(lty))
                                                     lty = 1
                                                   X = self$cube
                                                   Y = as.data.frame(self$response[1:nrow(X), ])
                                                   names(Y) = names(self$.response())
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
                                                         if (identical(" ", self$names()[[i]]))
                                                           xlab = factors[i]
                                                         else xlab = paste(factors[i], ": ", self$names()[[i]], sep = "")
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
                                                 },

                                                 lm = function(formula){
                                                   invisible(lm(formula, data = self$as.data.frame()))

                                                 },

                                                 desires = function(value){
                                                   if (!any(value$response == names(self$.response())))
                                                     stop(paste(value$response, "is not a response!"))
                                                   listPos = length(self$desirability) + 1
                                                   yName = value$response
                                                   isIn = (yName == names(self$desirability))
                                                   if (any(isIn))
                                                     listPos = (1:length(names(self$desirability)))[isIn]
                                                   x$desirability[[listPos]] = value
                                                   names(x$desirability)[listPos] = yName
                                                   invisible(x)
                                                 }

                                                 )
                      )



##Necesito clase doeFactor####################
doeFactor <- R6Class('doeFactor', public = list(low = -1,
                                                 high = 1,
                                                 name = "",
                                                 unit = "",
                                                 type = "numeric",

                                                 attributes = function(){
                                                   v <- c(self$low, self$high, self$name, self$unit, self$type)
                                                 },

                                                 .low = function(value){
                                                   if (missing(value)) {
                                                     return(unlist(self$low))
                                                   }
                                                   else{
                                                     boolOld = is.numeric(unlist(self$low))
                                                     self$low <- value
                                                     boolNew = is.numeric(self$low)
                                                     if (boolNew)
                                                       self$type = "numeric"
                                                     else self$type = "factor"
                                                     if (boolOld != boolNew)
                                                       print("Note: The types of the factors were changed!")
                                                     invisible(self)
                                                   }
                                                 },

                                                 .high = function(value){
                                                   if (missing(value)) {
                                                     return(unlist(self$high))
                                                   }
                                                   else{
                                                     boolOld = is.numeric(unlist(self$high))
                                                     self$high <- value
                                                     boolNew = is.numeric(self$high)
                                                     if (boolNew)
                                                       self$type = "numeric"
                                                     else self$type = "factor"
                                                     if (boolOld != boolNew)
                                                       print("Note: The types of the factors were changed!")
                                                     invisible(self)
                                                   }
                                                 }

                                                 )
                      )




##Necesito funcion randomize####
randomize <- function (fdo, random.seed, so = FALSE)
{
  if (missing(random.seed))
    set.seed(93275938)
  else set.seed(random.seed)
  j = 1
  temp = fdo$runOrder
  for (i in sort(unique(fdo$block[, 1]))) {
    pos = !is.na(match(fdo$block[, 1], i))
    count = sum(as.numeric(pos))
    if (so) {
      temp[pos, 1] = j:(j + (count - 1))
    }
    else {
      temp[pos, 1] = sample(j:(j + (count - 1)), count)
    }
    j = j + count
  }
  fdo$runOrder = temp
  return(fdo)
}


##Necesito .rsm#################
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


##Necesito .numFac#############
.numFac = function(fdo) {
  return(length(fdo$names()))
}


##Necesito .confoundings####
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


## Necesito .lociv####
.lociv = function(charVec) {
  lenVec = numeric(length = length(charVec))
  for (i in seq(along = charVec)) {
    lenVec[i] = length(strsplit(charVec[i], split = "")[[1]])
  }
  return(lenVec)
}

##Necesito . blockInteractions######
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
    if (!(nrow(unique(fdo$cube)) >= 2^.numFac(fdo)))
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
  lSet = fdo$names()
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



##Necesito .blockGenCol##########
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
      genTemp = .fdo$get(,gen[j])
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


####Necesito .blockCol##############
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


##Necesito funcion blocking######
blocking <- function (fdo, blocks, BoR = FALSE, random.seed, useTable = "rsm",
                      gen)
{
  override = FALSE
  Block = data.frame(Block = rep(1, fdo$nrow()))
  fdo$block = Block
  fdo = randomize(fdo, so = TRUE)
  if (missing(random.seed)) {
    runif(1)
    random.seed = .Random.seed[sample(1:626, 1)]
  }
  if (missing(gen))
    gen = NULL
  if (blocks <= 1) {
    Block = data.frame(Block = rep(1, fdo$nrow()))
    fdo$block = Block
    fdo = randomize(fdo, random.seed = random.seed)
    return(fdo)
  }
  if (nrow(fdo$star) > 0 | nrow(fdo$centerStar) > 0) {
    if (blocks == 2) {
      override = TRUE
      fdo = randomize(fdo, so = TRUE)
      numB1 = nrow(fdo$cube) + nrow(fdo$centerCube)
      numB2 = fdo$nrow() - numB1
      fdo$block = data.frame(Block = c(rep(1, numB1),
                                        rep(2, numB2)))
      fdo$blockGen = data.frame(B1 = rep(NA, fdo$nrow()))
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
    Block = .blockCol
    BlockGenCol = .blockGenCol
    fdo$block = Block
    fdo$blockGen = BlockGenCol
  }
  numCC = nrow(fdo$centerCube)
  if (numCC > 0) {
    ccFrame = as.data.frame(matrix(0, nrow = numCC, ncol = ncol(fdo$cube)))
    names(ccFrame) = names(fdo)
    fdo$centerCube = ccFrame
  }
  fdo = randomize(fdo, random.seed = random.seed)
  return(fdo)
}


###Necesito funcion "fracDesign"###################

fracDesign <- function (k = 3, p = 0, gen = NULL, replicates = 1, blocks = 1,
  centerCube = 0, random.seed = 1234)
{
  DB = FALSE
  STDfdo = FALSE
  if (p < 0 || p > 7)
    stop("p needs to be an integer between 0 and 7!")
  if (abs(p - round(p)) > .Machine$double.eps^0.5) {
    warning(paste("p needs to be an integer but is real! p was rounded to",
      round(p)))
    p = round(p)
  }
  if (p != 0) {
    gen = NULL
    for (i in 1:length(.fdoOrth)) {
      if (k == .fdoOrth[[i]]$k && p == .fdoOrth[[i]]$p) {
        STDfdo = TRUE
        return(fracDesign(k = .fdoOrth[[i]]$k, gen = .fdoOrth[[i]]$gen,
          replicates = replicates, blocks = blocks,
          centerCube = centerCube, random.seed = random.seed))
      }
    }
    if (STDfdo == FALSE)
      stop("No standard Design for the choosen combination of k and p (see: fracChoose())!")
  }
  if (!is.numeric(random.seed))
    stop("random.seed needs to be numeric")
  if (!is.numeric(blocks))
    stop("blocks needs to be numeric!")
  if (!is.numeric(replicates)){
    stop("replicates needs to be numeric!")
  } else {
    if (replicates < 0){
      stop("replicates need to >= 0")
    }
  }
  N <- 2^k
  X <- matrix(NA, nrow = N, ncol = k)
  for (j in 1:k) X[, j] <- rep(sort(rep(c(-1, 1), N/2^j)),
    2^(j - 1))
  X <- X[, ncol(X):1]
  if (is.null(gen)) {
    X = as.data.frame(X)
    names(X) = .NAMES[1:k]
  }
  origX = X
  if (replicates > 1) {
    for (i in 2:replicates) {
      X = rbind(X, origX)
    }
  }
  frameOut = data.frame(X)
  if (DB)
    print("juhu")
  if (!is.null(gen)) {
    listGen = vector("list", length(gen))
    .numFactors = numeric(0)
    charFactors = character(0)
    if (DB) {
      cat(paste("gen is: ", gen, "\n"))
      cat(paste("length of gen is: ", length(gen), "\n"))
      print(listGen)
    }
    temp = character(0)
    for (i in seq(along = gen)) {
      if (DB)
        cat("gen[", i, "] = ", gen[i], "\n")
      if (!is.character(gen[i]))
        stop("Defining Relations should contain characters only!")
      chars = strsplit(gen[i], split = character(0))[[1]]
      if (DB) {
        cat("\nchars: ")
        print(chars)
        cat("\n")
      }
      checkDupl = character(0)
      for (j in 1:length(chars)) {
        if (chars[j] %in% toupper(c(.NAMES[1:26], letters[1:26]))) {
          if (chars[j] %in% checkDupl)
            stop("Defining relations contain one or more duplicates!")
          checkDupl = c(checkDupl, chars[j])
          temp = c(temp, chars[j])
        }
      }
      if (DB) {
        cat("\ntemp: ")
        print(temp)
        cat("\n")
      }
    }
    temp = sort(unique(temp))
    numCharVec = 1:length(temp)
    names(numCharVec) = temp
    if (DB) {
      cat("Zuordnung Buchstabe und Spalte:\n")
      print(numCharVec)
      cat("\n")
    }
    for (i in seq(along = gen)) {
      if (DB)
        cat("gen[", i, "] = ", gen[i], "\n")
      if (!is.character(gen[i]))
        stop("Defining Relations should contain characters only!")
      chars = strsplit(gen[i], split = character(0))[[1]]
      numVec = numeric(0)
      charVec = character(0)
      allowedChars = c(.NAMES[1:26], letters[1:26], "=")
      for (j in 1:length(chars)) {
        if (chars[j] %in% allowedChars) {
          if ((chars[j] == "=") & (length(numVec) !=
            1))
            stop("check position of \"=\" in generators!")
          if (chars[j] != "=") {
            charVec = c(charVec, toupper(chars[j]))
            numVec = c(numVec, numCharVec[names(numCharVec) ==
              toupper(chars[j])])
          }
        }
      }
      if (DB) {
        cat("charVec for i = ", i, ": ", charVec, "\n")
        cat("numVec for i = ", i, ": ", numVec, "\n")
      }
      listGen[[i]] = numVec
      .numFactors = c(.numFactors, numVec)
      charFactors = c(charFactors, charVec)
    }
    if (DB)
      print("juhu")
    names(.numFactors) = charFactors
    if (length(unique(.numFactors)) > k)
      stop("number of distinct Factors in generators greater than k!")
    if (DB) {
      print(listGen)
      print(.numFactors)
      print(charFactors)
    }
    for (i in seq(along = listGen)) {
      ind <- trunc(listGen[[i]])
      if (any(abs(ind) > k))
        stop(paste("generator:", paste(ind[1], "=",
          paste(ind[-1], collapse = "*")), "includes undefined columns"))
      x <- rep(sign(ind[1]), N)
      for (j in ind[-1]) x <- x * X[, j]
      X[, abs(ind[1])] <- x
    }
    X <- unique(X)
    origX = X
    if (replicates > 1) {
      for (i in 2:replicates) {
        X = rbind(X, origX)
      }
    }
    frameOut = as.data.frame(X)
    names(frameOut) = names(numCharVec)
    if (k > length(temp)) {
      charsLeft = (.NAMES[1:26])[-match(charFactors, .NAMES[1:26])]
      naIndex = (1:k)[is.na(names(frameOut))]
      names(frameOut)[naIndex] = charsLeft[1:length(naIndex)]
    }
  }
  DesignOut <- facDesign.c$new()
  DesignOut$generator <-  gen
  DesignOut$cube <-  frameOut
  listFac <-  vector("list", ncol(frameOut))
  for (i in seq(along = listFac)){
    listFac[[i]] = doeFactor$new()
    listFac[[i]]$name = names(frameOut)[i]
  }


  DesignOut$factors = listFac
  if (DB)
    print(frameOut)
  if (DB)
    print("yes")
  if (DB)
    print("aha")
  if (centerCube >= 1) {
    temp = data.frame(matrix(rep(0, centerCube * k), ncol = k,
                             nrow = centerCube))
    names(temp) = names(frameOut)
    DesignOut$centerCube = temp
  }
  numRows = nrow(DesignOut$cube) + nrow(DesignOut$centerCube) + nrow(DesignOut$star) +
    nrow(DesignOut$centerStar)
  if (DB) {
    print(numRows)
    print("response")
  }
  DesignOut$response = data.frame(y = rep(NA, numRows))
  if (DB)
    print("response")
  standardOrder = data.frame(matrix(data = 1:numRows, nrow = numRows,
    ncol = 1))
  names(standardOrder) = "StandOrder"
  DesignOut$standardOrder <-  standardOrder
  if (DB)
    print("1")
  set.seed(random.seed)
  runOrder = as.data.frame(standardOrder[sample(1:numRows),])
  if (DB)
    print("2")
  names(runOrder) = "RunOrder"
  DesignOut$runOrder <- runOrder
  if (DB)
    print("3")
  temp = try(blocking(DesignOut, blocks = blocks)) #################Poner random.seed para variar con semilla
  if (inherits(temp, "try-error"))
    stop("Blocking not possible!")
  return(blocking(DesignOut, blocks = blocks))  #################Poner random.seed para variar con semilla
}

############## funcion facDesign########################
facDesign <- function (k = 3, p = 0, replicates = 1, blocks = 1, centerCube = 0, random.seed = 1234)
{
  frameOut = fracDesign(k = k, p = p, gen = NULL, replicates = replicates,
                        blocks = blocks, centerCube = centerCube, random.seed = random.seed)
  return(frameOut)
}


###############USO DE facDesign#######################################################
dfac <- facDesign(k = 3, centerCube = 4)
#dfac$names()
dfac$names(c('Factor 1', 'Factor 2', 'Factor 3'))
#dfac$names()
dfac$lows(c(80,120,1))
#dfac$lows()
dfac$highs(c(120,140,2))
#dfac$highs()
dfac$summary()

###############################################################################

####Necesito funcion .norm2d#####################################################
.norm2d <- function(x1, x2, mu1 = 160, mu2 = 165, rho = 0.7, sigma1 = 45, sigma2 = 22.5) {
  z = 1/(2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) * exp(-1/(2 * (1 - rho^2)) * (((x1 - mu1)/sigma1)^2 - 2 * rho * (x1 - mu1)/sigma1 * (x2 - mu2)/sigma2 +
                                                                                     ((x2 - mu2)/sigma2)^2))
  return(z)
}

####funcion simProc#####################################
simProc <- function(x1, x2, x3, noise = TRUE) {
  max_z = 0.0002200907
  min_z = 8.358082e-10
  yield = .norm2d(x1 = x1, x2 = x2)
  yield = yield - min_z
  yield = (yield/max_z) * 0.9
  if (noise)
    yield = yield + rnorm(length(yield), mean = 0, sd = 0.007)
  return(yield)
}



#####USO simProc######################################
#Primeros valores
rend <- simProc(x1=120,x2=140,x3=2)
#valores completos
rend <- c(simProc(120,140,1),simProc(80,140,1),simProc(120,140,2),simProc(120,120,1),simProc(90,130,1.5),simProc(90,130,1.5),simProc(80,120,2),simProc(90,130,1.5),simProc(90,130,1.5),simProc(120,120,2),simProc(80,140,2),simProc(80,120,1))

#Asignar rendimiento al diseño factorial
dfac$.response(rend)
dfac$.response()
######effectPlot###############################################################
dfac$effectPlot(classic = TRUE)


#####Necesito .letterPos .testFun###########################################

.letterPos <- function(LETTER) {
  if (!(nchar(LETTER) == 1))
    stop("factor names should be single characters only")
  return((1:26)[LETTERS[1:26] == LETTER])
}
.testFun <- function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = TRUE, trace.label = deparse(substitute(trace.factor)),
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

###### Necesito InteractionPlot#########################################################
interactionPlot <- function(fdo, y = NULL, response = NULL, fun = mean, main, col = 1:2, ...) { ###
  DB = FALSE
  mainmiss = FALSE
  if (missing(main))
    mainmiss = TRUE
  if (missing(fdo) || class(fdo)[1]!="facDesign")                                ###
    stop("fdo needs to be an object of class facDesign")                    ###
  parList = list(...)
  old.par <- par(no.readonly = TRUE)
  fdoName = deparse(substitute(fdo))                                          ###
  if(is.null(response)==FALSE)                                                ###
  {                                                                           ###
    temp=fdo$.response()[response]                                               ###
    fdo$.response(temp)                                                         ###
  }                                                                           ###
  diagNames = character(0)
  x = fdo$cube
  runIndex = order(fdo$runOrder[,1])
  x = x[runIndex[1:nrow(x)], ]
  y = fdo$.response()[1:nrow(x), ]
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
  for (r in 1:ncol(fdo$.response())) {
    if (r > 1)
      dev.new()
    y = fdo$.response()[1:nrow(x), r]
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
      main = paste("Interaction plot for", names(fdo$.response())[r], "in", fdoName)         ###
      title(main, outer = T, ...)
    }
    else title(main[r], outer = T, ...)
  }
  #    }
  par(old.par)
  invisible()
}



####Uso interactionPlot#########################################################
interactionPlot(dfac)

######lm##########################################################################
m1 <- dfac$lm(rend ~ A*B*C)
summary(m1)


#####Necesito .splotDev##########################################
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
##### funcion paretoPlot#################################
paretoPlot <- function(fdo, threeWay = FALSE, abs = TRUE, decreasing = TRUE, na.last = NA, alpha = 0.05, response = NULL, xlim, ylim, xlab, ylab, main, single = TRUE, ...) {  ###
  DB = FALSE
  if(single==FALSE)                                                           ###
    par(mfrow=.splitDev(length(fdo$.response()))[[2]])                           ###
  if(is.null(response)==FALSE)                                                ###
  {                                                                           ###
    temp=fdo$.response()[response]                                               ###
    fdo$.response(temp)                                                         ###
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
  for (j in 1:ncol(fdo$.response())) {
    par(mar = c(5.1, 4.1, 4.1, 4.1))
    if (j > 1 && single==TRUE) {
      dev.new()
      par(mar = c(5.1, 4.1, 4.1, 4.1))
    }
    if (!any(is.na(fdo$.response()[, j]))) {
      if (missing(ylab))
        ylabel = names(fdo$.response())[j]                                ###
      else                                                                ###
        ylabel = ylab                                                   ###
      form = paste("fdo$.response()[,", j, "]~")
      for (i in 1:ncol(fdo$cube)) {
        form = paste(form, names(fdo$cube)[i], sep = "")
        if (i < ncol(fdo$cube))
          form = paste(form, "*", sep = "")
      }
      if (DB == TRUE)
        print(form)
      lm.1 = lm(as.formula(form), data = fdo$as.data.frame())
      coefs = coef(lm.1)[-pmatch("(Intercept)", names(coef(lm.1)))]
      df.resid = df.residual(lm.1)
      num.c = nrow(fdo$centerCube)
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
      legend(location, legend = fdo$names(), pch = paste(names(fdo$cube), sep = ""), bg = "white", inset = 0.02)
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

########uso paretoPlot#################################
paretoPlot(dfac)



####necesito .lfkp####################
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
#####funcion normalPlot#####################
normalPlot <- function(fdo, threeWay = FALSE, na.last = NA, alpha = 0.05, response = NULL, sig.col = c("red1", "red2", "red3"), sig.pch = c(1,2,3), main, ylim, xlim, xlab, ylab, pch,  ###
                      col, border = "red", ...) {
  DB = FALSE
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  fdoName = deparse(substitute(fdo))                                          ###
  if(is.null(response)==FALSE)                                                ###
  {                                                                           ###
    temp=fdo$.response()[response]                                               ###
    fdo$.response(temp)                                                         ###
  }                                                                           ###
  parList = list(...)
  params = list()
  if (length(sig.col) < 3)
    sig.col = as.vector(matrix(sig.col, nrow = 1, ncol = 3))
  XLIM=FALSE;YLIM=FALSE                                                       ###
  if (!(class(fdo)[1] == "facDesign"))
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
  for (j in 1:ncol(fdo$.response())) {
    parList = list(...)                                                     ###
    params = list()                                                         ###
    leg.col = vector()                                                      ###
    p.col = vector()                                                        ###
    p.pch = vector()                                                        ###
    leg.txt = vector()                                                      ###
    main = paste("Normal plot for", names(fdo$.response())[j], "in", fdoName) ###
    if (j > 1)
      dev.new()
    form = paste("fdo$.response()[,", j, "]~")
    for (i in 1:ncol(fdo$cube)) {
      form = paste(form, names(fdo$cube)[i], sep = "")
      if (i < ncol(fdo$cube))
        form = paste(form, "*", sep = "")
    }
    if (DB == TRUE)
      print(paste("form:", form))
    lm.1 = lm(as.formula(form), data = fdo$as.data.frame())
    lm.1s = summary(lm.1)
    effect = coef(lm.1s)[row.names(coef(lm.1s)) != "(Intercept)", "t value"]
    print(effect)
    if (all(is.na(effect)))
      effect = 2 * coef(lm.1)[-pmatch("(Intercept)", names(coef(lm.1)))]      ###
    #            stop("effects could not be calculated")                            ###
    sig = summary(lm.1)$coefficients[, "Pr(>|t|)"][-pmatch("(Intercept)", names(coef(lm.1)))]
    df.resid = df.residual(lm.1)
    nc = nrow(fdo$centerCube)
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

##########uso normalPlot#################
normalPlot(dfac)

###Necesito clase desirability#################
###Necesito .desirFun###############
####funcion wirePlot###################
wirePlot <- function(x, y, z, data = NULL, xlim, ylim, zlim, main, xlab, ylab, border, sub, zlab, form = "fit", phi, theta, ticktype, col = 1, steps,
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
  if (class(data)[1] != "facDesign") {
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
    ylab = paste(y.c, ": ", fdo$names()[[y.c]])
  if (missing(xlab))
    xlab = paste(x.c, ": ", fdo$names()[[x.c]])
  if (missing(zlab))
    zlab = paste(x.c, ": ", names(fdo$.response()))
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
    xlim = c(min(fdo$get(, x.c)), max(fdo$get(, x.c)))
  if (missing(ylim))
    ylim = c(min(fdo$get(, y.c)), max(fdo$get(, y.c)))
  allVars = c(names(fdo$names()), names(fdo$.response))
  isct = intersect(c(x.c, y.c, z.c), c(fdo$names(), names(fdo$.response)))
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

###uso wirePlot###############################
wirePlot(A,B,rend,data=dfac)
names(dfac$.response())
