library(R6)

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
                                                     names(frameOut) = as.character(self$names())
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
m1 <- lm(formula(rend ~ A*B*C), data=dfac)
summary(m1)


library(qualityTools)
function (formula=rend ~ A*B*C, data=dfac, subset, weights, na.action, method = "qr",
          model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
          contrasts = NULL, offset, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame")
    return(mf)
  else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'",
                     method), domain = NA)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- model.offset(mf)
  mlm <- is.matrix(y)
  ny <- if (mlm)
    nrow(y)
  else length(y)
  if (!is.null(offset)) {
    if (!mlm)
      offset <- as.vector(offset)
    if (NROW(offset) != ny)
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    NROW(offset), ny), domain = NA)
  }
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (mlm) matrix(NA_real_, 0,
                                             ncol(y)) else numeric(), residuals = y, fitted.values = 0 *
                y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w !=
                                                                                0) else ny)
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if (is.null(w))
      lm.fit(x, y, offset = offset, singular.ok = singular.ok,
             ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,
                 ...)
  }
  class(z) <- c(if (mlm) "mlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr)
    z$qr <- NULL
  z
}
