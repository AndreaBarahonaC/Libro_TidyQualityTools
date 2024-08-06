library(R6)
library(ggplot2)
library(patchwork)
library(scales)
library(plotly)
library(gridExtra)
library(RColorBrewer)

### Funcion as.data.frame.facDesign####
as.data.frame.facDesign <- function(self, row.names = NULL, optional = FALSE, ...) {
  if (nrow(self$cube)>0) {
    frameOut = self$cube
    names(frameOut) = self$names()
  }
  else return(NULL)
  if (nrow(self$centerCube)>0){
    faux <- self$centerCube
    names(faux) <- self$names()
    frameOut = rbind(frameOut, faux)
  }
  if (nrow(self$star)>0)
    frameOut = rbind(frameOut, self$star)
  if (nrow(self$centerStar)>0)
    frameOut = rbind(frameOut, self$centerStar)
  aux <- list()
  for (i in 1:length(self$names())) {
    aux[[self$names()[i]]] <-.NAMES[i]
  }
  if (!is.null(self$factors) && length(self$factors) == dim(frameOut)[2]) {
    names(frameOut) = as.character(aux)
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
}

############DISEÑOS FACTORIALES 2^k#############
### necesito .helpAliasTable##################################################
### Funcion .replace2s####
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

### necesito aliasTable#######################################################
aliasTable <- function (fdo, degree, print = TRUE)
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
  if (print)
    print(round(alias.matrix, 2))
  invisible(alias.matrix)
}

### necesito .fdoOrth y .NAMES#######################################
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



### necesito .m.interaction.plot###########################
.m.interaction.plot <- function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = TRUE, trace.label = deparse(substitute(trace.factor)),
                                fixed = FALSE, xlab = deparse(substitute(x.factor)), ylab = ylabel, ytitle = TRUE, ylim = range(cells, na.rm = TRUE), lty = nc:1, col = 1, pch = c(1L:9, 0, letters), xpd = NULL,
                                leg.bg = par("bg"), leg.bty = "n", xtick = FALSE, xaxt = par("xaxt"), axes.x = TRUE, axes.y = TRUE, main, ...) {
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
  if (missing(main)) {
    main = paste("Effect Plot")
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

  df <- data.frame(x = xvals, y = c(cells))

  # PLOT
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(na.rm = TRUE) +
    ylim(ylim) + labs(x = xlab, y = ylab, title = main) + theme_bw() +

    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5))

  if (axes.x){
    p <- p +
      theme(axis.ticks.x = element_line(),
            axis.text.x = element_text()) +
      scale_x_continuous(breaks=c(-1, 1))
  }
  if (axes.y){
    p <- p +
      theme(axis.ticks.y = element_line(),
            axis.text.y = element_text())
  }
  if (ytitle == FALSE){
    p <-  p + theme(axis.title.y = element_blank())
  }
  if(main == ""){
    p <- p + labs(title = NULL)
  }
  if (deparse(substitute(mean)) == "mean"){
    p <- p + geom_hline(yintercept = median(df$y, na.rm = TRUE), linetype = "dashed", col = "#324B7A")
  }

  invisible(list(xVals = df$x, yVals = df$y, plot = p))
}

### necesito clase facDesign.c######################################################
facDesign.c <- R6Class("facDesign", public = list(name = NULL,
                                                 factors = NULL,
                                                 cube = data.frame(),
                                                 star = data.frame(),
                                                 centerCube = data.frame(),
                                                 centerStar = data.frame(),
                                                 generator = NULL,
                                                 response = data.frame(),
                                                 block = data.frame(),
                                                 blockGen = data.frame(),
                                                 runOrder = data.frame(),
                                                 standardOrder = data.frame(),
                                                 desireVal = NULL,
                                                 desirability = list(),
                                                 fits = NULL,

                                                 nrow = function(){
                                                   nrow(self$as.data.frame())
                                                 },

                                                 ncol = function(){
                                                   ncol(self$as.data.frame())
                                                 },

                                                 print = function(){
                                                   runIndex = order(self$runOrder[,1])
                                                   print(format(self$as.data.frame(), digits = 4))
                                                   invisible(self$as.data.frame())
                                                 },

                                                 .clear = function(){
                                                   self$standardOrder = data.frame()
                                                   self$runOrder = data.frame()
                                                   self$cube = data.frame()
                                                   self$centerStar = data.frame()
                                                   self$centerCube = data.frame()
                                                   self$star = data.frame()
                                                   self$block = data.frame()
                                                   self$blockGen = data.frame()
                                                   self$response = data.frame()
                                                   invisible(self)
                                                 },

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
                                                   if (nrow(self$cube)>0) {
                                                     frameOut = self$cube
                                                     names(frameOut) = self$names()
                                                   }
                                                   else return(NULL)
                                                   if (nrow(self$centerCube)>0){
                                                     faux <- self$centerCube
                                                     names(faux) <- self$names()
                                                     frameOut = rbind(frameOut, faux)
                                                   }
                                                   if (nrow(self$star)>0){
                                                     faux <- self$star
                                                     names(faux) <- self$names()
                                                     frameOut = rbind(frameOut, faux)
                                                   }
                                                   if (nrow(self$centerStar)>0){
                                                     faux <- self$centerStar
                                                     names(faux) <- self$names()
                                                     frameOut = rbind(frameOut, faux)
                                                   }
                                                   aux <- list()
                                                   for (i in 1:length(self$names())) {
                                                     aux[[self$names()[i]]] <-.NAMES[i]
                                                   }
                                                   if (!is.null(self$factors) && length(self$factors) == dim(frameOut)[2]) {
                                                     names(frameOut) = as.character(aux)
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
                                                   A = aliasTable(self, print = FALSE)
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
                                                   temp = aliasTable(self, print = FALSE)
                                                   if (ncol(temp) > 0) {
                                                     cat("\n---------\n\n")
                                                     self$identity()
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

                                                   if(is.null(response)==FALSE)
                                                   {                                                                           ###
                                                     temp=self$.response()[response]                                            ###
                                                     self$.response(temp)                                                      ###
                                                   }
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
                                                   if (missing(main)) {
                                                     main = paste("Effect Plot for", names(Y)[1])
                                                   }

                                                   list_plot <- list()
                                                   #############################################################################
                                                   for (j in 1:ncol(Y)){
                                                     counter = 0
                                                     cells = numeric(0)
                                                     for (i in 1:length(factors)){
                                                       cells = c(cells, as.vector(tapply(Y[, j], list(X[, factors[i]], rep(0, nrow(X))), fun)))
                                                       if (points)
                                                         cells = range(Y)
                                                     }
                                                     if (nextResponse & !single) {
                                                       dev.new()
                                                       par(mfrow = c(numRow, numCol))
                                                     }

                                                     # hacemos la primera afuera para ajustar los ejes
                                                     # 1. ------------------------------------
                                                     i <- 1
                                                     if ((counter != 0 & counter%%(numCol * numRow) == 0) & !single) {
                                                       dev.new()
                                                       par(mfrow = c(numRow, numCol))
                                                     }
                                                     if (missing(main)) {
                                                       main = paste("Effect Plot for", names(Y)[j])
                                                       mainmiss = TRUE
                                                     }
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
                                                     if (classic) {
                                                       p <- .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun,
                                                                                ylab = ylab, main = " ", ...)
                                                       list_plot[[paste0("p",j,i)]] <- p$plot
                                                     }
                                                     else {
                                                       p <- .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun,
                                                                                ylab = ylab, main = main, ...)
                                                       list_plot[[paste0("p",j,i)]] <- p$plot
                                                     }
                                                     counter = counter + 1
                                                     # 2. ------------------------------------
                                                     if(length(factors) >= 2){
                                                       for (i in 2:length(factors)) {
                                                         if ((counter != 0 & counter%%(numCol * numRow) == 0) & !single) {
                                                           dev.new()
                                                           par(mfrow = c(numRow, numCol))
                                                         }
                                                         if (missing(main)) {
                                                           main = paste("Effect Plot for", names(Y)[j])
                                                           mainmiss = TRUE
                                                         }
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
                                                         if (classic) {
                                                           aux <- .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun,
                                                                                      ylab = ylab, axes.y = FALSE, ytitle = FALSE , main = " ", ...)
                                                           list_plot[[paste0("p",j,i)]] <- aux$plot
                                                         }
                                                         else {
                                                           aux <- .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun,
                                                                                      ylab = ylab, main = main, ytitle = TRUE, ...)
                                                           list_plot[[paste0("p",j,i)]] <- aux$plot
                                                         }
                                                         counter = counter + 1
                                                       }
                                                     }
                                                     nextResponse = TRUE
                                                   }
                                                   # Obtener los nombres de todas las graficas que se crearon
                                                   grap <- c()
                                                   for(j in 1:ncol(Y)){
                                                     for(i in 1:length(factors)){
                                                       x <- paste0("p",j,i)
                                                       if(!x %in% grap){
                                                         grap <- c(grap, x)
                                                       }
                                                     }
                                                   }

                                                   p <- list_plot$p11
                                                   for(i in 2:length(grap)){
                                                     aux <- grap[i]
                                                     p <- p + list_plot[[aux]]
                                                   }

                                                   if(classic){
                                                     p <- p + plot_layout(ncol = numCol, nrow = numRow) +
                                                       plot_annotation(title = main,
                                                                       theme = theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")))
                                                     print(p)
                                                   }
                                                   else{print(p)}
                                                   par(mfcol=c(1,1))
                                                 },

                                                 lm = function(formula){
                                                   invisible(lm(formula, data = self$as.data.frame()))

                                                 },

                                                 desires = function(value){
                                                   if (missing(value)) {
                                                     return(self$desirability)
                                                   }
                                                   else{
                                                     if (!any(value$response == names(self$.response())))
                                                       stop(paste(value$response, "is not a response!"))
                                                     listPos = length(self$desirability) + 1
                                                     yName = value$response
                                                     isIn = (yName == names(self$desirability))
                                                     if (any(isIn))
                                                       listPos = (1:length(names(self$desirability)))[isIn]
                                                     self$desirability[[listPos]] = value
                                                     names(self$desirability)[listPos] = yName
                                                     invisible(self)
                                                   }

                                                 },

                                                 set.fits = function(value){
                                                   if (!identical(class(value), "lm"))
                                                     stop(paste(deparse(substitute(value)), "needs to an object of class lm"))
                                                   if (!any(names(value$model)[1] == names(self$.response())))
                                                     stop(paste("fitted response", names(value$model)[1], "could not be found in", deparse(substitute(x))))
                                                   listPos = length(self$fits) + 1
                                                   yName = names(value$model)[1]
                                                   isIn = (yName == names(self$fits))
                                                   if (any(isIn))
                                                     listPos = (1:length(names(self$fits)))[isIn]
                                                   self$fits[[listPos]] = value
                                                   names(self$fits)[listPos] = yName
                                                   invisible(self)

                                                 },

                                                 types = function(value){
                                                   if (missing(value)) {
                                                     v <- list()
                                                     for (i in 1:length(self$factors)) {
                                                       v[[self$names()[i]]] <- self$factors[[i]]$.type()
                                                     }
                                                     return(v)

                                                   }
                                                   else{
                                                     for (i in 1:length(self$factors)) {
                                                       if (!identical(value[i], "numeric") & !identical(value[i], "factor"))
                                                         stop(paste(value[i], "\ttype of factor needs to be 'numeric' or 'factor'"))
                                                       self$factors[[i]]$.type(as.character(value[i]))
                                                     }
                                                     invisible(self)
                                                   }
                                                 },

                                                 unit = function(value){
                                                   if (missing(value)) {
                                                     v <- list()
                                                     for (i in 1:length(self$factors)) {
                                                       v[[self$names()[i]]] <- self$factors[[i]]$.unit()
                                                     }
                                                     return(v)
                                                   }
                                                   else{
                                                     for (i in 1:length(self$factors)) if (length(value) > 1)
                                                       self$factors[[i]]$.unit(as.character(value[i]))
                                                     else self$factors[[i]]$.unit(as.character(value[1]))
                                                     invisible(self)
                                                   }

                                                 },

                                                 .star = function(value){
                                                   if (missing(value)) {
                                                     return(self$star)
                                                   }
                                                   else{
                                                     DB = FALSE
                                                     if (!is.data.frame(value))
                                                       stop("data.frame must be provided!")
                                                     if (.numFac(self) != ncol(value))
                                                       stop("number of columns not matching!")
                                                     if (nrow(value) == 0) {
                                                       return("TODO: remove star und Rest anpassen")
                                                     }
                                                     oldResponse = self$.response()
                                                     newDf = value
                                                     oldDf = self$star
                                                     numNewRow = nrow(newDf) - nrow(oldDf)
                                                     oldOrd = self$standardOrder
                                                     oldRunOrd = self$runOrder
                                                     len = nrow(oldOrd)
                                                     lenFirst = nrow(self$cube) + nrow(self$centerCube)
                                                     self$standardOrder = data.frame(StandOrd = 1:(len + numNewRow))
                                                     newRunOrd = data.frame()
                                                     if (numNewRow > 0) {
                                                       newNums = data.frame(newNums = seq(max(oldRunOrd) + 1, max(oldRunOrd) + numNewRow, by = 1))
                                                       if (DB)
                                                         print(newNums)
                                                       names(newNums) = names(oldRunOrd)
                                                       newRunOrd = data.frame(oldRunOrd[1:lenFirst, ])
                                                       if (DB)
                                                         print(newRunOrd)
                                                       names(newRunOrd) = names(oldRunOrd)
                                                       restFrame = data.frame(oldRunOrd[-c(1:lenFirst), ])
                                                       names(restFrame) = names(oldRunOrd)
                                                       newRunOrd = rbind(newRunOrd, newNums, restFrame)
                                                       if (DB)
                                                         print(newRunOrd)
                                                     }
                                                     else {
                                                       newRunOrd = data.frame(oldRunOrd[1:(lenFirst + nrow(newDf) + nrow(self$centerStar)), ])
                                                       names(newRunOrd) = names(oldRunOrd)
                                                     }
                                                     self$runOrder = newRunOrd
                                                     naFrame = as.data.frame(matrix(rep(NA, times = ncol(oldResponse) * nrow(newDf)), ncol = ncol(oldResponse)))
                                                     names(naFrame) = names(oldResponse)
                                                     newResponse = data.frame(oldResponse[1:lenFirst, ])
                                                     names(newResponse) = names(oldResponse)
                                                     restFrame = data.frame(oldResponse[-c(1:(lenFirst + nrow(oldDf))), ])
                                                     names(restFrame) = names(oldResponse)
                                                     newResponse = rbind(newResponse, naFrame, restFrame)
                                                     self$.response(newResponse)
                                                     if (DB) {
                                                       print(newResponse)
                                                       print("hinter response")
                                                     }
                                                     oldBlockGen = self$blockGen
                                                     if (ncol(oldBlockGen) > 0) {
                                                       if (DB)
                                                         print("TODO: BlockGen anpassen!")
                                                       newBlockGen = data.frame(oldBlockGen[1:lenFirst, ])
                                                       names(newBlockGen) = names(self$blockGen)
                                                       naFrameGen = as.data.frame(matrix(rep(NA, times = ncol(self$blockGen) * nrow(newDf)), ncol = ncol(self$blockGen)))
                                                       names(naFrameGen) = names(oldBlockGen)
                                                       restBlockGen = data.frame(oldBlockGen[-c(1:(lenFirst + nrow(oldDf))), ])
                                                       names(restBlockGen) = names(oldBlockGen)
                                                       newBlockGen = rbind(newBlockGen, naFrameGen, restBlockGen)
                                                       if (DB)
                                                         print(newBlockGen)
                                                       self$.blockGen(newBlockGen)
                                                     }
                                                     oldBlock = self$block
                                                     newBlock = data.frame(oldBlock[1:lenFirst, ])
                                                     names(newBlock) = names(oldBlock)
                                                     naFrame = as.data.frame(matrix(rep(max(newBlock) + 1, times = ncol(oldBlock) * nrow(newDf)),
                                                                                    ncol = ncol(oldBlock)))
                                                     names(naFrame) = names(oldBlock)
                                                     restBlock = data.frame(oldBlock[-c(1:(lenFirst + nrow(oldDf))), ])
                                                     names(restBlock) = names(oldBlock)
                                                     newBlock = rbind(newBlock, naFrame, restBlock)
                                                     self$.block(newBlock)
                                                     self$star <- newDf
                                                     invisible(self)
                                                   }
                                                 },

                                                 .blockGen = function(value){
                                                   if (missing(value)) {
                                                     return(self$blockGen)
                                                   }
                                                   else{
                                                     if (!is.vector(value) && !is.data.frame(value))
                                                       stop("vector or data.frame expected!")
                                                     if (is.vector(value) && (is.numeric(value) || is.na(value))) {
                                                       if (self$nrow() != length(value))
                                                         stop(paste("Number of rows for Design does not equal length of vector ", self$nrow(),
                                                                    " != ", length(value), " "))
                                                       self$blockGen <- as.data.frame(value)
                                                       names(self$blockGen) = deparse(substitute(value))

                                                     }
                                                     if (is.data.frame(value)) {
                                                       self$blockGen <- value

                                                     }
                                                     invisible(self)

                                                   }

                                                 },

                                                 .block = function(value){
                                                   if (missing(value)) {
                                                     return(self$block)
                                                   }
                                                   else{
                                                     if (!is.vector(value) && !is.data.frame(value))
                                                       stop("vector or data.frame expected!")
                                                     if (is.vector(value) && (is.numeric(value) || is.na(value))) {
                                                       if (self$nrow() != length(value))
                                                         stop(paste("Number of rows for Design does not equal length of vector ", nrow(object),
                                                                    " != ", length(value), " "))
                                                       self$block <- as.data.frame(value)
                                                       names(self$block) = deparse(substitute(value))

                                                     }
                                                     if (is.data.frame(value)) {
                                                       self$block <- value

                                                     }
                                                     invisible(self)

                                                   }
                                                 },

                                                 .centerCube = function(value){
                                                   if (missing(value)) {
                                                     return(self$centerCube)
                                                   }
                                                   else{
                                                     DB = FALSE
                                                     if (!is.data.frame(value))
                                                       stop("data.frame must be provided!")
                                                     if (.numFac(self) != ncol(value))
                                                       stop("number of columns not matching!")
                                                     if (nrow(value) == 0) {
                                                       return("TODO: remove CenterCube und Rest anpassen")
                                                     }
                                                     newDf = value
                                                     lenCube = nrow(self$cube)
                                                     oldDf = self$centerCube
                                                     oldRunOrd = self$runOrder
                                                     oldResponse = self$.response()
                                                     blockValues = unique(self$block[1:nrow(self$cube), ])
                                                     numBlocks = length(blockValues)
                                                     if (numBlocks > 1)
                                                       for (i in 1:(numBlocks - 1)) {
                                                         newDf = rbind(newDf, value)
                                                       }
                                                     if (DB)
                                                       print(newDf)
                                                     numNewRow = nrow(newDf) - nrow(oldDf)
                                                     oldOrd = self$standardOrder
                                                     len = nrow(oldOrd)
                                                     self$standardOrder = data.frame(StandOrd = 1:(len + numNewRow))
                                                     newRunOrd = data.frame()
                                                     if (numNewRow > 0) {
                                                       newNums = data.frame(newNums = seq(max(oldRunOrd) + 1, max(oldRunOrd) + numNewRow, by = 1))
                                                       names(newNums) = names(oldRunOrd)
                                                       if (DB) {
                                                         print("----")
                                                         print(newNums)
                                                       }
                                                       newRunOrd = data.frame(oldRunOrd[1:lenCube, ])
                                                       names(newRunOrd) = names(oldRunOrd)
                                                       restRunOrd = data.frame(oldRunOrd[-c(1:lenCube), ])
                                                       names(restRunOrd) = names(oldRunOrd)
                                                       newRunOrd = rbind(newRunOrd, newNums, restRunOrd)
                                                       if (DB) {
                                                         print("----")
                                                         print(oldRunOrd[-c(1:lenCube), ])
                                                         print("----")
                                                         print(newRunOrd)
                                                       }
                                                       self$runOrder = newRunOrd
                                                     }
                                                     else {
                                                       newRunOrd = data.frame(oldRunOrd[1:(lenCube + nrow(newDf)), ])
                                                       names(newRunOrd) = names(oldRunOrd)
                                                       restRunOrd = data.frame(oldRunOrd[-c(1:(lenCube + nrow(oldDf))), ])
                                                       names(restRunOrd) = names(oldRunOrd)
                                                       newRunOrd = rbind(newRunOrd, restRunOrd)
                                                       if (DB) {
                                                         print("----")
                                                         print(newRunOrd)
                                                       }
                                                       self$runOrder = newRunOrd
                                                     }
                                                     naFrame = as.data.frame(matrix(rep(NA, times = ncol(oldResponse) * nrow(newDf)), ncol = ncol(oldResponse)))
                                                     names(naFrame) = names(oldResponse)
                                                     newResponse = data.frame(oldResponse[1:lenCube, ])
                                                     names(newResponse) = names(oldResponse)
                                                     restResponse = data.frame(oldResponse[-c(1:(lenCube + nrow(oldDf))), ])
                                                     names(restResponse) = names(oldResponse)
                                                     newResponse = rbind(newResponse, naFrame, restResponse)
                                                     if (DB) {
                                                       print("newResponse_____")
                                                       print(newResponse)
                                                     }
                                                     self$.response(newResponse)
                                                     oldBlockGen = self$blockGen
                                                     if (ncol(oldBlockGen) > 0) {
                                                       if (DB)
                                                         print("TODO: BlockGen Spalte(n) anpassen")
                                                       newBlockGen = data.frame(oldBlockGen[1:lenCube, ])
                                                       names(newBlockGen) = names(self$blockGen)
                                                       naFrameGen = as.data.frame(matrix(rep(NA, times = ncol(self$blockGen) * nrow(newDf)), ncol = ncol(self$blockGen)))
                                                       names(naFrameGen) = names(oldBlockGen)
                                                       restFrame = data.frame(oldBlockGen[-c(1:(lenCube + nrow(oldDf))), ])
                                                       names(restFrame) = names(self$blockGen)
                                                       newBlockGen = rbind(newBlockGen, naFrameGen, restFrame)
                                                       self$.blockGen(newBlockGen)
                                                       if (DB)
                                                         print(newBlockGen)
                                                     }
                                                     oldBlock = self$block
                                                     newBlock = data.frame(oldBlock[1:lenCube, ])
                                                     names(newBlock) = names(self$block)
                                                     naFrame = as.data.frame(matrix(rep(blockValues, times = nrow(newDf)/numBlocks), ncol = 1))
                                                     restFrame = as.data.frame(oldBlock[-c(1:(lenCube + nrow(oldDf))), ])
                                                     names(restFrame) = names(self$block)
                                                     if (DB)
                                                       print(naFrame)
                                                     names(naFrame) = names(oldBlock)
                                                     newBlock = rbind(newBlock, naFrame, restFrame)
                                                     self$.block(newBlock)
                                                     self$centerCube <- newDf
                                                     invisible(self)
                                                   }

                                                 },

                                                 .centerStar = function(value){
                                                   if (missing(value)) {
                                                     return(self$centerStar)
                                                   }
                                                   else{
                                                     DB = FALSE
                                                     if (!is.data.frame(value))
                                                       stop("data.frame must be provided!")
                                                     if (.numFac(self) != ncol(value))
                                                       stop("number of columns not matching!")
                                                     if (nrow(value) == 0) {
                                                       return("TODO: remove CenterCube und Rest anpassen")
                                                     }
                                                     newDf = value
                                                     oldDf = self$centerStar
                                                     numNewRow = nrow(newDf) - nrow(oldDf)
                                                     oldResponse = self$.response()
                                                     lenRest = nrow(self$cube) + nrow(self$centerCube) + nrow(self$star)
                                                     oldRunOrd = self$runOrder
                                                     oldOrd = self$standardOrder
                                                     len = nrow(oldOrd)
                                                     self$standardOrder = data.frame(StandOrd = 1:(len + numNewRow))
                                                     newRunOrd = data.frame(oldRunOrd[1:lenRest, ])
                                                     names(newRunOrd) = names(oldRunOrd)
                                                     if (numNewRow > 0) {
                                                       newNums = data.frame(newNums = seq(max(oldRunOrd) + 1, max(oldRunOrd) + numNewRow, by = 1))
                                                       names(newNums) = names(oldRunOrd)
                                                       restFrame = data.frame(oldRunOrd[-c(1:lenRest), ])
                                                       names(restFrame) = names(oldRunOrd)
                                                       newRunOrd = rbind(newRunOrd, newNums, restFrame)
                                                       if (DB)
                                                         print(newRunOrd)
                                                       self$runOrder = newRunOrd
                                                     }
                                                     else {
                                                       newRunOrd = data.frame(oldRunOrd[1:(lenRest + nrow(newDf)), ])
                                                       names(newRunOrd) = names(oldRunOrd)
                                                       self$runOrder = newRunOrd
                                                     }
                                                     naFrame = as.data.frame(matrix(rep(NA, times = ncol(oldResponse) * nrow(newDf)), ncol = ncol(oldResponse)))
                                                     names(naFrame) = names(oldResponse)
                                                     newResponse = data.frame(oldResponse[1:lenRest, ])
                                                     names(newResponse) = names(self$.response())
                                                     newResponse = rbind(newResponse, naFrame)
                                                     if (DB)
                                                       print(" vor centerStar response")
                                                     self$.response(newResponse)
                                                     if (DB)
                                                       print("hinter centerStar response")
                                                     oldBlockGen = self$blockGen
                                                     if (ncol(oldBlockGen) > 0) {
                                                       print("TODO: BlockGen Spalte(n) anpassen")
                                                       newBlockGen = data.frame(oldBlockGen[1:lenRest, ])
                                                       names(newBlockGen) = names(self$blockGen)
                                                       naFrameGen = as.data.frame(matrix(rep(NA, times = ncol(self$blockGen) * nrow(newDf)), ncol = ncol(self$block)))
                                                       names(naFrameGen) = names(oldBlockGen)
                                                       restBlockGen = data.frame(oldBlockGen[-c(1:(lenRest + nrow(oldDf))), ])
                                                       names(restBlockGen) = names(oldBlockGen)
                                                       newBlockGen = rbind(newBlockGen, naFrameGen, restBlockGen)
                                                       if (DB)
                                                         print(newBlockGen)
                                                       self$.blockGen(newBlockGen)
                                                     }
                                                     oldBlock = self$block
                                                     newBlock = data.frame(oldBlock[1:lenRest, ])
                                                     names(newBlock) = names(self$block)
                                                     naFrame = as.data.frame(matrix(rep(max(self$block[1:nrow(self$cube), ]) + 1, times = ncol(self$block) *
                                                                                          nrow(newDf)), ncol = ncol(self$block)))
                                                     names(naFrame) = names(oldBlock)
                                                     restBlock = data.frame(oldBlock[-c(1:(lenRest + nrow(oldDf))), ])
                                                     names(restBlock) = names(oldBlock)
                                                     newBlock = rbind(newBlock, naFrame, restBlock)
                                                     self$.block(newBlock)
                                                     self$centerStar <- newDf
                                                     invisible(self)

                                                   }

                                                 },

                                                 .generators = function(value){
                                                   if (missing(value)) {
                                                     return(self$generator)

                                                   }
                                                   else{
                                                     self$generator <- value
                                                     invisible(self)
                                                   }
                                                 }




                                                 )
                      )

### necesito clase doeFactor####################
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
                                                 },

                                                .type = function(value){
                                                  if (missing(value)) {
                                                    return(self$type)
                                                  }
                                                  else{
                                                    self$type <- value
                                                    invisible(self)
                                                  }
                                                },

                                                .unit = function(value){
                                                  if (missing(value)){
                                                    return(self$unit)
                                                  }
                                                  else{
                                                    self$unit <- value
                                                    invisible(self)
                                                  }
                                                },

                                                names = function(value){
                                                  if (missing(value)) {
                                                    return(self$name)
                                                  }
                                                  else {
                                                    self$name <- value
                                                    invisible(self)
                                                  }
                                                },

                                                print = function(){
                                                  cat("Name: ", self$names(), "\n")
                                                  cat("low Setting: ", self$.low(), "\n")
                                                  cat("high setting: ", self$.high(), "\n")
                                                  cat("Unit: ", self$.unit(), "\n")
                                                  cat("type: ", self$.type(), "\n")
                                                  cat("\n")
                                                }

                                                )
                      )




### necesito funcion randomize####
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


### necesito .rsm#################
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


### necesito .numFac#############
.numFac = function(fdo) {
  return(length(fdo$names()))
}


### necesito .confoundings####
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


### necesito .lociv####
.lociv = function(charVec) {
  lenVec = numeric(length = length(charVec))
  for (i in seq(along = charVec)) {
    lenVec[i] = length(strsplit(charVec[i], split = "")[[1]])
  }
  return(lenVec)
}

### necesito .blockInteractions######
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



### necesito .blockGenCol##########
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


### necesito .blockCol##############
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


### necesito funcion blocking######
blocking <- function (fdo, blocks, BoR = FALSE, random.seed, useTable = "rsm",
                      gen)
{
  override = FALSE
  Block = data.frame(Block = rep(1, fdo$nrow()))
  fdo$.block(Block)
  fdo = randomize(fdo, so = TRUE)
  if (missing(random.seed)) {
    runif(1)
    random.seed = .Random.seed[sample(1:626, 1)]
  }
  if (missing(gen))
    gen = NULL
  if (blocks <= 1) {
    Block = data.frame(Block = rep(1, fdo$nrow()))
    fdo$.block(Block)
    fdo = randomize(fdo, random.seed = random.seed)
    return(fdo)
  }
  if (nrow(fdo$star) > 0 | nrow(fdo$centerStar) > 0) {
    if (blocks == 2) {
      override = TRUE
      fdo = randomize(fdo, so = TRUE)
      numB1 = nrow(fdo$cube) + nrow(fdo$centerCube)
      numB2 = fdo$nrow() - numB1
      fdo$.block(data.frame(Block = c(rep(1, numB1),
                                      rep(2, numB2))))
      fdo$.blockGen(data.frame(B1 = rep(NA, fdo$nrow())))
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
    fdo$.block(Block)
    fdo$.blockGen(BlockGenCol)
  }
  numCC = nrow(fdo$centerCube)
  if (numCC > 0) {
    ccFrame = as.data.frame(matrix(0, nrow = numCC, ncol = ncol(fdo$cube)))
    names(ccFrame) = names(fdo)
    fdo$.centerCube(ccFrame)
  }
  fdo = randomize(fdo, random.seed = random.seed)
  return(fdo)
}


### necesito funcion "fracDesign"###################

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

### Funcion facDesign########################
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

### necesito funcion .norm2d#####################################################
.norm2d <- function(x1, x2, mu1 = 160, mu2 = 165, rho = 0.7, sigma1 = 45, sigma2 = 22.5) {
  z = 1/(2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) * exp(-1/(2 * (1 - rho^2)) * (((x1 - mu1)/sigma1)^2 - 2 * rho * (x1 - mu1)/sigma1 * (x2 - mu2)/sigma2 +
                                                                                     ((x2 - mu2)/sigma2)^2))
  return(z)
}

### Funcion simProc#####################################
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



### USO simProc######################################
#Primeros valores
rend <- simProc(x1=120,x2=140,x3=2)
#valores completos
rend <- c(simProc(120,140,1),simProc(80,140,1),simProc(120,140,2),simProc(120,120,1),simProc(90,130,1.5),simProc(90,130,1.5),simProc(80,120,2),simProc(90,130,1.5),simProc(90,130,1.5),simProc(120,120,2),simProc(80,140,2),simProc(80,120,1))

#Asignar rendimiento al diseño factorial
dfac$.response(rend)
dfac$.response()


### effectPlot###############################################################
dfac$effectPlot(classic=TRUE)
dfac$effectPlot()


### necesito .letterPos .testFun###########################################
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


### necesito InteractionPlot#########################################################
interactionPlot <- function(fdo, y = NULL, response = NULL, fun = mean, main, col = 1:2, ...) {
  if (missing(main)) mainmiss = TRUE else mainmiss = FALSE
  if (missing(fdo) || class(fdo)[1] != "facDesign")
    stop("fdo needs to be an object of class facDesign")

  fdoName = deparse(substitute(fdo))
  if (!is.null(response)) {
    temp = fdo$.response()[response]
    fdo$.response(temp)
  }

  x <- fdo$cube
  runIndex <- order(fdo$runOrder[, 1])
  y <- fdo$.response()[1:nrow(x), ]
  numFac <- ncol(x)
  combMat <- combn(names(x), 2)

  plot_list <- list()

  for (r in 1:ncol(fdo$.response())) {
    y = fdo$.response()[1:nrow(x), r]

    for (i in 1:ncol(combMat)) {
      facName1 <- combMat[1, i]
      facName2 <- combMat[2, i]

      df = data.frame(fac1 = x[[facName1]], fac2 = x[[facName2]], response = y)

      p = ggplot(df, aes_string(x = "fac2", y = "response", color = "as.factor(fac1)")) +
        geom_line(stat = "summary", fun = fun,size=1.5) +
        labs(x = " ", y = " ", color = facName1) +
        theme_minimal()

      plot_list[[paste(facName1, facName2)]] = p
    }

    if (mainmiss) {
      main = paste("Interaction plot for", names(fdo$.response())[r], "in", fdoName)
    }

    plot_matrix <- vector("list", numFac * numFac)
    plot_idx <- 1

    for (j in 1:numFac) {
      for (i in 1:numFac) {
        if (i == j) {
          facName <- names(x)[i]
          p_diag <- ggplot() +
            labs(x = facName, y = "") +
            theme_void() +
            theme(
              plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
              axis.title.x = element_text(size = 20, face = "bold", margin = margin(0, 0, 10, 0)),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()
            )
          plot_matrix[[plot_idx]] <- p_diag
        } else if (i < j) {
          plot_matrix[[plot_idx]] <- plot_list[[paste(names(x)[i], names(x)[j])]]
        } else {
          plot_matrix[[plot_idx]] <- ggplot() + theme_void()  # Empty plot
        }
        plot_idx <- plot_idx + 1
      }
    }

    plot_matrix <- matrix(plot_matrix, ncol = numFac, byrow = TRUE)
    final_plot <- wrap_plots(plot_matrix, ncol = numFac)
    final_plot <- final_plot + plot_annotation(
      title = main[r],
      theme = theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 20)))
    )
    print(final_plot)
  }

  invisible()
}




### Uso interactionPlot#########################################################
interactionPlot(dfac)


### lm #########################################################################
m1 <- dfac$lm(rend ~ A*B*C)
summary(m1)


### necesito .splitDev ##########################################
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
### Funcion paretoPlot#################################
paretoPlot <- function(fdo, abs = TRUE, decreasing = TRUE, alpha = 0.05,
                       response = NULL, ylim, xlab, ylab, main, p.col, legend_left = TRUE) {
  # library(RColorBrewer)
  # Esta librería tiene los colores:
  # Set1, Set2, Set3, Pastel2, Pastel1,
  # Paired, Dark2, Accent

  if(is.null(response)==FALSE)
  {
    temp=fdo$.response()[response]
    fdo$.response(temp)
  }
  ylimMissing = FALSE
  if (missing(ylim)){
    ylimMissing = TRUE
  }
  if (missing(xlab))
    xlab = ""
  location = "topright"
  if (decreasing == F | abs == F | legend_left == T) {
    location = "topleft"
  }
  xVals = numeric(0)
  sig.neg = NULL
  sig.pos = NULL
  effect.list = vector("list")
  for (j in 1:ncol(fdo$.response())) {
    if (!any(is.na(fdo$.response()[, j]))) {
      if (missing(ylab))
        ylabel = names(fdo$.response())[j]                                ###
      else
        ylabel = ylab                                                   ###
      form = paste("fdo$.response()[,", j, "]~")
      for (i in 1:ncol(fdo$cube)) {
        form = paste(form, names(fdo$cube)[i], sep = "")
        if (i < ncol(fdo$cube))
          form = paste(form, "*", sep = "")
      }
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

        if(missing(p.col)){
          p.col = rep("lightblue", length(effect))
        }
        else{p.col = brewer.pal(length(effect), p.col)} #paste0("Set", p.col))

        if (ylimMissing)
          if (abs)
            ylim = (range(c(0, abs(effect), 1.3 * ME))) * 1.1
        else ylim = (range(c(effect, -1.3 * ME, 1.3 * ME))) * 1.1
        if (abs) {
          if (missing(ylabel))
            ylabel = ""

          p <- ggplot(data.frame(names = factor(names(effect), levels = names(effect)), effect_ = abs(as.vector(effect))),
                      aes(x = names, y = effect_, fill = names)) +
            geom_bar(stat = "identity", color = "black") +
            scale_fill_manual(values=c(p.col)) +
            theme_minimal()+
            geom_text(aes(label = round(effect,2)), vjust = -1, colour = "black") + # etiquetas sobre las barras
            labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  plot.title = element_text(hjust = 0.5),
                  legend.position = "none")

          if(ME >= ylim[1] & ME <= ylim[2] & SME >= ylim[1] & SME <= ylim[2]){
            p <- p +
              geom_hline(yintercept = ME, linetype = "dashed", color = "red") +
              geom_hline(yintercept = SME, linetype = "dashed", color = "red") +
              scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(abs(ME), abs(SME)),
                                                                                     labels = c(paste("ME = ", round(abs(ME), 2)), paste("SME = ", round(abs(SME), 2)))
              ))

          }
          else{
            if(SME >= ylim[1] & SME <= ylim[2]){
              p <- p + geom_hline(yintercept = SME, linetype = "dashed", color = "red") +
                scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(abs(SME)), labels = c(paste("SME = ", round(abs(SME), 2)))))
            }
            if(ME >= ylim[1] & ME <= ylim[2]){
              p <- p + geom_hline(yintercept = ME, linetype = "dashed", color = "red") +
                scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(abs(ME)), labels = c(paste("ME = ", round(abs(ME), 2)))))
            }
          }

        }
        else {
          if (missing(ylabel))
            ylabel = ""
          p <- ggplot(data.frame(names = names(effect), effect_ = abs(as.vector(effect))),
                      aes(x = names, y = effect_, fill = names)) +
            geom_bar(stat = "identity", color = "black") +
            scale_fill_manual(values=c(p.col)) +
            theme_minimal()+
            geom_text(aes(label = round(effect,2)), vjust = -1, colour = "black") + # etiquetas sobre las barras
            labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  plot.title = element_text(hjust = 0.5))

          if(ME >= ylim[1] & ME <= ylim[2] & SME >= ylim[1] & SME <= ylim[2]){
            p <- p +
              geom_hline(yintercept = ME, linetype = "dashed", color = "red") +
              geom_hline(yintercept = SME, linetype = "dashed", color = "red")
            scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(ME, SME), labels = c(paste("ME = ", round(abs(ME), 2)), paste("SME = ", round(abs(SME), 2)) )))

          }
          else{
            if(ME >= ylim[1] & ME <= ylim[2]){
              p <- p + geom_hline(yintercept = ME, linetype = "dashed", color = "red") +
                scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(ME), labels = c(paste("ME = ", round(ME, 2)))))

            }
            if(SME >= ylim[1] & SME <= ylim[2]){
              p <- p + geom_hline(yintercept = SME, linetype = "dashed", color = "red") +
                scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(SME), labels = c(paste("SME = ", round(SME, 2)))))
            }
          }



        }
      }

      else {
        if (missing(main))
          main = "Standardized main effects and interactions"
        effect = ((summary(lm.1)$coefficients[-pmatch("(Intercept)", names(coef(lm.1))), 1])/(summary(lm.1)$coefficients[-pmatch("(Intercept)", names(coef(lm.1))),2]))
        if (all(is.na(effect)))
          stop("effects could not be calculated")
        effect = effect[!is.na(effect)]
        effect.list[[j]] = effect
        if ((df.resid) > 0) {
          sig.pos = -qt(alpha/2, df.resid)
          sig.neg = +qt(alpha/2, df.resid)
        }
        # Ylimits ----

        if (ylimMissing)
          if (abs) {
            tempVec = c(effect, sig.pos)
            tempVec = tempVec[!is.na(tempVec)]
            ylim = c(0, 0.3 + max(tempVec))
          }
        else {
          tempVec1 = c(0, effect, sig.neg, sig.pos)
          tempVec1 = tempVec1[!is.na(tempVec1)]
          tempVec2 = c(abs(effect), sig.pos, sig.neg)
          tempVec2 = tempVec2[!is.na(tempVec2)]
          ylim = c(min(tempVec1)-0.3, max(tempVec2)+0.3)
        }

        if(missing(p.col)){
          p.col = rep("lightblue", length(effect))
        }
        else{p.col = brewer.pal(length(effect), p.col)} #paste0("Set", p.col))
        # Plot ---------
        if (abs) {
          effect = effect[order(abs(effect), na.last = TRUE, decreasing = decreasing)]
          effect = round(effect, 3)

          if (missing(ylabel))
            ylabel = ""

          # plot with abs
          df <- data.frame(Names = factor(names(effect), levels = names(effect)),
                           effect_ = abs(as.vector(effect)))

          p <- ggplot(data = df,
                      aes(x = Names, y = effect_, fill = Names)) +
            geom_bar(stat = "identity", color = "black") +
            scale_fill_manual(values=c(p.col)) +
            theme_minimal()+
            labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  plot.title = element_text(hjust = 0.5)) +
            geom_text(aes(label = round(effect,2)), vjust = -1, colour = "black") + # etiquetas sobre las barras
            geom_hline(yintercept = sig.pos, linetype = "dashed", color = "red") +
            annotate("text", x = max(as.numeric(df$Names)), y = sig.pos, label = paste(round(sig.pos, 2)), vjust = -0.5, color = "red")
        }
        else {
          effect = effect[order((effect), na.last = TRUE, decreasing = decreasing)]
          effect = round(effect, 3)

          if (missing(ylabel))
            ylabel = ""

          df <- data.frame(Names = factor(names(effect), levels = names(effect)),
                           effect_ = as.vector(effect))

          # Plot without abs
          p <- ggplot(df,
                      aes(x = Names, y = effect_, fill = Names)) +
            geom_bar(stat = "identity", color = "black") +
            scale_fill_manual(values=c(p.col)) +
            theme_minimal()+
            labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0)) +
            geom_text(aes(label = effect), vjust = ifelse(df$effect_ > 0, -0.5, 1.5) , colour = "black") + # etiquetas sobre las barras
            geom_hline(yintercept = sig.pos, linetype = "dashed", color = "red") +
            geom_hline(yintercept = sig.neg, linetype = "dashed", color = "red") +
            annotate("text", x = max(as.numeric(df$Names)), y = sig.pos, label = paste(round(sig.pos, 2)), vjust = -0.5, color = "red") +
            annotate("text", x = max(as.numeric(df$Names)), y = sig.neg, label = paste(round(sig.neg, 2)), vjust = 1.5, color = "red")
        }
        myDelta = diff(range(ylim)) * 0.02
      }
      # Legend ----
      titles <- data.frame(Name_title = paste0(names(fdo$cube),": ",fdo$names()))
      for(i in 1:dim(titles)[1]){
        titles$Pos_title[i] <- 0.95 - (0.05 * i)
      }
      caja <- ggplot(data.frame(x = 0,y = 0), aes(x = x, y = y)) +
        theme_bw() +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = -0.5,margin = margin(b = -12), size = 10)
        ) + xlim(c(0.24, 0.26)) + ylim(c(min(titles$Pos_title) - 0.05, max(titles$Pos_title) + 0.05))
      for(i in 1:dim(titles)[1]){
        caja <- caja + annotate("text", x = 0.25, y = titles$Pos_title[i], label = titles$Name_title[i], size = 3.5, hjust = 0.5)
      }
      # insert legend -----
      if(location == "topright"){
        p <- p + inset_element(caja, left = 0.75, right = 1, top = 1,  bottom = 0.80)
      }
      else{
        p <- p + inset_element(caja, left = 0.25, right = 0.05, top = 1,  bottom = 0.80)
      }

    }
  }

  print(p)
  invisible(list(effect.list, plot = p))
}
### Uso paretoPlot#################################
paretoPlot(dfac, decreasing = T, abs = F, p.col = "Pastel1")


#### necesito FitDistr ##############
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
    body(dens) <- parse(text = paste("densfun(x,", paste("parm[",1L:l, "]", collapse = ", "), ", ...)"))
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
#### necesito .lfkp####################
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
#### Funcion normalPlot#####################
normalPlot <- function(fdo, threeWay = FALSE, na.last = NA, alpha = 0.05, response = NULL, sig.col = c("red1", "red2", "red3"), sig.pch = c(1,2,3), main, ylim, xlim, xlab, ylab, pch,  ###
                       col, border = "red", ...){
  fdoName = deparse(substitute(fdo))
  if(is.null(response)==FALSE)
  {
    temp=fdo$.response()[response]
    fdo$.response(temp)
  }
  parList = list()
  parList = list(...)
  if (length(sig.col) < 3)
    sig.col = as.vector(matrix(sig.col, nrow = 1, ncol = 3))
  XLIM=FALSE;YLIM=FALSE
  if (!(class(fdo)[1] == "facDesign"))
    stop(paste(deparse(substitute(fdo)), "is not an object of class facDesign"))
  mainmiss = FALSE
  if (missing(main))
    mainmiss = TRUE
  if (missing(ylim))
    YLIM=TRUE
  if (missing(xlim))
    XLIM=TRUE
  if (missing(xlab))
    xlab = "Coefficients"
  if (missing(ylab))
    ylab = "Theoretical Quantiles"
  if (missing(pch))
    pch = 19
  if (missing(col))
    col = "black"

  list_plot = list()
  for(j in 1:ncol(fdo$.response())){
    parList = list(...)
    params = list()
    leg.col = vector()
    p.col = vector()
    p.pch = vector()
    leg.txt = vector()
    main = paste("Normal plot for", names(fdo$.response())[j], "in", fdoName)
    if (j > 1)
      dev.new()
    form = paste("fdo$.response()[,", j, "]~")

    for (i in 1:ncol(fdo$cube)) {
      form = paste(form, names(fdo$cube)[i], sep = "")
      if (i < ncol(fdo$cube))
        form = paste(form, "*", sep = "")
    }

    lm.1 = lm(as.formula(form), data = fdo$as.data.frame())
    lm.1s = summary(lm.1)
    effect = coef(lm.1s)[row.names(coef(lm.1s)) != "(Intercept)", "t value"]
    print(effect)
    if (all(is.na(effect)))
      effect = 2 * coef(lm.1)[-pmatch("(Intercept)", names(coef(lm.1)))]
    #            stop("effects could not be calculated")
    sig = summary(lm.1)$coefficients[, "Pr(>|t|)"][-pmatch("(Intercept)", names(coef(lm.1)))]
    df.resid = df.residual(lm.1)
    nc = nrow(fdo$centerCube)

    tQ = ppoints(effect)
    index = order(effect)
    sQ = effect[index]
    sig = sig[index]

    if (df.resid > 0) {
      # obtenemos el caracter de la cajita del p_value
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
    }else{p.col=col
    p.pch=pch}

    mid = round(length(tQ)/2)
    last = length(tQ)
    params$p = ppoints(effect)
    estimates = FitDistr(effect, "normal")   #estimates = MASS::fitdistr(effect, "normal")
    params$mean = estimates$estimate[["mean"]]
    params$sd = estimates$estimate[["sd"]]

    y = do.call(qnorm, params)

    if (XLIM)
      xlim = range(sQ)
    if (YLIM)
      ylim = range(y)

    # PLOT -----------------------
    df <- data.frame(sQ = names(sQ), value = as.numeric(sQ), y = y)

    p <- ggplot(df, aes(x = value, y = y, label = sQ)) +
      geom_point(col = p.col, pch = p.pch) +
      theme_classic() + lims(x = xlim, y = ylim) +
      labs(x = xlab, y = ylab, title = main) +
      geom_text(check_overlap = TRUE, vjust = 1) + theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))

    xp = c(qnorm(0.1), qnorm(0.99))
    yp = c(qnorm(0.1, mean = estimates$estimate[["mean"]], sd = estimates$estimate[["sd"]]), qnorm(0.99, mean = estimates$estimate[["mean"]], sd = estimates$estimate[["sd"]]))
    slope = (yp[2] - yp[1])/(xp[2] - xp[1])
    int = yp[1] - slope * xp[1]

    # line
    p <- p + geom_abline(intercept = int, slope = slope, col = border)

    # legend
    if (df.resid > 0){
      caja <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
        theme_bw() +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        xlim(c(0.25,0.30)) + ylim(c(0.24, 0.31))

      caja <- caja +
        annotate('text', x = 0.275, y = 0.28,
                 label = leg.txt, size = 3, hjust = 0.5, colour = leg.col)

      p <- p + inset_element(caja, left = 0.01, right = 0.2, top = 1, bottom = 0.85)
    }
    print(p)
    list_plot[[paste0("p",j)]] <- p
  }
  invisible(list(effect = effect, plots = list_plot))
}
#### Uso normalPlot#################
normalPlot(dfac)


#### necesito clase desirability.c#################
desirability.c <- R6Class("desirability", public = list(response = NULL,
                                                      low = NULL,
                                                      high = NULL,
                                                      target = NULL,
                                                      scale = NULL,
                                                      importance = NULL,
                                                      initialize = function(response=NULL, low=NULL, high=NULL, target=NULL, scale=NULL, importance=NULL) {
                                                        self$response <- response
                                                        self$low <- low
                                                        self$high <- high
                                                        self$target <- target
                                                        self$scale <- scale
                                                        self$importance <- importance
                                                      },
                                                      print = function(){
                                                        if (!is.numeric(self$target))
                                                          cat("Target is to", paste(self$target, "imize", sep = ""), self$response, "\n")
                                                        else cat("Target is ", self$target, " for", self$response, "\n")
                                                        cat("lower Bound: ", self$low, "\n")
                                                        cat("higher Bound: ", self$high, "\n")
                                                        if (is.numeric(self$target))
                                                          cat("Scale factor is: low =", self$scale[1], "and high =", self$scale[2], "\n")
                                                        else if (identical("min", self$target) | identical("max", self$target))
                                                          cat("Scale factor is: ", self$scale, "\n")
                                                        cat("importance: ", self$importance, "\n")
                                                        cat("\n")
                                                      },
                                                      plot = function(y, scale, main, xlab, ylab, type, col, numPoints = 500, ...){
                                                        xm1 = NULL
                                                        xm2 = NULL
                                                        ym = NULL
                                                        y = NULL
                                                        if (missing(main))
                                                          main = paste("Desirability function for", self$response)
                                                        if (missing(xlab))
                                                          xlab = self$response
                                                        if (missing(ylab))
                                                          ylab = "Desirability"
                                                        if (missing(type))
                                                          type = "l"
                                                        if (missing(scale))
                                                          scale = self$scale
                                                        if (missing(col))
                                                          col = 1:length(scale)
                                                        dFun = .desireFun(self$low, self$high, self$target, self$scale, self$importance)
                                                        xVals = seq(self$low - 0.04 * diff(range(self$low, self$high)), self$high + 0.04 * diff(range(self$low, self$high)), length = numPoints)
                                                        yVals = dFun(xVals)
                                                        plot(xVals, yVals, main = main, xlab = xlab, ylab = ylab, type = type, col = col, ...)
                                                        if (is.numeric(self$target)) {
                                                          xm1 = mean(c(par("usr")[1], self$target))
                                                          xm2 = mean(c(par("usr")[2], self$target))
                                                          ym1 = yVals[max((1:numPoints)[xVals <= xm1])]
                                                          ym2 = yVals[max((1:numPoints)[xVals <= xm2])]
                                                          text(xm1 + 0.025 * diff(range(par("usr")[1:2])), ym1, paste("scale =", scale[1]), adj = c(0, 0))
                                                          text(xm2 - 0.025 * diff(range(par("usr")[1:2])), ym2, paste("scale =", scale[2]), adj = c(1, 1))
                                                        }
                                                        else {
                                                          xm1 = mean(par("usr")[c(1, 2)])
                                                          ym1 = yVals[max((1:numPoints)[xVals <= xm1])]
                                                          if (identical(self$target, "max"))
                                                            text(xm1 + 0.025 * diff(range(par("usr")[1:2])), ym1 - 0.025 * diff(range(par("usr")[3:4])), paste("scale =", scale[1]), adj = c(0, 0))
                                                          else text(xm1 + 0.025 * diff(range(par("usr")[1:2])), ym1 + 0.025 * diff(range(par("usr")[3:4])), paste("scale =", scale[1]), adj = c(0, 1))
                                                        }
                                                        out = list(x = xVals, y = yVals)
                                                        names(out) = c(self$response, "desirability")
                                                        invisible(out)
                                                      }



                                                      )
                       )

#### necesito .desireFun###############
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


#### Funcion wirePlot###################
wirePlot <- function(x, y, z, data = NULL,
                     xlim, ylim, zlim, main,
                     xlab, ylab, sub, sub.a = TRUE, zlab,
                     form = "fit", col = "Rainbow", steps,
                     fun, plot = TRUE, show.scale = TRUE, n.scene = "scene") {
  form = form
  fact = NULL
  if (missing(steps))
    steps = 25
  fdo = data
  fit = NULL
  lm.1 = NULL

  # Col puede ser: "Rainbow", "Jet", "Earth", "Electric"
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

  if (missing(main))
    main = paste("Response Surface for", z.c)

  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[.NAMES[i]]] <-fdo$names()[i]
  }
  if (missing(ylab))
    ylab = paste(y.c, ": ", aux[[y.c]])
  if (missing(xlab))
    xlab = paste(x.c, ": ", aux[[x.c]])
  if (missing(zlab))
    zlab = paste(x.c, ": ", z.c)

  if (missing(xlim))
    xlim = c(min(fdo$get(, x.c)), max(fdo$get(, x.c)))
  if (missing(ylim))
    ylim = c(min(fdo$get(, y.c)), max(fdo$get(, y.c)))

  allVars = c(fdo$names(), names(fdo$.response()))
  isct = intersect(c(aux[[x.c]], aux[[y.c]], z.c), c(fdo$names(), names(fdo$.response())))

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
    obj = fdo$desires()[[z.c]]
    fun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
  }

  if (form %in% c("fit")) {
    lm.1 = fdo$fits[[z.c]]
    if (is.null(fit))
      form = "full"
  }

  if (form %in% c("quadratic", "full", "interaction", "linear")) {
    if (identical(form, "full")) {
      form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
      if (nrow(fdo$star) > 0)
        form = paste(form, "+ I(", x.c, "^2) + I(", y.c, "^2)")
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
  }

  if (is.null(form))
    stop(paste("invalid formula", form))
  if (is.null(lm.1))
    lm.1 = fdo$lm(form)
  if (missing(sub))
    sub = deparse(formula(lm.1))

  dcList = vector(mode = "list", length = length(fdo$names()))
  names(dcList) = names(aux)
  dcList[1:length(fdo$names())] = 0

  help.predict = function(x, y, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = x
    dcList[[y.c]] = y
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
  }

  xVec = seq(min(xlim), max(xlim), length = steps)
  yVec = seq(min(ylim), max(ylim), length = steps)

  mat = outer(xVec, yVec, help.predict, x.c, y.c, lm.1)

  if (is.function(fun))
    mat = try(apply(mat, c(1, 2), fun))
  if (identical(fun, "overall")) {
    main = "composed desirability"
    mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
    for (i in names(fdo$.response())) {
      obj = fdo$desires()[[i]]
      fun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
      temp = outer(xVec, yVec, help.predict, x.c, y.c, fits(fdo)[[i]])
      temp = try(apply(temp, c(1, 2), fun))
      mat = mat * temp
    }
    mat = mat^(1/length(names(fdo$response())))
  }

  if (missing(zlim))
    zlim = range(mat)

  p <- plot_ly(x = -yVec, y = xVec, z = mat, colorscale=col, scene = n.scene) %>%
    add_surface(showscale = show.scale) %>%
    layout(
      title = main,
      scene = list(
        xaxis = list(range = ylim, title = ylab, zeroline = FALSE),
        yaxis = list(range = xlim, title = xlab, zeroline = FALSE),
        zaxis = list(range = zlim, title = zlab, zeroline = FALSE),
        camera = list(eye = list(x=2, y=2, z=0.1))
      ),
      margin = list(l = 10, r = 15, t = 30, b = 20)
    )
  if(sub.a){
    p <- p %>%
      layout(
        annotations = list(
          list(
            text = sub,# Subtitulo
            x = 0.5,   # Posición x en la mitad de la gráfica
            y = -0.1,  # Posición y debajo de la gráfica
            showarrow = FALSE,
            font = list(size = 12)
          )
        )
      )
  }
  if (plot) {
    show(p)
  }
  invisible(list(x = xVec, y = yVec, z = mat, plot = p))
}
#### Uso wirePlot###############################
wirePlot(A,B,rend,data=dfac)

### Funcion contourPlot#####################
contourPlot <- function(x, y, z, data = NULL, xlim, ylim, main, xlab, ylab, border, sub, zlab, form = "fit", phi, theta, ticktype, col = 1, steps,
                       factors, fun, plot = TRUE, show.scale = TRUE) {
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

  if (is.null(data)) {
    cat("\n defaulting to filled.contour function\n")
    return("persp")
  }
  if (class(data)[1] != "facDesign") {
    cat("\n defaulting to filled.contour function using formula\n")
    return("persp")
  }

  x.c = deparse(substitute(x))
  y.c = deparse(substitute(y))
  z.c = deparse(substitute(z))

  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[.NAMES[i]]] <-fdo$names()[i]
  }
  if (missing(plot))
    plot = TRUE
  if (missing(main))
    main = paste("Filled Contour for", z.c)
  if (missing(ylab))
    ylab = paste(y.c, ": ", aux[[y.c]])
  if (missing(xlab))
    xlab = paste(x.c, ": ", aux[[x.c]])
  if (missing(factors))
    factors = NULL
  if (missing(xlim))
    xlim = c(min(fdo$get(, x.c)), max(fdo$get(, x.c)))
  if (missing(ylim))
    ylim = c(min(fdo$get(, y.c)), max(fdo$get(, y.c)))
  allVars = c(fdo$names(), names(fdo$.response()))
  isct = intersect(c(aux[[x.c]], aux[[y.c]], z.c), c(fdo$names(), names(fdo$.response())))

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
    obj = fdo$desires()[[z.c]]
    fun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
  }
  if (form %in% c("fit")) {
    lm.1 = fdo$fits[[z.c]]
    if (is.null(fit))
      form = "full"
  }
  if (form %in% c("quadratic", "full", "interaction", "linear")) {
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
      if (nrow(fdo$star) > 0)
        form = paste(form, "+ I(", x.c, "^2) + I(", y.c, "^2)")
    }
  }

  if (is.null(form))
    stop(paste("invalid formula", form))
  if (is.null(lm.1))
    lm.1 = fdo$lm(form)

  dcList = vector(mode = "list", length = length(fdo$names()))
  names(dcList) = names(aux)
  dcList[1:length(fdo$names())] = 0

  if (!is.null(factors)) {
    for (i in fdo$names()) dcList[[i]] = factors[[i]][1]
  }

  help.predict = function(x, y, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = x
    dcList[[y.c]] = y
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
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
      obj = fdo$desires()[[i]]
      fun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
      temp = outer(xVec, yVec, help.predict, x.c, y.c, fdo$fits[[i]])
      temp = try(apply(temp, c(1, 2), fun))
      mat = mat * temp
    }
    mat = mat^(1/length(names(fdo$.response())))
  }

  p <- plot_ly(x = xVec, y = yVec, z = mat, colors = color,
               type = "contour", contours = list(coloring = 'heatmap'), showscale = show.scale) %>%
    layout(
      title = main,
      xaxis = list(range = xlim, title = xlab, zeroline = FALSE),
      yaxis = list(range = ylim, title = ylab, zeroline = FALSE)
    )

  if (plot) {
    show(p)
  }

  invisible(list(x = xVec, y = yVec, z = mat, plot = p))
}
### Uso contourPlot########################
contourPlot(A,B,rend,data=dfac)




#########DISEÑOS FACTORIALES FRACCIONARIOS 2^k-p#################
###necesito confounds####
confounds <- function(x, depth = 2) {
  DB = FALSE
  varName = deparse(substitute(x))
  identityList = x$identity()
  x = x$cube
  if (length(identityList) < 1) {
    print(paste(varName, " contains no defining relations!"))
    invisible()
  }
  effect1 = numeric(0)
  effect2 = numeric(0)
  if (DB)
    identityList
  index = numeric(0)
  for (i in 1:(dim(x)[2])) {
    if (!(TRUE && all(unique(x[, i]) %in% c(-1, 1))))
      index = c(index, i)
  }
  if (length(index) > 0)
    x = x[, -index]
  if (DB)
    cat("Column(s) ", index, " are discarded for analysis\n")
  n = dim(x)[2]
  if (n <= 1)
    stop("Factorial Design contains only one row!")
  for (j in 1:length(identityList)) {
    ident = identityList[[j]]
    for (m in 1:n) {
      combMat = combn(1:n, m)
      for (i in 1:(dim(combMat)[2])) {
        isect = intersect(ident, combMat[, i])
        conf = setdiff(ident, isect)
        conf = sort(c(conf, setdiff(combMat[, i], isect)))
        effect1 = c(effect1, paste(sort(names(x)[as.numeric((combMat[, i]))]), sep = "",
                                   collapse = ""))
        effect2 = c(effect2, paste(sort(names(x)[conf]), sep = "", collapse = ""))
        if (DB) {
          cat("Effect(s) ", as.numeric((combMat[, i])), " aliased with Effect(s)", conf,
              "\n")
          cat("Effect(s)", names(x)[as.numeric((combMat[, i]))], " aliased with Effects ",
              names(x)[conf], "\n")
        }
      }
    }
  }
  if (DB) {
    cat(effect1, "\n")
    cat(effect2, "\n")
  }
  if (length(effect1) > 0)
    dupIndex = numeric(0)
  for (i in 1:length(effect1)) {
    if (DB) {
      cat("i: ", i, "\tlength(effect1): ", length(effect1), "\n")
      cat("effect 1 : ", effect1, "\n")
    }
    if (i > length(effect1))
      break
    index = (1:length(effect1))[effect2 == effect1[i]]
    if (DB)
      cat("index: ", index, "\n")
    dupIndex = numeric(0)
    for (j in index) {
      if (effect1[j] == effect2[i]) {
        if (i != j)
          dupIndex = c(dupIndex, j)
      }
    }
    if (length(dupIndex > 0)) {
      effect1 = effect1[-dupIndex]
      effect2 = effect2[-dupIndex]
    }
  }
  cat("\nAlias Structure:\n")
  for (i in 1:length(effect1)) {
    if ((length(strsplit(effect1[i], split = character(0))[[1]]) <= depth) && (length(strsplit(effect2[i],
                                                                                               split = character(0))[[1]]) <= depth))
      cat(effect1[i], "\tis confounded with\t", effect2[i], "\n")
    if (identical(depth, "all"))
      cat(effect1[i], "\tis confounded with\t", effect2[i], "\n")
  }
  invisible(effect1)
}

####Uso diseños factoriales fraccionarios####
dfacfrac <- fracDesign(k=3,gen='C=AB',centerCube = 4)
dfacfrac$summary()
aliasTable(dfacfrac)
confounds(dfacfrac)

####Funcion fracChoose####
fracChoose <- function() {
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
  if (class(fdo)[1] == "facDesign")
    return(fdo)
  else return(genList[[mat[y, x]]])
}

####Uso fracChoose()####
m1 <- fracChoose()
m1$summary()


#######DISEÑOS REPLICADOS Y PUNTOS CENTRALES###############
dfac2 <- facDesign(k = 3, centerCube = 2, replicates = 1)
dfac2$summary()

#######RESPUESTAS MÚLTIPLES###############################
set.seed(1234)
y2 <- rnorm(12,mean=120)
dfac$.response(data.frame(rend,y2))
dfac$.response()
#graficos
wirePlot(A, B, rend, data = dfac, form = "rend~A+B+C+A*B")
contourPlot(A, B, y2, data = dfac, form = "y2~A+B+C+A*B")

#AJUSTAR EL FACTOR
wirePlot(A, B, y2, data = dfac, factors = list(C=-1), form = "y2~A*B*C")
wirePlot(A, B, y2, data = dfac, factors = list(C=1), form = "y2~A*B*C")

#FITS
dfac$set.fits(dfac$lm(rend~A+B))
dfac$set.fits(dfac$lm(y2~A*B*C))
dfac$fits

#####Pasar a un entorno de proceso con un mayor rendimiento esperado######
#####Necesito clase steepAscent.c##########
steepAscent.c <- R6Class("facDesign", public = list(name = NULL,
                                                    X = data.frame(),
                                                    response = data.frame(),
                                                    .response = function(value){
                                                      if (missing(value)) {
                                                        return(self$response)
                                                      }
                                                      else{
                                                        if (is.vector(value)) {
                                                          temp = data.frame(value)
                                                          names(temp) = deparse(substitute(value))
                                                          if (nrow(self$X) == nrow(temp)) {
                                                            self$response = temp
                                                            invisible(self)
                                                          }
                                                          else{
                                                            stop("number of rows differ!")
                                                          }
                                                        }
                                                        else if (is.data.frame(value)) {
                                                          if (nrow(self$X) == nrow(value)) {
                                                            self$response = value
                                                            invisible(self)
                                                          }
                                                          else{
                                                            stop("number of rows differ!")
                                                          }
                                                        }
                                                        else{
                                                          stop(paste(deparse(substitute(value)), " needs to be a vector or data.frame"))
                                                        }
                                                      }
                                                    },
                                                    get = function(i, j){
                                                      bound = ncol(self$X)
                                                      if (j <= bound)
                                                        self$X[i, j]
                                                      else self$response[i, j - bound]
                                                    },
                                                    as.data.frame = function(row.names = NULL, optional = FALSE, ...){
                                                      return(cbind(self$X, self$response))
                                                    },
                                                    print = function(){
                                                      print(self$as.data.frame())
                                                    },
                                                    plot = function(y, ...){
                                                      Delta = (self$X)$Delta
                                                      frame = cbind(Delta, self$.response())
                                                      names(frame) = c("Delta", names(self$.response()))
                                                      plot(frame, ...)

                                                    }


                                                   )
                        )

#####Necesito code2real####
code2real = function(low, high, codedValue) {
  return((diff(c(low, high))/2) * codedValue + mean(c(low, high)))
}
######Necesito funcion steepAscent############
steepAscent <- function(factors, response, size = 0.2, steps = 5, data) {
  DB = FALSE
  if (missing(data))
    return("missing an object of class 'facDesign'")
  else fdo = data
  if (missing(factors) | length(factors) < 1)
    return("missing factors")
  if (!is.character(factors))
    return("factors needs to be a character")
  #names(names(fdo))
  model = data$fits[[response]]
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
    model = fdo$lm(form)
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
  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[.NAMES[i]]] <-fdo$names()[i]
  }
  for (i in seq(along = factors)) {
    frameOut[, i + initial] = code2real(fdo$lows()[[aux[i][[1]]]], fdo$highs()[[aux[i][[1]]]], x[i] * 0:steps)
    names(frameOut)[i + initial] = paste(factors[i], ".real", collapse = "", sep = "")
    if (DB)
      print(factors[i])
  }
  soa = steepAscent.c$new()
  soa$X = frameOut
  soa$.response(rep(NA, times = nrow(frameOut)))
  names(soa$response) = deparse(substitute(response))
  cat("\n")
  cat(paste("Steepest Ascent for", deparse(substitute(data)), "\n"))
  cat("\n")
  print(format(frameOut, digits = 3))
  invisible(soa)
}

######USO steepAscent################################
sao <- steepAscent(factors = c("A", "B"), response = "rend", data = dfac, steps = 20)
sao$print()
predicted <- simProc(sao$get(,5), sao$get(,6))
sao$.response(predicted)
sao$plot(type='b', col=2)





########DISEÑOS DE SUPERFICIE DE RESPUESTA#################
set.seed(1234)
fdo2 <- facDesign(k = 2, centerCube = 3)
fdo2$names(c("Factor1","Factor2"))
fdo2$lows(c(134,155))
fdo2$highs(c(155,175))
rend <- c(simProc(134,175),simProc(144.5,165.5),simProc(155,155),simProc(144.5,165.5),simProc(155,175),simProc(144.5,165.5),simProc(134,155))
fdo2$.response(rend)
###Necesito .nblock####
.nblock <- function(fdo) {
  if (class(fdo)[1] != "facDesign")
    stop(paste(deparse(substitute(fdo)), "needs to be an object of class 'facDesign'"))
  return(length(unique(fdo$block[[1]])))
}
###Necesito .starFrame####
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
###Necesito .alphaOrth, .alphaRot, .rsmOrth####
.alphaOrth <- function(k, p = 0, cc, cs) {
  alpha = sqrt((2^(k - p) * (2 * k + cs))/(2 * (2^(k - p) + cc)))
  return(alpha)
}
.alphaRot <- function(k, p = 0) {
  alpha = (2^(k - p))^0.25
  return(alpha)
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
###Necesito starDesign####
starDesign <- function(k, p = 0, alpha = c("both", "rotatable", "orthogonal"), cs, cc, data) {
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
    fdo = data$clone()
    k = ncol(fdo$cube)
    if (class(fdo)[1] != "facDesign") {
      stop(paste(deparse(substitute(data)), "needs to be an object of class 'facDesign'"))
    }
    if (nrow(fdo$star) > 0)
      stop(paste("star portion of", deparse(substitute(data)), "not empty"))
    k = length(fdo$names())
    nameVec = fdo$names()
    cc = nrow(fdo$centerCube)
    p = ncol(fdo$cube) - log(nrow(unique(fdo$cube)), 2)
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
    fdo$.star(starFrame)
  if (DB)
    print("starFrame added")
  if (cs > 0) {
    csFrame = as.data.frame(matrix(0, nrow = cs, ncol = k))
    names(csFrame) = nameVec
    if (!missing(data)) {
      fdo$.centerStar(csFrame)
      if (DB)
        print("csFrame added")
    }
  }
  if (cc > 0) {
    ccFrame = as.data.frame(matrix(0, nrow = cc, ncol = k))
    names(ccFrame) = nameVec
    if (!missing(data)) {
      fdo$.centerCube(ccFrame)
      if (DB)
        print("ccFrame added")
    }
  }
  if (!missing(data))
    return(fdo)
  else return(list(star = starFrame, centerStar = csFrame, centerCube = ccFrame))
}

###Uso starDesign#####
rsdo <- starDesign(data=fdo2)
rsdo$print()
rend2 <- c(rend,
  simProc(130, 165),
  simProc(155, 165),
  simProc(144, 155),
  simProc(144, 179),
  simProc(144, 165),
  simProc(144, 165),
  simProc(144, 165)
)

rsdo$.response(rend2)
rsdo$.response()

lm.3 <- rsdo$lm(rend2 ~ A*B + I(A^2) + I(B^2))
summary(lm.3)


####Visualizacion####
wirePlot(A, B, rend2, form = "rend2 ~ A*B + I(A^2) + I(B^2)", data = rsdo, theta = -70)
contourPlot(A, B, rend2, form = "rend2 ~ A*B + I(A^2) + I(B^2)", data = rsdo)
####filled.contour####
A = seq(40, 210, length = 100)
B = seq(90, 190, length = 100)
C = seq(90, 190, length = 100)
filled.contour(A, B,outer(A,B, simProc, noise = FALSE), xlab = "Factor 1", ylab = "Factor 2", color.palette = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))

####Necesito rsmDesign####
rsmDesign <- function(k = 3, p = 0, alpha = "rotatable", blocks = 1, cc = 1, cs = 1, fp = 1,
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
    temp = as.data.frame(matrix(0, nrow = cc, ncol = ncol(fdo$cube)))
    names(temp) = names(fdo$cube)
    fdo$centerCube(temp)
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
  names(starportion) = names(fdo$cube)
  fdo$.star(starportion)
  if (DB)
    print("star added")
  if (cs > 0) {
    temp = as.data.frame(matrix(0, nrow = cs, ncol = ncol(fdo$cube)))
    names(temp) = names(fdo$cube)
    fdo$.centerStar(temp)
  }
  #    return(fdo)
  fdo = blocking(fdo, blocks)
  return(fdo)
}
####Uso rsmDesign####
fdo <- rsmDesign(k=3, alpha=1.633, cc=0, cs=6)

####poner en orden estandar####
fdo <- randomize(fdo,so = TRUE)
fdo$summary()
###necesito rsmChoose()####
rsmChoose <- function() {
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

###uso rsmChoose()####
rsdo <- rsmChoose()



###MONTAJE SECUENCIAL#######################
fdo3 <- facDesign(k = 6)
fdo3$summary()
rsdo <- starDesign(alpha = "orthogonal", data = fdo3)
rsdo$summary()
###ALEATORIZACIÓN##########################
randomize(fdo, random.seed = 123)
fdo$summary()
randomize(fdo, so = TRUE)
fdo$summary()
###BLOQUEO###################################
blocking(fdo,blocks = 1)
fdo$summary()
###DESEABILIDADES############################
###Funcion desirability#############
desirability = function(response, low, high, target = "max", scale = c(1, 1), importance = 1, constraints) {
  if (low >= high)
    stop("the lower bound must be greater than the high bound!")
  if (any(scale <= 0))
    stop("the scale parameter must be greater than zero!")
  if (!is.numeric(target) & !identical(tolower(target), "min") & !identical(tolower(target), "max"))
    stop("target needs to be \"min\", \"max\" or a numeric value")
  return(desirability.c$new(response = deparse(substitute(response)), low = low, high = high, target = target, scale = scale, importance = importance))
}

##EJEMPLO:####
y1 <- c(102, 120, 117, 198, 103, 132, 132, 139, 102, 154, 96, 163, 116,
       153, 133, 133, 140, 142, 145, 142)
y2 <- c(900, 860, 800, 2294, 490, 1289, 1270, 1090, 770, 1690, 700, 1540,
       2184, 1784, 1300, 1300, 1145, 1090, 1260, 1344)
y3 <- c(470, 410, 570, 240, 640, 270, 410, 380, 590, 260, 520, 380, 520,
       290, 380, 380, 430, 430, 390, 390)
y4 <- c(67.5, 65, 77.5, 74.5, 62.5, 67, 78, 70, 76, 70, 63, 75, 65, 71,
       70, 68.5, 68, 68, 69, 70)

d1 <- desirability(y1, 120, 170, scale = c(1, 1), target = "max")
d3 <- desirability(y3, 400, 600, target = 500)
d1$print()
d1$plot(col = 2)
d3$plot(col = 2)


d<-desirability$new()
desirability


####Utilizacion de deseabilidades junto con experimentos desiñados#####
ddo <- rsmDesign(k = 3, alpha = 1.633, cc = 0, cs = 6)
ddo <- randomize(ddo,so=TRUE)
ddo$summary()
ddo$names(c("silica", "silan", "sulfur"))
ddo$highs(c(1.7, 60, 2.8))
ddo$lows(c(0.7, 40, 1.8))
ddo$summary()
#asignar response
ddo$.response(data.frame(y1, y2, y3, y4)[c(5, 2, 3, 8, 1, 6, 7, 4, 9:20), ])
d2 <- desirability(y2, 1000, 1300, target = "max")
d4 <- desirability(y4, 60, 75, target = 67.5)
#asignar desires
ddo$desires(d1)
ddo$desires(d2)
ddo$desires(d3)
ddo$desires(d4)
ddo$desires()
#asignar fits
ddo$set.fits(ddo$lm(y1 ~ A + B + C + A:B + A:C + B:C + I(A^2) + I(B^2) + I(C^2)))
ddo$set.fits(ddo$lm(y2 ~ A + B + C + A:B + A:C + B:C + I(A^2) + I(B^2) + I(C^2)))
ddo$set.fits(ddo$lm(y3 ~ A + B + C + A:B + A:C + B:C + I(A^2) + I(B^2) + I(C^2)))
ddo$set.fits(ddo$lm(y4 ~ A + B + C + A:B + A:C + B:C + I(A^2) + I(B^2) + I(C^2)))
ddo$fits
###Necesito .validizeConstraints####
.validizeConstraints = function(fdo, constraints) {
  X = fdo$as.data.frame()
  csOut = vector(mode = "list")
  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[fdo$names()[i]]] <-.NAMES[i]
  }
  for (i in aux) {
    csOut[[i]] = c(min(X[, i]), max(X[, i]))
  }
  if (missing(constraints))
    return(csOut)
  cs2 = constraints[fdo$names()]
  cs2 = cs2[!unlist(lapply(cs2, is.null))]
  cs2 = cs2[(unlist(lapply(cs2, length)) == 2)]
  csOut[names(cs2)] = cs2[names(cs2)]
  return(csOut)
}
###Necesito funcion overall####
overall <- function(fdo, steps = 20, constraints, ...) {
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
###Necesito .desHelp####
.desHelp = function(fdo, factors, ...) {
  if (length(factors) != length(fdo$names()))
    stop("not enough factors specified in factors")
  if (any(is.na(factors)))
    stop("factors contain NA")
  yCharSet = intersect(names(fdo$desires()), names(fdo$fits))
  desList = fdo$desires()
  fitList = fdo$fits
  yDes = vector(mode = "list")
  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[fdo$names()[i]]] <-.NAMES[i]
  }
  names(factors)<-unlist(aux)
  for (y in yCharSet) {
    obj = desList[[y]]
    dFun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
    lm.y = fitList[[y]]
    yHat = predict(lm.y, newdata = data.frame(factors))
    yDes[[y]] = dFun(yHat)
  }
  return(yDes)
}
###Necesito clase desOpt#################
desOpt <- R6Class("desOpt", public = list(facCoded = list(),
                                          facReal = list(),
                                          responses = list(),
                                          desirabilities = list(),
                                          overall = NULL,
                                          all = data.frame(),
                                          fdo = NULL,
                                          as.data.frame = function(x, row.names = NULL, optional = FALSE, ...) {
                                            return(x$all)
                                          },
                                          print = function(){
                                            cat(paste("\ncomposite (overall) desirability:", format(self$overall, digits = 3)))
                                            cat("\n")
                                            cat("\n")
                                            temp1 = do.call(data.frame, self$facCoded)
                                            temp2 = do.call(data.frame, self$facReal)
                                            facFrame = rbind(temp1, temp2)
                                            row.names(facFrame) = c("coded", "real")
                                            show(format(facFrame, digits = 3))
                                            temp1 = do.call(data.frame, self$responses)
                                            temp2 = do.call(data.frame, self$desirabilities)
                                            respDesFrame = rbind(temp1, temp2)
                                            row.names(respDesFrame) = c("Responses", "Desirabilities")
                                            cat("\n")
                                            show(format(respDesFrame, digits = 3))
                                          }
                                         )
                      )
###Necesito .dHelp#####################
.dHelp = function(model, dFun) {
  lm1 = model
  d1 = dFun
  out = function(newdata) {
    return(d1(predict(lm1, newdata = newdata)))
  }
  return(out)
}
###funcion optimum####
optimum <- function(fdo, constraints, steps = 25, type = "grid", start, ...) {
  DB = FALSE
  if (missing(fdo))
    stop("missing fdo!")
  X = fdo$as.data.frame()
  numFac = length(fdo$names())
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
  desOpt = desOpt$new()
  desOpt$fdo = fdo
  facCoded = NA
  desirabilities = NA
  overall = NA
  setList = list()
  dList = list()
  importances = list()
  yCharSet = intersect(names(fdo$desires()), names(fdo$fits))
  for (y in yCharSet) {
    obj = fdo$desires()[[y]]
    dFun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
    lm.y = fdo$fits[[y]]
    importances[[y]] = fdo$desires()[[y]]$importance
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
    names(facCoded) = fdo$names()
    desOpt$facCoded = facCoded
    overall = temp$value
    desirabilities = .desHelp(fdo, desOpt$facCoded)
  }
  if (type == "gosolnp") {
    #if (!require(Rsolnp, quietly = TRUE))
    #    stop("Package Rsolnp needs to be installed!")
    temp = Rsolnp::gosolnp(fun = dAllRsolnp, LB = lower, UB = upper)
    facCoded = as.list(temp$pars)
    names(facCoded) = fdo$names()
    desOpt$facCoded = facCoded
    overall = -rev(temp$values)[1]
    desirabilities = .desHelp(fdo, desOpt$facCoded)
  }
  if (type == "grid") {
    dVals = overall(fdo = fdo, constraints = constraints, steps = steps)
    index = order(dVals[, "overall"], decreasing = TRUE)[1]
    desOpt$all = dVals
    desOpt$facCoded = as.list(dVals[index, fdo$names()])
    desirabilities = as.list(dVals[index, names(fdo$.response())])
    names(desirabilities) = names(fdo$.response()) #fix for the case of having just one response
    overall = dVals[index, "overall"]
  }
  for (i in names(desOpt$facCoded)) {
    desOpt$facReal[[i]] = code2real(fdo$lows()[[i]], fdo$highs()[[i]], desOpt$facCoded[[i]])
  }
  desOpt$desirabilities = desirabilities
  desOpt$overall = overall
  newData = do.call(data.frame, desOpt$facCoded)
  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[fdo$names()[i]]] <-.NAMES[i]
  }
  names(newData)<-unlist(aux)
  for (i in names(desOpt$desirabilities)) {
    desOpt$responses[[i]] = predict(fdo$fits[[i]], newData)
  }
  return(desOpt)
}
###uso optimum####
optimum(ddo,type='optim')


