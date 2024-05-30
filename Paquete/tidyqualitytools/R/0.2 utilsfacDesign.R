
### .replace2s####
.replace2s = function(x) {
  if (!is.data.frame(x))
    stop(paste(deparse(substitute(x)), "needs to be a data.frame"))
  for (i in 1:ncol(x)) x[x[, i] == 2, i] = -1
  return(x)
}

### .helpAliasTable##################################################
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

### aliasTable#######################################################
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

### .fdoOrth y .NAMES#######################################
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



### .m.interaction.plot###########################
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
  p <- ggplot2::ggplot(df, aes(x = x, y = y)) +
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

### clase facDesign.c######################################################
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

### clase doeFactor####################
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




### randomize####
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


### .rsm#################
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


### .numFac#############
.numFac = function(fdo) {
  return(length(fdo$names()))
}


### .confoundings####
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


### .lociv####
.lociv = function(charVec) {
  lenVec = numeric(length = length(charVec))
  for (i in seq(along = charVec)) {
    lenVec[i] = length(strsplit(charVec[i], split = "")[[1]])
  }
  return(lenVec)
}

### .blockInteractions######
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



### .blockGenCol##########
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


### .blockCol##############
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


### blocking######
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


### .norm2d#####################################################
.norm2d <- function(x1, x2, mu1 = 160, mu2 = 165, rho = 0.7, sigma1 = 45, sigma2 = 22.5) {
  z = 1/(2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) * exp(-1/(2 * (1 - rho^2)) * (((x1 - mu1)/sigma1)^2 - 2 * rho * (x1 - mu1)/sigma1 * (x2 - mu2)/sigma2 +
                                                                                     ((x2 - mu2)/sigma2)^2))
  return(z)
}

### .letterPos .testFun###########################################
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


### .splitDev ##########################################
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
### FitDistr ##############
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
### .lfkp####################
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
### clase desirability.c#################
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

### .desireFun###############
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


### confounds####
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

### clase steepAscent.c##########
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

### code2real####
code2real = function(low, high, codedValue) {
  return((diff(c(low, high))/2) * codedValue + mean(c(low, high)))
}
### .nblock####
.nblock <- function(fdo) {
  if (class(fdo)[1] != "facDesign")
    stop(paste(deparse(substitute(fdo)), "needs to be an object of class 'facDesign'"))
  return(length(unique(fdo$block[[1]])))
}
### .starFrame####
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
### .alphaOrth, .alphaRot, .rsmOrth####
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
