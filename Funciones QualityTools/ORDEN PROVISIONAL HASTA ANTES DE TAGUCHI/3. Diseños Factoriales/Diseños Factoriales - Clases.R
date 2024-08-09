##########################################################################
###################### DISEÑOS FACTORIALES - CLASES ######################
##########################################################################

# Clase facDesign.c ----
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

                                                        names(newNums) = names(oldRunOrd)
                                                        newRunOrd = data.frame(oldRunOrd[1:lenFirst, ])

                                                        names(newRunOrd) = names(oldRunOrd)
                                                        restFrame = data.frame(oldRunOrd[-c(1:lenFirst), ])
                                                        names(restFrame) = names(oldRunOrd)
                                                        newRunOrd = rbind(newRunOrd, newNums, restFrame)

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

                                                      oldBlockGen = self$blockGen
                                                      if (ncol(oldBlockGen) > 0) {

                                                        newBlockGen = data.frame(oldBlockGen[1:lenFirst, ])
                                                        names(newBlockGen) = names(self$blockGen)
                                                        naFrameGen = as.data.frame(matrix(rep(NA, times = ncol(self$blockGen) * nrow(newDf)), ncol = ncol(self$blockGen)))
                                                        names(naFrameGen) = names(oldBlockGen)
                                                        restBlockGen = data.frame(oldBlockGen[-c(1:(lenFirst + nrow(oldDf))), ])
                                                        names(restBlockGen) = names(oldBlockGen)
                                                        newBlockGen = rbind(newBlockGen, naFrameGen, restBlockGen)

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

                                                      numNewRow = nrow(newDf) - nrow(oldDf)
                                                      oldOrd = self$standardOrder
                                                      len = nrow(oldOrd)
                                                      self$standardOrder = data.frame(StandOrd = 1:(len + numNewRow))
                                                      newRunOrd = data.frame()
                                                      if (numNewRow > 0) {
                                                        newNums = data.frame(newNums = seq(max(oldRunOrd) + 1, max(oldRunOrd) + numNewRow, by = 1))
                                                        names(newNums) = names(oldRunOrd)

                                                        newRunOrd = data.frame(oldRunOrd[1:lenCube, ])
                                                        names(newRunOrd) = names(oldRunOrd)
                                                        restRunOrd = data.frame(oldRunOrd[-c(1:lenCube), ])
                                                        names(restRunOrd) = names(oldRunOrd)
                                                        newRunOrd = rbind(newRunOrd, newNums, restRunOrd)

                                                        self$runOrder = newRunOrd
                                                      }
                                                      else {
                                                        newRunOrd = data.frame(oldRunOrd[1:(lenCube + nrow(newDf)), ])
                                                        names(newRunOrd) = names(oldRunOrd)
                                                        restRunOrd = data.frame(oldRunOrd[-c(1:(lenCube + nrow(oldDf))), ])
                                                        names(restRunOrd) = names(oldRunOrd)
                                                        newRunOrd = rbind(newRunOrd, restRunOrd)

                                                        self$runOrder = newRunOrd
                                                      }
                                                      naFrame = as.data.frame(matrix(rep(NA, times = ncol(oldResponse) * nrow(newDf)), ncol = ncol(oldResponse)))
                                                      names(naFrame) = names(oldResponse)
                                                      newResponse = data.frame(oldResponse[1:lenCube, ])
                                                      names(newResponse) = names(oldResponse)
                                                      restResponse = data.frame(oldResponse[-c(1:(lenCube + nrow(oldDf))), ])
                                                      names(restResponse) = names(oldResponse)
                                                      newResponse = rbind(newResponse, naFrame, restResponse)

                                                      self$.response(newResponse)
                                                      oldBlockGen = self$blockGen
                                                      if (ncol(oldBlockGen) > 0) {

                                                        newBlockGen = data.frame(oldBlockGen[1:lenCube, ])
                                                        names(newBlockGen) = names(self$blockGen)
                                                        naFrameGen = as.data.frame(matrix(rep(NA, times = ncol(self$blockGen) * nrow(newDf)), ncol = ncol(self$blockGen)))
                                                        names(naFrameGen) = names(oldBlockGen)
                                                        restFrame = data.frame(oldBlockGen[-c(1:(lenCube + nrow(oldDf))), ])
                                                        names(restFrame) = names(self$blockGen)
                                                        newBlockGen = rbind(newBlockGen, naFrameGen, restFrame)
                                                        self$.blockGen(newBlockGen)

                                                      }
                                                      oldBlock = self$block
                                                      newBlock = data.frame(oldBlock[1:lenCube, ])
                                                      names(newBlock) = names(self$block)
                                                      naFrame = as.data.frame(matrix(rep(blockValues, times = nrow(newDf)/numBlocks), ncol = 1))
                                                      restFrame = as.data.frame(oldBlock[-c(1:(lenCube + nrow(oldDf))), ])
                                                      names(restFrame) = names(self$block)

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

                                                      self$.response(newResponse)

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

# Clase doeFactor ----
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




# Clase desirability.c ----
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
                                                        plot = function(y, scale, main, xlab, ylab, line.width, col, numPoints = 500, ...){
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
                                                          if (missing(line.width))
                                                            line.width = 0.75
                                                          if (missing(scale))
                                                            scale = self$scale
                                                          if (missing(col))
                                                            col = "red"
                                                          dFun = .desireFun(self$low, self$high, self$target, self$scale, self$importance)
                                                          xVals = seq(self$low - 0.04 * diff(range(self$low, self$high)), self$high + 0.04 * diff(range(self$low, self$high)), length = numPoints)
                                                          yVals = dFun(xVals)

                                                          df <- data.frame(X = xVals, Y = yVals)

                                                          p <- ggplot(df, aes(x = X, y = Y)) +
                                                            geom_line(color = col, size = line.width) +
                                                            labs(title = main, x = xlab, y = ylab) +
                                                            theme_minimal() +
                                                            theme(plot.title = element_text(hjust = 0.5))
                                                          # annotate("text", x = (mean(range(df$X)) * 1.05), y = mean(range(df$Y)), label = "scale = 1", size = 5) # Añadir la etiqueta en el centro

                                                          if (is.numeric(self$target)) {
                                                            xm1 <- mean(c(min(xVals), self$target))
                                                            xm2 <- mean(c(max(xVals), self$target))
                                                            ym1 <- df$Y[which.min(abs(df$X - xm1))]
                                                            ym2 <- df$Y[which.min(abs(df$X - xm2))]
                                                            p <- p +
                                                              annotate("text", x = xm1*1.05, y = ym1, label = paste("scale =", scale[1])) +
                                                              annotate("text", x = xm2*0.95, y = ym2, label = paste("scale =", scale[2]))
                                                          } else {
                                                            xm1 <- mean(range(xVals))
                                                            ym1 <- df$Y[which.min(abs(df$X - xm1))]
                                                            if (identical(self$target, "max")) {
                                                              p <- p +
                                                                annotate("text", x = xm1*1.05, y = ym1, label = paste("scale =", scale[1]))
                                                            } else {
                                                              p <- p +
                                                                annotate("text", x = xm1*1.05, y = ym1, label = paste("scale =", scale[1]))
                                                            }
                                                          }
                                                          print(p)
                                                          out = list(x = xVals, y = yVals)
                                                          names(out) = c(self$response, "desirability")
                                                          invisible(out)
                                                        }
)
)

# Clase steepAscent.c ----
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
                                                    plot = function(y, main, xlab, ylab, l.col, p.col,
                                                                    line.type, point.shape,...){
                                                      Delta = (self$X)$Delta
                                                      frame = cbind(Delta, self$.response())
                                                      names(frame) = c("Delta", names(self$.response()))
                                                      if (missing(main))
                                                        main = ""
                                                      if (missing(xlab))
                                                        xlab = "Delta"
                                                      if (missing(ylab))
                                                        ylab = "predicted"
                                                      if (missing(l.col))
                                                        l.col = "red"
                                                      if (missing(p.col))
                                                        p.col = "red"
                                                      if (missing(line.type))
                                                        line.type = "dashed"
                                                      if (missing(line.type))
                                                        line.type = "dashed"
                                                      if (missing(point.shape))
                                                        point.shape = 16

                                                      p <- ggplot(frame, aes(x = Delta, y = predicted)) +
                                                        geom_line(color = l.col, linetype = line.type) +
                                                        geom_point(color = p.col, shape = point.shape) +
                                                        labs(title = main, x = xlab, y = ylab) +
                                                        theme_minimal()
                                                      print(p)
                                                    }


)
)

# Clase desOpt ----
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

