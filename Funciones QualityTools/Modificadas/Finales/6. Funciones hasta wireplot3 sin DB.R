library(R6)
library(ggplot2)
library(patchwork)
library(scales)
library(plotly)
library(gridExtra)
library(RColorBrewer)

############DISEÑOS FACTORIALES 2^k#############




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

### lm #########################################################################
m1 <- dfac$lm(rend ~ A*B*C)
summary(m1)


#########DISEÑOS FACTORIALES FRACCIONARIOS 2^k-p#################

####Uso diseños factoriales fraccionarios####
dfacfrac <- fracDesign(k=3,gen='C=AB',centerCube = 4)
dfacfrac$summary()
aliasTable(dfacfrac)
confounds(dfacfrac)

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

########DISEÑOS DE SUPERFICIE DE RESPUESTA#################
set.seed(1234)
fdo2 <- facDesign(k = 2, centerCube = 3)
fdo2$names(c("Factor1","Factor2"))
fdo2$lows(c(134,155))
fdo2$highs(c(155,175))
rend <- c(simProc(134,175),simProc(144.5,165.5),simProc(155,155),simProc(144.5,165.5),simProc(155,175),simProc(144.5,165.5),simProc(134,155))
fdo2$.response(rend)


####Visualizacion####
wirePlot(A, B, rend2, form = "rend2 ~ A*B + I(A^2) + I(B^2)", data = rsdo)
contourPlot(A, B, rend2, form = "rend2 ~ A*B + I(A^2) + I(B^2)", data = rsdo)
####filled.contour####
A = seq(40, 210, length = 100)
B = seq(90, 190, length = 100)
C = seq(90, 190, length = 100)
# filled.contour(A, B,outer(A,B, simProc, noise = FALSE), xlab = "Factor 1", ylab = "Factor 2", color.palette = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))
contourPlot(x = A, y = B, z = outer(A,B, simProc, noise = FALSE), xlab = "Factor 1", ylab = "Factor 2")


####poner en orden estandar####
fdo <- randomize(fdo,so = TRUE)
fdo$summary()



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


##EJEMPLO:####


# d<-desirability$new()
# desirability


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



################################################################################################33



###DISEÑOS DE MEZCLAS####
###Necesito .npp####
.npp = function(mdo) {
  pseudo = mdo$pseudo
  Type = mdo$Type
  temp = as.character(mdo$as.data.frame()$Type)
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
###Necesito .permut####
.permut <- function(x) {
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
###Necesito .simplexCentroid####
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
###Necesito clase mixDesign####
mixDesign.c <- R6Class("mixDesign", public = list(name = NULL,
                                                  factors =list(),
                                                  total = NULL,
                                                  lower = NULL,
                                                  design = data.frame(),
                                                  designType = NULL,
                                                  pseudo = data.frame(),
                                                  response = data.frame(),
                                                  Type = data.frame(),
                                                  block = data.frame(),
                                                  runOrder = data.frame(),
                                                  standardOrder = data.frame(),
                                                  desireVal = list(),
                                                  desirability = list(),
                                                  fits = data.frame(),

                                                  .factors = function(value){
                                                    if (missing(value)) {
                                                      return(self$factors)
                                                    }
                                                    else{
                                                      if (length(value) != ncol(self$pseudo))
                                                        stop("\nNumber of factors doesn't match with number of columns for factorial Design\n")
                                                      self$factors <- value
                                                      invisible(self)
                                                    }
                                                  },

                                                  names = function(value){
                                                    if(missing(value)){
                                                      aux <- list()
                                                      for (i in 1:length(self$factors)) {
                                                        aux[[.NAMES[i]]] <-self$factors[[i]]$name
                                                      }
                                                      return(aux)
                                                    }
                                                    else {
                                                      for (i in 1:length(self$factors)){
                                                        self$factors[[i]]$name = as.character(value[i])

                                                      }
                                                      invisible(self)
                                                    }

                                                  },

                                                  as.data.frame = function(row.names = NULL, optional = FALSE){
                                                    frameOut = cbind(self$standardOrder, self$runOrder, self$Type, self$pseudo, self$response)
                                                    return(frameOut)
                                                  },

                                                  print = function(){
                                                    print(format(self$as.data.frame(), digits = 4))
                                                    invisible(self$as.data.frame())
                                                  },

                                                  .response = function(value){
                                                    if (missing(value)) {
                                                      return(self$response)
                                                    }
                                                    else{
                                                      #print(deparse(substitute(value)))
                                                      if (!is.numeric(value) & !is.data.frame(value))
                                                        stop("vector or data.frame must be given")
                                                      if (is.numeric(value)) {
                                                        if (length(value) != nrow(self$pseudo))
                                                          stop("differing lengths")
                                                        temp = data.frame(value)
                                                        names(temp) = deparse(substitute(value))[1]
                                                        value = temp
                                                      }
                                                      if (is.data.frame(value)) {
                                                        if (nrow(value) != nrow(self$pseudo))
                                                          stop("differing number of rows")
                                                      }
                                                      self$response = value
                                                      invisible(self)
                                                    }

                                                  },

                                                  .nfp = function(){
                                                    x = self$.factors()
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

                                                  summary = function(){
                                                    cat(paste("Simplex", toupper(self$designType), "Design"))
                                                    cat("\n")
                                                    cat("Information about the factors:\n\n")
                                                    self$.nfp()
                                                    cat("\n-----------\n")
                                                    cat("\n")
                                                    .npp(self)
                                                    cat("\n-----------\n")
                                                    cat("\n")
                                                    cat("Information about the constraints:\n\n")
                                                    lower = object$lower
                                                    temp = character(0)
                                                    for (i in seq(along = lower)) temp = c(temp, paste(LETTERS[i], ">=", lower[i]))
                                                    cat(temp)
                                                    cat("\n")
                                                    cat("\n-----------\n")
                                                    cat("\n")
                                                    times = nrow(self$pseudo)
                                                    pseudo = format(self$pseudo, digits = 2)
                                                    design = format(self$design, digits = 2)
                                                    amount = design
                                                    if (self$total[2] != 1)
                                                      amount = format(self$design * self$total[2], digits = 2)
                                                    temp = c("                             ", "PseudoComponent", "_|_", "Proportion", "_|_", "Amount")
                                                    cat(temp)
                                                    cat("\n")
                                                    cat("\n")
                                                    temp = cbind(pseudo, `_` = rep(" ", times = times), `|` = rep("|", times = times), `_` = rep(" ", times = times), design)
                                                    temp = cbind(temp, `_` = rep(" ", times = times), `|` = rep("|", times = times), `_` = rep(" ", times = times), amount)
                                                    temp = cbind(self$standardOrder, self$runOrder, self$Type, `|` = rep("|", times = times), temp, `|` = rep("|", times = times), self$response)
                                                    show(temp)
                                                    cat("\n-----------\n")
                                                    cat("\n")
                                                    cat(paste("Mixture Total:", self$total[1], "equals", self$total[2]))
                                                    cat("\n")
                                                    cat("\n")
                                                    invisible(self$as.data.frame())
                                                  },

                                                  units = function(value){
                                                    if (missing(value)) {
                                                      v <- list()
                                                      for (i in 1:length(self$factors)) {
                                                        v[[unlist(self$names()[i])]] <- self$factors[[i]]$.unit()
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
                                                  }






)
)
###Funcion mixDesign####
mixDesign <- function(p, n = 3, type = "lattice", center = TRUE, axial = FALSE, delta, replicates = 1, lower, total = 1, randomize, seed) {

  frameOut = NA
  out = mixDesign.c$new()
  if (missing(p))
    stop("the number of factors p must be given")
  if (p <= 1 | !is.numeric(p))
    stop("invalid value for p")
  if (!(type %in% c("lattice", "centroid")))
    stop("type needs to be \"lattice\" or \"centroid\"")
  out$designType = type
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
  out$total = total
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
  out$lower = lower
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

  keepIndex = (1:nrow(frameOut))[!apply(Type, 1, "==", "center")]
  Type = data.frame(Type = Type[keepIndex, ])
  frameOut = frameOut[keepIndex, ]
  if (center) {
    center = data.frame(matrix(1/p, nrow = 1, ncol = p))

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

  }
  StandOrder = 1:nrow(frameOut)
  RunOrder = StandOrder
  if (randomize) {
    RunOrder = sample(1:nrow(frameOut), nrow(frameOut), replace = FALSE, prob = NULL)
  }
  frameOut = frameOut[order(RunOrder), ]
  row.names(frameOut) = frameOut$RunOrder
  out$pseudo = frameOut
  out$runOrder = data.frame(RunOrder = RunOrder)
  out$standardOrder = data.frame(StandOrder = StandOrder)
  out$Type = data.frame(Type = Type[order(RunOrder), ])
  out$response = data.frame(y = rep(NA, nrow(out$pseudo)))
  design = frameOut
  design[, ] = NA
  for (i in 1:ncol(frameOut)) {
    design[, i] = frameOut[, i] * (total[1] - sum(lower)) + lower[i]
  }
  out$design = design

  listFac = vector("list", p)
  for (i in seq(along = listFac)) listFac[[i]] = doeFactor$new()
  names(listFac) = LETTERS[1:p]
  out$.factors(listFac)
  if (out$total[2] != 1) {
    out$lows(lower * out$total[2])
    out$highs(1 * out$total[2])
    out$units("NA")
  }
  else if (any(out$lower != 0)) {
    out$lows(out$lower)
    out$highs(1 * out$total[1])
    outs$units("%")
  }
  else {
    out$lows(0)
    out$highs(1 * out$total[1])
    out$units("%")
  }
  return(out)
}
###uso mixDesign####
mdo <- mixDesign(3, 2, center = FALSE, axial = FALSE, randomize = FALSE, replicates = c(1, 1, 2, 3))
mdo$names(c("polyethylene", "polystyrene", "polypropylene"))
elongation <- c(11.0, 12.4, 15.0, 14.8, 16.1, 17.7, 16.4, 16.6, 8.8, 10.0, 10.0, 9.7, 11.8, 16.8, 16.0)
mdo$.response(elongation)
mdo$units()
###necesito .mfc####
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
###contourPlot3####
contourPlot3 = function(x, y, z, response, data = NULL, main, xlab, ylab, zlab, border, form = "linear", col = 1, col.text, cex.axis, axes = TRUE,
                        steps, factors) {
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
  nameVec = names(mdo$names())
  linStrings = "-1"
  for (i in seq(along = nameVec)) linStrings = paste(linStrings, "+", nameVec[i])

  combList = combn(nameVec, 2, simplify = FALSE)
  quadStrings = character(length = length(combList))
  for (i in seq(along = combList)){
    if (i == 1){
      quadStrings[i] = paste(combList[[i]][1], ":", combList[[i]][2])
    }
    else{
      quadStrings[i] = paste("+", combList[[i]][1], ":", combList[[i]][2])
    }
  }
  quadStrings = paste(quadStrings, collapse = "")

  if (identical(form, "linear")) {
    form = paste(r.c, "~", linStrings)

  }
  if (identical(form, "quadratic")) {
    form = paste(r.c, "~", linStrings, "+", quadStrings)

  }
  lm.1 = lm(formula = form, data = mdo$as.data.frame())

  dcList = vector(mode = "list", length = length(mdo$names()))
  names(dcList) = names(mdo$names())
  dcList[1:length(mdo$names())] = 0
  if (!is.null(factors)) {
    for (i in names(factors)) dcList[[i]] = factors[[i]][1]
  }

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
###uso contourPlot3####
contourPlot3(A, B, C, elongation, data = mdo, form = "quadratic")
###wirePlot3####
wirePlot3 = function(x, y, z, response, data = NULL, main, xlab, ylab, zlab, form = "linear", phi, theta, col = 1, steps, factors) {
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
  nameVec = names(mdo$names())
  linStrings = "-1"
  for (i in seq(along = nameVec)) linStrings = paste(linStrings, "+", nameVec[i])

  combList = combn(nameVec, 2, simplify = FALSE)
  quadStrings = character(length = length(combList))
  for (i in seq(along = combList)) if (i == 1)
    quadStrings[i] = paste(combList[[i]][1], ":", combList[[i]][2])
  else quadStrings[i] = paste("+", combList[[i]][1], ":", combList[[i]][2])
  quadStrings = paste(quadStrings, collapse = "")

  if (identical(form, "linear")) {
    form = paste(r.c, "~", linStrings)
  }
  if (identical(form, "quadratic")) {
    form = paste(r.c, "~", linStrings, "+", quadStrings)

  }
  lm.1 = lm(formula = form, data = mdo$as.data.frame())

  dcList = vector(mode = "list", length = length(mdo$names()))
  names(dcList) = names(mdo$names())
  dcList[1:length(mdo$names())] = 0
  if (!is.null(factors)) {
    for (i in names(factors)) dcList[[i]] = factors[[i]][1]
  }

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
###Uso wirePlot3####
wirePlot3(A, B, C, elongation, data = mdo, form = "quadratic", theta = -170)


