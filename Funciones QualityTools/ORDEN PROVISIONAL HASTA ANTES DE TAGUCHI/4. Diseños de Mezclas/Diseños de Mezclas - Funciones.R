##########################################################################
################ DISEÑOS DE MEZCLAS - FUNCIONES ##########################
##########################################################################

###Funcion mixDesign####
mixDesign <- function(p, n = 3, type = "lattice", center = TRUE, axial = FALSE, delta, replicates = 1, lower, total = 1, randomize, seed) {
  DB = FALSE
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
# Uso mixDesign
mdo <- mixDesign(3, 2, center = FALSE, axial = FALSE, randomize = FALSE, replicates = c(1, 1, 2, 3))

mdo$names(c("polyethylene", "polystyrene", "polypropylene"))
elongation <- c(11.0, 12.4, 15.0, 14.8, 16.1, 17.7, 16.4, 16.6, 8.8, 10.0, 10.0, 9.7, 11.8, 16.8, 16.0)
mdo$.response(elongation)

mdo$units()
mdo$summary()
