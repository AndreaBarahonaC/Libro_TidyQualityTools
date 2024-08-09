#################################################################################################
############################# DISEÑOS FACTORIALES - FUNCIONES ###################################
#################################################################################################

### CREAN LAS CLASES

# Función fracDesign ----
fracDesign <- function (k = 3, p = 0, gen = NULL, replicates = 1, blocks = 1,
                        centerCube = 0, random.seed = 1234)
{

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

  if (!is.null(gen)) {
    listGen = vector("list", length(gen))
    .numFactors = numeric(0)
    charFactors = character(0)

    temp = character(0)
    for (i in seq(along = gen)) {

      if (!is.character(gen[i]))
        stop("Defining Relations should contain characters only!")
      chars = strsplit(gen[i], split = character(0))[[1]]

      checkDupl = character(0)
      for (j in 1:length(chars)) {
        if (chars[j] %in% toupper(c(.NAMES[1:26], letters[1:26]))) {
          if (chars[j] %in% checkDupl)
            stop("Defining relations contain one or more duplicates!")
          checkDupl = c(checkDupl, chars[j])
          temp = c(temp, chars[j])
        }
      }
    }
    temp = sort(unique(temp))
    numCharVec = 1:length(temp)
    names(numCharVec) = temp

    for (i in seq(along = gen)) {

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

      listGen[[i]] = numVec
      .numFactors = c(.numFactors, numVec)
      charFactors = c(charFactors, charVec)
    }

    names(.numFactors) = charFactors
    if (length(unique(.numFactors)) > k)
      stop("number of distinct Factors in generators greater than k!")

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

  if (centerCube >= 1) {
    temp = data.frame(matrix(rep(0, centerCube * k), ncol = k,
                             nrow = centerCube))
    names(temp) = names(frameOut)
    DesignOut$centerCube = temp
  }
  numRows = nrow(DesignOut$cube) + nrow(DesignOut$centerCube) + nrow(DesignOut$star) +
    nrow(DesignOut$centerStar)

  DesignOut$response = data.frame(y = rep(NA, numRows))

  standardOrder = data.frame(matrix(data = 1:numRows, nrow = numRows,
                                    ncol = 1))
  names(standardOrder) = "StandOrder"
  DesignOut$standardOrder <-  standardOrder

  set.seed(random.seed)
  runOrder = as.data.frame(standardOrder[sample(1:numRows),])

  names(runOrder) = "RunOrder"
  DesignOut$runOrder <- runOrder

  temp = try(blocking(DesignOut, blocks = blocks)) #################Poner random.seed para variar con semilla
  if (inherits(temp, "try-error"))
    stop("Blocking not possible!")
  return(blocking(DesignOut, blocks = blocks))  #################Poner random.seed para variar con semilla
}

# Función facDesign ----
facDesign <- function (k = 3, p = 0, replicates = 1, blocks = 1, centerCube = 0, random.seed = 1234)
{
  frameOut = fracDesign(k = k, p = p, gen = NULL, replicates = replicates,
                        blocks = blocks, centerCube = centerCube, random.seed = random.seed)
  return(frameOut)
}


# USO DE facDesign
dfac <- facDesign(k = 3, centerCube = 4)
#dfac$names()
dfac$names(c('Factor 1', 'Factor 2', 'Factor 3'))
#dfac$names()
dfac$lows(c(80,120,1))
#dfac$lows()
dfac$highs(c(120,140,2))
#dfac$highs()
dfac$summary()
# effectPlot
dfac$effectPlot(classic=TRUE)
dfac$effectPlot()


# Función steepAscent ----
steepAscent <- function(factors, response, size = 0.2, steps = 5, data) {
  if (missing(data))
    return("missing an object of class 'facDesign'")
  else fdo = data
  if (missing(factors) | length(factors) < 1)
    return("missing factors")
  if (!is.character(factors))
    return("factors needs to be a character")
  #names(names(fdo))
  model = data$fits[[response]]
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
  b = numeric(length = length(factors))
  x = numeric(length = length(factors))
  names(x) = factors
  for (i in seq(along = factors)) {
    b[i] = coef(model)[factors[i]]
    if (i == 1) {
      x[i] = size * sign(b[i])
    }
    else {
      x[i] = (x[1]/b[1]) * b[i]
    }
  }

  Run = 1:(steps + 1)
  Delta = 0:steps
  frameOut = data.frame(Run, Delta)
  initial = ncol(frameOut)
  for (i in seq(along = factors)) {
    frameOut[, i + initial] = x[i] * 0:steps
    names(frameOut)[i + initial] = paste(factors[i], ".coded", collapse = "", sep = "")
  }
  initial = ncol(frameOut)
  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[.NAMES[i]]] <-fdo$names()[i]
  }
  for (i in seq(along = factors)) {
    frameOut[, i + initial] = code2real(fdo$lows()[[aux[i][[1]]]], fdo$highs()[[aux[i][[1]]]], x[i] * 0:steps)
    names(frameOut)[i + initial] = paste(factors[i], ".real", collapse = "", sep = "")
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
# USO steepAscent
sao <- steepAscent(factors = c("A", "B"), response = "rend", data = dfac, steps = 20)
sao$print()
predicted <- simProc(sao$get(j = 5), sao$get(j = 6))
sao$.response(predicted)
sao$plot()





# Función starDesign ----
starDesign <- function(k, p = 0, alpha = c("both", "rotatable", "orthogonal"), cs, cc, data) {
  fdo = NULL
  csFrame = NULL
  ccFrame = NULL
  starFrame = NULL
  blocks = 1
  alpha = alpha[1]
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
  if (!missing(data))
    fdo$.star(starFrame)
  if (cs > 0) {
    csFrame = as.data.frame(matrix(0, nrow = cs, ncol = k))
    names(csFrame) = nameVec
    if (!missing(data)) {
      fdo$.centerStar(csFrame)
    }
  }
  if (cc > 0) {
    ccFrame = as.data.frame(matrix(0, nrow = cc, ncol = k))
    names(ccFrame) = nameVec
    if (!missing(data)) {
      fdo$.centerCube(ccFrame)
    }
  }
  if (!missing(data))
    return(fdo)
  else return(list(star = starFrame, centerStar = csFrame, centerCube = ccFrame))
}

# Uso starDesign
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


# Función rsmDesign ----
rsmDesign <- function(k = 3, p = 0, alpha = "rotatable", blocks = 1, cc = 1, cs = 1, fp = 1,
                      sp = 1, faceCentered = FALSE) {
  if (blocks > 2^(k - 1) + 1)
    stop("Blocking not possible")
  if (alpha == "rotatable")
    alpha = .alphaRot(k, p)
  if (alpha == "orthogonal")
    alpha = .alphaOrth(k, p, cc = cc, cs = cs)
  if (alpha == "both") {
    found = FALSE
    for (i in seq(along = .rsmOrth)) {
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
  fdo = facDesign(k = k, p = p, replicates = fp)                              ###
  if (cc > 0) {
    temp = as.data.frame(matrix(0, nrow = cc, ncol = ncol(fdo$cube)))
    names(temp) = names(fdo$cube)
    fdo$centerCube(temp)
  }

  temp = .starFrame(k, alpha)
  starportion = data.frame()
  for (i in 1:sp) {
    starportion = rbind(temp, starportion)
  }
  names(starportion) = names(fdo$cube)
  fdo$.star(starportion)
  if (cs > 0) {
    temp = as.data.frame(matrix(0, nrow = cs, ncol = ncol(fdo$cube)))
    names(temp) = names(fdo$cube)
    fdo$.centerStar(temp)
  }
  #    return(fdo)
  fdo = blocking(fdo, blocks)
  return(fdo)
}
# Uso rsmDesign
fdo <- rsmDesign(k=3, alpha=1.633, cc=0, cs=6)

# Funcion desirability ----
desirability = function(response, low, high, target = "max", scale = c(1, 1), importance = 1, constraints) {
  if (low >= high)
    stop("the lower bound must be greater than the high bound!")
  if (any(scale <= 0))
    stop("the scale parameter must be greater than zero!")
  if (!is.numeric(target) & !identical(tolower(target), "min") & !identical(tolower(target), "max"))
    stop("target needs to be \"min\", \"max\" or a numeric value")
  return(desirability.c$new(response = deparse(substitute(response)), low = low, high = high, target = target, scale = scale, importance = importance))
}

# Uso
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
d1$plot()
d3$plot()


### OTRAS ----
# randomize ----
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
# as.data.frame.facDesign ----
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
# aliasTable ----
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

# blocking ----
blocking <- function (fdo, blocks, BoR = FALSE, random.seed, useTable = "rsm",
                      gen){
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

# simProc ----
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

# USO simProc
#Primeros valores
rend <- simProc(x1=120,x2=140,x3=2)
#valores completos
rend <- c(simProc(120,140,1),simProc(80,140,1),simProc(120,140,2),simProc(120,120,1),simProc(90,130,1.5),simProc(90,130,1.5),simProc(80,120,2),simProc(90,130,1.5),simProc(90,130,1.5),simProc(120,120,2),simProc(80,140,2),simProc(80,120,1))

#Asignar rendimiento al diseño factorial
dfac$.response(rend)
dfac$.response()


# confounds ----
confounds <- function(x, depth = 2) {

  varName = deparse(substitute(x))
  identityList = x$identity()
  x = x$cube
  if (length(identityList) < 1) {
    print(paste(varName, " contains no defining relations!"))
    invisible()
  }
  effect1 = numeric(0)
  effect2 = numeric(0)

  index = numeric(0)
  for (i in 1:(dim(x)[2])) {
    if (!(TRUE && all(unique(x[, i]) %in% c(-1, 1))))
      index = c(index, i)
  }
  if (length(index) > 0)
    x = x[, -index]

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

      }
    }
  }

  if (length(effect1) > 0)
    dupIndex = numeric(0)
  for (i in 1:length(effect1)) {
    if (i > length(effect1))
      break
    index = (1:length(effect1))[effect2 == effect1[i]]
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
# code2real ----
code2real = function(low, high, codedValue) {
  return((diff(c(low, high))/2) * codedValue + mean(c(low, high)))
}
# overall ----
overall <- function(fdo, steps = 20, constraints, ...) {
  importances = list()
  cs = list()
  if (!missing(constraints))
    cs = constraints
  l = vector(mode = "list", length = 0)
  fitList = fdo$fits
  if (length(fitList) < 1)
    stop(paste("no fits found in fits(", deparse(substitute(fdo)), ")"), sep = "")
  desList = fdo$desires()
  if (length(desList) < 1)
    stop(paste("no desirabilities found in desires(", deparse(substitute(fdo)), ")"), sep = "")
  X = fdo$cube
  newdata = NULL
  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[fdo$names()[i]]] <-.NAMES[i]
  }
  aux<-unlist(aux)
  for (i in aux) {
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

  newdata = expand.grid(l)
  names(newdata) = names(X)
  out = newdata
  yCharSet = intersect(names(fdo$desires()), names(fdo$fits))
  dFrame = data.frame(matrix(NA, nrow = nrow(newdata), ncol = length(yCharSet) + 1))
  names(dFrame) = c(yCharSet, "overall")
  dFrame[, "overall"] = 1
  for (y in yCharSet) {
    obj = desList[[y]]
    dFun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
    lm.y = fitList[[y]]
    importances[[y]] = fdo$desires()[[y]]$importance
    yHat = predict(lm.y, newdata = newdata, ...)
    yDes = dFun(yHat)
    dFrame[, y] = yDes
  }
  geomFac = 1/sum(unlist(importances))
  overall = apply(dFrame, 1, prod)^geomFac
  dFrame[, "overall"] = overall
  dFrame = cbind(out, dFrame)
  invisible(dFrame)
}

# optimum ----
optimum <- function(fdo, constraints, steps = 25, type = "grid", start, ...) {

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
    aux <- list()
    for (i in 1:length(fdo$names())) {
      aux[[fdo$names()[i]]] <-.NAMES[i]
    }
    aux<-unlist(aux)
    names(factors)<-unlist(aux)
    desOpt$facCoded = as.list(dVals[index, aux])
    names(desOpt$facCoded)<-fdo$names()
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
# Uso optimum
optimum(ddo,type='optim')



