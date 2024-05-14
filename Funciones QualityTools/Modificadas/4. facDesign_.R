library(R6)
library(MASS)

###R/COM_P.R#######################################
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

##Necesito clase facDesign######################################################
facDesign_ <- R6Class("facDesign", public = list(name = NULL,
                                                 factors = NULL,
                                                 cube = NULL,
                                                 star = NULL,
                                                 centerCube = NULL,
                                                 centerStar = NULL,
                                                 generator = NULL,
                                                 response = NULL,
                                                 block = NULL,
                                                 blockGen = NULL,
                                                 runOrder = NULL,
                                                 standardOrder = NULL,
                                                 desireVal = NULL,
                                                 desirability = NULL,
                                                 fits = NULL,

                                                 cube = function(x){
                                                   out <- self$cube
                                                   out
                                                 }
                                                 )
                      )

###Necesito funcion "fracDesign"

fracDesign_ <- function (k = 3, p = 0, gen = NULL, replicates = 1, blocks = 1,
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
  if (!is.numeric(replicates))
    stop("replicates needs to be numeric!")
  else if (replicates < 0)
    stop("replicates need to >= 0")
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
  DesignOut <-  facDesign_$new()
  DesignOut$generator <-  gen
  DesignOut$cube = frameOut
  listFac = vector("list", ncol(frameOut))
  for (i in seq(along = listFac)) listFac[i] = new("doeFactor")
  names(listFac) = names(frameOut)
  factors(DesignOut) = listFac
  if (DB)
    print(frameOut)
  if (DB)
    print("yes")
  if (DB)
    print("aha")
  numRows = nrow(cube(DesignOut)) + nrow(star(DesignOut)) +
    nrow(centerStar(DesignOut)) + nrow(centerCube(DesignOut))
  if (DB) {
    print(numRows)
    print("response")
  }
  DesignOut@response = data.frame(y = rep(NA, numRows))
  if (DB)
    print("response")
  standardOrder = data.frame(matrix(data = 1:numRows, nrow = numRows,
    ncol = 1))
  names(standardOrder) = "StandOrder"
  standardOrder
  standOrd(DesignOut) = standardOrder
  if (DB)
    print("1")
  set.seed(random.seed)
  runOrder = as.data.frame(standardOrder[sample(1:numRows),
    ])
  if (DB)
    print("2")
  names(runOrder) = "RunOrder"
  runOrd(DesignOut) = runOrder
  if (DB)
    print("3")
  if (centerCube >= 1) {
    temp = data.frame(matrix(rep(0, centerCube * k), ncol = k,
      nrow = centerCube))
    names(temp) = names(frameOut)
    centerCube(DesignOut) = temp
  }
  temp = try(blocking(DesignOut, blocks = blocks))
  if (inherits(temp, "try-error"))
    stop("Blocking not possible!")
  return(blocking(DesignOut, blocks = blocks))
}


facDesign_ <- function (k = 3, p = 0, replicates = 1, blocks = 1, centerCube = 0)
{
  frameOut = fracDesign(k = k, p = p, gen = NULL, replicates = replicates,
                        blocks = blocks, centerCube = centerCube)
  return(frameOut)
}



