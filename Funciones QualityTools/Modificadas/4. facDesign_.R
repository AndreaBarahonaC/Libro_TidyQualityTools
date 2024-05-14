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
        index = fdo$names() %in% temp[, j]
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



##Necesito clase facDesign######################################################
facDesign_c <- R6Class("facDesign", public = list(name = NULL,
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
                                                 }
                                                 )
                      )



##Necesito clase doeFactor####################
doeFactor_ <- R6Class('doeFactor', public = list(low = -1,
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
  DesignOut <- facDesign_c$new()
  DesignOut$generator <-  gen
  DesignOut$cube <-  frameOut
  listFac <-  vector("list", ncol(frameOut))
  for (i in seq(along = listFac)){
    listFac[[i]] = doeFactor_$new()
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
  numRows = nrow(DesignOut$cube) + nrow(DesignOut$centerCube) #+ nrow(DesignOut$star) +
    #nrow(DesignOut$centerStar)  ####################################################
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
  runOrder = as.data.frame(standardOrder[sample(1:numRows),
    ])
  if (DB)
    print("2")
  names(runOrder) = "RunOrder"
  DesignOut$runOrder <- runOrder
  if (DB)
    print("3")
  temp = try(blocking(DesignOut, blocks = blocks))
  if (inherits(temp, "try-error"))
    stop("Blocking not possible!")
  return(blocking(DesignOut, blocks = blocks))
}

############## funcion facDesign########################
facDesign_ <- function (k = 3, p = 0, replicates = 1, blocks = 1, centerCube = 0)
{
  frameOut = fracDesign_(k = k, p = p, gen = NULL, replicates = replicates,
                        blocks = blocks, centerCube = centerCube)
  return(frameOut)
}


###############USO DE facDesign#######################################################
set.seed(1234)
dfac <- facDesign_(k = 3, centerCube = 4)
dfac$names()
dfac$names(c('Factor 1', 'Factor 2', 'Factor 3'))
dfac$names()
dfac$lows(c(80,120,1))
dfac$lows()
dfac$highs(c(120,140,2))
dfac$highs()
dfac$summary()
