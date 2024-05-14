setClass(Class = "facDesign", representation = representation(name = "character", factors = "list",
                                                              cube = "data.frame", star = "data.frame", centerCube = "data.frame", centerStar = "data.frame",
                                                              generator = "ANY", response = "data.frame", block = "data.frame", blockGen = "data.frame", runOrder = "data.frame",
                                                              standardOrder = "data.frame", desireVal = "list", desirability = "list", fits = "list"))
setClass("doeFactor", representation = representation(low = "ANY", high = "ANY", name = "character",
                                                      unit = "character", type = "character"), prototype = prototype(low = -1, high = 1, name = "",
                                                                                                                     unit = "", type = "numeric"))



fracDesign = function(k = 3, p = 0, gen = NULL, replicates = 1, blocks = 1, centerCube = 0, random.seed = 1234) {    ###
  DB = FALSE
  STDfdo=FALSE
  if (p<0 || p>7)                                                                                                  ###
    stop("p needs to be an integer between 0 and 7!")                                                               ###
  if (abs(p - round(p)) > .Machine$double.eps^0.5)                                                                 ###
  {                                                                                                                ###
    warning(paste("p needs to be an integer but is real! p was rounded to", round(p)))                              ###
    p=round(p)                                                                                                      ###
  }                                                                                                                ###
  if (p != 0)                                                                                                      ###
  {                                                                                                                ###
    gen = NULL                                                                                                      ###
    for (i in 1:length(.fdoOrth))                                                                                   ###
    {                                                                                                              ###
      if (k==.fdoOrth[[i]]$k && p==.fdoOrth[[i]]$p)                                                                 ###
      {                                                                                                             ###
        STDfdo=TRUE                                                                                                 ###
        return(fracDesign(k=.fdoOrth[[i]]$k, gen=.fdoOrth[[i]]$gen,                                                 ###
                          replicates = replicates, blocks = blocks,                                                            ###
                          centerCube = centerCube, random.seed = random.seed))                                                 ###
      }                                                                                                             ###
    }                                                                                                              ###
    if(STDfdo==FALSE)                                                                                              ###
      stop("No standard Design for the choosen combination of k and p (see: fracChoose())!")                        ###
  }                                                                                                                ###
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
  for (j in 1:k) X[, j] <- rep(sort(rep(c(-1, 1), N/2^j)), 2^(j - 1))
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
          if ((chars[j] == "=") & (length(numVec) != 1))
            stop("check position of \"=\" in generators!")
          if (chars[j] != "=") {
            charVec = c(charVec, toupper(chars[j]))
            numVec = c(numVec, numCharVec[names(numCharVec) == toupper(chars[j])])
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
        stop(paste("generator:", paste(ind[1], "=", paste(ind[-1], collapse = "*")), "includes undefined columns"))
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
  DesignOut = new("facDesign")
  DesignOut@generator = gen
  cube(DesignOut) = frameOut
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
  numRows = nrow(cube(DesignOut)) + nrow(star(DesignOut)) + nrow(centerStar(DesignOut)) + nrow(centerCube(DesignOut))
  if (DB) {
    print(numRows)
    print("response")
  }
  DesignOut@response = data.frame(y = rep(NA, numRows))
  if (DB)
    print("response")
  standardOrder = data.frame(matrix(data = 1:numRows, nrow = numRows, ncol = 1))
  names(standardOrder) = "StandOrder"
  standardOrder
  standOrd(DesignOut) = standardOrder
  if (DB)
    print("1")
  set.seed(random.seed)
  runOrder = as.data.frame(standardOrder[sample(1:numRows), ])
  if (DB)
    print("2")
  names(runOrder) = "RunOrder"
  runOrd(DesignOut) = runOrder
  if (DB)
    print("3")
  if (centerCube >= 1) {
    temp = data.frame(matrix(rep(0, centerCube * k), ncol = k, nrow = centerCube))
    names(temp) = names(frameOut)
    centerCube(DesignOut) = temp
  }
  temp = try(blocking(DesignOut, blocks = blocks))
  if (inherits(temp, "try-error"))
    stop("Blocking not possible!")      #return(DesignOut)                  ###
  return(blocking(DesignOut, blocks = blocks))
}
facDesign = function(k = 3, p = 0, replicates = 1, blocks = 1, centerCube = 0) {  ###
  frameOut = fracDesign(k = k, p = p, gen = NULL, replicates = replicates, blocks = blocks, centerCube = centerCube)  ###
  return(frameOut)
}

