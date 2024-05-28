
### funcion fracDesign###################

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

### funcion facDesign########################
facDesign <- function (k = 3, p = 0, replicates = 1, blocks = 1, centerCube = 0, random.seed = 1234)
{
#' @title facDesign
#' @description Generates a 2^k full factorial design.
#' @param k numeric value giving the number of factors. By default k is set to ‘3’.
#' @param p numeric integer between ‘0’ and ‘7’. p is giving the number of additional factors in the response surface design by aliasing effects.
#' For further information see fracDesign and fracChoose.
#' By default p is set to ‘0’.
#' @param replicates numeric value giving the number of replicates per factor combination. By default replicates is set to ‘1’.
#' @param blocks numeric value giving the number of blocks. By default blocks is set to ‘1’. Blocking is only performed for k greater 2.
#' @param centerCube numeric value giving the number of centerpoints within the 2^k design. By default centerCube is set to ‘0’.
#' @return The function facDesign returns an object of class facDesign.

  frameOut = fracDesign(k = k, p = p, gen = NULL, replicates = replicates,
                        blocks = blocks, centerCube = centerCube, random.seed = random.seed)
  return(frameOut)
}


### funcion simProc#####################################
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



### funcion InteractionPlot#########################################################
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



### funcion paretoPlot#################################
paretoPlot <- function(fdo, threeWay = FALSE, abs = TRUE,
                        decreasing = TRUE, na.last = NA, alpha = 0.05,
                        response = NULL, ylim, xlab, ylab, main, single = TRUE, p.col, ...) {  ###
  if(single==FALSE)                                                           ###
    par(mfrow=.splitDev(length(fdo$.response()))[[2]])                           ###
  if(is.null(response)==FALSE)                                                ###
  {                                                                           ###
    temp=fdo$.response()[response]                                               ###
    fdo$.response(temp)                                                         ###
  }                                                                           ###
  ylimMissing = FALSE
  if (missing(ylim)){
    ylimMissing = TRUE
  }
  if (missing(xlab))
    xlab = ""
  location = "topright"
  if (decreasing == F) {
    location = "topleft"
  }
  xVals = numeric(0)
  sig.neg = NULL
  sig.pos = NULL
  effect.list = vector("list")
  suppressMessages(
    for (j in 1:ncol(fdo$.response())) {
      if (j > 1 && single==TRUE) {
        dev.new()
        par(mar = c(5.1, 4.1, 4.1, 4.1))
      }
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
          else{p.col = brewer.pal(length(effect), paste0("Set", p.col))}

          if (ylimMissing)
            if (abs)
              ylim = (range(c(0, abs(effect), 1.3 * ME))) * 1.1
          else ylim = (range(c(effect, -1.3 * ME, 1.3 * ME))) * 1.1
          if (abs) {
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
          effect = ((summary(lm.1)$coefficients[-pmatch("(Intercept)", names(coef(lm.1))), 1])/(summary(lm.1)$coefficients[-pmatch("(Intercept)", names(coef(lm.1))),
                                                                                                                           2]))
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
              ylim = c(0, 1.3 * max(tempVec))
            }
          else {
            tempVec1 = c(0, effect, sig.neg, sig.pos)
            tempVec1 = tempVec1[!is.na(tempVec1)]
            tempVec2 = c(abs(effect), sig.pos, sig.neg)
            tempVec2 = tempVec2[!is.na(tempVec2)]
            ylim = c(1.3 * min(tempVec1), 1.3 * max(tempVec2))
          }

          effect = effect[order(abs(effect), na.last = TRUE, decreasing = decreasing)]
          effect = round(effect, 3)

          if(missing(p.col)){
            p.col = rep("lightblue", length(effect))
          }
          else{
            p.col = brewer.pal(length(effect), paste0("Set", p.col))}
          # Plot ---------
          if (abs) {
            if (missing(ylabel))
              ylabel = ""
            # plot with abs
            p <- ggplot(data.frame(names = names(effect), effect_ = abs(as.vector(effect))),
                        aes(x = names, y = effect_, fill = names)) +
              geom_bar(stat = "identity", color = "black") +
              scale_fill_manual(values=c(p.col)) +
              theme_minimal()+
              labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
              theme(axis.text.x = element_text(angle = 90, hjust = 1),
                    plot.title = element_text(hjust = 0.5)) +
              geom_text(aes(label = round(effect,2)), vjust = -1, colour = "black") + # etiquetas sobre las barras
              geom_hline(yintercept = sig.pos, linetype = "dashed", color = "red") +
              scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(sig.pos), labels = c(round(sig.pos, 2))))
          }
          else {
            if (missing(ylabel))
              ylabel = ""
            # Plot without abs
            p <- ggplot(data.frame(names = names(effect), effect_ = abs(as.vector(effect))),
                        aes(x = names, y = effect_, fill = names)) +
              geom_bar(stat = "identity", color = "black") +
              scale_fill_manual(values=c(p.col)) +
              theme_minimal()+
              labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
              theme(axis.text.x = element_text(angle = 90, hjust = 1),
                    plot.title = element_text(hjust = 0.5)) +
              geom_text(aes(label = round(effect)), vjust = -1, colour = "black") + # etiquetas sobre las barras
              geom_hline(yintercept = sig.pos, linetype = "dashed", color = "red") +
              geom_hline(yintercept = sig.neg, linetype = "dashed", color = "red") +
              scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(sig.pos, sig.neg), labels = c(round(sig.pos, 2), round(sig.neg, 2))))
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
  )

  print(p)
  invisible(list(effect.list, plot = p))
  par(mfcol=c(1,1))
}
### funcion normalPlot#####################
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
### funcion wirePlot###################
wirePlot <- function(x, y, z, data = NULL, xlim, ylim, zlim, main, xlab, ylab, border, sub, zlab, form = "fit", phi, theta, ticktype, col = 1, steps,
                     factors, fun, plot) {
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
  }
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

  if (missing(plot))
    plot = TRUE
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

  if (missing(ticktype))
    ticktype = "detailed"
  if (missing(border))
    border = NULL
  if (missing(phi))
    phi = 30
  if (missing(theta))
    theta = -30
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

  if (!is.null(factors)) {
    for (i in names(factors)) dcList[[i]] = factors[[i]][1]
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

  if (is.function(col)) {
    nrMat <- nrow(mat)
    ncMat <- ncol(mat)
    jet.colors <- colorRampPalette(c("blue", "green"))
    nbcol <- 100
    color <- col(nbcol)
    matFacet <- mat[-1, -1] + mat[-1, -ncMat] + mat[-nrMat, -1] + mat[-nrMat, -ncMat]
    facetcol <- cut(matFacet, nbcol)
  }else {
    color = col
    facetcol = 1
  }

  if (missing(zlim))
    zlim = range(mat)

  p <- plot_ly(x = xVec, y = yVec, z = mat, colors = color) %>%
    add_surface() %>%
    layout(
      title = main,
      annotations = list(
        list(
          text = sub,# Subtitulo
          x = 0.5,   # Posición x en la mitad de la gráfica
          y = -0.1,  # Posición y debajo de la gráfica
          printarrow = FALSE,
          font = list(size = 12)
        )
      ),
      scene = list(
        xaxis = list(range = xlim, title = xlab, zeroline = FALSE),
        yaxis = list(range = ylim, title = ylab, zeroline = FALSE),
        zaxis = list(range = zlim, title = zlab, zeroline = FALSE),
        camera = list(eye = list(x=2, y=2, z=0.1))
      )
    )

  if (plot) {
    show(p)
  }
  invisible(list(x = xVec, y = yVec, z = mat, plot = p))
}

### funcion contourPlot#####################
contourPlot = function(x, y, z, data = NULL, xlim, ylim, main, xlab, ylab, border, sub, zlab, form = "fit", phi, theta, ticktype, col = 1, steps,
                       factors, fun, plot = TRUE) {
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

  p <- plot_ly(x = xVec, y = yVec, z = mat, colors = color, type = "contour", contours = list(coloring = 'heatmap')) %>%
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
### funcion fracChoose####
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

### funcion steepAscent############
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

### starDesign####
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

### rsmDesign####
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
### rsmChoose()####
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

### desirability#############
desirability = function(response, low, high, target = "max", scale = c(1, 1), importance = 1, constraints) {
  if (low >= high)
    stop("the lower bound must be greater than the high bound!")
  if (any(scale <= 0))
    stop("the scale parameter must be greater than zero!")
  if (!is.numeric(target) & !identical(tolower(target), "min") & !identical(tolower(target), "max"))
    stop("target needs to be \"min\", \"max\" or a numeric value")
  return(desirability.c$new(response = deparse(substitute(response)), low = low, high = high, target = target, scale = scale, importance = importance))
}
