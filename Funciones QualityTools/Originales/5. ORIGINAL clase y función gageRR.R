

setMethod("plot", signature(x = "gageRR"), function(x, y, main, xlab, ylab, col, lwd, fun = mean, ...) {
  horiz = FALSE
  parList = list(...)
  gdo = x
  yName = names(gdo)[5]
  aName = names(gdo)[3]
  bName = names(gdo)[4]
  abName = paste(aName, ":", bName, sep = "")
  if (missing(col))
    col = 2:(length(unique(gdo[, 3])) + 1)
  if (missing(lwd))
    lwd = 1
  par(mfrow = c(3, 2))
  temp = NULL
  Source = names(gdo@Varcomp)
  VarComp = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))]), 3)
  Contribution = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))])/as.numeric(gdo@Varcomp[length(gdo@Varcomp)]), 3)
  VarComp = t(data.frame(gdo@Varcomp))
  VarCompContrib = VarComp/gdo@Varcomp$totalVar
  Stdev = sqrt(VarComp)
  StudyVar = Stdev * gdo@Sigma
  StudyVarContrib = StudyVar/StudyVar["totalVar", ]
  if ((length(gdo@GageTolerance) > 0) && (gdo@GageTolerance > 0)) {
    ptRatio = StudyVar/gdo@GageTolerance
    temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib, ptRatio)
    contribFrame = data.frame(VarCompContrib, StudyVarContrib, ptRatio)
    names(temp)[6] = c("P/T Ratio")
    row.names(temp) = c(Source)
    SNR = sqrt(2 * (temp["bTob", "VarComp"]/temp["totalRR", "VarComp"]))
  }
  else {
    temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib)
    contribFrame = data.frame(VarCompContrib, StudyVarContrib)
  }
  bTob = paste(bName, "To", bName, sep = "")
  Source[Source == "bTob"] = bTob
  row.names(contribFrame) = Source
  if (gdo@method == "crossed")
    contribFrame = contribFrame[-match(c("totalVar", "a", "a_b"), row.names(temp)), ]
  else contribFrame = contribFrame[-match(c("totalVar"), row.names(temp)), ]


  numBars = ncol(contribFrame)
  ymax = max(max(contribFrame))
  # grafica 1
  main1 = NA
  if (missing(main) || is.na(main[1]))
    main1 = "Components of Variation"
  else main1 = main[1]
  xlab1 = NA
  if (missing(xlab) || is.na(xlab[1]))
    xlab1 = "component"
  else xlab1 = xlab[1]
  ylab1 = NA
  if (missing(ylab) || is.na(ylab[1]))
    ylab1 = ""
  else ylab1 = ylab[1]
  argList = list(...)
  redList = argList[names(argList) != "cex"]
  mybp = do.call(barplot, c(list(t(contribFrame), xlab = xlab1, ylab = ylab1, main = main1, names.arg = rep("", 4), axes = F, beside = T, ylim = c(0, 1.3 *
                                                                                                                                                     ymax), col = col[1:numBars]), redList))
  axis(1, at = colMeans(mybp), labels = names(as.data.frame(t(contribFrame))), ...)
  axis(2, ...)
  box()
  legend("topright", names(contribFrame), col = col[1:numBars], pch = c(15, 15), horiz = horiz, inset = 0.02)
  if (gdo@method == "crossed") {
    # 2. --------------------------------
    main2 = NA
    if (missing(main) || is.na(main[2]))
      main2 = paste(yName, "by", bName)
    else main2 = main[2]
    xlab2 = NA
    if (missing(xlab) || is.na(xlab[2]))
      xlab2 = bName
    else xlab2 = xlab[2]
    ylab2 = NA
    if (missing(ylab) || is.na(ylab[2]))
      ylab2 = yName
    else ylab2 = ylab[2]
    boxplot(split(gdo[, yName], gdo[, bName]), xlab = xlab2, ylab = ylab2, main = main2, ...)
    mByPa = split(gdo[, 5], as.numeric(gdo[, 4]))
    lines(sort(as.numeric(gdo[, 4])), lapply(mByPa, median)[sort(as.numeric(gdo[, 4]))], lwd = lwd)
    points(sort(as.numeric(gdo[, 4])), lapply(mByPa, median)[sort(as.numeric(gdo[, 4]))], lwd = lwd, pch = 13, cex = 2)

    # 3. -------------------------------
    main3 = NA
    if (missing(main) || is.na(main[3]))
      main3 = paste(yName, "by", aName)
    else main3 = main[3]
    xlab3 = NA
    if (missing(xlab) || is.na(xlab[3]))
      xlab3 = aName
    else xlab3 = xlab[3]
    ylab3 = NA
    if (missing(ylab) || is.na(ylab[3]))
      ylab3 = yName
    else ylab3 = ylab[3]
    colVec = .mapping(gdo[, 3], sort(unique(gdo[, 3])), col[1:length(unique(gdo[, 3]))])
    boxplot(split(gdo[, yName], gdo[, aName]), col = colVec, xlab = xlab3, ylab = ylab3, main = main3, ...)
    mByOp = split(gdo[, 5], as.numeric(gdo[, 3]))
    lines(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, mean)[sort(names(mByOp))], lwd = lwd)
    points(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, median)[sort(names(mByOp))], lwd = lwd, pch = 13, cex = 2)

    # 4. X_mean chart -------------
    agg = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = mean)
    tab = table(agg[, 2])
    sgSize = tab[1]
    aggSd = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = sd)
    tab = table(aggSd[, 2])
    sm = mean(aggSd[, 3])
    aggMean = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = mean)
    xm = mean(agg[, 3])
    UCL = xm + ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
    LCL = xm - ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
    values = c(UCL, xm, LCL)
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))
    plot(agg[, 3], type = "n", axes = FALSE, xlab = aName, ylab = expression(bar(x)), main = expression(paste(bar(x), " Chart")))
    box()
    abline(h = xm, col = 3)
    abline(h = UCL, col = 2)
    abline(h = LCL, col = 2)
    axis(2)
    axis(4, at = c(xm, UCL, LCL), labels = c("", "", ""))
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, xm, substitute(bar(x) == xm, list(xm = round(xm, 2))), adj = 0, srt = 0, xpd = TRUE)
    par(mar = old.par)
    j = 1
    for (i in 1:length(tab)) {
      lines(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
      points(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
      if (i < length(tab))
        abline(v = j + tab[i] - 1 + 0.5, lty = 2)
      axis(1, at = j, labels = names(tab[i]))
      j = j + tab[i]
    }
    # 5. ---------------------------------------------------------
    main4 = NA
    if (missing(main) || is.na(main[4]))
      main4 = paste("Interaction", abName)
    else main4 = main[4]
    xlab4 = NA
    if (missing(xlab) || is.na(xlab[4]))
      xlab4 = names(gdo)[4]
    else xlab4 = xlab[4]
    ylab4 = NA
    if (missing(ylab) || is.na(ylab[4]))
      ylab4 = paste(as.character(body(match.fun(fun)))[2], "of", names(gdo)[5])
    else ylab4 = ylab[4]
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))
    .aip(gdo[, 4], gdo[, 3], response = gdo[, 5], xlab = xlab4, ylab = ylab4, main = main4, col = col, type = "b", title = names(gdo)[3], ...)
    par(mar = old.par)


    # 6. R chart --------------------------------------------

    D3 = c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.284, 0.308, 0.329, 0.348)
    D4 = c(0, 3.267, 2.574, 2.282, 2.115, 2.004, 1.924, 1.864, 1.816, 1.777, 1.744, 1.716, 1.692, 1.671, 1.652)
    helpRange = function(x) {
      return(diff(range(x)))
    }
    aggForLimits = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = helpRange)
    Rm = mean(aggForLimits[, 3])
    UCL = D4[sgSize] * Rm
    LCL = D3[sgSize] * Rm
    agg = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = helpRange)
    tab = table(agg[, 2])
    sgSize = tab[1]
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))

    ###############
    plot(agg[, 3], ylim = c(0, max(max(agg[, 3]), UCL)), type = "n", xlab = aName, ylab = "R", axes = FALSE, main = "R Chart")
    axis(2)
    axis(4, at = c(Rm, UCL, LCL), labels = c("", "", ""))
    box()
    abline(h = Rm, col = 3)
    abline(h = UCL, col = 2)
    abline(h = LCL, col = 2)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, Rm, substitute(bar(R) == Rm, list(Rm = round(Rm, 2))), adj = 0, srt = 0, xpd = TRUE)
    j = 1
    for (i in 1:length(tab)) {
      lines(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
      points(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
      if (i < length(tab))
        abline(v = j + tab[i] - 1 + 0.5, lty = 2)
      axis(1, at = j, labels = names(tab[i]))
      j = j + tab[i]
    }
    par(mar = old.par)
  }
  if(gdo@method == "nested")
  {
    # 2.  **************************
    main2 = NA
    if (missing(main) || is.na(main[2]))
      main2 = paste(yName, "By", bName, "Within", aName)
    else main2 = main[2]
    xlab2 = NA
    if (missing(xlab) || is.na(xlab[2]))
      xlab2 = NA
    else xlab2 = xlab[2]
    ylab2 = NA
    if (missing(ylab) || is.na(ylab[2]))
      ylab2 = yName
    else ylab2 = ylab[2]
    agg = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = function(x) {
      return(x)
    })
    plot(1:nrow(agg), main = main2, xlab = xlab2, ylab = ylab2, ylim = range(agg[, 3]), axes = FALSE)
    axis(2)
    box()
    label2 = ""
    for (i in 1:nrow(agg)) {
      points(rep(i, length(agg[i, 3])), agg[i, 3])
      axis(1, at = i, labels = agg[i, 1])
      if (agg[i, 2] != label2) {
        axis(1, at = i, labels = agg[i, 2], line = 1, tick = FALSE)
        label2 = agg[i, 2]
      }
    }
    aggm = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = mean)
    lines(aggm[, 3])
    points(aggm[, 3], pch = 13, cex = 2)
    main3 = NA
    if (missing(main) || is.na(main[3]))
      main3 = paste(yName, "by", aName)
    else main3 = main[3]
    xlab3 = NA
    if (missing(xlab) || is.na(xlab[3]))
      xlab3 = aName
    else xlab3 = xlab[3]
    ylab3 = NA
    if (missing(ylab) || is.na(ylab[3]))
      ylab3 = yName
    else ylab3 = ylab[3]
    colVec = .mapping(gdo[, 3], sort(unique(gdo[, 3])), col[1:length(unique(gdo[, 3]))])
    boxplot(split(gdo[, yName], gdo[, aName]), col = colVec, xlab = xlab3, ylab = ylab3, main = main3, ...)
    mByOp = split(gdo[, 5], as.numeric(gdo[, 3]))
    lines(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, mean)[sort(names(mByOp))], lwd = lwd)
    points(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, mean)[sort(names(mByOp))], lwd = lwd, pch = 13, cex = 2)

    agg = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = mean)
    tab = table(agg[, 2])
    sgSize = tab[1]
    aggSd = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = sd)
    tab = table(aggSd[, 2])
    sm = mean(aggSd[, 3])
    aggMean = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = mean)
    xm = mean(agg[, 3])
    UCL = xm + ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
    LCL = xm - ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
    values = c(UCL, xm, LCL)
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))
    plot(agg[, 3], type = "n", axes = FALSE, xlab = aName, ylab = expression(bar(x)), main = expression(paste(bar(x), " Chart")))
    box()
    abline(h = xm, col = 3)
    abline(h = UCL, col = 2)
    abline(h = LCL, col = 2)
    axis(2)
    axis(4, at = c(xm, UCL, LCL), labels = c("", "", ""))
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, xm, substitute(bar(x) == xm, list(xm = round(xm, 2))), adj = 0, srt = 0, xpd = TRUE)
    par(mar = old.par)
    j = 1
    for (i in 1:length(tab)) {
      lines(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
      points(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
      if (i < length(tab))
        abline(v = j + tab[i] - 1 + 0.5, lty = 2)
      axis(1, at = j, labels = names(tab[i]))
      j = j + tab[i]
    }

    par(mar = old.par)
    D3 = c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.284, 0.308, 0.329, 0.348)
    D4 = c(0, 3.267, 2.574, 2.282, 2.115, 2.004, 1.924, 1.864, 1.816, 1.777, 1.744, 1.716, 1.692, 1.671, 1.652)
    helpRange = function(x) {
      return(diff(range(x)))
    }
    aggForLimits = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = helpRange)
    Rm = mean(aggForLimits[, 3])
    UCL = D4[sgSize] * Rm
    LCL = D3[sgSize] * Rm
    agg = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = helpRange)
    tab = table(agg[, 2])
    sgSize = tab[1]
    old.par = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 10.1))
    plot(agg[, 3], ylim = c(0, max(max(agg[, 3]), UCL)), type = "n", xlab = aName, ylab = "R", axes = FALSE, main = "R Chart")
    axis(2)
    axis(4, at = c(Rm, UCL, LCL), labels = c("", "", ""))
    box()
    abline(h = Rm, col = 3)
    abline(h = UCL, col = 2)
    abline(h = LCL, col = 2)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
    text(length(agg[, 3]) + length(agg[, 3]) * 0.075, Rm, substitute(bar(R) == Rm, list(Rm = round(Rm, 2))), adj = 0, srt = 0, xpd = TRUE)
    j = 1
    for (i in 1:length(tab)) {
      lines(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
      points(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
      if (i < length(tab))
        abline(v = j + tab[i] - 1 + 0.5, lty = 2)
      axis(1, at = j, labels = names(tab[i]))
      j = j + tab[i]

    }



    #        plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
    #        plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
    plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
  }
})
