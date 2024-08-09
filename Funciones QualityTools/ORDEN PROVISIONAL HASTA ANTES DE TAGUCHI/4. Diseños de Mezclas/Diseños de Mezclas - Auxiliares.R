#########################################################################
################ DISEÑOS DE MEZCLAS - AUXILIARES ########################
#########################################################################

# .npp ----
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
# .permut ----
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
# .simplexCentroid ----
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
# .mfc ----
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
