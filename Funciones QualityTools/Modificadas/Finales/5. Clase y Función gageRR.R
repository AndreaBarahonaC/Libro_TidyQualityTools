library(R6)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)
library(patchwork)

# Funciones auxiliares --------------
.aip <- function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = TRUE, trace.label = NULL,
                 fixed = FALSE, xlab = deparse(substitute(x.factor)), ylab ="Measurement" , ylim = NULL, lty = 1:length(unique(trace.factor)),
                 col = 1, pch = c(1L:9, 0, letters), xpd = NULL, leg.bg = par("bg"), leg.bty = "o", xtick = FALSE, xaxt = par("xaxt"), axes = TRUE, title = "", plot = TRUE, ...) {
  ylabel <- paste(ylab )
  type <- match.arg(type)

  # Asegurarse de que los factores son realmente factores
  x.factor <- factor(x.factor)
  trace.factor <- factor(trace.factor)

  # Calcular los valores de celda
  cellNew <- tapply(response, list(x.factor, trace.factor), fun)
  cellNew <- as.data.frame(as.table(cellNew))
  colnames(cellNew) <- c("x.factor", "trace.factor", "response")

  # Convertir x.factor a numérico sólo para la visualización, pero mantenerlo como factor en los datos
  cellNew$x.numeric <- as.numeric(cellNew$x.factor)

  # Crear el gráfico
  p <- ggplot(cellNew, aes(x = x.numeric, y = response, group = trace.factor, color = trace.factor, shape = trace.factor, linetype = trace.factor)) +
    geom_line() +
    geom_point(size = 2) +
    # geom_text(aes(label = round(response, 2)), vjust = -0.5) +  # Comentar o eliminar esta línea para ocultar números
    scale_x_continuous(breaks = unique(cellNew$x.numeric), labels = levels(cellNew$x.factor)) +  # Etiquetas de x.factor
    labs(x = xlab, y = ylabel, title = title, color = trace.label, shape = trace.label, linetype = trace.label) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  if (!is.null(ylim)) {
    p <- p + ylim(ylim)
  }

  if (legend) {
    p <- p + theme(legend.position="right",
                   legend.background = element_rect(fill = "white", colour = "black"))
  } else {
    p <- p + theme(legend.position = "none")
  }
  if(plot){
    print(p)
  }
  invisible(list(plot = p))
}

.c4 = function(n) {
  if (n > 1 && n < 342)
    sqrt(2/(n - 1)) * (factorial(n/2 - 1)/factorial((n - 1)/2 - 1))
  else stop("n needs to be bigger than 1 and smaller than 342")
}

.mapping = function(x, oldValues, newValues) {
  if (length(oldValues) != length(newValues)) {
    print("old and new")
    print(oldValues)
    print(newValues)
    warning(paste("unequal length of", deparse(substitute(oldValues)), "and", deparse(substitute(newValues))))
  }
  out = numeric(length(x))
  for (i in seq(along = newValues)) {
    index = (x == oldValues[i])
    out[index] = newValues[i]
  }
  return(out)
}

# Definición de la clase gageRR ------------------
gageRR <- R6Class("gageRR",
                  public = list(
                    X = NULL,
                    ANOVA = NULL,
                    RedANOVA = NULL,
                    method = NULL,
                    Estimates = NULL,
                    Varcomp = NULL,
                    Sigma = NULL,
                    GageName = NULL,
                    GageTolerance = NULL,
                    DateOfStudy = NULL,
                    PersonResponsible = NULL,
                    Comments = NULL,
                    b = NULL,
                    a = NULL,
                    y = NULL,
                    facNames = NULL,
                    numO = NULL,
                    numP = NULL,
                    numM = NULL,
                    initialize = function(X, ANOVA = NULL, RedANOVA = NULL, method = NULL, Estimates = NULL, Varcomp = NULL,
                                          Sigma = NULL, GageName = NULL, GageTolerance = NULL, DateOfStudy = NULL,
                                          PersonResponsible = NULL, Comments = NULL, b = NULL, a = NULL, y = NULL,
                                          facNames = NULL, numO = NULL, numP = NULL, numM = NULL) {
                      self$X <- X
                      self$ANOVA <- ANOVA
                      self$RedANOVA <- RedANOVA
                      self$method <- method
                      self$Estimates <- Estimates
                      self$Varcomp <- Varcomp
                      self$Sigma <- Sigma
                      self$GageName <- GageName
                      self$GageTolerance <- GageTolerance
                      self$DateOfStudy <- DateOfStudy
                      self$PersonResponsible <- PersonResponsible
                      self$Comments <- Comments
                      self$b <- b
                      self$a <- a
                      self$y <- y
                      self$facNames <- facNames
                      self$numO <- numO
                      self$numP <- numP
                      self$numM <- numM
                    },
                    show = function() {
                      print(as.data.frame(self$X))
                    },
                    subset = function(i, j) {
                      return(self$X[i, j])
                    },
                    summary = function() {
                      if (all(is.na(self$X$Measurement))) {
                        cat("Gage R&R Summary\n")
                        cat("-----------------\n")
                        cat("Method: ", self$method, "\n")
                        cat("Sigma: ", self$Sigma, "\n")
                        cat("Gage Name: ", self$GageName, "\n")
                        cat("Gage Tolerance: ", self$GageTolerance, "\n")
                        cat("Date of Study: ", self$DateOfStudy, "\n")
                        cat("Person Responsible: ", self$PersonResponsible, "\n")
                        cat("Comments: ", self$Comments, "\n")
                        cat("Operators: ", self$numO, "\n")
                        cat("Parts: ", self$numP, "\n")
                        cat("Measurements per Part: ", self$numM, "\n")
                      } else {
                        cat("\n")
                        cat("Operators:\t", self$numO, "\tParts:\t", self$numP, "\n")
                        cat("Measurements:\t", self$numM, "\tTotal:\t", nrow(self$X), "\n")
                        cat("----------\n")
                      }
                      return(invisible(self))
                    },
                    response = function() {
                      return(self$X$Measurement)
                    },
                    replace_response = function(value) {
                      self$X$Measurement = value
                      return(self)
                    },
                    names = function() {
                      return(names(as.data.frame(self)))
                    },
                    as_data_frame = function() {
                      return(as.data.frame(self))
                    },
                    get_tolerance = function() {
                      return(unlist(self$GageTolerance))
                    },
                    set_tolerance = function(value) {
                      if (!is.numeric(value))
                        stop("GageTolerance needs to be numeric")
                      self$GageTolerance = value
                      return(self)
                    },
                    get_sigma = function() {
                      return(unlist(self$Sigma))
                    },
                    set_sigma = function(value) {
                      if (!is.numeric(value))
                        stop("Sigma needs to be numeric")
                      self$Sigma = value
                      return(self)
                    },
                    plot = function(x, y, main=NULL, xlab=NULL, ylab=NULL, col, lwd, fun = mean, ...){
                      x = self
                      gdo <- x
                      yName <- x$facNames[1]
                      aName <- x$facNames[2]
                      bName <- x$facNames[3]
                      abName <- paste(aName, ":", bName, sep = "")
                      if (missing(col))
                        col <- 2:(length(unique(gdo$b)) + 1)
                      if (missing(lwd))
                        lwd <- 1
                      temp <- NULL
                      Source <- names(gdo$Varcomp)
                      VarComp <- round(as.numeric(gdo$Varcomp), 3)
                      Contribution <- round(as.numeric(gdo$Varcomp) / as.numeric(gdo$Varcomp[length(gdo$Varcomp)]), 3)
                      VarComp <- t(data.frame(gdo$Varcomp))
                      VarCompContrib <- VarComp / gdo$Varcomp$totalVar
                      Stdev <- sqrt(VarComp)
                      StudyVar <- Stdev * gdo$Sigma
                      StudyVarContrib <- StudyVar / StudyVar["totalVar", ]
                      if ((length(gdo$GageTolerance) > 0) && (gdo$GageTolerance > 0)) {
                        ptRatio <- StudyVar / gdo$GageTolerance
                        temp <- data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib, ptRatio)
                        contribFrame <- data.frame(VarCompContrib, StudyVarContrib, ptRatio)
                        names(temp)[6] <- c("P/T Ratio")
                        row.names(temp) <- c(Source)
                        SNR <- sqrt(2 * (temp["bTob", "VarComp"] / temp["totalRR", "VarComp"]))
                      } else {
                        temp <- data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib)
                        contribFrame <- data.frame(VarCompContrib, StudyVarContrib)
                      }
                      bTob <- paste(bName, "To", bName, sep = "")
                      Source[Source == "bTob"] <- bTob
                      row.names(contribFrame) <- Source
                      if (gdo$method == "crossed")
                        contribFrame <- contribFrame[-match(c("totalVar", "a", "a_b"), row.names(temp)), ]
                      else contribFrame <- contribFrame[-match(c("totalVar"), row.names(temp)), ]
                      numBars <-ncol(contribFrame)
                      
                      contrib_df <- as.data.frame(contribFrame)
                      contrib_df$component <- rownames(contribFrame)
                      contrib_df <- contrib_df %>% rownames_to_column(var = "Source") %>% filter(Source != "totalVar")
                      ymax <- max(max(contribFrame))
                      main1 <- NA
                      
                      contribFrame_long <- as.data.frame(contribFrame)
                      contribFrame_long$Component <- rownames(contribFrame_long)
                      contribFrame_long <- tidyr::gather(contribFrame_long, key = "Metric", value = "Value", -Component)
                      # 1. Components of Variation --------------------------------------------------
                      p1 <- ggplot(contribFrame_long, aes(x = Component, y = Value, fill = Metric)) +
                        geom_bar(stat = "identity", position = "dodge") +
                        labs(title = ifelse(is.null(main[1]), "Components of Variation", main[1]),
                             x = ifelse(is.null(xlab[1]), "Component", xlab[1]),
                             y = ifelse(is.null(ylab[1]), "", ylab[1])) +
                        theme_bw() +
                        scale_fill_manual(values = col[1:nlevels(factor(contribFrame_long$Metric))]) +
                        theme(legend.position = "top",
                              legend.background = element_rect(fill = "white", colour = "black"),
                              legend.key.size = unit(0.5, "cm"), # reduce the size of the legend key
                              legend.text = element_text(size = 8), # reduce the size of the legend text
                              legend.box.spacing = unit(0.1, 'cm'),
                              plot.title = element_text(hjust = 0.5)
                        )
                      
                      # ----------------------
                      if (gdo$method == "crossed") {
                        # 2. Measurement by part ------------------------------------------------------
                        main2 <- NA
                        if (missing(main) || is.na(main[2]))
                          main2 <- paste(yName, "by", bName)
                        else main2 <- main[2]
                        xlab2 <- NA
                        if (missing(xlab) || is.na(xlab[2]))
                          xlab2 <- bName
                        else xlab2 <- xlab[2]
                        ylab2 <- NA
                        if (missing(ylab) || is.na(ylab[2]))
                          ylab2 <- yName
                        else ylab2 <- ylab[2]
                        
                        p2 <- ggplot(gdo$X, aes_string(x = bName, y = yName)) +
                          geom_boxplot() +
                          stat_summary(fun = median, geom = "line", aes(group = 1), color = "red", linewidth = lwd) +
                          stat_summary(fun = median, geom = "point", color = "red") +
                          labs(title = ifelse(is.null(main2), paste(yName, "by", bName), main2),
                               x = ifelse(is.null(xlab2), bName, xlab2),
                               y = ifelse(is.null(ylab2), yName, ylab2)) +
                          theme_bw() +
                          theme(plot.title = element_text(hjust = 0.5))
                        
                        # 3. Measurement by operator --------------------------------------------------
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
                        
                        p3 <- ggplot(gdo$X, aes_string(x = aName, y = yName)) +
                          geom_boxplot(aes(fill = factor(gdo$X[, 3]))) +
                          stat_summary(fun = median, geom = "line", aes(group = 1), color = "red", linewidth = lwd) +
                          stat_summary(fun = median, geom = "point", color = "red", size = 3) +
                          labs(title = ifelse(is.null(main[3]), paste(yName, "by", aName), main[3]),
                               x = ifelse(is.null(xlab[3]), aName, xlab[3]),
                               y = ifelse(is.null(ylab[3]), yName, ylab[3]),
                               fill = "Factor") +
                          theme_bw() +
                          scale_fill_manual(values = col) +
                          theme(plot.title = element_text(hjust = 0.5), legend.position='none')
                        
                        # 4. X_mean Chart -------------------------------------------------------------------------------
                        agg = aggregate(gdo$X[, yName], list(gdo$X[, aName], gdo$X[, bName]), FUN = mean)
                        tab = table(agg[, 2])
                        sgSize = tab[1]
                        aggSd = aggregate(gdo$X[, yName], list(gdo$X[, bName], gdo$X[, aName]), FUN = sd)
                        tab = table(aggSd[, 2])
                        sm = mean(aggSd[, 3])
                        aggMean = aggregate(gdo$X[, yName], list(gdo$X[, bName], gdo$X[, aName]), FUN = mean)
                        xm = mean(agg[, 3])
                        UCL = xm + ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
                        LCL = xm - ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
                        values = c(UCL, xm, LCL)
                        
                        p4 <- ggplot(data.frame(x = 1:length(aggMean[, 3]), y = aggMean[, 3]), aes(x, y)) +
                          geom_point(shape = 1) +
                          geom_line() +
                          geom_hline(yintercept = xm, color = 3) +                          # línea xm
                          geom_hline(yintercept = UCL, color = "red") +                     # línea UCL
                          geom_hline(yintercept = LCL, color = "red") +                     # línea LCL
                          labs(x = aName, y = expression(bar(x))) +                         # títulos ejes
                          ggtitle(expression(paste(bar(x), " Chart"))) +                    # título
                          geom_vline(xintercept = cumsum(tab) - 0.5, linetype = "dashed") + # divide categorías A,B C
                          scale_x_continuous(breaks = cumsum(tab), labels = names(tab)) +   # y pone los nombres de las secciones
                          theme_minimal() +
                          theme(plot.title = element_text(hjust = 0.5)) +
                          scale_y_continuous(expand = c(0, 0),
                                             sec.axis = sec_axis(~ ., breaks = c(UCL, xm, LCL),
                                                                 labels= c(paste("UCL =", round(UCL, 2)),
                                                                           paste("X_mean =", round(xm,2)),
                                                                           paste("UCL =", round(UCL, 2))
                                                                 ) )) +
                          theme(axis.text.y.right = element_text(size = 8))
                        
                        # 5. Interaction Operator  ------------------------------------------------------------------
                        main4 <- NA
                        if (missing(main) || is.na(main[4]))
                          main4 <- paste("Interaction", abName)
                        else main4 <- main[4]
                        xlab4 <- NA
                        if (missing(xlab) || is.na(xlab[4]))
                          xlab4 <- colnames(gdo$X)[4]
                        else xlab4 <- xlab[4]
                        ylab4 <- NA
                        if (missing(ylab) || is.na(ylab[4]))
                          ylab4 <- paste(as.character(body(match.fun(fun)))[2], "of", colnames(gdo$X)[5])
                        else ylab4 <- ylab[4]
                        p5 <- .aip(gdo$X[, 4], gdo$X[, 3], response = gdo$X[, 5], xlab = xlab4, ylab = ylab4, title = "Interaction Operator: Part", legend = TRUE,col = col, type = "b", plot = FALSE)
                        
                        # 6. R chart -------------------------------------------------
                        D3 <- c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.284, 0.308, 0.329, 0.348)
                        D4 <- c(0, 3.267, 2.574, 2.282, 2.115, 2.004, 1.924, 1.864, 1.816, 1.777, 1.744, 1.716, 1.692, 1.671, 1.652)
                        helpRange = function(x) {
                          return(diff(range(x)))
                        }
                        aggForLimits <- aggregate(gdo$X[, yName], list(gdo$X[, aName], gdo$X[, bName]), FUN = helpRange)
                        Rm <- mean(aggForLimits[, 3])
                        sgSize <- length(unique(gdo$X[, bName]))
                        UCL <- D4[sgSize] * Rm
                        LCL <- D3[sgSize] * Rm
                        agg <- aggregate(gdo$X[, yName], list(gdo$X[, bName], gdo$X[, aName]), FUN = helpRange)
                        agg$Group.1 <- factor(agg$Group.1, levels = unique(agg$Group.1))
                        tab = table(agg[, 2])
                        sgSize = tab[1]
                        
                        p6 <-  ggplot(data.frame(x = 1:length(agg[, 3]), y = agg[, 3]), aes(x, y)) +
                          geom_point(shape = 1) +
                          geom_line() +
                          geom_hline(yintercept = Rm, color = 3) +                          # línea xm
                          geom_hline(yintercept = UCL, color = "red") +                     # línea UCL
                          geom_hline(yintercept = LCL, color = "red") +                     # línea LCL
                          labs(x = aName, y = "R") +                                        # títulos ejes
                          ggtitle("R Chart") +                                              # título
                          geom_vline(xintercept = cumsum(tab) - 0.5, linetype = "dashed") + # divide categorías A,B C
                          scale_x_continuous(breaks = cumsum(tab), labels = names(tab)) +   # y pone los nombres de las secciones
                          theme_bw() +
                          theme(plot.title = element_text(hjust = 0.5)) +
                          scale_y_continuous(expand = c(0, 0),
                                             sec.axis = sec_axis(~ ., breaks = c(UCL, Rm, LCL),
                                                                 labels= c(paste("UCL =", round(UCL, 2)),
                                                                           paste("R_mean =", round(Rm,2)),
                                                                           paste("UCL =", round(UCL, 2))
                                                                 ) )) +
                          theme(axis.text.y.right = element_text(size = 8))
                        
                        # UNION PLOTS -------------
                        p <- p1 + p3 + p5$plot + p2 + p4 + p6 + plot_layout(nrow = 3, byrow = FALSE)
                        
                      }
                      
                      if(gdo$method == "nested"){
                        # 2. Measurement by Part within operator --------------
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
                        agg = aggregate(gdo$X[, yName], list(gdo$X[, bName], gdo$X[, aName]), FUN = function(x){return(x)})
                        
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
                        aggm = aggregate(gdo$X[, yName], list(gdo$X[, bName], gdo$X[, aName]), FUN = mean)
                        lines(aggm[, 3])
                        points(aggm[, 3], pch = 13, cex = 2)
                        
                        # 3. Box Blot Measurement by Operator -----------------
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
                        
                        p3 <- ggplot(gdo$X, aes_string(x = aName, y = yName)) +
                          geom_boxplot(aes(fill = factor(gdo$X[, 3]))) +
                          stat_summary(fun = median, geom = "line", aes(group = 1), color = "red", linewidth = lwd) +
                          stat_summary(fun = median, geom = "point", color = "red", size = 3) +
                          labs(title = ifelse(is.null(main[3]), paste(yName, "by", aName), main[3]),
                               x = ifelse(is.null(xlab[3]), aName, xlab[3]),
                               y = ifelse(is.null(ylab[3]), yName, ylab[3]),
                               fill = "Factor") +
                          theme_bw() +
                          scale_fill_manual(values = col) +
                          theme(plot.title = element_text(hjust = 0.5), legend.position='none')
                        
                        
                        
                        
                        # 4. X_mean Chart ----------------------
                        agg = aggregate(gdo$X[, yName], list(gdo$X[, aName], gdo$X[, bName]), FUN = mean)
                        tab = table(agg[, 2])
                        sgSize = tab[1]
                        aggSd = aggregate(gdo$X[, yName], list(gdo$X[, bName], gdo$X[, aName]), FUN = sd)
                        tab = table(aggSd[, 2])
                        sm = mean(aggSd[, 3])
                        aggMean = aggregate(gdo$X[, yName], list(gdo$X[, bName], gdo$X[, aName]), FUN = mean)
                        xm = mean(agg[, 3])
                        UCL = xm + ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
                        LCL = xm - ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
                        values = c(UCL, xm, LCL)
                        
                        p4 <- ggplot(data.frame(x = 1:length(aggMean[, 3]), y = aggMean[, 3]), aes(x, y)) +
                          geom_point(shape = 1) +
                          geom_line() +
                          geom_hline(yintercept = xm, color = 3) +                          # línea xm
                          geom_hline(yintercept = UCL, color = "red") +                     # línea UCL
                          geom_hline(yintercept = LCL, color = "red") +                     # línea LCL
                          labs(x = aName, y = expression(bar(x))) +                         # títulos ejes
                          ggtitle(expression(paste(bar(x), " Chart"))) +                    # título
                          geom_vline(xintercept = cumsum(tab) - 0.5, linetype = "dashed") + # divide categorías A,B C
                          scale_x_continuous(breaks = cumsum(tab), labels = names(tab)) +   # y pone los nombres de las secciones
                          theme_minimal() +
                          theme(plot.title = element_text(hjust = 0.5)) +
                          scale_y_continuous(expand = c(0, 0),
                                             sec.axis = sec_axis(~ ., breaks = c(UCL, xm, LCL),
                                                                 labels= c(paste("UCL =", round(UCL, 2)),
                                                                           paste("X_mean =", round(xm,2)),
                                                                           paste("UCL =", round(UCL, 2))
                                                                 ) )) +
                          theme(axis.text.y.right = element_text(size = 8))
                        # 5. R Chart --------------------------------------
                        D3 = c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.284, 0.308, 0.329, 0.348)
                        D4 = c(0, 3.267, 2.574, 2.282, 2.115, 2.004, 1.924, 1.864, 1.816, 1.777, 1.744, 1.716, 1.692, 1.671, 1.652)
                        helpRange = function(x) {
                          return(diff(range(x)))
                        }
                        aggForLimits = aggregate(gdo$X[, yName], list(gdo$X[, aName], gdo$X[, bName]), FUN = helpRange)
                        Rm = mean(aggForLimits[, 3])
                        UCL = D4[sgSize] * Rm
                        LCL = D3[sgSize] * Rm
                        agg = aggregate(gdo$X[, yName], list(gdo$X[, bName], gdo$X[, aName]), FUN = helpRange)
                        tab = table(agg[, 2])
                        sgSize = tab[1]
                        
                        p5 <-  ggplot(data.frame(x = 1:length(agg[, 3]), y = agg[, 3]), aes(x, y)) +
                          geom_point(shape = 1) +
                          geom_line() +
                          geom_hline(yintercept = Rm, color = 3) +                          # línea xm
                          geom_hline(yintercept = UCL, color = "red") +                     # línea UCL
                          geom_hline(yintercept = LCL, color = "red") +                     # línea LCL
                          labs(x = aName, y = "R") +                                        # títulos ejes
                          ggtitle("R Chart") +                                              # título
                          geom_vline(xintercept = cumsum(tab) - 0.5, linetype = "dashed") + # divide categorías A,B C
                          scale_x_continuous(breaks = cumsum(tab), labels = names(tab)) +   # y pone los nombres de las secciones
                          theme_minimal() +
                          theme(plot.title = element_text(hjust = 0.5)) +
                          scale_y_continuous(expand = c(0, 0),
                                             sec.axis = sec_axis(~ ., breaks = c(UCL, Rm, LCL),
                                                                 labels= c(paste("UCL =", round(UCL, 2)),
                                                                           paste("R_mean =", round(Rm,2)),
                                                                           paste("UCL =", round(UCL, 2))
                                                                 ) )) +
                          theme(axis.text.y.right = element_text(size = 8))
                        
                        # UNION PLOTS -----------------------
                        p <- p1 / (p2 + p3) / (p4 + p5) # p2 falta cambiar a ggplot2
                      }
                      
                      show(p)
                      invisible(list(plot = p))
                    }
                  )
)

# función gageRRdesign modificada -------------
gageRRDesign = function(Operators = 3, Parts = 10, Measurements = 3, method = "crossed", sigma = 6, randomize = TRUE) {
  # Validación de argumentos
  if (!is.numeric(sigma))
    stop("sigma needs to be numeric")
  if (method != "nested" && method != "crossed")
    stop("Unknown method specified. Use 'method = nested' or 'method = crossed'.")
  Measurements <- as.integer(Measurements)
  if (!is.numeric(Measurements) || Measurements <= 0)
    stop("Number of Measurements per Part must be a positive integer.")


  opvec <- factor()
  partvec <- factor()

  yName <- "Measurement"
  aName <- "Operator"
  bName <- "Part"
  abName <- "Operator:Part"

  Operators <- unique(Operators)
  Parts <- unique(Parts)

  if (is.numeric(Operators))
    opvec <- factor(LETTERS[1:Operators[1]])
  if (is.character(Operators))
    opvec <- factor(Operators)

  if (length(unique(opvec)) > 26)
    stop("Too many Operators!")
  if (length(unique(opvec)) < 2)
    stop("Not enough Operators")

  if (is.numeric(Parts))
    partvec <- factor(LETTERS[1:Parts[1]])
  if (is.character(Parts))
    partvec <- factor(Parts)

  if (length(unique(partvec)) > 26)
    stop("Too many Parts!")
  if (length(unique(partvec)) < 2)
    stop("Too few Parts")

  Measurement <- rep(NA, (length(opvec) * length(partvec) * Measurements))
  outFrame <- data.frame()

  if (method == "crossed") {
    temp <- expand.grid(opvec, partvec)
    o <- rep(temp[, 1], Measurements)
    p <- rep(temp[, 2], Measurements)
  } else {
    p <- rep(sort(rep(partvec, length(opvec))), Measurements)
    o <- (rep(opvec, length(Measurement) / length(opvec)))
    p <- p[order(o,p)]
    o <- o[order(o,p)]
  }

  if (randomize)
    outFrame <- data.frame(StandardOrder = 1:length(Measurement), RunOrder = sample(1:length(Measurement), length(Measurement)), Operator = factor(o), Part = factor(p), Measurement)
  else
    outFrame <- data.frame(StandardOrder = 1:length(Measurement), RunOrder = 1:length(Measurement), Operator = factor(o), Part = factor(p), Measurement)

  outFrame <- outFrame[order(outFrame$RunOrder), ]
  # Valores predeterminados
  gageRRObj <- gageRR$new(
    X = outFrame,
    ANOVA = NULL,
    RedANOVA = NULL,
    method = method,
    Estimates = NULL,
    Varcomp = NULL,
    Sigma = sigma,
    GageName = NULL,
    GageTolerance = NULL,
    DateOfStudy = Sys.Date(),
    PersonResponsible = NULL,
    Comments = NULL,
    b = factor(p),
    a = factor(o),
    y = as.numeric(Measurement),
    facNames = c(yName, aName, bName, abName),
    numO = length(unique(opvec)),  # Número de operadores
    numP = length(unique(partvec)),  # Número de partes
    numM = Measurements  # Número de mediciones
  )

  return(gageRRObj)
}

gageRR_ = function(gdo, method = "crossed", sigma = 6, alpha = 0.25, DM = NULL, HM = NULL, tolerance = NULL, dig = 3, ...) {
  method <- method

  yName <- "Measurement"
  aName <- "Operator"
  bName <- "Part"

  abName <- if(method == "crossed") paste(aName, ":", bName, sep = "")
  else if(method == "nested") paste(bName, "(", aName, ")", sep = "")
  else NA

  bTobName <- paste(bName, "to", bName, sep = " ")

  if (is.null(tolerance)) tolerance <- gdo$get_tolerance()

  y <- gdo$X[[yName]]
  a <- gdo$X[[aName]]
  b <- gdo$X[[bName]]

  nestedFormula <- as.formula(paste(yName, "~", aName, "/", bName))
  crossedFormula <- as.formula(paste(yName, "~", aName, "*", bName))
  reducedFormula <- as.formula(paste(yName, "~", aName, "+", bName))

  if (method == "nested") {
    numA <- nlevels(a)
    numB <- nlevels(b)
    numMPP <- length(y) / (numB * numA)

    gdo$numO <- numA
    gdo$numP <- numB
    gdo$numM <- numMPP

    fit <- aov(nestedFormula, data = gdo$X)
    meanSq <- anova(fit)[, 3]

    gdo$ANOVA <- fit
    gdo$method <- "nested"

    MSa <- meanSq[1]
    MSab <- meanSq[2]
    MSe <- meanSq[3]

    Cerror <- MSe
    Cb <- (MSab - MSe) / numMPP
    Ca <- (MSa - MSab) / (numB * numMPP)

    if (Ca <= 0) Ca <- 0
    if (Cb <= 0) Cb <- 0

    Cab <- 0
    totalRR <- Ca + Cab + Cerror
    repeatability <- Cerror
    reproducibility <- Ca
    bTob <- Cb
    totalVar <- Cb + Ca + Cab + Cerror

    estimates <- list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
    varcomp <- list(totalRR = totalRR, repeatability = repeatability, reproducibility = reproducibility, bTob = bTob, totalVar = totalVar)

    gdo$Estimates <- estimates
    gdo$Varcomp <- varcomp
  }

  if (method == "crossed") {
    numA <- nlevels(a)
    numB <- nlevels(b)
    numMPP <- length(a) / (numA * numB)

    gdo$numO <- numA
    gdo$numP <- numB
    gdo$numM <- numMPP

    fit <- aov(crossedFormula, data = gdo$X)
    model <- anova(fit)

    gdo$ANOVA <- fit
    gdo$method <- "crossed"

    MSb <- MSa <- MSab <- MSe <- 0

    if (bName %in% row.names(model)) MSb <- model[bName, "Mean Sq"]
    else warning(paste("missing factor", bName, "in model"))

    if (aName %in% row.names(model)) MSa <- model[aName, "Mean Sq"]
    else warning(paste("missing factor", aName, "in model"))

    if (abName %in% row.names(model)) MSab <- model[abName, "Mean Sq"]
    else warning(paste("missing interaction", abName, "in model"))

    if ("Residuals" %in% row.names(model)) MSe <- model["Residuals", "Mean Sq"]
    else warning("missing Residuals in model")

    Cb <- Ca <- Cab <- Cerror <- 0

    Cb <- (MSb - MSab) / (numA * numMPP)
    Ca <- (MSa - MSab) / (numB * numMPP)
    Cab <- (MSab - MSe) / numMPP
    Cerror <- (MSe)

    gdo$RedANOVA <- gdo$ANOVA

    if ((Cab < 0) || (model[abName, "Pr(>F)"] >= alpha)) {
      redFit <- aov(reducedFormula, data = gdo$X)
      model <- anova(redFit)

      MSb <- MSa <- MSab <- MSe <- 0

      if (bName %in% row.names(model)) MSb <- model[bName, "Mean Sq"]
      else warning(paste("missing factor", bName, "in model"))

      if (aName %in% row.names(model)) MSa <- model[aName, "Mean Sq"]
      else warning(paste("missing factor", aName, "in model"))

      if ("Residuals" %in% row.names(model)) MSe <- model["Residuals", "Mean Sq"]
      else warning("missing Residuals in model")

      Cb <- Ca <- Cab <- Cerror <- 0

      Cb <- (MSb - MSe) / (numA * numMPP)
      Ca <- (MSa - MSe) / (numB * numMPP)
      Cab <- 0
      Cerror <- (MSe)

      gdo$RedANOVA <- redFit
    }

    gdo$method <- "crossed"
    Ca <- max(0, Ca)
    Cb <- max(0, Cb)
    Cab <- max(0, Cab)

    totalRR <- Ca + Cab + Cerror
    repeatability <- Cerror
    reproducibility <- Ca + Cab
    bTob <- max(0, Cb)
    totalVar <- Cb + Ca + Cab + Cerror

    estimates <- list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
    varcomp <- list(totalRR = totalRR, repeatability = repeatability, reproducibility = reproducibility, a = Ca, a_b = Cab, bTob = bTob, totalVar = totalVar)

    gdo$Estimates <- estimates
    gdo$Varcomp <- varcomp
  }

  cat("\n")
  cat(paste("AnOVa Table - ", gdo$method, "Design\n"))
  print(summary(gdo$ANOVA))
  cat("\n")
  cat("----------\n")

  if (!identical(gdo$RedANOVA, gdo$ANOVA) && gdo$method == "crossed") {
    cat(paste("AnOVa Table Without Interaction - ", gdo$method, "Design\n"))
    print(summary(gdo$RedANOVA))
    cat("\n")
    cat("----------\n")
  }

  Source <- names(gdo$Varcomp)
  Source[Source == "repeatability"] <- " repeatability"
  Source[Source == "reproducibility"] <- " reproducibility"
  Source[Source == "a_b"] <- paste("  ", abName)
  Source[Source == "a"] <- paste("  ", aName)
  Source[Source == "bTob"] <- bTobName

  VarComp <- round(as.numeric(gdo$Varcomp[c(1:length(gdo$Varcomp))]), 3)
  Contribution <- round(as.numeric(gdo$Varcomp[c(1:length(gdo$Varcomp))]) / as.numeric(gdo$Varcomp[length(gdo$Varcomp)]), 3)
  VarComp <- t(data.frame(gdo$Varcomp))
  VarCompContrib <- VarComp / gdo$Varcomp$totalVar
  Stdev <- sqrt(VarComp)
  StudyVar <- Stdev * gdo$Sigma
  StudyVarContrib <- StudyVar / StudyVar["totalVar", ]
  SNR <- 1
  ptRatio <- NULL
  temp <- NULL

  if ((length(gdo$GageTolerance) > 0) && (gdo$GageTolerance > 0)) {
    ptRatio <- StudyVar / gdo$GageTolerance
    temp <- data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib, ptRatio)
    names(temp)[6] <- c("P/T Ratio")
    row.names(temp) <- c(Source)
  } else {
    temp <- data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib)
    row.names(temp) <- c(Source)
  }

  cat("\n")
  cat("Gage R&R\n")
  tempout <- temp
  print(format(tempout, digits = dig))
  cat("\n")
  cat("---\n")
  cat(" * Contrib equals Contribution in %\n")

  SNRTemp <- sqrt(2) * (temp[bTobName, "Stdev"] / temp["totalRR", "Stdev"])
  if (SNRTemp > 1) SNR <- SNRTemp

  cat(paste(" **Number of Distinct Categories (truncated signal-to-noise-ratio) =", floor(SNR), "\n"))
  cat("\n")
  invisible(gdo)
}

# Ejemplo de uso: ------------------------
# crear un objeto de la clase 'gageRR'
mi_gageRR <- gageRR$new(
  X = data.frame(
    Operator = factor(c("A", "B", "C", "A", "B")),
    Part = factor(c("P1", "P1", "P2", "P2", "P3")),
    Measurement = c(10, 12, 11, 13, 9)
  ),
  ANOVA = NULL,  # Esto puede ser NULL inicialmente y luego calcularlo
  RedANOVA = NULL,  # Igual que ANOVA, puede ser NULL inicialmente
  method = "crossed",
  Estimates = list(),
  Varcomp = list(),
  Sigma = 0.5,
  GageName = "Gage1",
  GageTolerance = 0.1,
  DateOfStudy = "2024-05-15",
  PersonResponsible = "John Doe",
  Comments = "Sample gage R&R study",
  b = factor(c("A", "A", "B", "B", "C")),
  a = factor(c("P1", "P1", "P2", "P2", "P3")),
  y = c(10, 12, 11, 13, 9),
  facNames = c("Measurement", "Operator", "Part"),
  numO = 3,
  numP = 3,
  numM = 2
)

# Crear el objeto gageRRObj
design_example <- gageRRDesign(
  Operators = 3,
  Parts = 10,
  Measurements = 3,
  method = "crossed",
  sigma = 6,
  randomize = TRUE
)

# Crear un diseño para el estudio de Gage
design <- gageRRDesign(Operators = 3, Parts = 10, Measurements = 3, method = "crossed", sigma = 6, randomize = TRUE)
design$X$Measurement <- rnorm(nrow(design$X), mean = 10, sd = 2)

# Ejecutar la función gageRR_
result <- gageRR_(
  gdo = design,
  method = "crossed",   #  método "crossed"
  sigma = 6,            #  sigma
  alpha = 0.25,         # Nivel de significancia
  tolerance = NULL,     # Tolerancia
  dig = 3               # Número de dígitos a mostrar en los resultados
)
class(result)
result$plot()


