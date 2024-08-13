##############################################################
##################### gageRR - Clase #########################
##############################################################

# gageRR.c ----
gageRR.c <- R6Class("gageRR",
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
                      print = function() {
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
                      get.response = function() {
                        return(self$X$Measurement)
                      },
                      response = function(value) {
                        self$X$Measurement = value
                      },
                      names = function() {
                        return(names(as.data.frame(self$X)))
                      },
                      as.data.frame = function() {
                        return(as.data.frame(self$X))
                      },
                      get.tolerance = function() {
                        return(unlist(self$GageTolerance))
                      },
                      set.tolerance = function(value) {
                        if (!is.numeric(value))
                          stop("GageTolerance needs to be numeric")
                        self$GageTolerance = value
                        return(self)
                      },
                      get.sigma = function() {
                        return(unlist(self$Sigma))
                      },
                      set.sigma = function(value) {
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

                                legend.key.size = unit(0.3, "cm"), # reduce the size of the legend key
                                legend.text = element_text(size = 6), # reduce the size of the legend text
                                legend.box.spacing = unit(0.1, 'cm'),
                                legend.title = element_blank(),
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

                          p2 <- suppressWarnings(
                            ggplot(gdo$X, aes_string(x = bName, y = yName)) +
                              geom_boxplot() +
                              stat_summary(fun = median, geom = "line", aes(group = 1), color = "red", linewidth = lwd) +
                              stat_summary(fun = median, geom = "point", color = "red") +
                              labs(title = ifelse(is.null(main2), paste(yName, "by", bName), main2),
                                   x = ifelse(is.null(xlab2), bName, xlab2),
                                   y = ifelse(is.null(ylab2), yName, ylab2)) +
                              theme_bw() +
                              theme(plot.title = element_text(hjust = 0.5))
                          )


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
                            geom_hline(yintercept = xm, color = 3) +                          # lC-nea xm
                            geom_hline(yintercept = UCL, color = "red") +                     # lC-nea UCL
                            geom_hline(yintercept = LCL, color = "red") +                     # lC-nea LCL
                            labs(x = aName, y = expression(bar(x))) +                         # tC-tulos ejes
                            ggtitle(expression(paste(bar(x), " Chart"))) +                    # tC-tulo
                            geom_vline(xintercept = cumsum(tab) - 0.5, linetype = "dashed") + # divide categorC-as A,B C
                            scale_x_continuous(breaks = cumsum(tab), labels = names(tab)) +   # y pone los nombres de las secciones
                            theme_bw() +
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
                            geom_hline(yintercept = Rm, color = 3) +                          # lC-nea xm
                            geom_hline(yintercept = UCL, color = "red") +                     # lC-nea UCL
                            geom_hline(yintercept = LCL, color = "red") +                     # lC-nea LCL
                            labs(x = aName, y = "R") +                                        # tC-tulos ejes
                            ggtitle("R Chart") +                                              # tC-tulo
                            geom_vline(xintercept = cumsum(tab) - 0.5, linetype = "dashed") + # divide categorC-as A,B C
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
                            geom_hline(yintercept = xm, color = 3) +                          # lC-nea xm
                            geom_hline(yintercept = UCL, color = "red") +                     # lC-nea UCL
                            geom_hline(yintercept = LCL, color = "red") +                     # lC-nea LCL
                            labs(x = aName, y = expression(bar(x))) +                         # tC-tulos ejes
                            ggtitle(expression(paste(bar(x), " Chart"))) +                    # tC-tulo
                            geom_vline(xintercept = cumsum(tab) - 0.5, linetype = "dashed") + # divide categorC-as A,B C
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
                            geom_hline(yintercept = Rm, color = 3) +                          # lC-nea xm
                            geom_hline(yintercept = UCL, color = "red") +                     # lC-nea UCL
                            geom_hline(yintercept = LCL, color = "red") +                     # lC-nea LCL
                            labs(x = aName, y = "R") +                                        # tC-tulos ejes
                            ggtitle("R Chart") +                                              # tC-tulo
                            geom_vline(xintercept = cumsum(tab) - 0.5, linetype = "dashed") + # divide categorC-as A,B C
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

