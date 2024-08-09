#########################################################################
######################### DISTRIBUTION - CLASES #########################
#########################################################################

# Class Distr ----
Distr <- R6Class("Distr",
                 public = list(
                   x = NULL,
                   name = NULL,
                   parameters = NULL,
                   sd = NULL,
                   n = NULL,
                   loglik = NULL,

                   initialize = function(x, name, parameters, sd, n, loglik) {
                     self$x <- x
                     self$name <- name
                     self$parameters <- parameters
                     self$sd <- sd
                     self$n <- n
                     self$loglik <- loglik
                   },

                   plot = function(main = NULL, xlab = NULL, xlim = NULL, ylim = NULL, ylab = NULL, line.col = "red", box=TRUE,line.width = 1, ...)
                   {
                     object <- self
                     xVals <- object$x
                     parameters <- object$parameters
                     lq <- NULL
                     uq <- NULL
                     y <- NULL

                     if (missing(line.col)) {
                       line.col <- "red"
                     }
                     if (missing(line.width)) {
                       line.width <- 1
                     }
                     if (missing(main)) {
                       main <- object$name
                     }
                     if (missing(xlab)) {
                       xlab <- "x"
                     }
                     if (missing(ylab)) {
                       ylab <- "Density"
                     }

                     distr <- object$name
                     qFun <- .charToDistFunc(distr, type = "q")
                     dFun <- .charToDistFunc(distr, type = "d")
                     adTestStats <- .myADTest(xVals, distr)

                     if (class(adTestStats) == "adtest") {
                       A <- adTestStats$statistic
                       p <- adTestStats$p.value
                     } else {
                       A <- NA
                       p <- NA
                     }

                     histObj <- hist(xVals, plot = FALSE)
                     df <- data.frame(
                       mid = histObj$mids,
                       density = histObj$density
                     )
                     width <- diff(df$mid)[1]

                     if (missing(xlim)) {
                       lq <- do.call(qFun, c(list(1e-04), as.list(parameters)))
                       uq <- do.call(qFun, c(list(0.9999), as.list(parameters)))
                       xlim <- range(lq, uq, xVals)
                     }

                     xPoints <- seq(xlim[1], xlim[2], length = 200)
                     yPoints <- do.call(dFun, c(list(xPoints), as.list(parameters)))

                     if (missing(ylim)) {
                       ylim <- range(0, histObj$density, yPoints)
                     }

                     # Histograma
                     p1 <- ggplot(df, aes(x = mid, y = density)) +
                       geom_bar(stat = "identity", width = width, fill = "lightblue", color = "black", alpha = 0.5) +
                       labs(y = ylab, x = xlab, title = main) + xlim(xlim) + ylim(ylim) +
                       theme_minimal() + theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
                       guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
                       geom_line(data = data.frame(x = xPoints, y = yPoints), aes(x = x, y = y), color = line.col, linewidth = line.width) + # densidad
                       theme(legend.position = "none")

                     # Caja de Info
                     if (box==FALSE) {
                       p1
                     }
                     else {
                       p2 <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
                         theme_bw() +
                         theme(
                           axis.text = element_blank(),
                           axis.ticks = element_blank(),
                           axis.title = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()
                         ) +
                         xlim(c(0.25,0.26)) + ylim(c(0.19, 0.36))
                       {
                         # n y A
                         p2 <- p2 +
                           annotate('text', x = 0.25, y = 0.25,
                                    label = paste("A==", round(as.numeric(A), digits = 3)),
                                    parse = TRUE, size = 3, hjust = 0)
                         # p
                         if (!is.null(adTestStats$smaller) && adTestStats$smaller){
                           p2 <- p2 +
                             annotate('text',x = 0.25,y = 0.20,
                                      label = paste("p<", round(as.numeric(p), digits =3)),
                                      parse = TRUE,size = 3,hjust = 0)
                         }
                         if (!is.null(adTestStats$smaller) && !adTestStats$smaller){
                           p2 <- p2 +
                             annotate('text',x = 0.25, y = 0.20,
                                      label = paste("p>=", round(as.numeric(p),digits = 3)),
                                      parse = TRUE,size = 3,hjust = 0)
                         }
                         if (is.null(adTestStats$smaller)){
                           p2 <- p2 +
                             annotate('text',x = 0.25,y = 0.20,
                                      label = paste("p==", round(as.numeric(p), digits = 3)),
                                      parse = TRUE,size = 3,hjust = 0)
                         }

                         # mean y sd
                         p2 <- p2 + annotate('text', x = 0.25, y = 0.35,
                                             label = paste("mean==", round(self$parameters[[1]], digits = 3)),
                                             parse = TRUE, size = 3, hjust = 0) +
                           annotate('text', x = 0.25, y = 0.30,
                                    label = paste("sd==", round(self$parameters[[2]], digits = 3)),
                                    parse = TRUE, size = 3, hjust = 0)
                         }

                       p1 + inset_element(p2, left = 0.7, right = 1, top = 1, bottom = 0.60)
                     }
                   }
                 )
)


# Class DistrCollection ----
DistrCollection <- R6::R6Class("DistrCollection",
                               public = list(
                                 distr = NULL,
                                 initialize = function() {
                                   self$distr <- list()
                                 },
                                 add = function(distr) {
                                   self$distr <- append(self$distr, list(distr))
                                 },
                                 get = function(i) {
                                   self$distr[[i]]
                                 },
                                 print = function() {
                                   cat("\n")
                                   for (i in seq_along(self$distr)) {
                                     temp <- self$distr[[i]]
                                     cat("\n")
                                     cat("fitted distribution is", temp$name, ":\n")
                                     print(temp$parameters)
                                     cat("\n")
                                   }
                                 },
                                 summary = function() {
                                   numDist <- length(self$distr)
                                   gofMatrix <- data.frame(matrix(nrow = numDist, ncol = 3))
                                   names(gofMatrix) <- c("Distribution", "A", "p.value")
                                   cat("\n------ Fitted Distribution and estimated parameters ------\n")
                                   for (i in seq_along(self$distr)) {
                                     distrObj <- self$distr[[i]]
                                     x <- distrObj$x
                                     distribution <- distrObj$name
                                     parameters <- distrObj$parameters
                                     statistic <- NA
                                     p.value <- NA
                                     temp <- .myADTest(x, distribution)
                                     try(statistic <- as.numeric(temp$statistic), silent = TRUE)
                                     try(p.value <- as.numeric(temp$p.value), silent = TRUE)
                                     gofMatrix[i, ] <- c(distribution, as.numeric(statistic), as.numeric(p.value))
                                     cat("\n")
                                     cat("fitted distribution is", distribution, ":\n")
                                     print(parameters)
                                   }
                                   cat("\n")
                                   cat("\n------ Goodness of Fit - Anderson Darling Test ------\n")
                                   cat("\n")
                                   gofMatrixPrint <- gofMatrix
                                   gofMatrixPrint[, 2] <- signif(as.numeric(gofMatrixPrint[, 2]), 4)
                                   gofMatrixPrint[, 3] <- signif(as.numeric(gofMatrixPrint[, 3]), 4)
                                   print(gofMatrixPrint)
                                 },
                                 plot = function(xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, line.col = "red", line.width = 1, box = TRUE , ...) {
                                   distrList <- self$distr
                                   numDist <- length(self$distr)
                                   numColWin <- ceiling(numDist/2)
                                   if (missing(xlim)) {
                                     xlim <- .xyLimits(self)$xlim
                                   }
                                   if (missing(ylim)) {
                                     ylim <- .xyLimits(self)$ylim
                                   }
                                   if (missing(line.col)) {
                                     line.col <- "red"
                                   }
                                   if (missing(line.width)) {
                                     line.width <- 1
                                   }
                                   p <- distrList[[1]]$plot(xlab = xlab, ylab = ylab, line.col = line.col, line.width = line.width, box = box)
                                   for (i in 2:length(distrList)) {
                                     p <- p+distrList[[i]]$plot(xlab = xlab, ylab = ylab, line.col = line.col, line.width = line.width, box = box)
                                   }
                                   p
                                 }
                               )
)
