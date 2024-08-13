###Clase MSALinearity####
MSALinearity <-R6Class("MSALinearity", public = list(X = data.frame(),
                                                     Y = data.frame(),
                                                     model = NULL,
                                                     conf.level = NULL,
                                                     Linearity = NULL,
                                                     GageName = NULL,
                                                     GageTolerance = NULL,
                                                     DateOfStudy = NULL,
                                                     PersonResponsible = NULL,
                                                     Comments = NULL,
                                                     facNames = NULL,
                                                     response = function(value){
                                                       if (missing(value)) {
                                                         out <- self$Y
                                                         return(out)
                                                       }
                                                       else{
                                                         if (is.vector(value) == TRUE)
                                                           value = data.frame(matrix(value, ncol = ncol(self$Y)))
                                                         if (is.matrix(value) == TRUE)
                                                           value = data.frame(value)
                                                         self$Y = value
                                                         invisible(self)
                                                       }

                                                     },
                                                     summary = function(){
                                                       cat("----------------------", fill = TRUE)
                                                       print(self)
                                                       cat("----------------------", fill = TRUE)
                                                       if(!is.null(self$model)){
                                                         print(summary(self$model))
                                                         cat("----------------------", fill = TRUE)
                                                       }
                                                       if(!is.null(self$Linearity)){
                                                         names(self$Linearity) = "Linearity:"
                                                         print(self$Linearity)
                                                       }
                                                     },
                                                     plot = function(ylim, col, pch, lty = c(1, 2), ...){
                                                       conf.level = self$conf.level
                                                       g = nrow(self$X[2])
                                                       m = ncol(self$Y)
                                                       if (missing(col))
                                                         col = c(1, 2, 1, 4)
                                                       if (missing(pch))
                                                         pch = c(20, 18)
                                                       bias = self$Y
                                                       mbias = numeric(g)
                                                       for (i in 1:g) bias[i, ] = self$Y[i, ] - self$X[i, 2]
                                                       for (i in 1:g) mbias[i] = mean(as.numeric(bias[i, ]))
                                                       if (missing(ylim))
                                                         ylim = c(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE))
                                                       plot(x = self$X$Ref, y = mbias, ylim = ylim, col = col[2], pch = pch[2], ylab = "Bias", xlab = "Reference Values", ...)
                                                       for (i in 1:g) points(x = rep(self$X$Ref[i], length = m), y = bias[i, ], col = col[1], pch = pch[1])
                                                       abline(h = 0, lty = 3, col = "gray")
                                                       BIAS = numeric()
                                                       ref = numeric()
                                                       for (i in 1:g) {
                                                         BIAS = c(BIAS, as.numeric(bias[i, ]))
                                                         ref = c(ref, rep(self$X$Ref[i], length = m))
                                                       }
                                                       lm.1 = lm(formula = BIAS ~ ref)
                                                       a = lm.1[[1]][2]
                                                       names(a) = "slope"
                                                       b = lm.1[[1]][1]
                                                       names(b) = "intercept"
                                                       y.vec = numeric()
                                                       for (i in 1:g) y.vec = c(y.vec, self$Y[i, ])
                                                       pre = predict.lm(lm.1, interval = "confidence", level = conf.level)
                                                       lines(ref, pre[, 1], col = col[3], lty = lty[1])
                                                       lines(ref, pre[, 2], col = col[4], lty = lty[2])
                                                       lines(ref, pre[, 3], col = col[4], lty = lty[2])
                                                       legend("topright", legend = c("Single Bias", "Mean Bias", "Regression", paste(conf.level * 100, "% conf.level")), pch = c(pch, -1, -1), col = col, lty = c(-1,-1, lty), inset = 0.04)
                                                     },
                                                     print = function(){
                                                       print(self$as.data.frame())
                                                     },
                                                     as.data.frame = function(row.names = NULL, optional = FALSE, ...){
                                                       return(cbind(self$X,self$Y))
                                                     }
                                                    )

                      )

###funcion gageLin####
gageLin <- function(object, conf.level = 0.95, ylim, col, pch, lty = c(1, 2), stats = TRUE, plot = TRUE){
  if (class(object)[1]!="MSALinearity")
    stop("object needs to be from class 'MSALinearity'")
  object$conf.level = conf.level
  g = nrow(object$X[2])
  m = ncol(object$Y)
  if (missing(col))
    col = c(1, 2, 1, 4)
  if (missing(pch))
    pch = c(20, 18)
  bias = object$Y
  mbias = numeric(g)
  for (i in 1:g) bias[i, ] = object$Y[i, ] - object$X[i, 2]
  for (i in 1:g) mbias[i] = mean(as.numeric(bias[i, ]))
  if (missing(ylim))
    ylim = c(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE))
  if (plot == TRUE)
    plot(x = object$X$Ref, y = mbias, ylim = ylim, col = col[2], pch = pch[2], ylab = "Bias", xlab = "Reference Values")
  if (plot == TRUE) {
    for (i in 1:g) points(x = rep(object$X$Ref[i], length = m), y = bias[i, ], col = col[1], pch = pch[1])
    abline(h = 0, lty = 3, col = "gray")
  }
  BIAS = numeric()
  ref = numeric()
  for (i in 1:g) {
    BIAS = c(BIAS, as.numeric(bias[i, ]))
    ref = c(ref, rep(object$X$Ref[i], length = m))
  }
  lm.1 = lm(formula = BIAS ~ ref)
  object$model <- lm.1
  a = lm.1[[1]][2]
  names(a) = "slope"
  b = lm.1[[1]][1]
  names(b) = "intercept"
  y.vec = numeric()
  for (i in 1:g) y.vec = c(y.vec, object$Y[i, ])
  if (plot == TRUE) {
    pre = predict.lm(lm.1, interval = "confidence", level = conf.level)
    lines(ref, pre[, 1], col = col[3], lty = lty[1])
    lines(ref, pre[, 2], col = col[4], lty = lty[2])
    lines(ref, pre[, 3], col = col[4], lty = lty[2])
  }
  test = numeric(g + 1)
  for (i in 1:g) test[i] = t.test(bias[, i], mu = 0, conf.level = conf.level)["p.value"]
  test[g + 1] = t.test(BIAS, mu = 0, conf.level = conf.level)["p.value"]
  Linearity = abs(lm.1[[1]][2]) * 100
  object$Linearity = Linearity
  names(Linearity) = "LINEARITY:"
  if (plot == TRUE) {
    legend("topright", legend = c("Single Bias", "Mean Bias", "Regression", paste(conf.level * 100, "% conf.level")), pch = c(pch, -1, -1), col = col, lty = c(-1,-1, lty), inset = 0.04)
  }
  if (stats == TRUE) {
    cat("----------------------", fill = TRUE)
    cat("BIAS:", fill = TRUE)
    print(bias)
    cat("----------------------", fill = TRUE)
    cat("MEAN OF BIAS:", fill = TRUE)
    temp = mbias
    names(temp) = rownames(object$Y)
    print(temp)
    cat("----------------------", fill = TRUE)
    cat("LINEAR MODEL:", fill = TRUE)
    print(summary(lm.1))
    cat("----------------------", fill = TRUE)
    print(Linearity)
    return(object)
  }
  else return(object)
}
###funcion gageLinDesign####
gageLinDesign <- function(ref, n = 5) {
  numColY = n
  numRowY = length(ref)
  X = data.frame(cbind(Part = 1:length(ref), Ref = ref))
  Y = data.frame(matrix(NA, nrow = numRowY, ncol = numColY))
  for (i in 1:n) names(Y)[i] = paste("n", i, sep = "")
  gageLinDesign = MSALinearity$new()
  gageLinDesign$X = X
  gageLinDesign$Y = Y
  return(gageLinDesign)
}

#ejemplo####
# results of run A-E
A=c(2.7,2.5,2.4,2.5,2.7,2.3,2.5,2.5,2.4,2.4,2.6,2.4)
B=c(5.1,3.9,4.2,5,3.8,3.9,3.9,3.9,3.9,4,4.1,3.8)
C=c(5.8,5.7,5.9,5.9,6,6.1,6,6.1,6.4,6.3,6,6.1)
D=c(7.6,7.7,7.8,7.7,7.8,7.8,7.8,7.7,7.8,7.5,7.6,7.7)
E=c(9.1,9.3,9.5,9.3,9.4,9.5,9.5,9.5,9.6,9.2,9.3,9.4)

# create Design
test=gageLinDesign(ref=c(2,4,6,8,10),n=12)
# create data.frame for results
Messungen=data.frame(rbind(A,B,C,D,E))
# enter results in Design
test$response(Messungen)
test$summary()


# no plot and no return
MSALin=gageLin(test,stats=FALSE,plot=FALSE)
# plot only
plot(MSALin)
MSALin$plot()
# summary
test$summary()

