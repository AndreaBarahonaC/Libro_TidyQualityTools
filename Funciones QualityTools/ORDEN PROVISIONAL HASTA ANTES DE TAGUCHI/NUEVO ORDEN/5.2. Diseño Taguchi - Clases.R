############################################################################
####################### DISEÑO TAGUCHI - CLASES ############################
############################################################################

# Clase taguchiFactor ----
taguchiFactor <- R6Class("taguchiFactor", public = list(values = NA,
                                                        name = "",
                                                        unit = "",
                                                        type = "numeric",
                                                        attributes = function(){
                                                          v <- c(self$values,self$name, self$unit, self$type)
                                                        },
                                                        .values = function(value){
                                                          if (missing(value)) {
                                                            return(self$values)
                                                          }
                                                          else{
                                                            self$values <- value
                                                            invisible(self)
                                                          }
                                                        },
                                                        .unit = function(value){
                                                          if(missing(value)){
                                                            return(self$unit)
                                                          }
                                                          else{
                                                            self$unit <- value
                                                            invisible(self)
                                                          }
                                                        },
                                                        names = function(value){
                                                          if(missing(value)){
                                                            return(self$name)
                                                          }
                                                          else{
                                                            self$name <- value
                                                            invisible(self)
                                                          }
                                                        }
)
)
# Clase taguchiDesign.c ----
taguchiDesign.c <- R6Class("taguchiDesign", public = list(name = NULL,
                                                          factors = list(),
                                                          design = data.frame(),
                                                          designType = NULL,
                                                          replic = data.frame(),
                                                          response = data.frame(),
                                                          Type = data.frame(),
                                                          block = data.frame(),
                                                          runOrder = data.frame(),
                                                          standardOrder = data.frame(),
                                                          desireVal = list(),
                                                          desirability = list(),
                                                          fits = data.frame(),
                                                          values = function(value){
                                                            if(missing(value)){
                                                              listOut = vector(mode = "list")
                                                              for (i in names(self$design)) {
                                                                listOut[[i]] = self$factors[[i]]$.values()
                                                              }
                                                              return(listOut)
                                                            }
                                                            else{
                                                              for (i in names(value)) {
                                                                if (i %in% names(self$design))
                                                                  if (length(value[[i]]) == length(unique(self$design[, i])))
                                                                    self$factors[[i]]$.values(value[[i]])
                                                                else stop("Number of values greater or less than number of factor settings!")
                                                              }
                                                              invisible(self)
                                                            }

                                                          },
                                                          units = function(value){
                                                            if (missing(value)) {
                                                              v <- list()
                                                              for (i in 1:length(self$factors)) {
                                                                v[[unlist(self$names()[i])]] <- self$factors[[i]]$.unit()
                                                              }
                                                              return(v)
                                                            }
                                                            else{
                                                              for (i in 1:length(self$factors)) if (length(value) > 1)
                                                                self$factors[[i]]$.unit(as.character(value[i]))
                                                              else self$factors[[i]]$.unit(as.character(value[1]))
                                                              invisible(self)
                                                            }
                                                          },
                                                          .factors = function(value){
                                                            if (missing(value)) {
                                                              return(self$factors)
                                                            }
                                                            else{
                                                              if (length(value) != ncol(self$design))
                                                                stop("\nNumber of factors doesn't match with number of columns for factorial Design\n")
                                                              self$factors <- value
                                                              invisible(self)
                                                            }
                                                          },
                                                          names = function(value){
                                                            if(missing(value)){
                                                              aux <- list()
                                                              for (i in 1:length(self$factors)) {
                                                                aux[[.NAMES[i]]] <-self$factors[[i]]$name
                                                              }
                                                              return(aux)
                                                            }
                                                            else {
                                                              for (i in 1:length(self$factors)){
                                                                self$factors[[i]]$name = as.character(value[i])

                                                              }
                                                              invisible(self)
                                                            }
                                                          },
                                                          as.data.frame = function(row.names = NULL, optional = FALSE, ...){
                                                            frameOut = cbind(self$standardOrder, self$runOrder, self$replic, self$design, self$response)
                                                            return(frameOut)
                                                          },
                                                          print = function(){
                                                            print(format(as.data.frame(self), digits = 4))
                                                          },
                                                          .response = function(value){
                                                            if(missing(value)){
                                                              return(self$response)
                                                            }
                                                            else{
                                                              if (!is.numeric(value) & !is.data.frame(value))
                                                                stop("vector or data.frame must be given")
                                                              if (is.numeric(value)) {
                                                                if (length(value) != nrow(self$design))
                                                                  stop("differing lengths")
                                                                temp = data.frame(value)
                                                                names(temp) = deparse(substitute(value))[1]
                                                                value = temp
                                                              }
                                                              if (is.data.frame(value)) {
                                                                if (nrow(value) != nrow(self$design))
                                                                  stop("differing number of rows")
                                                              }
                                                              self$response = value
                                                              invisible(self)
                                                            }
                                                          },
                                                          .nfp = function(){
                                                            x = self$.factors()
                                                            DB = FALSE
                                                            if (is.list(x) && length(x[[1]]) > 0) {
                                                              numAttr = length(x[[1]]$attributes())
                                                              .numFac = length(x)
                                                              #len = 0
                                                              # for (i in names(x)) if (length(x[[i]]$values) > len)
                                                              #   len = length(x[[i]]$values)
                                                              #numAttr = numAttr + len
                                                              numrows = numAttr #- 1
                                                              frameOut = data.frame(matrix(NA, ncol = .numFac, nrow = numrows))
                                                              names(frameOut) = names(x)
                                                              rownames(frameOut) = c(paste("value", 1:len), "name", "unit", "type")
                                                              for (i in names(x)) {
                                                                vin = 1:length(x[[i]]$values)
                                                                frameOut[vin, i] = x[[i]]$values
                                                                frameOut[numrows - 2, i] = x[[i]]$name
                                                                frameOut[numrows - 1, i] = x[[i]]$unit
                                                                frameOut[numrows, i] = x[[i]]$type
                                                              }
                                                              print(frameOut)
                                                            }

                                                          },
                                                          summary = function(){
                                                            cat(paste("Taguchi", toupper(self$designType), "Design"))
                                                            cat("\n")
                                                            cat("Information about the factors:\n\n")
                                                            self$.nfp()
                                                            cat("\n")
                                                            cat("-----------\n")
                                                            cat("\n")
                                                            print(self$as.data.frame())
                                                            cat("\n")
                                                            cat("-----------\n")
                                                            cat("\n")
                                                          },
                                                          effectPlot = function(factors, fun = mean, response = NULL, single = FALSE, points = FALSE, classic = FALSE,  ###
                                                                                axes = TRUE, lty, xlab, ylab, main, ylim, ...){
                                                            oldMar = par("mar")
                                                            oldOma = par("oma")
                                                            oldMfrow = par("mfrow")
                                                            oldMfcol = par("mfcol")
                                                            on.exit(par(mar = oldMar, oma = oldOma, mfrow = oldMfrow, mfcol = oldMfcol))
                                                            if(is.null(response)==FALSE)                                                ###
                                                            {                                                                           ###
                                                              temp=self$.response()[response]                                            ###
                                                              self$.response(temp)                                                      ###
                                                            }                                                                           ###
                                                            ylabmiss = FALSE
                                                            xlabmiss = FALSE
                                                            mainmiss = FALSE
                                                            ylimmiss = FALSE
                                                            if (missing(ylim))
                                                              ylimmiss = TRUE
                                                            if (missing(lty))
                                                              lty = 1
                                                            X = self$design
                                                            Y = self$.response()
                                                            if (!missing(factors))
                                                              k = length(factors)
                                                            else #(missing(factors))                                                    ###
                                                            {
                                                              k = ncol(X)
                                                              factors = names(X)
                                                            }
                                                            numCol = 1
                                                            numRow = 1
                                                            if (!single && missing(factors)) {                                          ###
                                                              if (ncol(X) == 2) {
                                                                numCol = 2
                                                                numRow = 1
                                                              }
                                                              if (ncol(X) > 2) {
                                                                numCol = 2
                                                                numRow = 2
                                                              }
                                                            }
                                                            if (!single && !missing(factors)) {                                         ###
                                                              if (length(factors) == 2) {                                             ###
                                                                numCol = 2                                                          ###
                                                                numRow = 1                                                          ###
                                                              }                                                                       ###
                                                              if (length(factors) == 3) {                                             ###
                                                                numCol = 3                                                          ###
                                                                numRow = 1                                                          ###
                                                              }                                                                       ###
                                                              if (length(factors) == 4) {                                             ###
                                                                numCol = 2                                                          ###
                                                                numRow = 2                                                          ###
                                                              }                                                                       ###
                                                              if (length(factors) == 5) {                                             ###
                                                                numCol = 3                                                          ###
                                                                numRow = 2                                                          ###
                                                              }                                                                       ###
                                                              if (length(factors) == 6) {                                             ###
                                                                numCol = 3                                                          ###
                                                                numRow = 2                                                          ###
                                                              }                                                                       ###
                                                              if (length(factors) > 6) {                                              ###
                                                                numRow = ceiling(sqrt(length(factors)))                             ###
                                                                numCol = ceiling(sqrt(length(factors)))                             ###
                                                              }                                                                       ###
                                                            }                                                                           ###
                                                            if (classic) {
                                                              numCol = ncol(X)
                                                              numRow = 1
                                                            }
                                                            if (!single)
                                                              par(mfrow = c(numRow, numCol))
                                                            nextResponse = FALSE
                                                            for (j in 1:ncol(Y)) {
                                                              counter = 0
                                                              cells = numeric(0)
                                                              for (i in 1:length(factors)) {
                                                                cells = c(cells, as.vector(tapply(Y[, j], list(X[, factors[i]], rep(0, nrow(X))), fun)))
                                                                if (points)
                                                                  cells = range(Y)
                                                              }
                                                              if (nextResponse & !single) {
                                                                dev.new()
                                                                par(mfrow = c(numRow, numCol))
                                                              }
                                                              for (i in 1:length(factors)) {
                                                                if ((counter != 0 & counter%%(numCol * numRow) == 0) & !single) {
                                                                  dev.new()
                                                                  par(mfrow = c(numRow, numCol))
                                                                }
                                                                if (missing(main)) {
                                                                  main = paste("Effect Plot for", names(Y)[j])
                                                                  mainmiss = TRUE
                                                                }
                                                                if (mainmiss)
                                                                  main = paste("Effect Plot for", names(Y)[j])
                                                                if (missing(xlab)) {
                                                                  xlab = factors[i]
                                                                  xlabmiss = TRUE
                                                                }
                                                                if (xlabmiss) {
                                                                  if (identical(" ", names(self)[[i]]))
                                                                    xlab = factors[i]
                                                                  else xlab = paste(factors[i], ": ", names(self)[[i]], sep = "")
                                                                }
                                                                if (missing(ylab)) {
                                                                  ylab = paste(deparse(substitute(fun)), "of ", names(Y)[j])
                                                                  ylabmiss = TRUE
                                                                }
                                                                if (ylabmiss)
                                                                  ylab = paste(deparse(substitute(fun)), "of ", names(Y)[j])
                                                                if (ylimmiss)
                                                                  ylim = range(cells, na.rm = TRUE)
                                                                if (classic & i == 1) {
                                                                  par(mar = c(5, 0, 0, 0) + 0.1)
                                                                  par(oma = c(-0.1, 4, 4, 1) + 0.1)
                                                                }
                                                                if (classic) {
                                                                  .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun,
                                                                                      ylab = ylab, legend = FALSE, axes = FALSE, main = " ", ...)
                                                                  grid(NA, 2)
                                                                  axis(1, at = X[, factors[i]])
                                                                  if (i == 1)
                                                                    axis(2)
                                                                  box()
                                                                  title(main, outer = TRUE)
                                                                }
                                                                else {
                                                                  .m.interaction.plot(x.factor = X[, factors[i]], trace.factor = rep(0, nrow(X)), response = Y[, j], lty = lty, ylim = ylim, xlab = xlab, fun = fun,
                                                                                      ylab = ylab, legend = FALSE, axes = axes, main = main, ...)
                                                                  grid(NA, 2)
                                                                }
                                                                if (points)
                                                                  points(X[, factors[i]], Y[, j], ...)
                                                                counter = counter + 1
                                                              }
                                                              nextResponse = TRUE
                                                            }

                                                          },
                                                          identity = function(){
                                                            identity = character(0)
                                                            identityList = vector(mode = "list", length = 0)
                                                            resolution = numeric(0)
                                                            temp = NULL
                                                            A = aliasTable(self)
                                                            if (any(dim(A) == 0))
                                                              return(identityList)
                                                            temp = as.matrix(A["Identity", ])
                                                            boolTemp = apply(temp, 2, as.logical)
                                                            identity = row.names(temp)[boolTemp[, 1]]
                                                            if (length(identity) > 0) {
                                                              charList = strsplit(toupper(identity), split = "")
                                                              identityList = lapply(charList, match, LETTERS[1:26])
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
                                                          }

)
)


