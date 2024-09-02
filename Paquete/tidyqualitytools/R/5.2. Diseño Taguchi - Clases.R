############################################################################
####################### DISEÑO TAGUCHI - CLASES ############################
############################################################################

# Clase taguchiFactor ----
#' @title taguchiFactor
#' @description An R6 class representing a factor in a Taguchi design.
#' @field values A vector containing the levels or values associated with the factor. Default is \code{NA}.
#' @field name A character string specifying the name of the factor. Default is an empty string \code{""}.
#' @field unit A character string specifying the unit of measurement for the factor. Default is an empty string \code{""}.
#' @field type A character string specifying the type of the factor, which can be either \code{"numeric"} or \code{"categorical"}. Default is \code{"numeric"}.
taguchiFactor <- R6Class("taguchiFactor", public = list(values = NA,
                                                        name = "",
                                                        unit = "",
                                                        type = "numeric",

                                                        #' @description Get the attributes of the factor.
                                                        attributes = function(){
                                                          v <- c(self$values, self$name, self$unit, self$type)
                                                        },

                                                        #' @description Get and set the \code{values} for the factors in an object of class \code{taguchiFactor}.
                                                        #' @param value New values, If missing value get the \code{values}.
                                                        .values = function(value){
                                                          if (missing(value)) {
                                                            return(self$values)
                                                          }
                                                          else{
                                                            self$values <- value
                                                            invisible(self)
                                                          }
                                                        },

                                                        #' @description Get and set the \code{units} for the factors in an object of class \code{taguchiFactor}.
                                                        #' @param value New unit, If missing value get the \code{units}.
                                                        .unit = function(value){
                                                          if(missing(value)){
                                                            return(self$unit)
                                                          }
                                                          else{
                                                            self$unit <- value
                                                            invisible(self)
                                                          }
                                                        },

                                                        #' @description Get and set the \code{names} in an object of class \code{taguchiFactor}.
                                                        #' @param value New names, If missing value get the \code{names}.
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
#' @title taguchiDesign
#' @description An R6 class representing a Taguchi experimental design.
#' @field name A character string specifying the name of the design. Default is \code{NULL}.
#' @field factors A list of factors included in the Taguchi design. Each factor is typically an instance of the \code{taguchiFactor} class.
#' @field design A `data.frame` representing the design matrix of the experiment. This includes the levels of each factor for every run of the experiment. Default is an empty \code{data.frame}.
#' @field designType A character string specifying the type of Taguchi design used. Default is \code{NULL}.
#' @field replic A `data.frame` containing the replication information for the design. Default is an empty \code{data.frame}.
#' @field response A `data.frame` storing the response values collected from the experiment. Default is an empty \code{data.frame}.
#' @field Type A `data.frame` specifying the type of responses or factors involved in the design. Default is an empty \code{data.frame}.
#' @field block A `data.frame` indicating any blocking factors used in the design. Default is an empty \code{data.frame}.
#' @field runOrder A `data.frame` detailing the order in which the experimental runs were conducted. Default is an empty \code{data.frame}.
#' @field standardOrder A `data.frame` detailing the standard order of the experimental runs. Default is an empty \code{data.frame}.
#' @field desireVal A list storing desired values for responses in the experiment. Default is an empty list.
#' @field desirability A list storing desirability functions used to evaluate the outcomes of the experiment. Default is an empty list.
#' @field fits A `data.frame` containing model fits or other statistical summaries from the analysis of the experimental data. Default is an empty \code{data.frame}.
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

                                                          #' @description Get and set the \code{values} for an object of class \code{taguchiDesign}.
                                                          #' @param value New value, If missing value get the \code{values}.
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

                                                          #' @description Get and set the \code{units} for an object of class \code{taguchiDesign}.
                                                          #' @param value New units, If missing value get the \code{units}.
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

                                                          #' @description Get and set the \code{factors} in an object of class \code{taguchiDesign}.
                                                          #' @param value New factors, If missing value get the \code{factors}.
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

                                                          #' @description Get and set the \code{names} in an object of class \code{taguchiDesign}.
                                                          #' @param value New names, If missing value get the \code{names}.
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

                                                          #' @description Methods for function \code{print} in Package \code{base}.
                                                          print = function(){
                                                            print(format(as.data.frame(self), digits = 4))
                                                          },

                                                          #' @description Get and set the the \code{response} in an object of class \code{taguchiDesign}.
                                                          #' @param value New response, If missing value get the \code{response}.
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

                                                          #' @description Prints a summary of the factors attributes including their low, high, name, unit, and type.
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

                                                          #' @description Methods for function \code{summary} in Package \code{base}.
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

                                                          #' @description Plots the effects of factors on the response variables.
                                                          #' @param factors Factors to be plotted.
                                                          #' @param fun Function applied to the response variables (e.g., mean).
                                                          #' @param response Optional; specifies which response variables to plot.
                                                          #' @param single Logical; if TRUE, plots effects for single factor; otherwise, for combinations of factors.
                                                          #' @param points Logical; if TRUE, plots data points.
                                                          #' @param classic Logical; if TRUE, uses classic plotting style.
                                                          #' @param axes Logical; if TRUE, includes axes in the plot.
                                                          #' @param lty Line type for plotting.
                                                          #' @param xlab Label for the x-axis.
                                                          #' @param ylab Label for the y-axis.
                                                          #' @param main Main title for the plot.
                                                          #' @param ylim Limits for the y-axis.
                                                          #' @param ... Additional plotting parameters.
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

                                                          #' @description Calculates the alias table for a fractional factorial design and prints an easy to read summary of the defining relations such as 'I = ABCD' for a standard 2^(4-1) factorial design.
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


