####################################################################################
####################### DISEÑO PLACKETT BURMAN - CLASES ############################
####################################################################################

# Clase pbFactor ----
pbFactor <- R6Class("pbFactor", public = list(values = NA,
                                              name = "",
                                              unit = "",
                                              type = "numeric",
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
                                              },
                                              attributes = function(){
                                                v <- c(self$values,self$name, self$unit, self$type)
                                              }




)
)

# Clase pbDesign ----
pbDesign.c <- R6Class("pbDesign", public = list(name = NULL,
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
                                                  print(format(self$as.data.frame(), digits = 4))
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
                                                  cat(paste("Plackett-Burman", toupper(self$designType), "Design"))
                                                  cat("\n")
                                                  cat("Information about the factors:\n\n")
                                                  self$.nfp
                                                  cat("\n")
                                                  cat("-----------\n")
                                                  cat("\n")
                                                  print(self$as.data.frame())
                                                  cat("\n")
                                                  cat("-----------\n")
                                                  cat("\n")
                                                }


)
)


