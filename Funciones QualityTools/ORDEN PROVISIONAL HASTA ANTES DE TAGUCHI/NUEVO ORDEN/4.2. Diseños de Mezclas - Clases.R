######################################################################
################ DISEÑOS DE MEZCLAS - CLASES #########################
######################################################################

# mixDesign.c ----
mixDesign.c <- R6Class("mixDesign", public = list(name = NULL,
                                                  factors =list(),
                                                  total = NULL,
                                                  lower = NULL,
                                                  design = data.frame(),
                                                  designType = NULL,
                                                  pseudo = data.frame(),
                                                  response = data.frame(),
                                                  Type = data.frame(),
                                                  block = data.frame(),
                                                  runOrder = data.frame(),
                                                  standardOrder = data.frame(),
                                                  desireVal = list(),
                                                  desirability = list(),
                                                  fits = data.frame(),

                                                  .factors = function(value){
                                                    if (missing(value)) {
                                                      return(self$factors)
                                                    }
                                                    else{
                                                      if (length(value) != ncol(self$pseudo))
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

                                                  as.data.frame = function(row.names = NULL, optional = FALSE){
                                                    frameOut = cbind(self$standardOrder, self$runOrder, self$Type, self$pseudo, self$response)
                                                    return(frameOut)
                                                  },

                                                  print = function(){
                                                    print(format(self$as.data.frame(), digits = 4))
                                                    invisible(self$as.data.frame())
                                                  },

                                                  .response = function(value){
                                                    if (missing(value)) {
                                                      return(self$response)
                                                    }
                                                    else{
                                                      #print(deparse(substitute(value)))
                                                      if (!is.numeric(value) & !is.data.frame(value))
                                                        stop("vector or data.frame must be given")
                                                      if (is.numeric(value)) {
                                                        if (length(value) != nrow(self$pseudo))
                                                          stop("differing lengths")
                                                        temp = data.frame(value)
                                                        names(temp) = deparse(substitute(value))[1]
                                                        value = temp
                                                      }
                                                      if (is.data.frame(value)) {
                                                        if (nrow(value) != nrow(self$pseudo))
                                                          stop("differing number of rows")
                                                      }
                                                      self$response = value
                                                      invisible(self)
                                                    }

                                                  },

                                                  .nfp = function(){
                                                    x = self$.factors()
                                                    atr <- c('low','high','name','unit','type')
                                                    if (is.list(x) && length(x[[1]]) > 0) {
                                                      numAttr = length(x[[1]]$attributes())
                                                      .numFac = length(x)
                                                      frameOut = data.frame(matrix(ncol = .numFac, nrow = numAttr ))
                                                      for (i in 1:numAttr ) {
                                                        charVec = character(0)
                                                        for (j in 1:.numFac) {
                                                          charVec = c(charVec, atr[i], "\t\t")
                                                          frameOut[i, j] = x[[j]]$attributes()[i]
                                                        }
                                                      }
                                                      names(frameOut) = self$names()
                                                      rownames(frameOut) = atr[1:numAttr ]
                                                    }
                                                    else {
                                                      stop("no list given or length of list < 1")
                                                    }
                                                    print(frameOut)
                                                  },

                                                  summary = function(){
                                                    cat(paste("Simplex", toupper(self$designType), "Design"))
                                                    cat("\n")
                                                    cat("Information about the factors:\n\n")
                                                    self$.nfp()
                                                    cat("\n-----------\n")
                                                    cat("\n")
                                                    .npp(self)
                                                    cat("\n-----------\n")
                                                    cat("\n")
                                                    cat("Information about the constraints:\n\n")
                                                    lower = object$lower
                                                    temp = character(0)
                                                    for (i in seq(along = lower)) temp = c(temp, paste(LETTERS[i], ">=", lower[i]))
                                                    cat(temp)
                                                    cat("\n")
                                                    cat("\n-----------\n")
                                                    cat("\n")
                                                    times = nrow(self$pseudo)
                                                    pseudo = format(self$pseudo, digits = 2)
                                                    design = format(self$design, digits = 2)
                                                    amount = design
                                                    if (self$total[2] != 1)
                                                      amount = format(self$design * self$total[2], digits = 2)
                                                    temp = c("                             ", "PseudoComponent", "_|_", "Proportion", "_|_", "Amount")
                                                    cat(temp)
                                                    cat("\n")
                                                    cat("\n")
                                                    temp = cbind(pseudo, `_` = rep(" ", times = times), `|` = rep("|", times = times), `_` = rep(" ", times = times), design)
                                                    temp = cbind(temp, `_` = rep(" ", times = times), `|` = rep("|", times = times), `_` = rep(" ", times = times), amount)
                                                    temp = cbind(self$standardOrder, self$runOrder, self$Type, `|` = rep("|", times = times), temp, `|` = rep("|", times = times), self$response)
                                                    show(temp)
                                                    cat("\n-----------\n")
                                                    cat("\n")
                                                    cat(paste("Mixture Total:", self$total[1], "equals", self$total[2]))
                                                    cat("\n")
                                                    cat("\n")
                                                    invisible(self$as.data.frame())
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

                                                  lows = function(value){
                                                    if (missing(value)) {
                                                      listOut = vector(mode = "list")
                                                      for (i in seq(along = self$factors)) {
                                                        listOut[self$factors[[i]]$name] = self$factors[[i]]$.low()
                                                      }
                                                      return(listOut)
                                                    }
                                                    else {
                                                      for (i in seq(along = self$factors)) {
                                                        self$factors[[i]]$.low(value[i])
                                                      }
                                                      invisible(self)
                                                    }
                                                  },

                                                  highs = function(value){
                                                    if (missing(value)) {
                                                      listOut = vector(mode = "list")
                                                      for (i in seq(along = self$factors)) {
                                                        listOut[self$factors[[i]]$name] = self$factors[[i]]$.high()
                                                      }
                                                      return(listOut)
                                                    }
                                                    else {
                                                      for (i in seq(along = self$factors)) {
                                                        self$factors[[i]]$.high(value[i])
                                                      }
                                                      invisible(self)
                                                    }
                                                  }






)
)
