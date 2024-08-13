facDesign <- function (k = 3, p = 0, replicates = 1, blocks = 1, centerCube = 0, random.seed = 1234)
{
  #' @title facDesign
  #' @description Generates a 2^k full factorial design.
  #' @param k numeric value giving the number of factors. By default k is set to ‘3’.
  #' @param p numeric integer between ‘0’ and ‘7’. p is giving the number of additional factors in the response surface design by aliasing effects.
  #' For further information see fracDesign and fracChoose.
  #' By default p is set to ‘0’.
  #' @param replicates numeric value giving the number of replicates per factor combination. By default replicates is set to ‘1’.
  #' @param blocks numeric value giving the number of blocks. By default blocks is set to ‘1’. Blocking is only performed for k greater 2.
  #' @param centerCube numeric value giving the number of centerpoints within the 2^k design. By default centerCube is set to ‘0’.
  #' @return The function facDesign returns an object of class facDesign.
  
  frameOut = fracDesign(k = k, p = p, gen = NULL, replicates = replicates,
                        blocks = blocks, centerCube = centerCube, random.seed = random.seed)
  return(frameOut)
}

### clase facDesign.c######################################################

#' @title facDesign Class
#' @description facDesign Class
#' @field name name of de facDesign
#' @field factors description
#' @field cube description
#' @field star description
#' @field centerCube description
#' @field centerStar description
#' @field generator description
#' @field response description
#' @field block description
#' @field blockGen description
#' @field runOrder description
#' @field standardOrder description
#' @field desireVal description
#' @field desirability description
#' @field fits description
facDesign.c <- R6Class("facDesign", public = list(name = NULL,
                                                  factors = NULL,
                                                  cube = data.frame(),
                                                  star = data.frame(),
                                                  centerCube = data.frame(),
                                                  centerStar = data.frame(),
                                                  generator = NULL,
                                                  response = data.frame(),
                                                  block = data.frame(),
                                                  blockGen = data.frame(),
                                                  runOrder = data.frame(),
                                                  standardOrder = data.frame(),
                                                  desireVal = NULL,
                                                  desirability = list(),
                                                  fits = NULL,
                                                  
                                                  #' @description Get the number of row Design
                                                  nrow = function(){
                                                    nrow(self$as.data.frame())
                                                  },
                                                  
                                                  ncol = function(){
                                                    ncol(self$as.data.frame())
                                                  },
                                                  
                                                  print = function(){
                                                    runIndex = order(self$runOrder[,1])
                                                    print(format(self$as.data.frame(), digits = 4))
                                                    invisible(self$as.data.frame())
                                                  },
                                                  
                                                  .clear = function(){
                                                    self$standardOrder = data.frame()
                                                    self$runOrder = data.frame()
                                                    self$cube = data.frame()
                                                    self$centerStar = data.frame()
                                                    self$centerCube = data.frame()
                                                    self$star = data.frame()
                                                    self$block = data.frame()
                                                    self$blockGen = data.frame()
                                                    self$response = data.frame()
                                                    invisible(self)
                                                  },
                                                  
                                                  #' @description Set the factors names
                                                  #' @param value factors names
                                                  names = function(value){
                                                    if(missing(value)){
                                                      n <- c()
                                                      for (i in 1:length(self$factors)) {
                                                        n[i] <- self$factors[[i]]$name
                                                      }
                                                      return(n)
                                                    }
                                                    else {
                                                      for (i in 1:length(self$factors)){
                                                        self$factors[[i]]$name = as.character(value[i])
                                                        
                                                      }
                                                      invisible(self)
                                                    }
                                                    
                                                  },
                                                  
                                                  