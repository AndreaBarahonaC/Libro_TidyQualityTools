###############################################################################
####################### DISEÑO TAGUCHI - FUNCIONES ############################
###############################################################################

# Funcion taguchiDesign ----
taguchiDesign <- function(design, randomize = TRUE, replicates = 1) {

  #' @title taguchiDesign: Taguchi Designs
  #' @description Function to create a taguchi design.
  #' @param design A character string specifying the orthogonal array of the Taguchi design. The available options are:
  #' \itemize{
  #'   \item {"L4_2" for three two-level factors.}
  #'   \item {"L8_2" for seven two-level factors.}
  #'   \item {"L9_3" for four three-level factors.}
  #'   \item {"L12_2" for 11 two-level factors.}
  #'   \item {"L16_2" for 16 two-level factors}
  #'   \item {"L16_4" for 16 four-level factors.}
  #'   \item {"L18_2_3" for one two-level and seven three-level factors.}
  #'   \item {"L25_5" for six five-level factors.}
  #'   \item {"L27_3" for 13 three-level factors.}
  #'   \item {"L32_2" for 32 two-level factors.}
  #'   \item {"L32_2_4" for one two-level factor and nine four-level factors.}
  #'   \item {"L36_2_3_a" for 11 two-level factors and 12 three-level factors.}
  #'   \item {"L36_2_3_b" for three two-level factors and 13 three-level factors.}
  #'   \item {"L50_2_5" for one two-level factor and eleven five-level factors.}
  #'   \item {"L8_4_2" for one four-level factor and four two-level factors.}
  #'   \item {"L16_4_2_a" for one four-level factor and 12 two-level factors.}
  #'   \item {"L16_4_2_b" for two four-level factors and nine two-level factors.}
  #'   \item {"L16_4_2_c" for three four-level factors and six two-level factors.}
  #'   \item {"L16_4_2_d" for five four-level factors and two two-level factors.}
  #'   \item {"L18_6_3" for one six-level factor and six three-level factors.}
  #' }
  #' @param randomize A logical value (`TRUE`/`FALSE`) that specifies whether to randomize the RunOrder of the design. By default, `randomize` is set to `TRUE`.
  #' @param replicates An integer specifying the number of replicates for each run in the design.
  #' @return A `taguchiDesign` returns an object of class `taguchiDesign`.
  #' @details An overview of possible taguchi designs is possible with `taguchiChoose`.
  #' @seealso
  #'  \item{\code{\link{facDesig}}}{for 2^k factorial designs.}
  #' \item{\code{\link{rsmDesign}}}{for response surface designs.}
  #' \item{\code{\link{fracDesig}}}{for fractional factorial design.}
  #' \item{\code{\link{pbDesig}}}{for response surface designs.}
  #' \item{\code{\link{gageRRDesig}}}{for gage designs.}
  #' @examples
  #' set.seed(1234)
  #' tdo <- taguchiDesign("L9_3")
  #' tdo$values(list(A = c(20, 40, 60), B = c("material 1", "material 2", "material 3"), C = c(1, 2, 3)))
  #' tdo$names(c("Factor 1", "Factor 2", "Factor 3", "Factor 4"))
  #' tdo$summary()
  #' tdo$.response(rnorm(9))
  #'  tdo$effectPlot(col = 2)

  odo = NA
  type = "single"
  for (i in seq(along = .oaList)) {
    pmatch(design, .oaList[[i]]$id)
    if (!is.na(pmatch(design, .oaList[[i]]$id))) {

      temp = .oaList[[i]]
      design = temp$design
      repVec = rep(1, nrow(design))
      if (replicates > 1) {
        X = temp$design
        for (i in 1:(replicates - 1)) {
          design = rbind(design, X)
          repVec = c(repVec, rep(i + 1, times = nrow(X)))
        }
      }
      Replicate = data.frame(Replicate = as.numeric(repVec))

      odo = taguchiDesign.c$new()
      odo$design = design
      names(odo$design) = .NAMES[1:ncol(design)]
      odo$name = temp$id
      odo$designType = temp$type
      odo$replic = Replicate
      StandOrder = 1:nrow(odo$design)
      RunOrder = StandOrder
      if (randomize) {
        RunOrder = sample(1:nrow(odo$design), nrow(odo$design), replace = FALSE, prob = NULL)
      }
      odo$design = odo$design[order(RunOrder), ]
      odo$replic = data.frame(Replicate = odo$replic[order(RunOrder), 1])
      row.names(odo$design) = odo$design$RunOrder
      odo$runOrder = data.frame(RunOrder = data.frame(RunOrder = RunOrder)[order(RunOrder), ])
      odo$standardOrder = data.frame(StandOrder = data.frame(StandOrder = StandOrder)[order(RunOrder), ])
      odo$response = data.frame(y = rep(NA, nrow(odo$design)))
      tfList = vector("list", ncol(design))
      for (i in seq(along = tfList)) tfList[[i]] = taguchiFactor$new()
      names(tfList) = names(odo$design)
      odo$.factors(tfList)
      valList = list(length = length(odo$names()))
      for (i in names(odo$names())) valList[[i]] = sort(unique(odo$design[, i]))
      odo$values(valList)
      return(odo)
    }
  }
  return(NA)
}

# Funcion oaChoose ----
oaChoose <- function(factors1, factors2, level1, level2, ia) {
  #' @title oaChoose: Taguchi Designs
  #' @description Shows a matrix of possible taguchi designs.
  #' @param factors1 Number of factors on level1.
  #' @param factors2 Number of factors on level2.
  #' @param level1 Number of levels on level1.
  #' @param level2 Number of levels on level2.
  #' @param ia Number of interactions.
  #' @details `oaChoose` returns possible taguchi designs. Specifying the number of factor1 factors with level1 levels (factors1 = 2, level1 = 3 means 2 factors with 3 factor levels) and factor2 factors with level2 levels and desired interactions one or more taguchi designs are suggested.
  #' If all parameters are set to ‘0’, a matrix of possible taguchi designs is shown.
  #' @return `oaChoose` returns an object of class `taguchiDesign`.
  #' @seealso
  #'  \item{\code{\link{facDesig}}}{for 2^k factorial designs.}
  #' \item{\code{\link{rsmDesign}}}{for response surface designs.}
  #' \item{\code{\link{fracDesig}}}{for fractional factorial design.}
  #' \item{\code{\link{gageRRDesig}}}{for gage designs.}

  params = list(factors1 = 0, factors2 = 0, level1 = 0, level2 = 0, ia = 0)
  if (!missing(ia))
    params$ia = ia
  if (!missing(factors2))
    params$factors2 = factors2
  if (!missing(level2))
    params$level2 = level2
  do.call(taguchiChoose, params)
}

# Funcion taguchiChoose ----
taguchiChoose <- function(factors1 = 0, factors2 = 0, level1 = 0, level2 = 0, ia = 0) {
  #' @title taguchiChoose: Taguchi Designs
  #' @description Shows a matrix of possible taguchi designs
  #' @param factors1 Integer number of factors on level1. By default set to ‘0’.
  #' @param factors2 Integer number of factors on level2. By default set to ‘0’.
  #' @param level1 Integer number of levels on level1.
  #' @param level2 Integer number of levels on level2. By default set to ‘0’.
  #' @param ia Integer number of interactions. By default set to ‘0’.
  #' @details `taguchiChoose` returns possible taguchi designs.
  #' Specifying the number of factor1 factors with level1 levels (factors1 = 2, level1 = 3 means 2 factors with 3 factor levels) and factor2 factors with level2 levels and desired interactions one or more taguchi designs are suggested.
  #' If all parameters are set to 0, a matrix of possible taguchi designs is shown.
  #' @return `taguchiChoose` returns an object of class `taguchiDesign`.
  #' @seealso
  #'  \item{\code{\link{facDesig}}}{for 2^k factorial designs.}
  #' \item{\code{\link{rsmDesign}}}{for response surface designs.}
  #' \item{\code{\link{fracDesig}}}{for fractional factorial design.}
  #' \item{\code{\link{gageRRDesig}}}{for gage designs.}

  if (factors1 == 0 & factors2 == 0 & level1 == 0 & level2 == 0 & ia == 0) {
    temp = vector(mode = "character", length = length(.oaList))
    for (i in 1:length(.oaList)) temp[i] = .oaList[[i]]$id
    temp = c(temp, rep(" ", (length(temp)%/%6 + 1) * 6 - length(temp)))
    mat = data.frame(matrix(temp, ncol = 6, byrow = TRUE))
    names(mat) = rep(" ", ncol(mat))
    print(mat)
    cat("\n")
    cat("Choose a design using e.g. taguchiDesign(\"L4_2\")")
    cat("\n")
  }
  else {

    if (factors2 <= 0)
      level2 = 0

    Anzahl_Spalten = factors1 + factors2 + ia
    ss = list()
    for (i in seq(along = .oaList)) {
      li = .oaList[[i]]
      if (li$factors1 >= factors1 & li$factors2 >= factors2 & (li$levels1 == level1 | li$levels1 == level2) & (li$levels2 == level2 | li$levels2 == level1) &
          li$anzahl_spalten >= Anzahl_Spalten)
        ss[i] = li$id
    }
    out = as.character(ss)
    out = out[out != "NULL"]
    if (length(out) > 0) {
      cat(paste(factors1, "factors on", level1, "levels and", factors2, "factors on", level2, "levels with", ia, "desired interactions to be estimated\n"))
      cat("\n")
      cat("Possible Designs:\n")
      cat("\n")
      cat(paste(out, sep = " | "))
      cat("\n")
      cat("\n")
      cat(paste("Use taguchiDesign(\"", out[1], "\") or different to create a taguchi design object\n", sep = ""))
    }
    else {
      cat("No Design Found\n")
      cat("\n")
      out = NA
    }
    invisible(out)
  }
}

#Arreglar taguchiChoose para que sirva oaChoose####
# taguchiChoose()
# oaChoose()
