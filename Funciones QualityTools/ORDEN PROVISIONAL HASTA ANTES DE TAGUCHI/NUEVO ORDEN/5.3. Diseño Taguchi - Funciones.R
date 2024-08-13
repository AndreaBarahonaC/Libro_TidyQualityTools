###############################################################################
####################### DISEÑO TAGUCHI - FUNCIONES ############################
###############################################################################

# Funcion taguchiDesign ----
taguchiDesign <- function(design, randomize = TRUE, replicates = 1) {
  DB = FALSE
  odo = NA
  type = "single"
  for (i in seq(along = .oaList)) {
    pmatch(design, .oaList[[i]]$id)
    if (!is.na(pmatch(design, .oaList[[i]]$id))) {
      if (DB)
        print(design)
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
      if (DB)
        print(Replicate)
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
###ejemplo uso taguchiDesign####
set.seed(1234)
tdo <- taguchiDesign("L9_3")
tdo$values(list(A = c(20, 40, 60), B = c("material 1", "material 2", "material 3"), C = c(1, 2, 3)))
tdo$names(c("Factor 1", "Factor 2", "Factor 3", "Factor 4"))
tdo$summary()
tdo$.response(rnorm(9))
tdo$effectPlot(col = 2)
# Funcion oaChoose ----
oaChoose <- function(factors1, factors2, level1, level2, ia) {
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
    DB = FALSE
    if (factors2 <= 0)
      level2 = 0
    if (DB)
      print(level2)
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
taguchiChoose()
oaChoose()
