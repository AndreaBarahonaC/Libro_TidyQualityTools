##################################################################
##################### gageRR - FUNCIONES #########################
##################################################################

# gageRRDesign ----
gageRRDesign = function(Operators = 3, Parts = 10, Measurements = 3, method = "crossed", sigma = 6, randomize = TRUE) {
  # ValidaciC3n de argumentos
  if (!is.numeric(sigma))
    stop("sigma needs to be numeric")
  if (method != "nested" && method != "crossed")
    stop("Unknown method specified. Use 'method = nested' or 'method = crossed'.")
  Measurements <- as.integer(Measurements)
  if (!is.numeric(Measurements) || Measurements <= 0)
    stop("Number of Measurements per Part must be a positive integer.")


  opvec <- factor()
  partvec <- factor()

  yName <- "Measurement"
  aName <- "Operator"
  bName <- "Part"
  abName <- "Operator:Part"

  Operators <- unique(Operators)
  Parts <- unique(Parts)

  if (is.numeric(Operators))
    opvec <- factor(LETTERS[1:Operators[1]])
  if (is.character(Operators))
    opvec <- factor(Operators)

  if (length(unique(opvec)) > 26)
    stop("Too many Operators!")
  if (length(unique(opvec)) < 2)
    stop("Not enough Operators")

  if (is.numeric(Parts))
    partvec <- factor(LETTERS[1:Parts[1]])
  if (is.character(Parts))
    partvec <- factor(Parts)

  if (length(unique(partvec)) > 26)
    stop("Too many Parts!")
  if (length(unique(partvec)) < 2)
    stop("Too few Parts")

  Measurement <- rep(NA, (length(opvec) * length(partvec) * Measurements))
  outFrame <- data.frame()

  if (method == "crossed") {
    temp <- expand.grid(opvec, partvec)
    o <- rep(temp[, 1], Measurements)
    p <- rep(temp[, 2], Measurements)
  } else {
    p <- rep(sort(rep(partvec, length(opvec))), Measurements)
    o <- (rep(opvec, length(Measurement) / length(opvec)))
    p <- p[order(o,p)]
    o <- o[order(o,p)]
  }

  if (randomize)
    outFrame <- data.frame(StandardOrder = 1:length(Measurement), RunOrder = sample(1:length(Measurement), length(Measurement)), Operator = factor(o), Part = factor(p), Measurement)
  else
    outFrame <- data.frame(StandardOrder = 1:length(Measurement), RunOrder = 1:length(Measurement), Operator = factor(o), Part = factor(p), Measurement)

  outFrame <- outFrame[order(outFrame$RunOrder), ]
  # Valores predeterminados
  gageRRObj <- gageRR$new(
    X = outFrame,
    ANOVA = NULL,
    RedANOVA = NULL,
    method = method,
    Estimates = NULL,
    Varcomp = NULL,
    Sigma = sigma,
    GageName = NULL,
    GageTolerance = NULL,
    DateOfStudy = Sys.Date(),
    PersonResponsible = NULL,
    Comments = NULL,
    b = factor(p),
    a = factor(o),
    y = as.numeric(Measurement),
    facNames = c(yName, aName, bName, abName),
    numO = length(unique(opvec)),  # NC:mero de operadores
    numP = length(unique(partvec)),  # NC:mero de partes
    numM = Measurements  # NC:mero de mediciones
  )

  return(gageRRObj)
}

# gageRR ----
gageRR = function(gdo, method = "crossed", sigma = 6, alpha = 0.25, DM = NULL, HM = NULL, tolerance = NULL, dig = 3, ...) {
  method <- method

  yName <- "Measurement"
  aName <- "Operator"
  bName <- "Part"

  abName <- if(method == "crossed") paste(aName, ":", bName, sep = "")
  else if(method == "nested") paste(bName, "(", aName, ")", sep = "")
  else NA

  bTobName <- paste(bName, "to", bName, sep = " ")

  if (is.null(tolerance)) tolerance <- gdo$get.tolerance()

  y <- gdo$X[[yName]]
  a <- gdo$X[[aName]]
  b <- gdo$X[[bName]]

  nestedFormula <- as.formula(paste(yName, "~", aName, "/", bName))
  crossedFormula <- as.formula(paste(yName, "~", aName, "*", bName))
  reducedFormula <- as.formula(paste(yName, "~", aName, "+", bName))

  if (method == "nested") {
    numA <- nlevels(a)
    numB <- nlevels(b)
    numMPP <- length(y) / (numB * numA)

    gdo$numO <- numA
    gdo$numP <- numB
    gdo$numM <- numMPP

    fit <- aov(nestedFormula, data = gdo$X)
    meanSq <- anova(fit)[, 3]

    gdo$ANOVA <- fit
    gdo$method <- "nested"

    MSa <- meanSq[1]
    MSab <- meanSq[2]
    MSe <- meanSq[3]

    Cerror <- MSe
    Cb <- (MSab - MSe) / numMPP
    Ca <- (MSa - MSab) / (numB * numMPP)

    if (Ca <= 0) Ca <- 0
    if (Cb <= 0) Cb <- 0

    Cab <- 0
    totalRR <- Ca + Cab + Cerror
    repeatability <- Cerror
    reproducibility <- Ca
    bTob <- Cb
    totalVar <- Cb + Ca + Cab + Cerror

    estimates <- list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
    varcomp <- list(totalRR = totalRR, repeatability = repeatability, reproducibility = reproducibility, bTob = bTob, totalVar = totalVar)

    gdo$Estimates <- estimates
    gdo$Varcomp <- varcomp
  }

  if (method == "crossed") {
    numA <- nlevels(a)
    numB <- nlevels(b)
    numMPP <- length(a) / (numA * numB)

    gdo$numO <- numA
    gdo$numP <- numB
    gdo$numM <- numMPP

    fit <- aov(crossedFormula, data = gdo$X)
    model <- anova(fit)

    gdo$ANOVA <- fit
    gdo$method <- "crossed"

    MSb <- MSa <- MSab <- MSe <- 0

    if (bName %in% row.names(model)) MSb <- model[bName, "Mean Sq"]
    else warning(paste("missing factor", bName, "in model"))

    if (aName %in% row.names(model)) MSa <- model[aName, "Mean Sq"]
    else warning(paste("missing factor", aName, "in model"))

    if (abName %in% row.names(model)) MSab <- model[abName, "Mean Sq"]
    else warning(paste("missing interaction", abName, "in model"))

    if ("Residuals" %in% row.names(model)) MSe <- model["Residuals", "Mean Sq"]
    else warning("missing Residuals in model")

    Cb <- Ca <- Cab <- Cerror <- 0

    Cb <- (MSb - MSab) / (numA * numMPP)
    Ca <- (MSa - MSab) / (numB * numMPP)
    Cab <- (MSab - MSe) / numMPP
    Cerror <- (MSe)

    gdo$RedANOVA <- gdo$ANOVA

    if ((Cab < 0) || (model[abName, "Pr(>F)"] >= alpha)) {
      redFit <- aov(reducedFormula, data = gdo$X)
      model <- anova(redFit)

      MSb <- MSa <- MSab <- MSe <- 0

      if (bName %in% row.names(model)) MSb <- model[bName, "Mean Sq"]
      else warning(paste("missing factor", bName, "in model"))

      if (aName %in% row.names(model)) MSa <- model[aName, "Mean Sq"]
      else warning(paste("missing factor", aName, "in model"))

      if ("Residuals" %in% row.names(model)) MSe <- model["Residuals", "Mean Sq"]
      else warning("missing Residuals in model")

      Cb <- Ca <- Cab <- Cerror <- 0

      Cb <- (MSb - MSe) / (numA * numMPP)
      Ca <- (MSa - MSe) / (numB * numMPP)
      Cab <- 0
      Cerror <- (MSe)

      gdo$RedANOVA <- redFit
    }

    gdo$method <- "crossed"
    Ca <- max(0, Ca)
    Cb <- max(0, Cb)
    Cab <- max(0, Cab)

    totalRR <- Ca + Cab + Cerror
    repeatability <- Cerror
    reproducibility <- Ca + Cab
    bTob <- max(0, Cb)
    totalVar <- Cb + Ca + Cab + Cerror

    estimates <- list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
    varcomp <- list(totalRR = totalRR, repeatability = repeatability, reproducibility = reproducibility, a = Ca, a_b = Cab, bTob = bTob, totalVar = totalVar)

    gdo$Estimates <- estimates
    gdo$Varcomp <- varcomp
  }

  cat("\n")
  cat(paste("AnOVa Table - ", gdo$method, "Design\n"))
  print(summary(gdo$ANOVA))
  cat("\n")
  cat("----------\n")

  if (!identical(gdo$RedANOVA, gdo$ANOVA) && gdo$method == "crossed") {
    cat(paste("AnOVa Table Without Interaction - ", gdo$method, "Design\n"))
    print(summary(gdo$RedANOVA))
    cat("\n")
    cat("----------\n")
  }

  Source <- names(gdo$Varcomp)
  Source[Source == "repeatability"] <- " repeatability"
  Source[Source == "reproducibility"] <- " reproducibility"
  Source[Source == "a_b"] <- paste("  ", abName)
  Source[Source == "a"] <- paste("  ", aName)
  Source[Source == "bTob"] <- bTobName

  VarComp <- round(as.numeric(gdo$Varcomp[c(1:length(gdo$Varcomp))]), 3)
  Contribution <- round(as.numeric(gdo$Varcomp[c(1:length(gdo$Varcomp))]) / as.numeric(gdo$Varcomp[length(gdo$Varcomp)]), 3)
  VarComp <- t(data.frame(gdo$Varcomp))
  VarCompContrib <- VarComp / gdo$Varcomp$totalVar
  Stdev <- sqrt(VarComp)
  StudyVar <- Stdev * gdo$Sigma
  StudyVarContrib <- StudyVar / StudyVar["totalVar", ]
  SNR <- 1
  ptRatio <- NULL
  temp <- NULL

  if ((length(gdo$GageTolerance) > 0) && (gdo$GageTolerance > 0)) {
    ptRatio <- StudyVar / gdo$GageTolerance
    temp <- data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib, ptRatio)
    names(temp)[6] <- c("P/T Ratio")
    row.names(temp) <- c(Source)
  } else {
    temp <- data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib)
    row.names(temp) <- c(Source)
  }

  cat("\n")
  cat("Gage R&R\n")
  tempout <- temp
  print(format(tempout, digits = dig))
  cat("\n")
  cat("---\n")
  cat(" * Contrib equals Contribution in %\n")

  SNRTemp <- sqrt(2) * (temp[bTobName, "Stdev"] / temp["totalRR", "Stdev"])
  if (SNRTemp > 1) SNR <- SNRTemp

  cat(paste(" **Number of Distinct Categories (truncated signal-to-noise-ratio) =", floor(SNR), "\n"))
  cat("\n")
  invisible(gdo)
}

# Ejemplo de uso: ----
# crear un objeto de la clase 'gageRR'
mi_gageRR <- gageRR.c$new(
  X = data.frame(
    Operator = factor(c("A", "B", "C", "A", "B")),
    Part = factor(c("P1", "P1", "P2", "P2", "P3")),
    Measurement = c(10, 12, 11, 13, 9)
  ),
  ANOVA = NULL,  # Esto puede ser NULL inicialmente y luego calcularlo
  RedANOVA = NULL,  # Igual que ANOVA, puede ser NULL inicialmente
  method = "crossed",
  Estimates = list(),
  Varcomp = list(),
  Sigma = 0.5,
  GageName = "Gage1",
  GageTolerance = 0.1,
  DateOfStudy = "2024-05-15",
  PersonResponsible = "John Doe",
  Comments = "Sample gage R&R study",
  b = factor(c("A", "A", "B", "B", "C")),
  a = factor(c("P1", "P1", "P2", "P2", "P3")),
  y = c(10, 12, 11, 13, 9),
  facNames = c("Measurement", "Operator", "Part"),
  numO = 3,
  numP = 3,
  numM = 2
)

# Crear el objeto gageRRObj
design_example <- gageRRDesign(
  Operators = 3,
  Parts = 10,
  Measurements = 3,
  method = "crossed",
  sigma = 6,
  randomize = TRUE
)

# Crear un diseño para el estudio de Gage
design <- gageRRDesign(Operators = 3, Parts = 10, Measurements = 3, method = "crossed", sigma = 6, randomize = TRUE)
design$X$Measurement <- rnorm(nrow(design$X), mean = 10, sd = 2)

# Ejecutar la función gageRR_
result <- gageRR(
  gdo = design,
  method = "crossed",   #  método "crossed"
  sigma = 6,            #  sigma
  alpha = 0.25,         # Nivel de significancia
  tolerance = NULL,     # Tolerancia
  dig = 3               # Número de dígitos a mostrar en los resultados
)
class(result)
result$plot()
