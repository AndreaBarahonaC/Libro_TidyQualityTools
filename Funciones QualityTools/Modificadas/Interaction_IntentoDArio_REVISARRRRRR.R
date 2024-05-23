###### Necesito InteractionPlot#########################################################
#INTENTO DARIO
library(ggplot2)
library(gridExtra)

#ORIGINAL
dfac <- facDesign(k = 3, centerCube = 4)
  interactionPlot <- function(fdo, y = NULL, response = NULL, fun = mean, main, col = 1:2, ...) { ###
    DB = FALSE
    mainmiss = FALSE
    if (missing(main))
      mainmiss = TRUE
    if (missing(fdo) || class(fdo)[1]!="facDesign")                                ###
      stop("fdo needs to be an object of class facDesign")                    ###
    parList = list(...)
    old.par <- par(no.readonly = TRUE)
    fdoName = deparse(substitute(fdo))                                          ###
    if(is.null(response)==FALSE)                                                ###
    {                                                                           ###
      temp=fdo$.response()[response]                                               ###
      fdo$.response(temp)                                                         ###
    }                                                                           ###
    diagNames = character(0)
    x = fdo$cube
    runIndex = order(fdo$runOrder[,1])
    x = x[runIndex[1:nrow(x)], ]
    y = fdo$.response()[1:nrow(x), ]
    numFac = ncol(x)
    combMat = combn(names(x), 2)
    if (numFac == 2) {
      facName2 = combMat[1, 1]
      facName1 = combMat[2, 1]
      temp = with(cbind(y, x), .testFun(eval(parse(text = facName2)), eval(parse(text = facName1)), xlab = facName1, response = y, trace.label = facName1,
                                        ylim = range(y), axes = F, fun, title = facName1, col = col, ...))
      tempList = parList
      tempList$col = 1
      tempList$lwd = 1
      tempList$side = c(2)
      do.call(axis, tempList)
      box()
      tempList$at = temp[[1]]
      tempList$labels = temp[[2]]
      tempList$side = c(1)
      do.call(axis, tempList)
      invisible()
    }

    for (r in 1:ncol(fdo$.response())) {
      if (r > 1)
        dev.new()
      y = fdo$.response()[1:nrow(x), r]
      par(mfrow = c(numFac, numFac))
      par(mar = c(0, 0, 0, 0))
      par(oma = c(0, 0, 8, 8))
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
      for (i in 1:ncol(combMat)) {
        facName1 = combMat[1, i]
        facName2 = combMat[2, i]
        rowNum = .letterPos(facName1)
        colNum = .letterPos(facName2)
        if (DB) {
          print(numFac)
          print(rowNum)
          print(colNum)
          print(facName1)
          print(facName2)
          print(y)
          print(fun)
          cat(paste(i, "\t", c(rowNum, colNum)))
        }
        par(mfg = c(rowNum, colNum))
        temp = with(cbind(x, y), .testFun(eval(parse(text = facName2)), eval(parse(text = facName1)), response = y, trace.label = facName1, ylim = range(y),
                                          axes = F, fun = fun, title = facName1, col = col, ...))
        if (colNum == numFac) {
          tempList = parList
          tempList$col = 1
          tempList$lwd = 1
          tempList$side = 4
          do.call(axis, tempList)
        }
        if (rowNum == 1) {
          tempList = parList
          tempList$col = 1
          tempList$lwd = 1
          tempList$side = 3
          tempList$at = temp[[1]]
          tempList$labels = temp[[2]]
          do.call(axis, tempList)
        }
        box(which = "plot")
        par(mfg = c(rowNum, rowNum))
        plot(c(-1, 1), c(-1, 1), type = "n", axes = F, xlab = "", ylab = "", main = "")
        text(0, 0, facName1, cex = 4)
        diagNames = c(diagNames, facName1)
      }
      par(mfg = c(numFac, numFac))
      plot(c(-1, 1), c(-1, 1), type = "n", axes = F, xlab = "", ylab = "", main = "")
      text(0, 0, setdiff(names(x), diagNames), cex = 4)
      if (mainmiss) {
        main = paste("Interaction plot for", names(fdo$.response())[r], "in", fdoName)         ###
        title(main, outer = T, ...)
      }
      else title(main[r], outer = T, ...)
    }
    #    }
    par(old.par)
    invisible()
  }

#INTENTO 1
interactionPlot <- function(fdo, y = NULL, response = NULL, fun = mean, main, col = 1:2, ...) {
  if (missing(main)) mainmiss = TRUE else mainmiss = FALSE
  if (missing(fdo) || class(fdo)[1] != "facDesign")
    stop("fdo needs to be an object of class facDesign")

  fdoName = deparse(substitute(fdo))
  if (!is.null(response)) {
    temp = fdo$.response()[response]
    fdo$.response(temp)
  }

  x <- fdo$cube
  runIndex <- order(fdo$runOrder[, 1])
  x <- x[runIndex[1:nrow(x)], ]
  y <- fdo$.response()[1:nrow(x), ]
  numFac <- ncol(x)
  combMat <- combn(names(x), 2)

  plot_list <- list()

  for (r in 1:ncol(fdo$.response())) {
    y = fdo$.response()[1:nrow(x), r]

    for (i in 1:ncol(combMat)) {
      facName1<- combMat[1, i]
      facName2 <- combMat[2, i]

      df = data.frame(fac1 = x[[facName1]], fac2 = x[[facName2]], response = y)

      p = ggplot(df, aes_string(x = "fac1", y = "response", color = "as.factor(fac2)")) +
        geom_line(stat = "summary", fun = fun) +
        labs(x = facName1, y = "Response", color = facName2) +
        theme_minimal()

      plot_list[[paste(facName1, facName2)]] = p
    }

    if (mainmiss) {
      main = paste("Interaction plot for", names(fdo$.response())[r], "in", fdoName)
    }

    plot_matrix <- vector("list", numFac * numFac)
    plot_idx <- 1

    for (j in 1:numFac) {
      for (i in 1:numFac) {
        if (i == j) {
          p_diag <- ggplot() +
            geom_text(aes(0.5, 0.5, label = names(x)[i]), size = 8) +
            theme_void()
          plot_matrix[[plot_idx]] <- p_diag
        } else if (i < j) {
          plot_matrix[[plot_idx]] <- plot_list[[paste(names(x)[i], names(x)[j])]]
        } else {
          plot_matrix[[plot_idx]] <- ggplot() + theme_void()  # Empty plot
        }
        plot_idx <- plot_idx + 1
      }
    }

    plot_matrix <- matrix(plot_matrix, ncol = numFac, byrow = TRUE)
    final_plot <- wrap_plots(plot_matrix, ncol = numFac)
    final_plot <- final_plot + plot_annotation(title = main[r])
    print(final_plot)
  }

  invisible()
}



#INTENTO 2

interactionPlotGG <- function(fdo, y = NULL, response = NULL, fun = mean, main = NULL, col = c("red", "blue"), ...) {
  if (missing(fdo) || class(fdo)[1] != "facDesign") {
    stop("fdo needs to be an object of class facDesign")
  }

  if (!is.null(response)) {
    temp <- fdo$.response()[response]
    fdo$.response(temp)
  }

  x <- fdo$cube
  runIndex <- order(fdo$runOrder[, 1])
  x <- x[runIndex[1:nrow(x)], ]
  y <- fdo$.response()[1:nrow(x), ]

  # Convertir todas las variables en factores
  x <- as.data.frame(lapply(x, factor))

  numFac <- ncol(x)
  combMat <- combn(names(x), 2)

  plot_list <- matrix(list(), nrow = numFac, ncol = numFac)

  # Crear gráficos de interacción y ubicarlos en la matriz
  for (i in 1:ncol(combMat)) {
    facName1 <- combMat[1, i]
    facName2 <- combMat[2, i]

    data <- cbind(y, x)

    # Remover valores no finitos
    data <- data[complete.cases(data), ]

    p <- ggplot(data, aes_string(x = facName2, y = "y", group = facName1, color = facName1)) +
      stat_summary(fun = fun, geom = "line") +
      labs(title = NULL, x = NULL, y = NULL) +
      theme_minimal() +
      scale_color_manual(values = col) +
      theme(legend.position = "none")

    # Ubicar en la posición correcta
    pos1 <- which(names(x) == facName1)
    pos2 <- which(names(x) == facName2)
    plot_list[[pos1, pos2]] <- p
  }

  # Llenar la diagonal con los nombres de los factores
  for (i in 1:numFac) {
    plot_list[[i, i]] <- ggplot() +
      geom_text(aes(0.5, 0.5, label = names(x)[i]), size = 10) +
      theme_void()
  }

  # Crear gráficos vacíos en la parte inferior de la diagonal
  for (i in 1:numFac) {
    for (j in 1:numFac) {
      if (i > j) {
        plot_list[[i, j]] <- ggplot() + theme_void()
      } else if (i < j) {
        plot_list[[i, j]] <- plot_list[[i, j]] + theme(legend.position = "right")
      }
    }
  }

  plot_grid <- wrap_plots(plot_list, ncol = numFac)

  if (is.null(main)) {
    main <- paste("Interaction plot for", names(fdo$.response())[1], "in", deparse(substitute(fdo)))
  }

  plot_grid <- plot_grid + plot_annotation(title = main)

  print(plot_grid)
}


interactionPlotGG(dfac)

interactionPlot(dfac)
class(dfac)
