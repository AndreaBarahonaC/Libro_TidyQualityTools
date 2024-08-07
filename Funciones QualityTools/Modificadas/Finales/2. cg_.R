library(ggplot2)
library(patchwork)
library(gridExtra)


cg_RunChart <- function (x, target, tolerance, ref.interval, facCg, facCgk,
                         n = 0.2, col = "black", pch = 19,
                         xlim = NULL, ylim = NULL, main = "Run Chart",
                         conf.level = 0.95, cgOut = TRUE)
{
  if (missing(x)) {
    stop("x must be given as a vector")
  }

  if (missing(target)) {
    target <- mean(x)
    targetmissing <- FALSE
  } else {
    targetmissing <- TRUE
  }

  if (missing(ref.interval)) {
    ref.interval <- qnorm(0.99865) - qnorm(0.00135)
  }

  sd <- sd(x)
  mean <- mean(x)
  ref.ar <- qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval, mean, sd)

  if (missing(facCg)) {
    facCg <- 0.2
  }

  if (missing(facCgk)) {
    facCgk <- 0.1
  }

  if (missing(tolerance)) {
    width <- ref.ar/facCg
    tolerance <- numeric(2)
    tolerance[1] <- mean - width/2
    tolerance[2] <- mean + width/2
  }

  quant1 <- qnorm((1 - ref.interval)/2, mean, sd)
  quant2 <- qnorm(ref.interval + (1 - ref.interval)/2, mean, sd)

  if (length(tolerance) != 2) {
    stop("tolerance has wrong length")
  }

  if (missing(xlim)) {
    xlim <- c(0, length(x))
  }

  if (missing(ylim)) {
    ylim <- c(min(x, target - n/2 * (abs(diff(tolerance))), quant1, quant2),
              max(x, target + n/2 * (abs(diff(tolerance))), quant1, quant2))
  }

  if (missing(main)) {
    main <- "Run Chart"
  }

  Cg <- (facCg * tolerance[2]-tolerance[1])/ref.interval
  Cgk <- (facCgk * abs(target-mean(x))/(ref.interval/2))

  # Create a data frame for plotting
  df <- data.frame(x = x, y = x)

  # Add target line
  df$y_target <- target

  # Calculate the upper and lower control limits
  df$y_lower <- quant1
  df$y_upper <- quant2
  # Calculate the tolerance limits
  df$y_tolerance_lower <- tolerance[1]
  df$y_tolerance_upper <- tolerance[2]

  # Calculate the mean
  df$y_mean <- mean

  # Calculate the Cg and Cgk values
  df$Cg <- Cg
  df$Cgk <- Cgk

  # Create the ggplot
  # 1. Principal plot and target line
  p <- ggplot(df, aes(x = seq_along(x), y = x)) +
    geom_point(color = col, shape = pch) +
    geom_line(color = col, linetype = "solid") +
    scale_x_continuous(limits = c(xlim[1] - 0.05 * xlim[2], xlim[2]), expand = c(0, 0)) +
    labs(title = main, x = "Index", y = "x") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
    geom_hline(aes(yintercept = target)) # Linea target

  # 2. Red Plot (Lowess)
  p <- p + geom_smooth(method = "loess", color = "red", se = FALSE, span = 1.25, size = 0.25,)

  # 3. Green lines
  p <- p + geom_hline(aes(yintercept = mean), linetype = "dashed", color = "seagreen")+  #center line
    geom_hline(aes(yintercept = quant1), linetype = "dashed", color = "seagreen") + # Bottom line
    geom_hline(aes(yintercept = quant2), linetype = "dashed", color = "seagreen")   # Top line
  # 4. Xtar +- 0.1
  p <- p + geom_hline(yintercept = c(target + n/2 * (abs(diff(tolerance))), target - n/2 * (abs(diff(tolerance)))), color = "#012B78", linetype = "solid") # Agregar líneas
  #Label
  p <- p + scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis =
                                sec_axis(~ .,breaks = c(target,mean,quant1,quant2,target + n/2 * (abs(diff(tolerance))),
                                                        target - n/2 * (abs(diff(tolerance)))),
                                         labels=c("target",
                                                  expression(bar(x)),
                                                  substitute(x[a * b], list(a = round(((1 - ref.interval)/2) * 100, 3), b = "%")),
                                                  substitute(x[a * b], list(a = round(((ref.interval + (1 - ref.interval)/2)) * 100,3), b = "%")),
                                                  substitute(x[tar] + a, list(a = round(n/2, 4))),
                                                  substitute(x[tar] - a, list(a = round(n/2, 4)))))) +
    theme(axis.text.y.right = element_text(size = 15))

  # Label  Cg and Cgk
  if (cgOut == TRUE) {
    caja <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
      theme_bw() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank()
      ) +
      xlim(c(0.25,0.26)) + ylim(c(0.24, 0.31))

    caja <- caja +
      annotate('text', x = 0.25, y = 0.30,
               label = paste("Cg: ", round(Cg,digits = 6)),
               parse = TRUE, size = 3, hjust = 0) +
      annotate('text', x = 0.25, y = 0.25,
               label = paste("Cgk:", round(Cgk,digits = 6)),
               parse = TRUE, size = 3, hjust = 0)

    p <- p + inset_element(caja, left = 0.7, right = 1, top = 1, bottom = 0.75)
    suppressMessages(show(p))
  }
  else{
    suppressMessages(show(p))
  }
  invisible(list(Cg, Cgk))
}


cg_HistChart <- function (x, target, tolerance, ref.interval, facCg, facCgk,
                          n = 0.2, col, xlim, ylim, main, conf.level = 0.95, cgOut = TRUE)
{
  if (missing(x))
    stop("x must be given as a vector")
  if (missing(target)) {
    target = mean(x)
    targetmissing = FALSE
  }
  else targetmissing = TRUE
  if (missing(ref.interval))
    ref.interval = pnorm(3) - pnorm(-3)
  sd = sd(x)
  mean = mean(x)
  ref.ar = qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval,
                                                 mean, sd)
  if (missing(facCg))
    facCg = 0.2
  if (missing(facCgk))
    facCgk = 0.1
  if (missing(tolerance))
    warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
  if (missing(tolerance)) {
    width = ref.ar/facCg
    tolerance = numeric(2)
    tolerance[1] = mean(x) - width/2
    tolerance[2] = mean(x) + width/2
  }
  quant1 = qnorm((1 - ref.interval)/2, mean, sd)
  quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean,
                 sd)
  if (length(tolerance) != 2)
    stop("tolerance has wrong length")
  if (missing(col))
    col = "lightblue"
  if (missing(xlim))
    xlim = c(0, length(x))
  if (missing(ylim))
    ylim = c(min(x, target - n/2 * (abs(diff(tolerance))),
                 quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))),
                                      quant1, quant2))
  if (missing(main))
    main = paste("Histogram of", deparse(substitute(x)),
                 "- target")

  Cg <- (facCg * tolerance[2]-tolerance[1])/ref.interval
  Cgk <- (facCgk * abs(target-mean(x))/(ref.interval/2))

  # Calculos previos
  x.c <- x - target
  temp <- hist(x.c, plot = FALSE)
  # Obtenemos la información para el histograma
  df <- data.frame(
    mid = temp$mids,
    density = temp$density
  )
  width <- diff(df$mid)[1] # Ancho de cada barra
  # Histograma
  p <- ggplot(df, aes(x = mid, y = density)) +
    geom_bar(stat = "identity", width = width, fill = "lightblue", color = "black", alpha = 0.5) +
    labs(y = "Density", x = "x.c", title = main) +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
    guides(color = guide_legend(title.position = "top", title.hjust = 0.5))
  # Linea x=0
  p <- p + geom_vline(xintercept = 0, color = "red")
  # Intervalos de Confianza - Azules
  test = t.test(x.c, mu = 0, conf.level = conf.level)
  p <- p + geom_vline(aes(xintercept = test$conf.int[1], color = "Confidence interval"), linetype = "dashed", col = "blue") +
    geom_vline(aes(xintercept = test$conf.int[2], color = "Confidence interval"), linetype = "dashed", col = "blue") +
    theme(legend.position = "none")
  # Curva de densidad
  p <- p + geom_line(data = data.frame(x = density(x.c)$x, y = density(x.c)$y), aes(x = x, y = y), color = "black", linewidth = 0.5)

  # Hipotesis, P_val y t_val
  caja <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.border = element_blank()
    ) +
    xlim(c(0.25,0.26)) + ylim(c(0.24, 0.36))

  caja <- caja +
    annotate('text', x = 0.25, y = 0.35,
             label = paste("H[0]==", "0"),
             parse = TRUE, size = 3, hjust = 0) +
    annotate('text', x = 0.25, y = 0.30,
             label = paste("t-value: ", round(test$statistic, digits = 3)),
             parse = TRUE, size = 3, hjust = 0) +
    annotate('text', x = 0.25, y = 0.25,
             label = paste("p-value: ", round(test$p.value, 3)),
             parse = TRUE, size = 3, hjust = 0)

  p <- p + inset_element(caja, left = 0, right = 0.35, top = 1, bottom = 0.65)

  # Añadir label del Cg y Cgk
  if (cgOut == TRUE) {
    caja <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
      theme_bw() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank()
      ) +
      xlim(c(0.25,0.26)) + ylim(c(0.24, 0.31))

    caja <- caja +
      annotate('text', x = 0.25, y = 0.30,
               label = paste("Cg: ", round(Cg,digits = 6)),
               parse = TRUE, size = 3, hjust = 0) +
      annotate('text', x = 0.25, y = 0.25,
               label = paste("Cgk:", round(Cgk,digits = 6)),
               parse = TRUE, size = 3, hjust = 0)

    p <- p + inset_element(caja, left = 0.7, right = 1, top = 1, bottom = 0.75)
    show(p)
  }
  else{
    caja <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
      theme_bw() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank()
      ) +
      xlim(c(0.25,0.26)) + ylim(c(0.24, 0.31))

    caja <- caja +
      annotate('text', x = 0.25, y = 0.30,
               label = "----- conf. int", size = 3, hjust = 0, colour = "blue")

    p <- p + inset_element(caja, left = 0.75, right = 1, top = 1, bottom = 0.75)
    show(p)
  }
  invisible(list(Cg, Cgk))
}

cg_ToleranceChart <- function (x, target, tolerance, ref.interval, facCg, facCgk,
                               n = 0.2, col, pch, xlim, ylim, main, conf.level = 0.95,
                               cgOut = TRUE)
{
  if (missing(x))
    stop("x must be given as a vector")
  if (missing(target)) {
    target = mean(x)
    targetmissing = FALSE
  }
  else targetmissing = TRUE
  if (missing(ref.interval))
    ref.interval = pnorm(3) - pnorm(-3)
  sd = sd(x)
  mean = mean(x)
  ref.ar = qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval,
                                                 mean, sd)
  if (missing(facCg))
    facCg = 0.2
  if (missing(facCgk))
    facCgk = 0.1
  if (missing(tolerance))
    warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
  if (missing(tolerance)) {
    width = ref.ar/facCg
    tolerance = numeric(2)
    tolerance[1] = mean(x) - width/2
    tolerance[2] = mean(x) + width/2
  }
  quant1 = qnorm((1 - ref.interval)/2, mean, sd)
  quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean,
                 sd)
  if (length(tolerance) != 2)
    stop("tolerance has wrong length")
  if (missing(col))
    col = 1
  if (missing(pch))
    pch = 19
  if (missing(xlim))
    xlim = c(0, length(x))
  if (missing(ylim))
    ylim = c(min(x, target - n/2 * (abs(diff(tolerance))),
                 quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))),
                                      quant1, quant2))
  if (missing(main))
    main = "Tolerance View"
  Cg <- (facCg * tolerance[2]-tolerance[1])/ref.interval
  Cgk <- (facCgk * abs(target-mean(x))/(ref.interval/2))

  # Gráfica
  p <- ggplot(data.frame(x = 1:length(x), y = x), aes(x = x, y = y)) +
    geom_point(color = col, shape = pch) +
    geom_line(color = col, linetype = "solid")+
    geom_hline(aes(yintercept = target), linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = tolerance[1]), linetype = "dashed", color = "blue") +
    geom_hline(aes(yintercept = tolerance[2]), linetype = "dashed", color = "blue") +
    geom_hline(aes(yintercept = (target + n/2 * (tolerance[2] - tolerance[1]))), color = "black") +
    geom_hline(aes(yintercept = (target - n/2 * (tolerance[2] - tolerance[1]))), color = "black") +
    scale_color_manual(values = c("Data" = col)) +
    labs(x = "", y = "x", color = "Variable", title = main)+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
    theme(legend.position = "none")

  if (cgOut == TRUE) {
    caja <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
      theme_bw() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank()
      ) +
      xlim(c(0.25,0.26)) + ylim(c(0.24, 0.31))

    caja <- caja +
      annotate('text', x = 0.25, y = 0.30,
               label = paste("Cg: ", round(Cg,digits = 6)),
               parse = TRUE, size = 3, hjust = 0) +
      annotate('text', x = 0.25, y = 0.25,
               label = paste("Cgk:", round(Cgk,digits = 6)),
               parse = TRUE, size = 3, hjust = 0)

    p <- p + inset_element(caja, left = 0.7, right = 1, top = 1, bottom = 0.75)
    show(p)
  }
  else{
    show(p)
  }
  invisible(list(Cg, Cgk))
}


cg <- function (x, target, tolerance, ref.interval, facCg, facCgk, n = 0.2,
                col, pch, xlim, ylim, conf.level = 0.95){
    old.par <- par(no.readonly = TRUE)
    if (missing(x))
      stop("x must be given as a vector")
    if (missing(target)) {
      target = mean(x)
      targetmissing = FALSE
    }
    else
      targetmissing = TRUE
    if (missing(ref.interval))
      ref.interval = pnorm(3) - pnorm(-3)
    sd = sd(x)
    mean = mean(x)
    ref.ar = qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval,
                                                   mean, sd)
    if (missing(facCg))
      facCg = 0.2
    if (missing(facCgk))
      facCgk = 0.1
    if (missing(tolerance))
      warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
    if (missing(tolerance)) {
      width = ref.ar / facCg
      tolerance = numeric(2)
      tolerance[1] = mean(x) - width / 2
      tolerance[2] = mean(x) + width / 2
    }
    quant1 = qnorm((1 - ref.interval) / 2, mean, sd)
    quant2 = qnorm(ref.interval + (1 - ref.interval) / 2, mean,
                   sd)
    if (length(tolerance) != 2)
      stop("tolerance has wrong length")
    if (missing(col))
      col = 1
    if (missing(pch))
      pch = 19
    if (missing(xlim))
      xlim = c(0, length(x))
    if (missing(ylim))
      ylim = c(min(x, target - n / 2 * (abs(diff(
        tolerance
      ))),
      quant1, quant2),
      max(x, target + n / 2 * (abs(diff(
        tolerance
      ))),
      quant1, quant2))
    Cg <- (facCg * tolerance[2] - tolerance[1]) / ref.interval
    Cgk <- (facCgk * abs(target - mean(x)) / (ref.interval / 2))

    # Plots

    # RunChart
    df1 <- data.frame(x = x, y = x)
    df1$y_target <- target
    df1$y_lower <- quant1
    df1$y_upper <- quant2
    df1$y_tolerance_lower <- tolerance[1]
    df1$y_tolerance_upper <- tolerance[2]
    df1$y_mean <- mean
    df1$Cg <- Cg
    df1$Cgk <- Cgk
    p1 <- ggplot(df1, aes(x = seq_along(x), y = x)) +
      geom_point(color = col, shape = pch) +
      geom_line(color = col, linetype = "solid") +
      scale_x_continuous(limits = c(xlim[1] - 0.05 * xlim[2], xlim[2]),
                         expand = c(0, 0)) +
      labs(title = "Run Chart", x = "Index", y = "x") +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      geom_hline(aes(yintercept = target)) +
      geom_smooth(
        method = "loess",
        color = "red",
        se = FALSE,
        span = 1.25,
        size = 0.25
      ) +
      geom_hline(aes(yintercept = mean),
                 linetype = "dashed",
                 color = "seagreen") +  #center line
      geom_hline(aes(yintercept = quant1),
                 linetype = "dashed",
                 color = "seagreen") + # Bottom line
      geom_hline(aes(yintercept = quant2),
                 linetype = "dashed",
                 color = "seagreen") +
      geom_hline(
        yintercept = c(target + n / 2 * (abs(diff(
          tolerance
        ))), target - n / 2 * (abs(diff(
          tolerance
        )))),
        color = "#012B78",
        linetype = "solid"
      ) +
      scale_y_continuous(
        limits = ylim,
        expand = c(0, 0),
        sec.axis =
          sec_axis(
            ~ .,
            breaks = c(
              target,
              mean,
              quant1,
              quant2,
              target + n / 2 * (abs(diff(tolerance))),
              target - n / 2 * (abs(diff(tolerance)))
            ),
            labels = c(
              "target",
              expression(bar(x)),
              substitute(x[a * b], list(a = round(((1 - ref.interval) / 2
              ) * 100, 3), b = "%")),
              substitute(x[a * b], list(a = round(((ref.interval + (1 - ref.interval) /
                                                      2)
              ) * 100,
              3), b = "%")),
              substitute(x[tar] + a, list(a = round(n / 2, 4))),
              substitute(x[tar] - a, list(a = round(n / 2, 4)))
            )
          )
      ) + theme(axis.text.y.right = element_text(size = 15))

    # HistChart
    x.c <- x - target
    temp <- hist(x.c, plot = FALSE)
    df3 <- data.frame(mid = temp$mids,
                      density = temp$density)
    width <- diff(df3$mid)[1] # Ancho de cada barra
    p3 <- ggplot(df3, aes(x = mid, y = density)) +
      geom_bar(
        stat = "identity",
        width = width,
        fill = "lightblue",
        color = "black",
        alpha = 0.5
      ) +
      labs(y = "Density", x = "x.c", title = paste("Histogram of",deparse(substitute(x)),"- target")) +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      guides(color = guide_legend(title.position = "top", title.hjust = 0.5)) +
      geom_vline(xintercept = 0, color = "red")

    test = t.test(x.c, mu = 0, conf.level = conf.level)
    p3 <-
      p3 + geom_vline(
        aes(xintercept = test$conf.int[1], color = "Confidence interval"),
        linetype = "dashed",
        col = "blue"
      ) +
      geom_vline(
        aes(xintercept = test$conf.int[2], color = "Confidence interval"),
        linetype = "dashed",
        col = "blue"
      ) +
      theme(legend.position = "none") +
      geom_line(
        data = data.frame(x = density(x.c)$x, y = density(x.c)$y),
        aes(x = x, y = y),
        color = "black",
        linewidth = 0.5
      )

    p3 <-
      p3 + annotation_custom(grob = grid::textGrob(
        label = c(expression(paste(H[0], " : Bias = 0"))),
        x = unit(0.05, "npc") + unit(0.05, "cm"),
        y = unit(1, "npc") - unit(0.05, "cm"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 6, fontface = "bold")
      )) +
      annotation_custom(grob = grid::textGrob(
        label = c(paste("t-value: ", round(test$statistic, 3))),
        x = unit(0.05, "npc") + unit(0.05, "cm"),
        y = unit(1, "npc") - unit(0.4, "cm"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 6)
      )) +
      annotation_custom(grob = grid::textGrob(
        label = c(paste("p-value: ", round(test$p.value, 3))),
        x = unit(0.05, "npc") + unit(0.05, "cm"),
        y = unit(1, "npc") - unit(0.65, "cm"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 6)
      ))

    # Tolerance View
    p4 <-
      ggplot(data.frame(x = 1:length(x), y = x), aes(x = x, y = y)) +
      geom_point(color = col, shape = pch) +
      geom_line(color = col, linetype = "solid") +
      geom_hline(aes(yintercept = target),
                 linetype = "dashed",
                 color = "red") +
      geom_hline(aes(yintercept = tolerance[1]),
                 linetype = "dashed",
                 color = "blue") +
      geom_hline(aes(yintercept = tolerance[2]),
                 linetype = "dashed",
                 color = "blue") +
      geom_hline(aes(yintercept = (target + n / 2 * (
        tolerance[2] - tolerance[1]
      ))), color = "black") +
      geom_hline(aes(yintercept = (target - n / 2 * (
        tolerance[2] - tolerance[1]
      ))), color = "black") +
      scale_color_manual(values = c("Data" = col)) +
      labs(
        x = "",
        y = "x",
        color = "Variable",
        title = "Tolerance View"
      ) +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      theme(legend.position = "none")

    # Data Box
    p2 <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
      theme_bw() +
      annotate(
        'text',
        x = 0.25,
        y = 0.40,
        label = paste("bar(x)==", round(mean, 2)),
        parse = TRUE,
        size = 3,
        hjust = 0
      ) +
      annotate(
        'text',
        x = 0.25,
        y = 0.35,
        label = paste("s==", round(sd, 2)),
        parse = TRUE,
        size = 3,
        hjust = 0
      ) +
      annotate(
        'text',
        x = 0.25,
        y = 0.3,
        label = paste("target==", round(target, 5)),
        parse = TRUE,
        size = 3,
        hjust = 0
      ) +
      annotate(
        'text',
        x = 0.25,
        y = 0.25,
        label = paste("C[g]==", round(Cg, 2)),
        parse = TRUE,
        size = 3,
        hjust = 0
      ) +
      annotate(
        'text',
        x = 0.25,
        y = 0.20,
        label = paste("C[gk]==", round(Cgk, 2)),
        parse = TRUE,
        size = 3,
        hjust = 0
      ) +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      xlim(c(0.15, 0.5)) + ylim(c(0.1, 0.5))

    design <- "
  112
  113
  114
  "
    p <- p1 + p2 + p3 + p4 + plot_layout(design = design)

    suppressMessages(show(p))

    invisible(list(Cg, Cgk))
  }


x <- c(9.991, 10.013, 10.001, 10.007, 10.010, 10.013, 10.008, 10.017, 10.005, 10.005, 10.002,
      10.017, 10.005, 10.002, 9.996, 10.011, 10.009 , 10.006, 10.008, 10.003, 10.002, 10.006,
      10.010, 9.992, 10.013)

cg(x, target = 10.003, tolerance = c(9.903, 10.103))

