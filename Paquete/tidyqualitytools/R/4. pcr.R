# Correr los archivo:
# "Class Distr, pplot y qqPlot"
# Funciones Auxilicares

pcr <- function (x, distribution = "normal", lsl, usl, target, boxcox = FALSE,
                  lambda = c(-5, 5), main, xlim, ylim, grouping = NULL, std.dev = NULL,
                  conf.level = 0.9973002, start, lineWidth = 1, lineCol = "red",
                  lineType = "solid", specCol = "red3", specWidth = 1, cex.text = 2,
                  cex.val = 1.5, cex.col = "darkgray", plot = TRUE, bounds.lty = 3,
                  bounds.col = "red", ...) {
  data.name = deparse(substitute(x))[1]

  parList = list(...)
  if (is.null(parList[["col"]]))
    parList$col = "lightblue"
  if (is.null(parList[["border"]]))
    parList$border = 1
  if (is.null(parList[["lwd"]]))
    parList$lwd = 1
  if (is.null(parList[["cex.axis"]]))
    parList$cex.axis = 1.5
  if (missing(lsl))
    lsl = NULL
  if (missing(usl))
    usl = NULL
  if (missing(target))
    target = NULL
  if (missing(lambda))
    lambda = c(-5, 5)
  if (!is.numeric(lambda))
    stop("lambda needs to be numeric")
  if (any(x < 0) && any(distribution == c("beta", "chi-squared",
                                          "exponential", "f", "geometric", "lognormal", "log-normal",
                                          "negative binomial", "poisson", "weibull", "gamma")))
    stop("choosen distribution needs all values in x to be > 0!")
  if (any(x > 1) && distribution == "beta")
    stop("choosen distribution needs all values in x to be between 0 and 1!")

  paramsList = vector(mode = "list", length = 0)
  estimates = vector(mode = "list", length = 0)
  varName = deparse(substitute(x))
  dFun = NULL
  pFun = NULL
  qFun = NULL
  cp = NULL
  cpu = NULL
  cpl = NULL
  cpk = NULL
  ppt = NULL
  ppl = NULL
  ppu = NULL
  xVec = numeric(0)
  yVec = numeric(0)

  if (is.vector(x))
    x = as.data.frame(x)
  any3distr = FALSE
  not3distr = FALSE

  if (distribution == "weibull3" || distribution == "lognormal3" ||
      distribution == "gamma3")
    any3distr = TRUE
  if (distribution != "weibull3" && distribution != "lognormal3" &&
      distribution != "gamma3")
    not3distr = TRUE
  if (boxcox) {
    distribution = "normal"
    if (length(lambda) >= 2) {
      temp = boxcox(x[, 1] ~ 1, lambda = seq(min(lambda),
                                             max(lambda), 1/10), plotit = FALSE)
      i = order(temp$y, decreasing = TRUE)[1]
      lambda = temp$x[i]
    }
    x = as.data.frame(x[, 1]^lambda)
  }
  numObs = nrow(x)
  if (!is.null(grouping))
    if (is.vector(grouping))
      grouping = as.data.frame(grouping)
  center = colMeans(x)
  if (!is.null(x) & !is.null(grouping)) {
    if (nrow(x) != nrow(grouping))
      stop(paste("length of ", deparse(substitute(grouping)),
                 " differs from length of ", varName))
  }
  if (missing(main))
    if (boxcox)
      main = paste("Process Capability using box cox transformation for",
                   varName)
  else main = paste("Process Capability using", as.character(distribution),
                    "distribution for", varName)
  if (is.null(std.dev)) {
    if (is.null(grouping))
      std.dev = .sdSg(x)
    else std.dev = .sdSg(x, grouping)
  }
  if (conf.level < 0 | conf.level > 1)
    stop("conf.level must be a value between 0 and 1")
  confHigh = conf.level + (1 - conf.level)/2
  confLow = 1 - conf.level - (1 - conf.level)/2
  distWhichNeedParameters = c("weibull", "logistic", "gamma",
                              "exponential", "f", "geometric", "chi-squared", "negative binomial",
                              "poisson")
  if (is.character(distribution)) {
    dis = distribution
    if (identical(distribution, "weibull3"))
      dis = "weibull3"
    if (identical(distribution, "gamma3"))
      dis = "gamma3"
    if (identical(distribution, "lognormal3"))
      dis = "lognormal3"
    qFun = .charToDistFunc(dis, type = "q")
    pFun = .charToDistFunc(dis, type = "p")
    dFun = .charToDistFunc(dis, type = "d")
    if (is.null(qFun) & is.null(pFun) & is.null(dFun))
      stop(paste(deparse(substitute(y)), "distribution could not be found!"))
  }
  if (TRUE) {
    fitList = vector(mode = "list", length = 0)
    fitList$x = x[, 1]
    fitList$densfun = dis
    if (!missing(start))
      fitList$start = start
    if (not3distr) {
      fittedDistr = do.call(FitDistr, fitList)
      estimates = as.list(fittedDistr$estimate)
      paramsList = estimates
    }
    if (distribution == "weibull3") {
      paramsList = .weibull3(x[, 1])
      estimates = paramsList
    }
    if (distribution == "lognormal3") {
      paramsList = .lognormal3(x[, 1])
      estimates = paramsList
    }
    if (distribution == "gamma3") {
      paramsList = .gamma3(x[, 1])
      estimates = paramsList
    }

  }
  paramsList = c(paramsList, .lfkp(parList, formals(qFun)))
  if (distribution == "normal") {
    paramsList$mean = center
    paramsList$sd = std.dev
    estimates = paramsList
  }
  if (boxcox) {
    if (!is.null(lsl))
      lsl = lsl^lambda
    if (!is.null(usl))
      usl = usl^lambda
    if (!is.null(target))
      target = target^lambda
  }
  if (is.null(lsl) && is.null(usl)) {
    paramsList$p = confLow
    lsl = do.call(qFun, paramsList)
    paramsList$p = confHigh
    usl = do.call(qFun, paramsList)
  }

  if (identical(lsl, usl))
    stop("lsl == usl")
  if (!is.null(lsl) && !is.null(target) && target < lsl)
    stop("target is less than lower specification limit")
  if (!is.null(usl) && !is.null(target) && target > usl)
    stop("target is greater than upper specification limit")
  if (!is.null(lsl) && !is.null(usl))
    if (lsl > usl) {
      temp = lsl
      lsl = usl
      usl = temp
    }
  paramsList$p = c(confLow, 0.5, confHigh)
  paramsListTemp = .lfkp(paramsList, formals(qFun))
  qs = do.call(qFun, paramsListTemp)
  paramsListTemp = .lfkp(paramsList, formals(pFun))
  if (!is.null(lsl) && !is.null(usl))
    cp = (usl - lsl)/(qs[3] - qs[1])
  if (!is.null(usl)) {
    cpu = (usl - qs[2])/(qs[3] - qs[2])
    paramsListTemp$q = usl
    ppu = 1 - do.call(pFun, paramsListTemp)
  }
  if (!is.null(lsl)) {
    cpl = (qs[2] - lsl)/(qs[2] - qs[1])
    paramsListTemp$q = lsl
    ppl = do.call(pFun, paramsListTemp)
  }
  cpk = min(cpu, cpl)
  ppt = sum(ppl, ppu)

  if(plot==TRUE){
    # ----------------------------- IF PLOT == TRUE -----------------------------------------------------------
    if (missing(xlim)) {
      xlim <- range(x[, 1], usl, lsl)
      xlim <- xlim + diff(xlim) * c(-0.2, 0.2)
    }
    xVec <- seq(min(xlim), max(xlim), length = 200)
    dParamsList = .lfkp(paramsList, formals(dFun))
    dParamsList$x = xVec

    yVec = do.call(dFun, dParamsList)
    histObj <- hist(x[, 1], plot = FALSE)
    if (missing(ylim)) {
      ylim <- range(histObj$density, yVec)
      ylim <- ylim + diff(ylim) * c(0, 0.05)
    }

    # 1. Histograma --------------------------------------------------------------------------
    # Calculos previos
    x.c <- x[, 1]
    temp <- hist(x.c, plot = FALSE)
    # Obtenemos la información para el histograma
    df <- data.frame(
      mid = temp$mids,
      density = temp$density
    )
    width <- diff(df$mid)[1] # Ancho de cada barra
    # Histograma
    p1 <- ggplot(df, aes(x = mid, y = density)) +
      geom_bar(stat = "identity", width = width, fill = "lightblue", color = "black", alpha = 0.5) +
      labs(y = "", x = "", title = "") +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
      guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
      geom_line(data = data.frame(x = xVec, y = yVec), aes(x = x, y = y), color = "red", linewidth = 0.5) + # densidad
      theme(legend.position = "none")

    #  etiquetas de los límites
    if (!is.null(lsl) & !is.null(usl)){
      p1 <- p1 +
        geom_vline(aes(xintercept = usl, color = "Confidence interval"), linetype = "dashed", col = "red") + # USL
        geom_vline(aes(xintercept = lsl, color = "Confidence interval"), linetype = "dashed", col = "red") + # LSL
        scale_x_continuous(limits = xlim, expand = c(0, 0),
                           sec.axis = sec_axis(~ ., breaks = c(lsl, usl),
                                               labels = c(paste("LSL =",format(lsl, digits = 3)), paste("USL =",format(usl, digits = 3)))
                           )) +
        theme(axis.text.y.right = element_text(size = 15))
    }else{
      if(!is.null(lsl)){
        p1 <- p1 +
          geom_vline(aes(xintercept = lsl, color = "Confidence interval"), linetype = "dashed", col = "red") + # LSL
          scale_x_continuous(limits = xlim, expand = c(0, 0),
                             sec.axis = sec_axis(~ ., breaks = lsl, labels = paste("LSL =",format(lsl, digits = 3)) )) +
          theme(axis.text.y.right = element_text(size = 15))
      }
      if(!is.null(usl)){
        p1 <- p1 +
          geom_vline(aes(xintercept = usl, color = "Confidence interval"), linetype = "dashed", col = "red") + # USL
          scale_x_continuous(limits = xlim, expand = c(0, 0),
                             sec.axis = sec_axis(~ ., breaks = usl,labels = paste("USL =",format(usl, digits = 3)))) +
          theme(axis.text.y.right = element_text(size = 15))
      }
    }

    # 2. Cajita de info Cp's --------------------------------------------------------------------------
    p2 <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
      theme_bw() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      xlim(c(0.24, 0.26)) + ylim(c(0.21, 0.43))
    {
      if(is.null(cpu))
        p2 <- p2 + annotate('text', x = 0.25, y = 0.40,
                            label = paste("C[pkL] == ", "NA"),
                            parse = TRUE, size = 3.5, hjust = 0.5)
      else p2 <- p2 + annotate('text', x = 0.25, y = 0.40,
                               label = paste("C[pkU]==", round(cpu, 2)),
                               parse = TRUE, size = 3.5, hjust = 0.5)
      if(is.null(cpl))
        p2 <- p2 + annotate('text', x = 0.25, y = 0.35,
                            label = paste("C[pkL] == ", "NA"),
                            parse = TRUE, size = 3.5, hjust = 0.5)
      else p2 <- p2 + annotate('text', x = 0.25, y = 0.35,
                               label = paste("C[pkL]==", round(cpl, 2)),
                               parse = TRUE, size = 3.5, hjust = 0.5)
      if(is.null(cpk))
        p2 <- p2 + annotate('text', x = 0.25, y = 0.30,
                            label = paste("C[pkL] == ", "NA"),
                            parse = TRUE, size = 3.5, hjust = 0.5)
      else p2 <- p2 + annotate('text', x = 0.25, y = 0.30,
                               label = paste("C[pk]==", round(cpk, 2)),
                               parse = TRUE, size = 3.5, hjust = 0.5)
      if(is.null(cp))
        p2 <- p2 + annotate('text', x = 0.25, y = 0.25,
                            label = paste("C[pkL] == ", "NA"),
                            parse = TRUE, size = 3.5, hjust = 0.5)
      else p2 <- p2 + annotate('text',x = 0.25,y = 0.25,
                               label = paste("C[p]==", round(cp, 2)),
                               parse = TRUE,size = 3.5,hjust = 0.5)
      }
    # 3. Cajita de info n, means, sd --------------------------------------------------------------------------
    index = 1:(length(estimates) + 3)
    names(x) = data.name
    if(not3distr){
      names(x) = data.name
      adTestStats = .myADTest(x, distribution)
      A = numeric()
      p = numeric()
      if (adTestStats$class == "adtest"){
        A = adTestStats$statistic$A
        p = adTestStats$p.value$p
      }

      # Caja de Info
      p3 <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
        theme_bw() +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        xlim(c(0.24, 0.26)) + ylim(c(0.19, 0.43))
      {
        # n y A
        p3 <- p3 + annotate('text', x = 0.25, y = 0.40,
                            label = paste("n==", numObs),
                            parse = TRUE, size = 3, hjust = 0.5) +
          annotate('text', x = 0.25, y = 0.35,
                   label = paste("A==", format(as.numeric(A), digits = 3)),
                   parse = TRUE, size = 3, hjust = 0.5)

        # p
        if (!is.null(adTestStats$smaller) && adTestStats$smaller){
          p3 <- p3 + annotate(
            'text',
            x = 0.25,
            y = 0.30,
            label = paste("p<", format(as.numeric(p), digits =3)),
            parse = TRUE,
            size = 3,
            hjust = 0.5
          )
        }
        if (!is.null(adTestStats$smaller) && !adTestStats$smaller){
          p3 <- p3 + annotate(
            'text',
            x = 0.25,
            y = 0.30,
            label = paste("p>=", format(as.numeric(p),digits = 3)),
            parse = TRUE,
            size = 3,
            hjust = 0.5
          )
        }
        if (is.null(adTestStats$smaller)){
          p3 <- p3 + annotate(
            'text',
            x = 0.25,
            y = 0.30,
            label = paste("p==", format(as.numeric(p), digits = 3)),
            parse = TRUE,
            size = 3,
            hjust = 0.5
          )
        }

        # mean y sd
        p3 <- p3 + annotate('text', x = 0.25, y = 0.25,
                            label = paste(names(estimates)[1], "==", format(estimates[[names(estimates)[1]]], digits = 3)),
                            parse = TRUE, size = 3, hjust = 0.5) +
          annotate('text', x = 0.25, y = 0.20,
                   label = paste(names(estimates)[2], "==", format(estimates[[names(estimates)[2]]], digits = 3)),
                   parse = TRUE, size = 3, hjust = 0.5)
        }
    }
    if(any3distr){
      p3 <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
        theme_bw() +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        xlim(c(0.24, 0.26)) + ylim(c(0.21, 0.43))
      {
        # n y A
        p3 <- p3 + annotate('text', x = 0.25, y = 0.40,
                            label = paste("n==", numObs),
                            parse = TRUE, size = 3, hjust = 0.5) +
          annotate('text', x = 0.25, y = 0.35,
                   label = paste("A = ", "*"),
                   parse = FALSE, size = 3, hjust = 0.5) +
          annotate('text', x = 0.25, y = 0.30,
                   label = paste("p = ", "*"),
                   parse = FALSE, size = 3, hjust = 0.5) +
          annotate('text', x = 0.25, y = 0.25,
                   label = paste("mean==", format(estimates[[1]], digits = 3)),
                   parse = TRUE, size = 3, hjust = 0.5) +
          annotate('text', x = 0.25, y = 0.20,
                   label = paste("sd==", format(estimates[[2]], digits = 3)),
                   parse = TRUE, size = 3, hjust = 0.5)
        }
    }

    # 4. qqPlot --------------------------------------------------------------------------
    p4 <- qqPlot(x[, 1], y = distribution, xlab = "", ylab = "", main = "",
                 axes = FALSE, bounds.lty = bounds.lty, bounds.col = bounds.col, grapic = FALSE, axis.y.right = TRUE, bw.theme = TRUE)

    # Unimos las 4 primeras gráficas
    main_plot <- p1 + (p2 / p3 / p4$plot) + plot_layout(widths = c(5, 1))

    # 5. Cajita de info 4 (Expected Fraction Nonconforming) --------------------------------------------------------------------------
    p5 <- ggplot(data.frame(x = c(-1, 1),y = c(0.5, 5)), aes(x = x, y = y)) +
      theme_bw() +
      ggtitle("Expected Fraction Nonconforming") +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5, vjust = -0.5,margin = margin(b = -12),size = 10)
      )

    p5_left <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
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
    {
      p5_left <- p5_left +
        annotate('text', x = 0.25, y = 0.35, label = "-----", size = 3, hjust = 0, colour = "white")
      # Pt
      p5_left <- p5_left +
        annotate("text", x = 0.25, y = 0.33, label = paste("p[t]==", format(ppt, digits = 6)),
                 parse = TRUE, size = 3.5, hjust = 0)
      # PL
      if(is.null(ppl)){
        p5_left <- p5_left +
          annotate("text", x = 0.25, y = 0.3, label = paste("p[L]==", "0"),
                   parse = TRUE, size = 3.5, hjust = 0)
      }else{
        p5_left <- p5_left +
          annotate("text", x = 0.25, y = 0.3, label = paste("p[L]==", format(ppl, digits = 6)),
                   parse = TRUE, size = 3.5, hjust = 0)
      }
      # PU
      if(is.null(ppu)){
        p5_left <- p5_left +
          annotate("text", x = 0.25, y = 0.27, label = paste("p[U]==", "0"),
                   parse = TRUE, size = 3.5, hjust = 0)
      }else{
        p5_left <- p5_left +
          annotate("text", x = 0.25, y = 0.27, label = paste("p[U]==", format(ppu, digits = 6)),
                   parse = TRUE, size = 3.5, hjust = 0)
      }
      }

    p5_right <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
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
    {
      p5_right <- p5_right +
        annotate('text', x = 0.25, y = 0.35, label = "-----", size = 3, hjust = 0, colour = "white")

      # ppm
      p5_right <- p5_right +
        annotate("text", x = 0.25, y = 0.33, label = paste("ppm==", format(ppt * 1e+06, digits = 6)),
                 parse = TRUE, size = 3.5, hjust = 0)
      if(is.null(ppl)){
        p5_right <- p5_right +
          annotate("text", x = 0.25, y = 0.3, label = paste("ppm==", "0"),
                   parse = TRUE, size = 3.5, hjust = 0)
      }else{
        p5_right <- p5_right +
          annotate("text", x = 0.25, y = 0.3, label = paste("ppm==", format(ppl * 1e+06, digits = 6)),
                   parse = TRUE, size = 3.5, hjust = 0)
      }
      if(is.null(ppu)){
        p5_right <- p5_right +
          annotate("text", x = 0.25, y = 0.27, label = paste("ppm==", "0"),
                   parse = TRUE, size = 3.5, hjust = 0)
      }else{
        p5_right <- p5_right +
          annotate("text", x = 0.25, y = 0.27, label = paste("ppm==", format(ppu * 1e+06, digits = 6)),
                   parse = TRUE, size = 3.5, hjust = 0)
      }
      }

    p5 <- p5 + inset_element(p5_left, left = 0, right = 0.5, top = 0,  bottom = 0.80)+
      inset_element(p5_right, left = 0.5, right = 1, top = 0,  bottom = 0.80)

    # Caja info 6. Observed --------------------------------------------------------------------------
    obsL = 0
    obsU = 0
    p6 <- ggplot(data.frame(x = 0,y = 0), aes(x = x, y = y)) +
      theme_bw() +
      ggtitle("Observed") +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = -0.5,margin = margin(b = -12), size = 10)
      ) +
      xlim(c(0.24, 0.26)) + ylim(c(0.27, 0.36))
    if (!is.null(lsl)){
      obsL = (sum(x < lsl)/length(x)) * 1e+06
      p6 <- p6 + annotate("text", x = 0.25, y = 0.31, label = paste("ppm==", format(obsL, digits = 6)),
                          parse = TRUE, size = 3.5, hjust = 0.5)
    } else{
      p6 <- p6 + annotate("text", x = 0.25, y = 0.31, label = paste("ppm==", 0),
                          parse = TRUE, size = 3.5, hjust = 0.5)
    }
    if (!is.null(usl)){
      obsU = (sum(x > usl)/length(x)) * 1e+06
      p6 <- p6 + annotate("text", x = 0.25, y = 0.28, label = paste("ppm==", format(obsU, digits = 6)),
                          parse = TRUE, size = 3.5, hjust = 0.5)
    } else{
      p6 <- p6 + annotate("text", x = 0.25, y = 0.28, label = paste("ppm==", 0),
                          parse = TRUE, size = 3.5, hjust = 0.5)
    }
    p6 <- p6 + annotate("text", x = 0.25, y = 0.34, label = paste("ppm==", format(obsL + obsU, digits = 6)),
                        parse = TRUE, size = 3.5, hjust = 0.5)

    # UNION --------------------------------------------------------------------------
    box_bottom <- p5 + p6 + plot_layout(widths = c(2, 1))
    main_plot <- p1 / box_bottom + plot_layout(heights = c(1, 0.5))
    box_right <- p2 / p3 / p4$plot
    box_right <- box_right/plot_spacer()
    main_plot <- (main_plot | box_right) + plot_layout(ncol = 2, widths = c(5, 1))
    main_plot <- main_plot + plot_annotation(
      title = main,
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
    if(not3distr){
      print.adtest(adTestStats)
      show(main_plot)
      invisible(list(lambda = lambda, cp = cp, cpk = cpk,
                     cpl = cpl, cpu = cpu, ppt = ppt, ppl = ppl, ppu = ppu,
                     A = A, usl = usl, lsl = lsl, target = target, plot = main_plot))
    }else{
      show(main_plot)
      invisible(list(lambda = lambda, cp = cp, cpk = cpk,
                     cpl = cpl, cpu = cpu, ppt = ppt, ppl = ppl, ppu = ppu,
                     usl = usl, lsl = lsl, target = target, plot = main_plot))
    }
  }
  ## ---------------- end if plot == true --------------------------------------------------------------------------
  invisible(list(lambda = lambda, cp = cp, cpk = cpk, cpl = cpl,
                 cpu = cpu, ppt = ppt, ppl = ppl, ppu = ppu, usl = usl,
                 lsl = lsl, target = target, plot = main_plot))
}

.pcr = function(x, distribution = "normal", lsl, usl, target, boxcox = FALSE, lambda = c(-5,5), main, xlim, ylim, grouping = NULL, std.dev = NULL, conf.level = 0.9973002, start, lineWidth = 1,
                lineCol = "red", lineType = "solid", specCol = "red3", specWidth = 1, cex.text = 2, cex.val = 1.5, cex.col = "darkgray", plot = TRUE, ...) {
  data.name = deparse(substitute(x))[1]                                 ####
  parList = list(...)
  if (is.null(parList[["col"]]))
    parList$col = "lightblue"
  if (is.null(parList[["border"]]))
    parList$border = 1
  if (is.null(parList[["lwd"]]))
    parList$lwd = 1
  if (is.null(parList[["cex.axis"]]))
    parList$cex.axis = 1.5
  if (missing(lsl))
    lsl = NULL
  if (missing(usl))
    usl = NULL
  if (missing(target))
    target = NULL
  if (missing(lambda))
    lambda = c(-5, 5)
  if (!is.numeric(lambda))
    stop("lambda needs to be numeric")
  paramsList = vector(mode = "list", length = 0)
  estimates = vector(mode = "list", length = 0)
  varName = deparse(substitute(x))
  dFun = NULL
  pFun = NULL
  qFun = NULL
  cp = NULL
  cpu = NULL
  cpl = NULL
  cpk = NULL
  ppt = NULL
  ppl = NULL
  ppu = NULL
  xVec = numeric(0)
  yVec = numeric(0)
  if (is.vector(x))
    x = as.data.frame(x)                                                                          ####
  any3distr=FALSE;not3distr=FALSE                                            ####
  if(distribution=="weibull3" || distribution=="lognormal3" || distribution=="gamma3")####
    any3distr=TRUE                                                             ####
  if (distribution!="weibull3" && distribution!="lognormal3" && distribution!="gamma3")####
    not3distr=TRUE                                                                  ####
  if (boxcox) {
    distribution = "normal"
    if (length(lambda) >= 2) {
      temp = boxcox(x[, 1] ~ 1, lambda = seq(min(lambda), max(lambda), 1/10), plotit = FALSE)
      i = order(temp$y, decreasing = TRUE)[1]
      lambda = temp$x[i]
    }
    x = as.data.frame(x[, 1]^lambda)
  }
  numObs = nrow(x)
  if (!is.null(grouping))
    if (is.vector(grouping))
      grouping = as.data.frame(grouping)
  center = colMeans(x)
  if (!is.null(x) & !is.null(grouping)) {
    if (nrow(x) != nrow(grouping))
      stop(paste("length of ", deparse(substitute(grouping)), " differs from length of ", varName))
  }
  if (missing(main))
    if (boxcox)
      main = paste("Process Capability using box cox transformation for", varName)
  else main = paste("Process Capability using", as.character(distribution), "distribution for",
                    varName)
  if (is.null(std.dev)) {
    if (is.null(grouping))
      std.dev = .sdSg(x)
    else std.dev = .sdSg(x, grouping)
  }
  if (conf.level < 0 | conf.level > 1)
    stop("conf.level must be a value between 0 and 1")
  confHigh = conf.level + (1 - conf.level)/2
  confLow = 1 - conf.level - (1 - conf.level)/2
  distWhichNeedParameters = c("weibull", "logistic", "gamma", "exponential", "f", "geometric",
                              "chi-squared", "negative binomial", "poisson")
  if (is.character(distribution)) {
    dis=distribution                                                           ####
    if (identical(distribution,"weibull3"))                                     ####
      dis="weibull3"                                                             ####
    if (identical(distribution,"gamma3"))                                       ####
      dis="gamma3"                                                               ####
    if (identical(distribution,"lognormal3"))                                   ####
      dis="lognormal3"                                                           ####
    qFun = .charToDistFunc(dis, type = "q")                                 ####
    pFun = .charToDistFunc(dis, type = "p")                                 ####
    dFun = .charToDistFunc(dis, type = "d")                                 ####
    if (is.null(qFun) & is.null(pFun) & is.null(dFun))
      stop(paste(deparse(substitute(y)), "distribution could not be found!"))
  }
  if (TRUE) {                                                                 #### distribution!="weibull3" && distribution!="lognormal3" && distribution!="gamma3"
    fitList = vector(mode = "list", length = 0)
    fitList$x = x[, 1]
    fitList$densfun = dis                                                   ####
    if (!missing(start))
      fitList$start = start
    if (not3distr)                                                          ####
    {                                                                       ####
      fittedDistr = do.call(MASS::fitdistr, fitList)
      estimates = as.list(fittedDistr$estimate)
      paramsList = estimates
    }                                                                       ####
    if (distribution=="weibull3")                                           ####
    {                                                                       ####
      paramsList= .weibull3(x[,1])                                           ####
      estimates = paramsList                                                 ####
    }                                                                       ####
    if (distribution=="lognormal3")                                         ####
    {                                                                       ####
      paramsList= .lognormal3(x[,1])                                         ####
      estimates = paramsList                                                 ####
    }                                                                       ####
    if (distribution=="gamma3")                                             ####
    {                                                                       ####
      paramsList= .gamma3(x[,1])                                             ####
      estimates = paramsList                                                 ####
    }                                                                       ####
  }
  paramsList = c(paramsList, .lfkp(parList, formals(qFun)))
  if (distribution == "normal") {
    paramsList$mean = center
    paramsList$sd = std.dev
    estimates = paramsList
  }
  if (boxcox) {
    if (!is.null(lsl))
      lsl = lsl^lambda
    if (!is.null(usl))
      usl = usl^lambda
    if (!is.null(target))
      target = target^lambda
  }
  if (is.null(lsl) && is.null(usl)) {
    paramsList$p = confLow
    lsl = do.call(qFun, paramsList)
    paramsList$p = confHigh
    usl = do.call(qFun, paramsList)
  }
  if (identical(lsl, usl))
    stop("lsl == usl")
  if (!is.null(lsl) && !is.null(target) && target < lsl)
    stop("target is less than lower specification limit")
  if (!is.null(usl) && !is.null(target) && target > usl)
    stop("target is greater than upper specification limit")
  if (!is.null(lsl) && !is.null(usl))
    if (lsl > usl) {
      temp = lsl
      lsl = usl
      usl = temp
    }
  paramsList$p = c(confLow, 0.5, confHigh)
  paramsListTemp = .lfkp(paramsList, formals(qFun))                           ####
  qs = do.call(qFun, paramsListTemp)                                         ####
  paramsListTemp = .lfkp(paramsList, formals(pFun))                           ####
  if (!is.null(lsl) && !is.null(usl))
    cp = (usl - lsl)/(qs[3] - qs[1])
  if (!is.null(usl)) {
    cpu = (usl - qs[2])/(qs[3] - qs[2])
    paramsListTemp$q = usl                                                  ####
    ppu = 1 - do.call(pFun, paramsListTemp)                                 ####
  }
  if (!is.null(lsl)) {
    cpl = (qs[2] - lsl)/(qs[2] - qs[1])
    paramsListTemp$q = lsl                                                  ####
    ppl = do.call(pFun, paramsListTemp)                                     ####
  }
  cpk = min(cpu, cpl)
  ppt = sum(ppl, ppu)

  if(plot==TRUE){
    if (missing(xlim)) {
      xlim <- range(x[, 1], usl, lsl)
      xlim <- xlim + diff(xlim) * c(-0.2, 0.2)
    }
    xVec <- seq(min(xlim), max(xlim), length = 200)
    dParamsList = .lfkp(paramsList, formals(dFun))
    dParamsList$x = xVec

    yVec = do.call(dFun, dParamsList)
    histObj <- hist(x[, 1], plot = FALSE)
    if (missing(ylim)) {
      ylim <- range(histObj$density, yVec)
      ylim <- ylim + diff(ylim) * c(0, 0.05)
    }

    # 1. Histograma --------------------------------------------------------------------------
    x.c <- x[, 1]
    temp <- hist(x.c, plot = FALSE)
    df <- data.frame(
      mid = temp$mids,
      density = temp$density
    )
    width <- diff(df$mid)[1] # Ancho de cada barra
    # Histograma
    p1 <- ggplot(df, aes(x = mid, y = density)) +
      geom_bar(stat = "identity", width = width, fill = "lightblue", color = "black", alpha = 0.5) +
      labs(y = "", x = "", title = main) +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
      guides(color = guide_legend(title.position = "top", title.hjust = 0.5))+
      geom_line(data = data.frame(x = xVec, y = yVec), aes(x = x, y = y), color = "red", linewidth = 0.5) + # densidad
      theme(legend.position = "none")

    #  etiquetas de los límites
    if (!is.null(lsl) & !is.null(usl)){
      p1 <- p1 +
        geom_vline(aes(xintercept = usl, color = "Confidence interval"), linetype = "dashed", col = "red") + # USL
        geom_vline(aes(xintercept = lsl, color = "Confidence interval"), linetype = "dashed", col = "red") + # LSL
        scale_x_continuous(limits = xlim, expand = c(0, 0),
                           sec.axis = sec_axis(~ ., breaks = c(lsl, usl),
                                               labels = c(paste("LSL =",format(lsl, digits = 3)), paste("USL =",format(usl, digits = 3)))
                           )) +
        theme(axis.text.y.right = element_text(size = 15))
    }else{
      if(!is.null(lsl)){
        p1 <- p1 +
          geom_vline(aes(xintercept = lsl, color = "Confidence interval"), linetype = "dashed", col = "red") + # LSL
          scale_x_continuous(limits = xlim, expand = c(0, 0),
                             sec.axis = sec_axis(~ ., breaks = lsl, labels = paste("LSL =",format(lsl, digits = 3)) )) +
          theme(axis.text.y.right = element_text(size = 15))
      }
      if(!is.null(usl)){
        p1 <- p1 +
          geom_vline(aes(xintercept = usl, color = "Confidence interval"), linetype = "dashed", col = "red") + # USL
          scale_x_continuous(limits = xlim, expand = c(0, 0),
                             sec.axis = sec_axis(~ ., breaks = usl,labels = paste("USL =",format(usl, digits = 3)))) +
          theme(axis.text.y.right = element_text(size = 15))
      }
    }
    return(list(lambda = lambda, cp = cp, cpk = cpk, cpl = cpl, cpu = cpu,        ####
                ppt = ppt, ppl = ppl, ppu = ppu, usl = usl,                  ####
                lsl = lsl, target = target, plot = p1))
  }

  return(list(lambda = lambda, cp = cp, cpk = cpk, cpl = cpl, cpu = cpu,        ####
              ppt = ppt, ppl = ppl, ppu = ppu, usl = usl,                  ####
              lsl = lsl, target = target))                                 ####
}

