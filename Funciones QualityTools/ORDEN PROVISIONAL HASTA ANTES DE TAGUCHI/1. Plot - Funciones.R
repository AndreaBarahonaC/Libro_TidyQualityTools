###########################################################################
######################### FUNCIONES PARA GRÁFICAS #########################
###########################################################################

# ParetoChart ----------
paretoChart <- function (x, weight, showTable = TRUE, showPlot = TRUE,
                         main, col, border, xlab, ylab = "Frequency", percentVec, ...){
  varName = deparse(substitute(x))[1]
  corp.col = "#C4B9FF"
  corp.border = "#9E0138"
  if (!is.vector(x) & !is.data.frame(x) & !is.table(x))
    stop("x should be a vector, dataframe or a table")
  if (is.table(x)) {
    xtable = x
  }
  if (is.vector(x)) {
    if (!is.null(names(x)))
      xtable = as.table(x)
    else xtable = table(x)
  }
  if (!missing(weight)) {
    if (!is.numeric(weight))
      stop("weight must be numeric!")
    if (is.null(names(weight)))
      stop("weight is missing names for matching!")
    else {
      if (FALSE %in% (sort(names(weight)) == sort(names(xtable))))
        stop("names of weight and table do not match!")
      else {
        for (i in 1:length(xtable)) {
          xtable[i] = weight[names(weight) == names(xtable)[i]] *
            xtable[i]
        }
      }
    }
  }
  else {
    weight = FALSE
  }
  if (missing(showTable))
    showTable = TRUE
  if (missing(xlab))
    xlab = ""
  if (missing(main))
    main = paste("Pareto Chart for", varName)
  if (missing(col))
    col = corp.col
  if (missing(border))
    border = corp.border
  if (missing(percentVec))
    percentVec = seq(0, 1, by = 0.25)
  call <- match.call(expand.dots = TRUE)
  # Plot
  if (length(xtable) > 1) {
    ylim = c(min(xtable), max(xtable) * 1.025)
    xtable = c(sort(xtable, decreasing = TRUE, na.last = TRUE))
    cumFreq = cumsum(xtable)
    sumFreq = sum(xtable)
    percentage = xtable/sum(xtable) * 100
    cumPerc = cumFreq/sumFreq * 100


    data <- data.frame(Frequency = xtable,
                       Cum.Frequency = cumFreq,
                       Percentage = round(percentage, digits = 2),
                       Cum.Percentage = round(cumPerc, digits = 2))
    tabla <- t(data)

    p <- ggplot(data, aes(x = reorder(names(xtable), -xtable), y = Frequency)) +
      geom_col(aes(fill = "Frequency"), width = 0.7) +
      geom_point(aes(y = Cum.Frequency, color = "Cumulative Percentage")) +
      geom_line(aes(y = Cum.Frequency, group = 1, color = "Cumulative Percentage")) +
      scale_y_continuous(name = ylab,
                         sec.axis = sec_axis(~ . / sum(xtable),
                                             name = "Cumulative Percentage",
                                             labels = percentVec)) +
      scale_x_discrete(name = xlab) +
      scale_color_manual(values = c(border, border)) +
      scale_fill_manual(values = col) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = main)+theme(plot.title = element_text(hjust = 0.5,face = "bold"))

  }
  else {
    warning("data should have at least two categories!")
  }
  if(showPlot == TRUE){
    if(showTable == TRUE){
      show(p/tableGrob(tabla))

    }
    else {
      show(p)
    }
  }
  else{
    show(tabla)
  }

  invisible(list(plot = p, table = tabla))
}

# qqPlot ---------------------
qqPlot <- function(x, y, confbounds = TRUE, alpha, main, xlab, ylab, xlim, ylim, border = "red",
                   bounds.col = "black", bounds.lty = 1, start, grapic = TRUE, axis.y.right = FALSE, bw.theme = FALSE,...){
  parList = list()
  if (is.null(parList[["col"]])){
    parList$col = 1:2
  }
  if (is.null(parList[["pch"]])){
    parList$pch = 19
  }
  if (is.null(parList[["lwd"]])){
    parList$lwd = 0.5
  }
  if (is.null(parList[["cex"]])){
    parList$cex = 1
  }
  if(inherits(x, "DistrCollection")){
    distList <- x$distr
    grap <- qqPlot(distList[[1]]$x, grapic = FALSE, ylab = "", xlab = "", main = paste(distList[[1]]$name,"distribution"))
    for (i in 2:length(distList)){
      aux <- qqPlot(distList[[i]]$x, grapic = FALSE, ylab = "", xlab = "", main = paste(distList[[i]]$name,"distribution"))
      grap$plot <-  grap$plot + aux$plot
    }
    show(grap$plot + plot_annotation(title = "QQ Plot for a Collection Distribution"))
    invisible()
  }
  else{
    if (missing(y))
      y = "normal"
    if(missing(alpha))
      alpha = 0.05
    if (alpha <=0 || alpha >=1)
      stop(paste("alpha should be between 0 and 1!"))
    if (missing(main))
      main = paste("QQ Plot for", deparse(substitute(y)), "distribution")
    if (missing(xlab))
      xlab = paste("Quantiles for", deparse(substitute(x)))
    if (missing(ylab))
      ylab = paste("Quantiles from", deparse(substitute(y)), "distribution")
    if (is.numeric(y)) {
      cat("\ncalling (original) qqplot from namespace stats!\n")
      return(stats::qqplot(x, y, ...))
    }
    qFun = NULL
    theoretical.quantiles = NULL
    xs = sort(x)
    distribution = tolower(y)
    distWhichNeedParameters = c("weibull", "logistic", "gamma","exponential", "f",
                                "geometric", "chi-squared", "negative binomial",
                                "poisson")

    threeParameterDistr = c("weibull3", "lognormal3", "gamma3")
    threeParameter = distribution %in% threeParameterDistr
    if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)
    if(is.character(distribution)){
      qFun = .charToDistFunc(distribution, type = "q")
      if (is.null(qFun))
        stop(paste(deparse(substitute(y)), "distribution could not be found!"))
    }
    # Puntos teoricos
    theoretical.probs = ppoints(xs)
    # Quantiles
    xq = NULL
    yq = quantile(xs, prob = c(0.25, 0.75))


    dots <- list(...)

    if(TRUE){
      fitList = .lfkp(parList, formals(qFun))
      fitList$x = xs
      fitList$densfun = distribution
      if(!missing(start))
        fitList$start = start
      if(!threeParameter){
        fittedDistr = do.call(FitDistr, fitList)
        parameter = fittedDistr$estimate

        #save the distribution parameter#
        thethas = fittedDistr$estimate
        # save the cariance-covariance matrix
        varmatrix = fittedDistr$vcov

      }else{
        parameter = do.call(paste(".",distribution, "3", sep = ""), list(xs) )    ####
        threshold = parameter$threshold
      }

      parameter = .lfkp(as.list(parameter), formals(qFun))
      params = .lfkp(parList, formals(qFun))
      parameter = .lfrm(as.list(parameter), params)
      parameter = c(parameter, params)
      theoretical.quantiles = do.call(qFun, c(list(c(theoretical.probs)), parameter))

      if(!threeParameter){
        confIntCapable = c("exponential", "log-normal", "logistic", "normal", "weibull", "gamma", "beta", "cauchy")
        getConfIntFun = .charToDistFunc(distribution, type = ".confint")
        if(confbounds == TRUE){
          if(distribution %in% confIntCapable){
            confInt = getConfIntFun(xs, thethas, varmatrix, alpha)
          }
        }
      }

      xq <- do.call(qFun, c(list(c(0.25, 0.75)), parameter))
    }
    else {
      params =.lfkp(parList, formals(qFun))
      params$p = theoretical.probs
      theoretical.quantiles = do.call(qFun, params)
      params$p = c(0.25, 0.75)
      xq = do.call(qFun, params)
    }

    params =.lfkp(parList, c(formals(plot.default), par()))

    if(!threeParameter){
      params$y = theoretical.quantiles
    }else{
      params$y = theoretical.quantiles+threshold
    }
    params$x = xs
    params$xlab = xlab
    params$ylab = ylab
    params$main = main
    if (!(is.null(params$col[1]) || is.na(params$col[1])))
      params$col = params$col[1]

    params$lwd = 1

    ############ PLOT ############
    p <- ggplot(data = data.frame(x=params$x, y=params$y), mapping=aes(x=x, y=y)) +
      geom_point() + labs(x = xlab, y = ylab, title = main) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))

    params =.lfkp(parList, c(formals(abline), par()))
    params$a = 0
    params$b = 1
    params$col = border
    p <- p + geom_abline(intercept = params$a, slope = params$b, col = params$col, lwd = params$lwd)

    if(!threeParameter){
      if(confbounds == TRUE){
        if(distribution %in% confIntCapable){
          params =.lfkp(parList, c(formals(lines), par()))
          params$x = confInt[[3]]
          params$y = confInt[[1]]
          params$col = bounds.col
          params$lty = bounds.lty
          # La curva de abajo
          p <- p + geom_line(data = data.frame(x=params$x, y=params$y), aes(x = x, y = y),
                             col = params$col, lty = params$lty, lwd = params$lwd)
          params$x = confInt[[3]]
          params$y = confInt[[2]]
          params$col = bounds.col
          params$lty = bounds.lty
          # curva de arriba
          p <- p + geom_line(data = data.frame(x=params$x, y=params$y), aes(x = x, y = y),
                             col = params$col, lty = params$lty, lwd = params$lwd)
        }
      }
    }
    if(axis.y.right){
      p <- p + scale_y_continuous(position = "right")
    }
    if(bw.theme){
      p <- p + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    }
    if(main == ""){
      p <- p + labs(title = NULL)
    }
    if(grapic){
      show(p)
      invisible(list(x = theoretical.quantiles, y = xs, int = params$a, slope = params$b, plot = p))
    }
    else{
      invisible(list(x = theoretical.quantiles, y = xs, int = params$a, slope = params$b, plot = p))
    }
  }
}

# ppPlot ---------------------
ppPlot <- function (x, distribution, confbounds = TRUE, alpha, probs, main, xlab, ylab, xlim, ylim, border = "red", bounds.col = "black", bounds.lty = 1, start, grapic = TRUE, axis.y.right = FALSE, bw.theme = FALSE,...)
{
  conf.level = 0.95
  conf.lines = TRUE
  if (!(is.numeric(x) | inherits(x, "DistrCollection")))
    stop(paste(deparse(substitute(x)), " needs to be numeric or an object of class distrCollection"))
  parList = list(...)
  if (is.null(parList[["col"]]))
    parList$col = c("black", "red", "gray")
  if (is.null(parList[["pch"]]))
    parList$pch = 19
  if(is.null(parList[["lwd"]]))
    parList$lwd = 1
  if (is.null(parList[["cex"]]))
    parList$cex = 1
  qFun = NULL
  xq = NULL
  yq = NULL
  x1 = NULL
  if(missing(alpha))
    alpha = 0.05
  if (alpha <=0 || alpha >=1)
    stop(paste("alpha should be between 0 and 1!"))
  if (missing(probs))
    probs = ppoints(11)
  else if (min(probs) <= 0 || max(probs) >= 1)
    stop("probs should be values within (0,1)!")
  probs = round(probs, 2)
  if (is.numeric(x)) {
    x1 <- sort(na.omit(x))
    if (missing(xlim))
      xlim = c(min(x1) - 0.1 * diff(range(x1)), max(x1) +0.1 * diff(range(x1)))
  }
  if (missing(distribution))
    distribution = "normal"
  if (missing(ylim))
    ylim = NULL
  if (missing(main))
    main = paste("Probability Plot for", deparse(substitute(distribution)),
                 "distribution")
  if (missing(xlab))
    xlab = deparse(substitute(x))
  if (missing(ylab))
    ylab = "Probability"
  if(inherits(x, "DistrCollection")){
    distList <- x$distr
    grap <- ppPlot(distList[[1]]$x, grapic = FALSE, ylab = "", xlab = "", main = paste(distList[[1]]$name,"distribution"))
    for (i in 2:length(distList)){
      aux <- ppPlot(distList[[i]]$x, grapic = FALSE, ylab = "", xlab = "", main = paste(distList[[i]]$name,"distribution"))
      grap$plot <-  grap$plot + aux$plot
    }
    show(grap$plot + plot_annotation(title = "QQ Plot for a Collection Distribution"))
    invisible()
  }
  distWhichNeedParameters = c("weibull", "gamma", "logistic","exponential","f",
                              "geometric", "chi-squared", "negative binomial",
                              "poisson")
  # new
  threeParameterDistr = c("weibull3", "lognormal3", "gamma3")
  threeParameter = distribution %in% threeParameterDistr
  if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)
  # end new
  if (is.character(distribution)) {
    qFun = .charToDistFunc(distribution, type = "q")
    pFun = .charToDistFunc(distribution, type = "p")
    dFun = .charToDistFunc(distribution, type = "d")
    if (is.null(qFun))
      stop(paste(deparse(substitute(y)), "distribution could not be found!"))
  }
  dots <- list(...)
  if (TRUE) {
    fitList = .lfkp(parList, formals(qFun))
    fitList$x = x1
    fitList$densfun = distribution
    if (!missing(start))
      fitList$start = start
    if(!threeParameter){
      fittedDistr = do.call(FitDistr, fitList)
      parameter = fittedDistr$estimate
      #save the distribution parameter#
      thethas = fittedDistr$estimate
      # save the cariance-covariance matrix
      varmatrix = fittedDistr$vcov
    }else{
      parameter = do.call(paste(".",distribution, "3", sep = ""), list(x1) )    ####
      print(parameter[3])
      threshold = parameter$threshold
    }
    parameter = .lfkp(as.list(parameter), formals(qFun))
    params = .lfkp(parList, formals(qFun))
    parameter = .lfrm(as.list(parameter), params)
    print(parameter)
    parameter = c(parameter, params)
    # new
    if(!threeParameter){
      # array containing names of the distributions, for which conf intervals can be computed
      confIntCapable = c("exponential", "log-normal", "logistic", "normal", "weibull", "gamma", "beta", "cauchy")
      getConfIntFun = .charToDistFunc(distribution, type = ".confint")
      # if possible, compute the conf intervals
      if(confbounds == TRUE){
        if(distribution %in% confIntCapable){
          confInt = getConfIntFun(x1, thethas, varmatrix, alpha)
        }
      }# end of my code
    }
    y = do.call(qFun, c(list(ppoints(x1)), as.list(parameter)))
    yc = do.call(qFun, c(list(ppoints(x1)), as.list(parameter)))
    cv = do.call(dFun, c(list(yc), as.list(parameter)))
    print(cv)
    axisAtY = do.call(qFun, c(list(probs), as.list(parameter)))
    yq = do.call(qFun, c(list(c(0.25, 0.75)), as.list(parameter)))
    xq = quantile(x1, probs = c(0.25, 0.75))
  }
  else {
    params = .lfkp(parList, formals(qFun))
    params$p = ppoints(x1)
    y = do.call(qFun, params)
    params$p = probs
    axisAtY = do.call(qFun, params)
    params$p = c(0.25, 0.75)
    yq = do.call(qFun, params)
    xq = quantile(x1, probs = c(0.25, 0.75))
  }
  params = .lfkp(parList, c(formals(plot.default), par()))
  params$x = x1
  params$y = y
  params$xlab = xlab
  params$ylab = ylab
  params$main = main
  params$xlim = xlim
  params$axes = FALSE
  params$lwd = 1
  if (!(is.null(params$col[1]) || is.na(params$col[1])))
    params$col = params$col[1]
  # PLOT
  p <- ggplot(data.frame(x = x1, y = y), aes(x = x, y = y)) +
    geom_point(size = 1, color = "black") +
    labs(x = xlab, y = ylab, title = main) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(scale = 1/max(x1), suffix = "", accuracy = 0.01))+
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  # Regression line
  params = .lfkp(parList, c(formals(abline), par()))
  if(!threeParameter){
    params$a = 0
  }else{
    params$a = -threshold
  }
  params$b = 1
  params$col = border
  p <- p + geom_abline(intercept = params$a, slope = params$b, color = border)
  if(!threeParameter){
    if(confbounds == TRUE){
      if(distribution %in% confIntCapable){
        params =.lfkp(parList, c(formals(lines), par()))
        params$x = confInt[[3]]
        params$y = confInt[[1]]
        params$col = bounds.col
        params$lty = bounds.lty
        # La curva de abajo
        p <- p + geom_line(data = data.frame(x=params$x, y=params$y), aes(x = x, y = y),
                           col = params$col, lty = params$lty, lwd = params$lwd)
        params$x = confInt[[3]]
        params$y = confInt[[2]]
        params$col = bounds.col
        params$lty = bounds.lty
        # curva de arriba
        p <- p + geom_line(data = data.frame(x=params$x, y=params$y), aes(x = x, y = y),
                           col = params$col, lty = params$lty, lwd = params$lwd)
      }
    }
  }
  if(axis.y.right){
    p <- p + scale_y_continuous(position = "right")
  }
  if(bw.theme){
    p <- p + theme_bw()
  }
  if(main == ""){
    p <- p + labs(title = NULL)
  }
  if(grapic){
    show(p)
    invisible(list(x = x, y = y, int = params$a, slope = params$b, plot = p))
  }
  else{
    invisible(list(x = x, y = y, int = params$a, slope = params$b, plot = p))
  }
}

set.seed(1234)
x <- rnorm(20, mean = 20)
distribution <- "normal"
ppPlot(x, distribution, bounds.lty = 3, bounds.col = "red")

# paretoPlot ----
paretoPlot <- function(fdo, abs = TRUE, decreasing = TRUE, alpha = 0.05,
                       response = NULL, ylim, xlab, ylab, main, p.col, legend_left = TRUE) {
  # library(RColorBrewer)
  # Esta librería tiene los colores:
  # Set1, Set2, Set3, Pastel2, Pastel1,
  # Paired, Dark2, Accent

  if(is.null(response)==FALSE)
  {
    temp=fdo$.response()[response]
    fdo$.response(temp)
  }
  ylimMissing = FALSE
  if (missing(ylim)){
    ylimMissing = TRUE
  }
  if (missing(xlab))
    xlab = ""
  location = "topright"
  if (decreasing == F | abs == F | legend_left == T) {
    location = "topleft"
  }
  xVals = numeric(0)
  sig.neg = NULL
  sig.pos = NULL
  effect.list = vector("list")
  for (j in 1:ncol(fdo$.response())) {
    if (!any(is.na(fdo$.response()[, j]))) {
      if (missing(ylab))
        ylabel = names(fdo$.response())[j]                                ###
      else
        ylabel = ylab                                                   ###
      form = paste("fdo$.response()[,", j, "]~")
      for (i in 1:ncol(fdo$cube)) {
        form = paste(form, names(fdo$cube)[i], sep = "")
        if (i < ncol(fdo$cube))
          form = paste(form, "*", sep = "")
      }
      lm.1 = lm(as.formula(form), data = fdo$as.data.frame())
      coefs = coef(lm.1)[-pmatch("(Intercept)", names(coef(lm.1)))]
      df.resid = df.residual(lm.1)
      num.c = nrow(fdo$centerCube)

      if (df.resid == 0) {
        effect = 2 * coefs
        effect = effect[!is.na(effect)]
        effect.list[[j]] = effect
        if (missing(main))
          main = "Lenth Plot of effects"
        limits = TRUE
        faclab = NULL
        m = length(effect)
        d = m/3
        s0 = 1.5 * median(abs(effect))
        rmedian = effect[abs(effect) < 2.5 * s0]
        PSE = 1.5 * median(abs(rmedian))
        ME = qt(1 - alpha/2, d) * PSE
        Gamma = (1 + (1 - alpha)^(1/m))/2
        SME = qt(Gamma, d) * PSE
        n = length(effect)

        if(missing(p.col)){
          p.col = rep("lightblue", length(effect))
        }
        else{p.col = brewer.pal(length(effect), p.col)} #paste0("Set", p.col))

        if (ylimMissing)
          if (abs)
            ylim = (range(c(0, abs(effect), 1.3 * ME))) * 1.1
        else ylim = (range(c(effect, -1.3 * ME, 1.3 * ME))) * 1.1
        if (abs) {
          if (missing(ylabel))
            ylabel = ""

          p <- ggplot(data.frame(names = factor(names(effect), levels = names(effect)), effect_ = abs(as.vector(effect))),
                      aes(x = names, y = effect_, fill = names)) +
            geom_bar(stat = "identity", color = "black") +
            scale_fill_manual(values=c(p.col)) +
            theme_minimal()+
            geom_text(aes(label = round(effect,2)), vjust = -1, colour = "black") + # etiquetas sobre las barras
            labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  plot.title = element_text(hjust = 0.5),
                  legend.position = "none")

          if(ME >= ylim[1] & ME <= ylim[2] & SME >= ylim[1] & SME <= ylim[2]){
            p <- p +
              geom_hline(yintercept = ME, linetype = "dashed", color = "red") +
              geom_hline(yintercept = SME, linetype = "dashed", color = "red") +
              scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(abs(ME), abs(SME)),
                                                                                     labels = c(paste("ME = ", round(abs(ME), 2)), paste("SME = ", round(abs(SME), 2)))
              ))

          }
          else{
            if(SME >= ylim[1] & SME <= ylim[2]){
              p <- p + geom_hline(yintercept = SME, linetype = "dashed", color = "red") +
                scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(abs(SME)), labels = c(paste("SME = ", round(abs(SME), 2)))))
            }
            if(ME >= ylim[1] & ME <= ylim[2]){
              p <- p + geom_hline(yintercept = ME, linetype = "dashed", color = "red") +
                scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(abs(ME)), labels = c(paste("ME = ", round(abs(ME), 2)))))
            }
          }

        }
        else {
          if (missing(ylabel))
            ylabel = ""
          p <- ggplot(data.frame(names = names(effect), effect_ = abs(as.vector(effect))),
                      aes(x = names, y = effect_, fill = names)) +
            geom_bar(stat = "identity", color = "black") +
            scale_fill_manual(values=c(p.col)) +
            theme_minimal()+
            geom_text(aes(label = round(effect,2)), vjust = -1, colour = "black") + # etiquetas sobre las barras
            labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  plot.title = element_text(hjust = 0.5))

          if(ME >= ylim[1] & ME <= ylim[2] & SME >= ylim[1] & SME <= ylim[2]){
            p <- p +
              geom_hline(yintercept = ME, linetype = "dashed", color = "red") +
              geom_hline(yintercept = SME, linetype = "dashed", color = "red")
            scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(ME, SME), labels = c(paste("ME = ", round(abs(ME), 2)), paste("SME = ", round(abs(SME), 2)) )))

          }
          else{
            if(ME >= ylim[1] & ME <= ylim[2]){
              p <- p + geom_hline(yintercept = ME, linetype = "dashed", color = "red") +
                scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(ME), labels = c(paste("ME = ", round(ME, 2)))))

            }
            if(SME >= ylim[1] & SME <= ylim[2]){
              p <- p + geom_hline(yintercept = SME, linetype = "dashed", color = "red") +
                scale_y_continuous(limits = ylim, expand = c(0, 0),sec.axis = sec_axis(~ ., breaks = c(SME), labels = c(paste("SME = ", round(SME, 2)))))
            }
          }



        }
      }

      else {
        if (missing(main))
          main = "Standardized main effects and interactions"
        effect = ((summary(lm.1)$coefficients[-pmatch("(Intercept)", names(coef(lm.1))), 1])/(summary(lm.1)$coefficients[-pmatch("(Intercept)", names(coef(lm.1))),2]))
        if (all(is.na(effect)))
          stop("effects could not be calculated")
        effect = effect[!is.na(effect)]
        effect.list[[j]] = effect
        if ((df.resid) > 0) {
          sig.pos = -qt(alpha/2, df.resid)
          sig.neg = +qt(alpha/2, df.resid)
        }
        # Ylimits ----

        if (ylimMissing)
          if (abs) {
            tempVec = c(effect, sig.pos)
            tempVec = tempVec[!is.na(tempVec)]
            ylim = c(0, 0.3 + max(tempVec))
          }
        else {
          tempVec1 = c(0, effect, sig.neg, sig.pos)
          tempVec1 = tempVec1[!is.na(tempVec1)]
          tempVec2 = c(abs(effect), sig.pos, sig.neg)
          tempVec2 = tempVec2[!is.na(tempVec2)]
          ylim = c(min(tempVec1)-0.3, max(tempVec2)+0.3)
        }

        if(missing(p.col)){
          p.col = rep("lightblue", length(effect))
        }
        else{p.col = brewer.pal(length(effect), p.col)} #paste0("Set", p.col))
        # Plot ---------
        if (abs) {
          effect = effect[order(abs(effect), na.last = TRUE, decreasing = decreasing)]
          effect = round(effect, 3)

          if (missing(ylabel))
            ylabel = ""

          # plot with abs
          df <- data.frame(Names = factor(names(effect), levels = names(effect)),
                           effect_ = abs(as.vector(effect)))

          p <- ggplot(data = df,
                      aes(x = Names, y = effect_, fill = Names)) +
            geom_bar(stat = "identity", color = "black") +
            scale_fill_manual(values=c(p.col)) +
            theme_minimal()+
            labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  plot.title = element_text(hjust = 0.5)) +
            geom_text(aes(label = round(effect,2)), vjust = -1, colour = "black") + # etiquetas sobre las barras
            geom_hline(yintercept = sig.pos, linetype = "dashed", color = "red") +
            annotate("text", x = max(as.numeric(df$Names)), y = sig.pos, label = paste(round(sig.pos, 2)), vjust = -0.5, color = "red")
        }
        else {
          effect = effect[order((effect), na.last = TRUE, decreasing = decreasing)]
          effect = round(effect, 3)

          if (missing(ylabel))
            ylabel = ""

          df <- data.frame(Names = factor(names(effect), levels = names(effect)),
                           effect_ = as.vector(effect))

          # Plot without abs
          p <- ggplot(df,
                      aes(x = Names, y = effect_, fill = Names)) +
            geom_bar(stat = "identity", color = "black") +
            scale_fill_manual(values=c(p.col)) +
            theme_minimal()+
            labs(title = main, x = xlab, y = ylabel) + ylim(c(ylim)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0)) +
            geom_text(aes(label = effect), vjust = ifelse(df$effect_ > 0, -0.5, 1.5) , colour = "black") + # etiquetas sobre las barras
            geom_hline(yintercept = sig.pos, linetype = "dashed", color = "red") +
            geom_hline(yintercept = sig.neg, linetype = "dashed", color = "red") +
            annotate("text", x = max(as.numeric(df$Names)), y = sig.pos, label = paste(round(sig.pos, 2)), vjust = -0.5, color = "red") +
            annotate("text", x = max(as.numeric(df$Names)), y = sig.neg, label = paste(round(sig.neg, 2)), vjust = 1.5, color = "red")
        }
        myDelta = diff(range(ylim)) * 0.02
      }
      # Legend ----
      titles <- data.frame(Name_title = paste0(names(fdo$cube),": ",fdo$names()))
      for(i in 1:dim(titles)[1]){
        titles$Pos_title[i] <- 0.95 - (0.05 * i)
      }
      caja <- ggplot(data.frame(x = 0,y = 0), aes(x = x, y = y)) +
        theme_bw() +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = -0.5,margin = margin(b = -12), size = 10)
        ) + xlim(c(0.24, 0.26)) + ylim(c(min(titles$Pos_title) - 0.05, max(titles$Pos_title) + 0.05))
      for(i in 1:dim(titles)[1]){
        caja <- caja + annotate("text", x = 0.25, y = titles$Pos_title[i], label = titles$Name_title[i], size = 3.5, hjust = 0.5)
      }
      # insert legend -----
      if(location == "topright"){
        p <- p + inset_element(caja, left = 0.75, right = 1, top = 1,  bottom = 0.80)
      }
      else{
        p <- p + inset_element(caja, left = 0.25, right = 0.05, top = 1,  bottom = 0.80)
      }

    }
  }

  print(p)
  invisible(list(effect.list, plot = p))
}
# Uso paretoPlot
paretoPlot(dfac, decreasing = T, abs = F, p.col = "Pastel1")

# normalPlot ----
normalPlot <- function(fdo, threeWay = FALSE, na.last = NA, alpha = 0.05, response = NULL, sig.col = c("red1", "red2", "red3"), sig.pch = c(1,2,3), main, ylim, xlim, xlab, ylab, pch,  ###
                       col, border = "red", ...){
  fdoName = deparse(substitute(fdo))
  if(is.null(response)==FALSE)
  {
    temp=fdo$.response()[response]
    fdo$.response(temp)
  }
  parList = list()
  parList = list(...)
  if (length(sig.col) < 3)
    sig.col = as.vector(matrix(sig.col, nrow = 1, ncol = 3))
  XLIM=FALSE;YLIM=FALSE
  if (!(class(fdo)[1] == "facDesign"))
    stop(paste(deparse(substitute(fdo)), "is not an object of class facDesign"))
  mainmiss = FALSE
  if (missing(main))
    mainmiss = TRUE
  if (missing(ylim))
    YLIM=TRUE
  if (missing(xlim))
    XLIM=TRUE
  if (missing(xlab))
    xlab = "Coefficients"
  if (missing(ylab))
    ylab = "Theoretical Quantiles"
  if (missing(pch))
    pch = 19
  if (missing(col))
    col = "black"

  list_plot = list()
  for(j in 1:ncol(fdo$.response())){
    parList = list(...)
    params = list()
    leg.col = vector()
    p.col = vector()
    p.pch = vector()
    leg.txt = vector()
    main = paste("Normal plot for", names(fdo$.response())[j], "in", fdoName)
    if (j > 1)
      dev.new()
    form = paste("fdo$.response()[,", j, "]~")

    for (i in 1:ncol(fdo$cube)) {
      form = paste(form, names(fdo$cube)[i], sep = "")
      if (i < ncol(fdo$cube))
        form = paste(form, "*", sep = "")
    }

    lm.1 = lm(as.formula(form), data = fdo$as.data.frame())
    lm.1s = summary(lm.1)
    effect = coef(lm.1s)[row.names(coef(lm.1s)) != "(Intercept)", "t value"]
    print(effect)
    if (all(is.na(effect)))
      effect = 2 * coef(lm.1)[-pmatch("(Intercept)", names(coef(lm.1)))]
    #            stop("effects could not be calculated")
    sig = summary(lm.1)$coefficients[, "Pr(>|t|)"][-pmatch("(Intercept)", names(coef(lm.1)))]
    df.resid = df.residual(lm.1)
    nc = nrow(fdo$centerCube)

    tQ = ppoints(effect)
    index = order(effect)
    sQ = effect[index]
    sig = sig[index]

    if (df.resid > 0) {
      # obtenemos el caracter de la cajita del p_value
      for (k in seq(along = sig)) {
        setted = FALSE
        if (abs(sig)[k] < 0.01) {
          if (!setted) {
            p.col[k] = sig.col[1]
            p.pch[k] = sig.pch[1]
            leg.txt = c(leg.txt, "p < 0.01")
            leg.col = c(leg.col, p.col)
            setted = TRUE
          }
        }
        if (abs(sig)[k] < 0.05) {
          if (!setted) {
            p.col[k] = sig.col[2]
            p.pch[k] = sig.pch[2]
            leg.txt = c(leg.txt, "p < 0.05")
            leg.col = c(leg.col, p.col)
            setted = TRUE
          }
        }
        if (abs(sig)[k] < 0.1) {
          if (!setted) {
            p.col[k] = sig.col[3]
            p.pch[k] = sig.pch[3]
            leg.txt = c(leg.txt, "p < 0.1")
            leg.col = c(leg.col, p.col)
            setted = TRUE
          }
        }
        if (abs(sig)[k] >= 0.1) {
          if (!setted) {
            p.col[k] = col
            p.pch[k] = pch
            leg.txt = c(leg.txt, "p >= 0.1")
            leg.col = c(leg.col, p.col)
            setted = TRUE
          }
        }
      }
      leg.txt = unique(leg.txt)
      leg.col = unique(leg.col)
    }else{p.col=col
    p.pch=pch}

    mid = round(length(tQ)/2)
    last = length(tQ)
    params$p = ppoints(effect)
    estimates = FitDistr(effect, "normal")   #estimates = MASS::fitdistr(effect, "normal")
    params$mean = estimates$estimate[["mean"]]
    params$sd = estimates$estimate[["sd"]]

    y = do.call(qnorm, params)

    if (XLIM)
      xlim = range(sQ)
    if (YLIM)
      ylim = range(y)

    # PLOT -----------------------
    df <- data.frame(sQ = names(sQ), value = as.numeric(sQ), y = y)

    p <- ggplot(df, aes(x = value, y = y, label = sQ)) +
      geom_point(col = p.col, pch = p.pch) +
      theme_classic() + lims(x = xlim, y = ylim) +
      labs(x = xlab, y = ylab, title = main) +
      geom_text(check_overlap = TRUE, vjust = 1) + theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))

    xp = c(qnorm(0.1), qnorm(0.99))
    yp = c(qnorm(0.1, mean = estimates$estimate[["mean"]], sd = estimates$estimate[["sd"]]), qnorm(0.99, mean = estimates$estimate[["mean"]], sd = estimates$estimate[["sd"]]))
    slope = (yp[2] - yp[1])/(xp[2] - xp[1])
    int = yp[1] - slope * xp[1]

    # line
    p <- p + geom_abline(intercept = int, slope = slope, col = border)

    # legend
    if (df.resid > 0){
      caja <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
        theme_bw() +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        xlim(c(0.25,0.30)) + ylim(c(0.24, 0.31))

      caja <- caja +
        annotate('text', x = 0.275, y = 0.28,
                 label = leg.txt, size = 3, hjust = 0.5, colour = leg.col)

      p <- p + inset_element(caja, left = 0.01, right = 0.2, top = 1, bottom = 0.85)
    }
    print(p)
    list_plot[[paste0("p",j)]] <- p
  }
  invisible(list(effect = effect, plots = list_plot))
}
# Uso normalPlot
normalPlot(dfac)


# InteractionPlot ----
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
  y <- fdo$.response()[1:nrow(x), ]
  numFac <- ncol(x)
  combMat <- combn(names(x), 2)

  plot_list <- list()

  for (r in 1:ncol(fdo$.response())) {
    y = fdo$.response()[1:nrow(x), r]

    for (i in 1:ncol(combMat)) {
      facName1 <- combMat[1, i]
      facName2 <- combMat[2, i]

      df = data.frame(fac1 = x[[facName1]], fac2 = x[[facName2]], response = y)

      p = ggplot(df, aes_string(x = "fac2", y = "response", color = "as.factor(fac1)")) +
        geom_line(stat = "summary", fun = fun,size=1.5) +
        labs(x = " ", y = " ", color = facName1) +
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
          facName <- names(x)[i]
          p_diag <- ggplot() +
            labs(x = facName, y = "") +
            theme_void() +
            theme(
              plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
              axis.title.x = element_text(size = 20, face = "bold", margin = margin(0, 0, 10, 0)),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()
            )
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
    final_plot <- final_plot + plot_annotation(
      title = main[r],
      theme = theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 20)))
    )
    print(final_plot)
  }

  invisible()
}
# Uso interactionPlot
interactionPlot(dfac)
# wirePlot ----
wirePlot <- function(x, y, z, data = NULL,
                     xlim, ylim, zlim, main,
                     xlab, ylab, sub, sub.a = TRUE, zlab,
                     form = "fit", col = "Rainbow", steps,
                     factors, fun, plot = TRUE, show.scale = TRUE,
                     n.scene = "scene") {
  form = form
  fact = NULL
  if (missing(steps))
    steps = 25
  fdo = data
  fit = NULL
  lm.1 = NULL

  # Col puede ser: "Rainbow", "Jet", "Earth", "Electric"
  if (is.null(data)) {
    cat("\n defaulting to persp function\n")
    return("persp")
  }
  if (class(data)[1] != "facDesign") {
    cat("\n defaulting to persp function using formula\n")
    return("persp")
  }

  x.c = deparse(substitute(x))
  y.c = deparse(substitute(y))
  z.c = deparse(substitute(z))

  if (missing(main))
    main = paste("Response Surface for", z.c)

  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[.NAMES[i]]] <-fdo$names()[i]
  }
  if (missing(ylab))
    ylab = paste(y.c, ": ", aux[[y.c]])
  if (missing(xlab))
    xlab = paste(x.c, ": ", aux[[x.c]])
  if (missing(zlab))
    zlab = paste(x.c, ": ", z.c)

  if (missing(xlim))
    xlim = c(min(fdo$get(, x.c)), max(fdo$get(, x.c)))
  if (missing(ylim))
    ylim = c(min(fdo$get(, y.c)), max(fdo$get(, y.c)))

  allVars = c(fdo$names(), names(fdo$.response()))
  isct = intersect(c(aux[[x.c]], aux[[y.c]], z.c), c(fdo$names(), names(fdo$.response())))

  if (length(isct) < length(c(x.c, y.c, z.c))) {
    d = setdiff(isct, allVars)
    stop(paste(d, "could not be found\n"))
  }

  if (missing(fun))
    fun = NULL
  if (!is.function(fun) & !is.null(fun))
    if (!(fun %in% c("overall", "desirability")))
      stop("fun should be a function, \"overall\" or \"desirability\"")
  if (identical(fun, "desirability")) {
    obj = fdo$desires()[[z.c]]
    fun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
  }

  if (form %in% c("fit")) {
    lm.1 = fdo$fits[[z.c]]
    if (is.null(fit))
      form = "full"
  }

  if (form %in% c("quadratic", "full", "interaction", "linear")) {
    if (identical(form, "full")) {
      form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
      if (nrow(fdo$star) > 0)
        form = paste(form, "+ I(", x.c, "^2) + I(", y.c, "^2)")
    }
    if (identical(form, "interaction")) {
      form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
    }
    if (identical(form, "linear")) {
      form = paste(z.c, "~", x.c, "+", y.c)
    }
    if (identical(form, "quadratic")) {
      form = paste(z.c, "~I(", x.c, "^2) + I(", y.c, "^2)")
    }
  }

  if (is.null(form))
    stop(paste("invalid formula", form))
  if (is.null(lm.1))
    lm.1 = fdo$lm(form)
  if (missing(sub))
    sub = deparse(formula(lm.1))

  dcList = vector(mode = "list", length = length(fdo$names()))
  names(dcList) = names(aux)
  dcList[1:length(fdo$names())] = 0

  help.predict = function(x, y, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = x
    dcList[[y.c]] = y
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
  }

  xVec = seq(min(xlim), max(xlim), length = steps)
  yVec = seq(min(ylim), max(ylim), length = steps)

  mat = outer(xVec, yVec, help.predict, x.c, y.c, lm.1)

  if (is.function(fun))
    mat = try(apply(mat, c(1, 2), fun))
  if (identical(fun, "overall")) {
    main = "composed desirability"
    mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
    for (i in names(fdo$.response())) {
      obj = fdo$desires()[[i]]
      fun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
      temp = outer(xVec, yVec, help.predict, x.c, y.c, fits(fdo)[[i]])
      temp = try(apply(temp, c(1, 2), fun))
      mat = mat * temp
    }
    mat = mat^(1/length(names(fdo$response())))
  }

  if (missing(zlim))
    zlim = range(mat)

  p <- plot_ly(x = -yVec, y = xVec, z = mat, colorscale=col, scene = n.scene) %>%
    add_surface(showscale = show.scale) %>%
    layout(
      title = main,
      scene = list(
        xaxis = list(range = ylim, title = ylab, zeroline = FALSE),
        yaxis = list(range = xlim, title = xlab, zeroline = FALSE),
        zaxis = list(range = zlim, title = zlab, zeroline = FALSE),
        camera = list(eye = list(x=2, y=2, z=0.1))
      ),
      margin = list(l = 10, r = 15, t = 30, b = 20)
    )
  if(sub.a){
    p <- p %>%
      layout(
        annotations = list(
          list(
            text = sub,# Subtitulo
            x = 0.5,   # Posición x en la mitad de la gráfica
            y = -0.1,  # Posición y debajo de la gráfica
            showarrow = FALSE,
            font = list(size = 12)
          )
        )
      )
  }
  if (plot) {
    show(p)
  }
  invisible(list(x = xVec, y = yVec, z = mat, plot = p))
}
# Uso wirePlot
wirePlot(A,B,rend,data=dfac)

# contourPlot ----
contourPlot <- function(x, y, z, data = NULL, xlim, ylim, main, xlab, ylab, zlab, form = "fit", col = 1, steps,
                        factors, fun, plot = TRUE, show.scale = TRUE) {
  form = form
  fact = NULL
  if (missing(steps))
    steps = 25
  fdo = data
  fit = NULL
  lm.1 = NULL
  if (!is.function(col)) {
    if (identical(col, 1))
      col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if (identical(col, 2))
      col = colorRampPalette(c("blue", "white", "red"), space = "Lab")
    if (identical(col, 3))
      col = colorRampPalette(c("blue", "white", "orange"))
    if (identical(col, 4))
      col = colorRampPalette(c("gold", "white", "firebrick"))
    if (identical(col, 5))
      col = colorRampPalette(c("blue4", "lightblue1", "lightgreen", "green4"))
  }

  if (is.null(data) | class(data)[1] != "facDesign") {
    if(length(x) == length(y)){
      if(dim(z)[1] == length(x) & dim(z)[2] == length(x)){
        if(missing(main))
          main = "Filled Contour"
        if(missing(xlab))
          xlab = ""
        if(missing(ylab))
          ylab = ""

        if (is.function(col)) {
          mat <- z
          nrMat <- nrow(mat)
          ncMat <- ncol(mat)
          nbcol <- 1000
          color <- col(nbcol)
          matFacet <- mat[-1, -1] + mat[-1, -ncMat] + mat[-nrMat, -1] + mat[-nrMat, -ncMat]
          facetcol <- cut(matFacet, nbcol)
        }
        else {
          color = col
          facetcol = 1
        }

        p <- plot_ly(x = x, y = y, type = "contour", z = z, autocontour = TRUE, colors = color,
                     contours = list(coloring = 'heatmap'), line = list(smoothing = 0),
                     showscale = show.scale) %>%
          layout(
            title = main,
            xaxis = list(title = xlab, zeroline = FALSE),
            yaxis = list(title = ylab, zeroline = FALSE)
          )

        if (plot) {
          show(p)
        }

        invisible(list(x = x, y = y, z = z, plot = p))
      }

    }
    else{
      cat("\n defaulting to filled.contour function\n")
      return("persp")
    }
  }
  else{x.c = deparse(substitute(x))
  y.c = deparse(substitute(y))
  z.c = deparse(substitute(z))

  aux <- list()
  for (i in 1:length(fdo$names())) {
    aux[[.NAMES[i]]] <-fdo$names()[i]
  }
  if (missing(plot))
    plot = TRUE
  if (missing(main))
    main = paste("Filled Contour for", z.c)
  if (missing(ylab))
    ylab = paste(y.c, ": ", aux[[y.c]])
  if (missing(xlab))
    xlab = paste(x.c, ": ", aux[[x.c]])
  if (missing(factors))
    factors = NULL
  if (missing(xlim))
    xlim = c(min(fdo$get(, x.c)), max(fdo$get(, x.c)))
  if (missing(ylim))
    ylim = c(min(fdo$get(, y.c)), max(fdo$get(, y.c)))
  allVars = c(fdo$names(), names(fdo$.response()))
  isct = intersect(c(aux[[x.c]], aux[[y.c]], z.c), c(fdo$names(), names(fdo$.response())))

  if (length(isct) < length(c(x.c, y.c, z.c))) {
    d = setdiff(isct, allVars)
    stop(paste(d, "could not be found\n"))
  }
  if (missing(fun))
    fun = NULL
  if (!is.function(fun) & !is.null(fun))
    if (!(fun %in% c("overall", "desirability")))
      stop("fun should be a function, \"overall\" or \"desirability\"")
  if (identical(fun, "desirability")) {
    obj = fdo$desires()[[z.c]]
    fun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
  }
  if (form %in% c("fit")) {
    lm.1 = fdo$fits[[z.c]]
    if (is.null(fit))
      form = "full"
  }
  if (form %in% c("quadratic", "full", "interaction", "linear")) {
    if (identical(form, "interaction")) {
      form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
    }
    if (identical(form, "linear")) {
      form = paste(z.c, "~", x.c, "+", y.c)
    }
    if (identical(form, "quadratic")) {
      form = paste(z.c, "~I(", x.c, "^2) + I(", y.c, "^2)")
    }
    if (identical(form, "full")) {
      form = paste(z.c, "~", x.c, "+", y.c, "+", x.c, ":", y.c)
      if (nrow(fdo$star) > 0)
        form = paste(form, "+ I(", x.c, "^2) + I(", y.c, "^2)")
    }
  }

  if (is.null(form))
    stop(paste("invalid formula", form))
  if (is.null(lm.1))
    lm.1 = fdo$lm(form)

  dcList = vector(mode = "list", length = length(fdo$names()))
  names(dcList) = names(aux)
  dcList[1:length(fdo$names())] = 0

  if (!is.null(factors)) {
    for (i in fdo$names()) dcList[[i]] = factors[[i]][1]
  }

  help.predict = function(x, y, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = x
    dcList[[y.c]] = y
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
  }

  xVec = seq(min(xlim), max(xlim), length = steps)
  yVec = seq(min(ylim), max(ylim), length = steps)
  mat = outer(xVec, yVec, help.predict, x.c, y.c, lm.1)

  if (is.function(col)) {
    nrMat <- nrow(mat)
    ncMat <- ncol(mat)
    nbcol <- 1000
    color <- col(nbcol)
    matFacet <- mat[-1, -1] + mat[-1, -ncMat] + mat[-nrMat, -1] + mat[-nrMat, -ncMat]
    facetcol <- cut(matFacet, nbcol)
  }
  else {
    color = col
    facetcol = 1
  }
  if (is.function(fun))
    mat = try(apply(mat, c(1, 2), fun))
  if (identical(fun, "overall")) {
    main = "composed desirability"
    mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
    for (i in names(response(fdo))) {
      obj = fdo$desires()[[i]]
      fun = .desireFun(obj$low, obj$high, obj$target, obj$scale, obj$importance)
      temp = outer(xVec, yVec, help.predict, x.c, y.c, fdo$fits[[i]])
      temp = try(apply(temp, c(1, 2), fun))
      mat = mat * temp
    }
    mat = mat^(1/length(names(fdo$.response())))
  }

  p <- plot_ly(x = xVec, y = yVec, z = mat, colors = color,
               type = "contour", autocontour = TRUE, line = list(smoothing = 0), contours = list(coloring = 'heatmap'),
               showscale = show.scale) %>%
    layout(
      title = main,
      xaxis = list(range = xlim, title = xlab, zeroline = FALSE),
      yaxis = list(range = ylim, title = ylab, zeroline = FALSE)
    )

  if (plot) {
    show(p)
  }

  invisible(list(x = xVec, y = yVec, z = mat, plot = p))
  }

}
#  Uso contourPlot
contourPlot(A,B,rend,data=dfac)


############### CHOOSE
# fracChoose ----
fracChoose <- function() {
  genList = list(6 * 9)
  genList = list(c("C = AB"), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c("D = ABC"), c("D = AB", "E = AC"), c("D = AB",
                                                                                                                                                      "E = AC", "F = BC"), c("D = AB", "E = AC", "F = BC", "G = ABC"), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c("E = ABCD"), c("E = ABC", "F = BCD"),
                 c("E = ABC", "F = BCD", "G = ACD"), c("E = BCD", "F = ACD", "G = ABC", "H = ABD"), c("E = ABC", "F = BCD", "G = ACD", "H = ABD", "J = ABCD"), c("E = ABC",
                                                                                                                                                                 "F = BCD", "G = ACD", "H = ABD", "J = ABCD", "K = AB"), c("E = ABC", "F = BCD", "G = ACD", "H = ABD", "J = ABCD", "K = AB", "L = AC"), c(NULL),
                 c(NULL), c(NULL), c("F = ABCDE"), c("F = ABCD", "G = ABDE"), c("F = ABC", "G = ABD", "H = BCDE"), c("F = BCDE", "G = ACDE", "H = ABDE", "J = ABCE"),
                 c("F = ABCD", "G = ABCE", "H = ABDE", "J = ACDE", "K = BCDE"), c("F = ABC", "G = BCD", "H = CDE", "J = ACD", "K = AEF", "L = ADEF"), c(NULL), c(NULL),
                 c(NULL), c(NULL), c("G = ABCDEF"), c("G = ABCD", "H = ABEF"), c("G = ABCD", "H = ACEF", "J = CDEF"), c("G = BCDF", "H = ACDF", "J = ABDE", "K = ABCE"),
                 c("G = CDE", "H = ABCD", "J = ABF", "K = BDEF", "L = ADEF"), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c("H = ABCDEFG"), c("H = ACDFG", "J = BCEFG"),
                 c("H = ABCG", "J = BCDE", "K = ACDF"), c("H = ABCG", "J = BCDE", "K = ACDF", "L = ABCDEFG"))
  resList = list(6 * 9)
  resList = list(c(3), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(4), c(3), c(3), c(3), c(NULL), c(NULL), c(NULL),
                 c(NULL), c(NULL), c(NULL), c(5), c(4), c(4), c(4), c(3), c(3), c(3), c(NULL), c(NULL), c(NULL), c(6), c(4), c(4), c(4), c(4), c(4), c(NULL), c(NULL),
                 c(NULL), c(NULL), c(7), c(5), c(4), c(4), c(4), c(NULL), c(NULL), c(NULL), c(NULL), c(NULL), c(8), c(6), c(5), c(5))
  facMat = matrix(rep(3:11, 6), ncol = 9, byrow = TRUE)
  runMat = matrix(c(rep(2^2, 9), rep(2^3, 9), rep(2^4, 9), rep(2^5, 9), rep(2^6, 9), rep(2^7, 9)), ncol = 9, byrow = TRUE)
  par(mfrow = c(6, 9))
  par(mar = c(0, 0, 0, 0))
  par(oma = c(4, 4, 4, 4))
  colList = vector(mode = "list")
  colList[3] = "red"
  colList[4] = "yellow"
  colList[5] = "green"
  colList[6] = "green"
  colList[7] = "green"
  colList[8] = "green"
  k = 3
  N = 2^2
  m = 0
  for (i in seq(along = genList)) {
    res = unlist(resList[[i]])
    plot(0, 0, xaxs = "i", yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "", bg = "red", fg = "green")
    box()
    if (!is.null(res))
      rect(0, 0, 1, 1, col = colList[[res]])
    yPos = 0.04
    xPos = 0.6
    item = rev(genList[[i]])
    for (j in seq(along = item)) {
      text(xPos, yPos + (j - 1) * 0.125, item[j], adj = c(0, 0), cex = 0.8)
    }
    p = log2((2^k)/N)
    if (!is.null(res)) {
      romNum = as.character(as.roman(res))
      text(0.1, 0.9, do.call("expression", list(substitute(2[romNum]^(k - p), list(k = k, p = p, romNum = romNum)))), adj = c(0, 1), cex = 1.5)
    }
    k = k + 1
    if ((i%%9) == 0) {
      N = N * 2
      k = 3
    }
  }
  cat("\nChoose a fractional factorial design by clicking into the appropriate field")
  cat("\nWaiting for your selection:")
  cat("\n\n")
  flush.console()
  mtext("number of runs N", side = 2, line = 2.5, outer = TRUE)
  mtext("number of variables k", side = 3, line = 2.5, outer = TRUE)
  for (numFac in 1:9) {
    mtext(numFac + 2, at = c(-0.05 + (1/9) * numFac), outer = TRUE, line = 0.5)
  }
  for (k in 1:6) {
    mtext(2^(k + 1), at = (7 - k)/6 - (0.5 * (1/6)), side = 2, outer = TRUE, line = 0.5)
  }

  xyList = NULL
  xyList = try(locator(1), silent = TRUE)
  x = 1
  y = 1
  if (!is.null(xyList)) {
    x = ceiling(xyList$x + 8)
    y = ceiling(6 - xyList$y)
  }
  mat = matrix(1:54, ncol = 9, byrow = TRUE)
  fdo = NULL
  if (!(x %in% 1:ncol(mat)) || !(y %in% 1:nrow(mat)))
    return(fracDesign(k = 3, gen = NULL, replicates = 1))
  else index = mat[y, x]
  k = facMat[y, x]
  generator = genList[[index]]
  N = runMat[y, x]
  if (!is.null(generator)) {
    fdo = try(do.call("fracDesign", list(k = k, gen = generator)), silent = TRUE)
  }
  if (N >= 2^k & is.null(generator)) {
    replicates = N/(2^k)
    fdo = try(fracDesign(k = k, gen = NULL, replicates = replicates), silent = TRUE)
  }
  if (class(fdo)[1] == "facDesign")
    return(fdo)
  else return(genList[[mat[y, x]]])
}
# Uso fracChoose()
m1 <- fracChoose()
m1$summary()

# rsmChoose() ----
rsmChoose <- function() {
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  colFun = colorRampPalette(c("yellow", "red"), space = "rgb")
  colPalette = colFun(169)
  numRows = 6
  numCol = 9
  blockVals = c(1, 2, 3, 5, 9, 17)
  factorVals = c(2, 3, 4, 5, 5, 6, 6, 7, 7)
  rsmList = .rsmOrth
  plot.new()
  par(mfrow = c(6, 9))
  par(mar = c(0, 0, 0, 0))
  par(oma = c(4, 4, 4, 4))
  for (i in 1:6) for (j in 1:9) {
    par(mfg = c(i, j))
    plot(0, 0, xaxs = "i", yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n",
         xlab = "", ylab = "", bg = "red", fg = "green")
    box()
  }
  for (i in seq(along = rsmList)) {
    temp = rsmList[[i]]
    par(mfg = c(temp$row, temp$col))
    par(mfg = c(temp$row, temp$col))
    plot(0, 0, xaxs = "i", yaxs = "i", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, type = "n",
         xlab = "", ylab = "", bg = "red", fg = "green")
    rect(0, 0, 1, 1, col = colPalette[2^((temp$k) - (temp$p))])
    text(0.1, 0.9, paste("N =", 2^((temp$k) - (temp$p)) + temp$cc * (temp$blocks - 1) + temp$cs +
                           (temp$k + temp$p) * 2), adj = c(0, 1), cex = 1.25)
    text(0.1, 0.75, paste("k =", temp$k), adj = c(0, 1), cex = 1.25)
    text(0.1, 0.6, paste("p =", temp$p), adj = c(0, 1), cex = 1.25)
    text(0.1, 0.45, ".centerPoints", adj = c(0, 1), cex = 1.25)
    text(0.1, 0.3, paste("Cube:", temp$cc), adj = c(0, 1), cex = 1.25)
    text(0.1, 0.15, paste("Axial:", temp$cs), adj = c(0, 1), cex = 1.25)
    box()
  }
  x = 1/18 + (0:8) * 2/18
  mtext(factorVals, at = x, side = 3, line = 0.5, outer = TRUE)
  mtext("number of factors k", at = 0.5, side = 3, line = 2.5, outer = TRUE)
  x = 1/12 + (5:0) * 2/12
  mtext(blockVals, at = x, side = 2, line = 0.5, outer = TRUE, las = 2)
  mtext("number of blocks", at = 0.5, side = 2, line = 2.5, outer = TRUE)
  cat("\nChoose a response surface design by clicking into the appropriate field")
  cat("\nWaiting for your selection:")
  cat("\n\n")
  flush.console()
  x = numeric(0)
  y = numeric(0)
  xyList = locator(1)                                                         ###
  print(xyList)
  x = ceiling(xyList$x + 8)
  y = ceiling(5 - xyList$y)

  if (length(x) < 1)
    return(rsmDesign(k = 2, p = 0, blocks = 2, alpha = "both"))
  if (length(y) < 1)
    return(rsmDesign(k = 2, p = 0, blocks = 2, alpha = "both"))
  #    if (!(x %in% factorVals) || !(y %in% blockVals))                           ###
  #        return(rsmDesign(k = 2, p = 0, blocks = 2, alpha = "both"))            ###
  blocks = blockVals[y]
  k = factorVals[x]
  if (x==5 || x==7 || x==9 )                                                  ###
    p = 1                                                                      ###
  else                                                                        ###
    p = 0                                                                      ###

  for (i in seq(along = rsmList)) {
    if (rsmList[[i]]$k == k)
      if (rsmList[[i]]$blocks == blocks)
        if (rsmList[[i]]$p == p)                                        ###
          return(rsmDesign(k = k, p = rsmList[[i]]$p, blocks = blocks, ###
                           alpha = "both", cc = rsmList[[i]]$cc,                 ###
                           cs = rsmList[[i]]$cs))                                ###
  }
  return(cat("\nno selection recognized\n"))
}

# Uso rsmChoose()
rsdo <- rsmChoose()



# contourPlot3 ----
contourPlot3 = function(x, y, z, response, data = NULL, main, xlab, ylab, zlab, border, form = "linear", col = 1, col.text, cex.axis, axes = TRUE,
                        steps, factors) {
  DB = FALSE
  out = list()
  mdo = data
  x.c = deparse(substitute(x))
  y.c = deparse(substitute(y))
  z.c = deparse(substitute(z))
  r.c = deparse(substitute(response))
  if (missing(col))
    col = 1
  if (missing(col.text))
    col.text = 1
  if (missing(cex.axis))
    cex.axis = 1
  if (missing(main))
    main = paste("Response Surface for", r.c)
  if (missing(ylab))
    ylab = y.c
  if (missing(xlab))
    xlab = x.c
  if (missing(zlab))
    zlab = z.c
  if (missing(border))
    border = "white"
  if (missing(factors))
    factors = NULL
  if (missing(steps))
    steps = 100
  col.axis = par("col.axis")
  if (!is.function(col)) {
    if (identical(col, 1))
      col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if (identical(col, 2))
      col = colorRampPalette(c("blue", "white", "red"), space = "Lab")
    if (identical(col, 3))
      col = colorRampPalette(c("blue", "white", "orange"))
    if (identical(col, 4))
      col = colorRampPalette(c("gold", "white", "firebrick"))
  }
  nameVec = names(mdo$names())
  linStrings = "-1"
  for (i in seq(along = nameVec)) linStrings = paste(linStrings, "+", nameVec[i])
  if (DB)
    print(linStrings)
  combList = combn(nameVec, 2, simplify = FALSE)
  quadStrings = character(length = length(combList))
  for (i in seq(along = combList)){
    if (i == 1){
      quadStrings[i] = paste(combList[[i]][1], ":", combList[[i]][2])
    }
    else{
      quadStrings[i] = paste("+", combList[[i]][1], ":", combList[[i]][2])
    }
  }
  quadStrings = paste(quadStrings, collapse = "")
  if (DB)
    print(quadStrings)
  if (identical(form, "linear")) {
    form = paste(r.c, "~", linStrings)
    if (DB)
      print(form)
  }
  if (identical(form, "quadratic")) {
    form = paste(r.c, "~", linStrings, "+", quadStrings)

  }
  lm.1 = lm(formula = form, data = mdo$as.data.frame())
  if (DB)
    print(lm.1)
  dcList = vector(mode = "list", length = length(mdo$names()))
  names(dcList) = names(mdo$names())
  dcList[1:length(mdo$names())] = 0
  if (!is.null(factors)) {
    for (i in names(factors)) dcList[[i]] = factors[[i]][1]
  }
  if (DB)
    print(dcList)
  help.predict = function(a, b, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = 2 * b/sqrt(3)
    dcList[[y.c]] = 1 - (2 * b/sqrt(3)) - (a - b/sqrt(3))
    dcList[[z.c]] = a - b/sqrt(3)
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
  }
  a = seq(0, 1, length = steps)
  b = seq(0, sqrt(3)/2, length = steps)
  mat = outer(a, b, help.predict, x.c, y.c, lm.1)
  acc = nrow(mat)
  v = seq(acc, 1, length = acc)
  w = c(seq(2, acc, length = acc/2), seq(acc, 2, length = acc/2))
  mat[outer(w, v, `+`) <= acc] = NA
  .mfc(mat, main = main, col = col, axes = FALSE, key.axes = axis(4))
  if (axes == TRUE) {
    segments(0.5, 0, 0.5, 1, col = col.axis)
    coox1 = rep(0.49, 11)
    coox2 = rep(0.51, 11)
    cooy = seq(0, 1, length = 11)
    for (i in 2:10) {
      segments(coox1[i], cooy[i], coox2[i], cooy[i], col = col.axis)
      text(coox2[i] + 0.01, cooy[i], labels = (i - 1)/10, cex = cex.axis, col = col.text)
    }
    segments(0, 0, 0.75, 0.5, col = col.axis)
    coox1 = seq(0.745, -0.005, length = 11)
    coox2 = seq(0.755, 0.005, length = 11)
    cooy1 = seq(0.51, 0.01, length = 11)
    cooy2 = seq(0.49, -0.01, length = 11)
    for (i in 2:10) {
      segments(coox1[i], cooy1[i], coox2[i], cooy2[i], col = col.axis)
      text(coox2[i], cooy1[i] + 0.01, labels = (i - 1)/10, cex = cex.axis, col = col.text)
    }
    segments(0.25, 0.5, 1, 0, col = col.axis)
    coox1 = seq(0.245, 0.995, length = 11)
    coox2 = seq(0.255, 1.005, length = 11)
    cooy1 = seq(0.49, -0.01, length = 11)
    cooy2 = seq(0.51, 0.01, length = 11)
    for (i in 2:10) {
      segments(coox1[i], cooy1[i], coox2[i], cooy2[i], col = col.axis)
      text(coox1[i] + 0.01, cooy2[i] + 0.01, labels = (i - 1)/10, cex = cex.axis, col = col.text)
    }
  }
  segments(-0.005, 0, 0.495, 1, lwd = 5, col = border)
  segments(0.505, 1, 1.005, 0, lwd = 5, col = border)
  segments(1, -0.005, 0, -0.005, lwd = 5, col = border)
  mtext(ylab, 1, at = -0.025, cex = 1.5)
  mtext(xlab, 3, at = 0.5, cex = 1.5, line = 0.1)
  mtext(zlab, 1, at = 1.025, cex = 1.5)
  invisible(mat)
}
# Uso contourPlot3
contourPlot3(A, B, C, elongation, data = mdo, form = "quadratic")

# wirePlot3 ----
wirePlot3 = function(x, y, z, response, data = NULL, main, xlab, ylab, zlab, form = "linear", phi, theta, col = 1, steps, factors) {
  DB = FALSE
  out = list()
  mdo = data
  x.c = deparse(substitute(x))
  y.c = deparse(substitute(y))
  z.c = deparse(substitute(z))
  r.c = deparse(substitute(response))
  if (missing(col))
    col = 1
  if (missing(main))
    main = paste("Response Surface for", r.c)
  if (missing(ylab))
    ylab = y.c
  if (missing(xlab))
    xlab = x.c
  if (missing(zlab))
    zlab = z.c
  if (missing(phi))
    phi = 30
  if (missing(theta))
    theta = 30
  if (missing(factors))
    factors = NULL
  if (missing(steps))
    steps = 100
  if (!is.function(col)) {
    if (identical(col, 1))
      col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if (identical(col, 2))
      col = colorRampPalette(c("blue", "white", "red"), space = "Lab")
    if (identical(col, 3))
      col = colorRampPalette(c("blue", "white", "orange"))
    if (identical(col, 4))
      col = colorRampPalette(c("gold", "white", "firebrick"))
  }
  phi = phi%%360
  .phi = phi
  theta = theta%%360
  .theta = theta
  nameVec = names(mdo$names())
  linStrings = "-1"
  for (i in seq(along = nameVec)) linStrings = paste(linStrings, "+", nameVec[i])
  if (DB)
    print(linStrings)
  combList = combn(nameVec, 2, simplify = FALSE)
  quadStrings = character(length = length(combList))
  for (i in seq(along = combList)) if (i == 1)
    quadStrings[i] = paste(combList[[i]][1], ":", combList[[i]][2])
  else quadStrings[i] = paste("+", combList[[i]][1], ":", combList[[i]][2])
  quadStrings = paste(quadStrings, collapse = "")
  if (DB)
    print(quadStrings)
  if (identical(form, "linear")) {
    form = paste(r.c, "~", linStrings)
    if (DB)
      print(form)
  }
  if (identical(form, "quadratic")) {
    form = paste(r.c, "~", linStrings, "+", quadStrings)
    if (DB)
      print(form)
  }
  lm.1 = lm(formula = form, data = mdo$as.data.frame())
  if (DB)
    print(lm.1)
  dcList = vector(mode = "list", length = length(mdo$names()))
  names(dcList) = names(mdo$names())
  dcList[1:length(mdo$names())] = 0
  if (!is.null(factors)) {
    for (i in names(factors)) dcList[[i]] = factors[[i]][1]
  }
  if (DB)
    print(dcList)
  help.predict = function(a, b, x.c, y.c, lm.1, ...) {
    dcList[[x.c]] = 2 * b/sqrt(3)
    dcList[[y.c]] = 1 - (2 * b/sqrt(3)) - (a - b/sqrt(3))
    dcList[[z.c]] = a - b/sqrt(3)
    temp = do.call(data.frame, dcList)
    invisible(predict(lm.1, temp))
  }
  a = seq(0, 1, length = steps)
  b = seq(0, sqrt(3)/2, length = steps)
  mat = outer(a, b, help.predict, x.c, y.c, lm.1)
  acc = nrow(mat)
  sca = sin(1/3 * pi)
  ncmat = ncol(mat)
  v = seq(acc, 1, length = acc)
  w = c(seq(2, acc, length = acc/2), seq(acc, 2, length = acc/2))
  mat[outer(w, v, `+`) <= acc] = NA
  if (is.function(col)) {
    nrMat <- nrow(mat)
    ncMat <- ncol(mat)
    nbcol <- 100
    color <- col(nbcol)
    matFacet = mat[-1, -1] + mat[-1, -ncmat] + mat[-acc, -1] + mat[-acc, -ncmat]
    facetcol <- cut(matFacet, nbcol)
  }
  else {
    color = col
    facetcol = 1
  }
  maxim = max(mat, na.rm = TRUE) * acc
  minim = min(mat, na.rm = TRUE) * acc
  per = persp(x = seq(0, acc, length = acc), y = seq(0, acc * sca, length = ncmat), mat * acc, phi = .phi, theta = .theta, scale = TRUE, col = "transparent",
              border = FALSE, box = FALSE, main = main, xlab = xlab, ylab = ylab)
  lineList = contourLines(x = seq(0, acc, length = acc), y = seq(0, acc * sca, length = ncmat), mat)
  for (i in seq(along = lineList)) lines(trans3d(lineList[[i]]$x, lineList[[i]]$y, z = minim, pmat = per))
  if (.phi < 90) {
    lines(trans3d(x = seq(0, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
    lines(trans3d(x = seq(acc, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
    lines(trans3d(x = 0:acc, y = 0, z = maxim, pmat = per), lty = 2)
  }
  if (.theta > 323 || .theta < 37) {
    lines(trans3d(x = acc/2, y = acc * sca, z = minim:maxim, pmat = per), lty = 2)
    lines(trans3d(x = 0, y = 0, z = minim:maxim, pmat = per), lty = 2)
    lines(trans3d(x = acc, y = 0, z = minim:maxim, pmat = per), lty = 2)
  }
  if (.theta > 37 && .theta < 156)
    lines(trans3d(x = 0, y = 0, z = minim:maxim, pmat = per), lty = 2)
  if (.theta > 156 && .theta < 323) {
    lines(trans3d(x = acc, y = 0, z = minim:maxim, pmat = per), lty = 2)
  }
  lines(trans3d(x = seq(0, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = minim, pmat = per), lty = 1, lwd = 2)
  lines(trans3d(x = seq(acc, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = minim, pmat = per), lty = 1, lwd = 2)
  lines(trans3d(x = 0:acc, y = 0, z = minim, pmat = per), lty = 1, lwd = 2)
  text(trans3d(x = acc/2 + acc/50, y = acc * sca + acc * sca/50, z = minim, pmat = per), labels = xlab, lwd = 2)
  text(trans3d(x = -acc/50, y = -acc * sca/50, z = minim, pmat = per), labels = ylab, lwd = 2)
  text(trans3d(x = acc + acc/50, 0, z = minim, pmat = per), labels = zlab, cex = 1, lwd = 2)
  par(new = TRUE)
  persp(x = seq(0, acc, length = acc), y = seq(0, acc * sca, length = ncmat), mat * acc, phi = .phi, theta = .theta, scale = TRUE, col = color[facetcol],
        border = FALSE, box = FALSE)
  if (.phi > 0) {
    lines(trans3d(x = seq(0, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
    lines(trans3d(x = seq(acc, acc/2, length = 10), y = seq(0, acc * sca, length = 10), z = maxim, pmat = per), lty = 2)
    lines(trans3d(x = 0:acc, y = 0, z = maxim, pmat = per), lty = 2)
  }
  if (.theta > 37 && .theta < 156) {
    lines(trans3d(x = acc/2, y = acc * sca, z = minim:maxim, pmat = per), lty = 2)
    lines(trans3d(x = acc, y = 0, z = minim:maxim, pmat = per), lty = 2)
  }
  if (.theta > 156 && .theta < 323) {
    lines(trans3d(x = acc/2, y = acc * sca, z = minim:maxim, pmat = per), lty = 2)
    lines(trans3d(x = 0, y = 0, z = minim:maxim, pmat = per), lty = 2)
  }
  if (TRUE) {
    zlim = range(mat, finite = TRUE, na.rm = TRUE)
    leglevel = pretty(zlim, 6)
    legcol = col(length(leglevel))
    legpretty = as.character(abs(leglevel))
    temp = character(length(leglevel))
    temp[leglevel > 0] = "+"
    temp[leglevel < 0] = "-"
    temp[leglevel == 0] = " "
    legpretty = paste(temp, legpretty, sep = "")
    if (.theta <= 180)
      legend("topright", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
    if (.theta > 180)
      legend("topleft", inset = 0.02, legend = paste(">", legpretty), col = legcol, bg = "white", pt.cex = 1.5, cex = 0.75, pch = 15)
  }
  invisible(mat)
}
# Uso wirePlot3
wirePlot3(A, B, C, elongation, data = mdo, form = "quadratic", theta = -170)
