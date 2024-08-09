######################################################################
###################### DISTRIBUTION - FUNCIONES ######################
######################################################################

# distribution ----
distribution <- function(x = NULL, distrib = "weibull", start, ...) {
  distr_coll <- DistrCollection$new()
  if (is.character(distrib))
    distrib = tolower(distrib)
  allDistr = c("beta", "cauchy", "chi-squared", "exponential", "f", "gamma", "geometric", "log-normal", "logistic", "negative binomial", "normal", "poisson",
               "t", "weibull")
  if (distrib %in% allDistr){
    distrVec = distrib
  }
  else{distrVec = c("normal")}
  if (identical(distrib, "all"))
    distrVec = allDistr
  if (identical(distrib, "quality"))
    distrVec = c("normal", "log-normal", "exponential", "weibull")
  for (i in seq(along = distrVec)) {
    temp <- suppressWarnings(FitDistr(x, densfun = distrVec[i]))
    fit <- Distr$new(x = x,
                     name = distrVec[i],
                     parameters = temp$estimate,
                     sd = temp$sd,
                     loglik = temp$loglik,
                     n = length(x))
    distr_coll$add(fit)
  }
  return(distr_coll)
}

# FitDistr ----
FitDistr <- function (x, densfun, start, ...){
  myfn <- function(parm, ...) -sum(log(dens(parm, ...)))
  mylogfn <- function(parm, ...) -sum(dens(parm, ..., log = TRUE))
  mydt <- function(x, m, s, df, log) dt((x - m)/s, df, log = TRUE) -
    log(s)
  Call <- match.call(expand.dots = TRUE)
  if (missing(start))
    start <- NULL
  dots <- names(list(...))
  dots <- dots[!is.element(dots, c("upper", "lower"))]
  if (missing(x) || length(x) == 0L || mode(x) != "numeric")
    stop("'x' must be a non-empty numeric vector")
  if (any(!is.finite(x)))
    stop("'x' contains missing or infinite values")
  if (missing(densfun) || !(is.function(densfun) || is.character(densfun)))
    stop("'densfun' must be supplied as a function or name")
  control <- list()
  n <- length(x)
  if (is.character(densfun)) {
    distname <- tolower(densfun)
    densfun <- switch(distname, beta = dbeta, cauchy = dcauchy,
                      `chi-squared` = dchisq, exponential = dexp, f = df,
                      gamma = dgamma, geometric = dgeom, `log-normal` = dlnorm,
                      lognormal = dlnorm, logistic = dlogis, `negative binomial` = dnbinom,
                      normal = dnorm, poisson = dpois, t = mydt, weibull = dweibull,
                      NULL)
    if (is.null(densfun))
      stop("unsupported distribution")
    if (distname %in% c("lognormal", "log-normal")) {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "log-Normal"), domain = NA)
      if (any(x <= 0))
        stop("need positive values to fit a log-Normal")
      lx <- log(x)
      sd0 <- sqrt((n - 1)/n) * sd(lx)
      mx <- mean(lx)
      estimate <- c(mx, sd0)
      sds <- c(sd0/sqrt(n), sd0/sqrt(2 * n))
      names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      vc <- matrix(c(sds[1]^2, 0, 0, sds[2]^2), ncol = 2,
                   dimnames = list(names(sds), names(sds)))
      names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dlnorm(x, mx, sd0, log = TRUE))), class = "FitDistr"))
    }
    if (distname == "normal") {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "Normal"), domain = NA)
      sd0 <- sqrt((n - 1)/n) * sd(x)
      mx <- mean(x)
      estimate <- c(mx, sd0)
      sds <- c(sd0/sqrt(n), sd0/sqrt(2 * n))
      names(estimate) <- names(sds) <- c("mean", "sd")
      vc <- matrix(c(sds[1]^2, 0, 0, sds[2]^2), ncol = 2,
                   dimnames = list(names(sds), names(sds)))
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dnorm(x, mx,
                                                                 sd0, log = TRUE))), class = "FitDistr"))
    }
    if (distname == "poisson") {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "Poisson"), domain = NA)
      estimate <- mean(x)
      sds <- sqrt(estimate/n)
      names(estimate) <- names(sds) <- "lambda"
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("lambda",
                                                              "lambda"))
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dpois(x, estimate,
                                                                 log = TRUE))), class = "FitDistr"))
    }
    if (distname == "exponential") {
      if (any(x < 0))
        stop("Exponential values must be >= 0")
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "exponential"), domain = NA)
      estimate <- 1/mean(x)
      sds <- estimate/sqrt(n)
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("rate",
                                                              "rate"))
      names(estimate) <- names(sds) <- "rate"
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dexp(x, estimate,
                                                                log = TRUE))), class = "FitDistr"))
    }
    if (distname == "geometric") {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "geometric"), domain = NA)
      estimate <- 1/(1 + mean(x))
      sds <- estimate * sqrt((1 - estimate)/n)
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("prob",
                                                              "prob"))
      names(estimate) <- names(sds) <- "prob"
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(dgeom(x, estimate,
                                                                 log = TRUE))), class = "FitDistr"))
    }
    if (distname == "weibull" && is.null(start)) {
      if (any(x <= 0))
        stop("Weibull values must be > 0")
      lx <- log(x)
      m <- mean(lx)
      v <- var(lx)
      shape <- 1.2/sqrt(v)
      scale <- exp(m + 0.572/shape)
      start <- list(shape = shape, scale = scale)
      start <- start[!is.element(names(start), dots)]
    }
    if (distname == "gamma" && is.null(start)) {
      if (any(x < 0))
        stop("gamma values must be >= 0")
      m <- mean(x)
      v <- var(x)
      start <- list(shape = m^2/v, rate = m/v)
      start <- start[!is.element(names(start), dots)]
      control <- list(parscale = c(1, start$rate))
    }
    if (distname == "negative binomial" && is.null(start)) {
      m <- mean(x)
      v <- var(x)
      size <- if (v > m)
        m^2/(v - m)
      else 100
      start <- list(size = size, mu = m)
      start <- start[!is.element(names(start), dots)]
    }
    if (is.element(distname, c("cauchy", "logistic")) &&
        is.null(start)) {
      start <- list(location = median(x), scale = IQR(x)/2)
      start <- start[!is.element(names(start), dots)]
    }
    if (distname == "t" && is.null(start)) {
      start <- list(m = median(x), s = IQR(x)/2, df = 10)
      start <- start[!is.element(names(start), dots)]
    }
  }
  if (is.null(start) || !is.list(start))
    stop("")
  nm <- names(start)
  f <- formals(densfun)
  args <- names(f)
  m <- match(nm, args)
  if (any(is.na(m)))
    stop("'start' specifies names which are not arguments to 'densfun'")
  formals(densfun) <- c(f[c(1, m)], f[-c(1, m)])
  dens <- function(parm, x, ...) densfun(x, parm, ...)
  if ((l <- length(nm)) > 1L)
    body(dens) <- parse(text = paste("densfun(x,", paste("parm[",
                                                         1L:l, "]", collapse = ", "), ", ...)"))
  Call[[1L]] <- quote(stats::optim)
  Call$densfun <- Call$start <- NULL
  Call$x <- x
  Call$par <- start
  Call$fn <- if ("log" %in% args)
    mylogfn
  else myfn
  Call$hessian <- TRUE
  if (length(control))
    Call$control <- control
  if (is.null(Call$method)) {
    if (any(c("lower", "upper") %in% names(Call)))
      Call$method <- "L-BFGS-B"
    else if (length(start) > 1L)
      Call$method <- "BFGS"
    else Call$method <- "Nelder-Mead"
  }
  res <- eval.parent(Call)
  if (res$convergence > 0L)
    stop("optimization failed")
  vc <- solve(res$hessian)
  sds <- sqrt(diag(vc))
  structure(list(estimate = res$par, sd = sds, vcov = vc,
                 loglik = -res$value, n = n), class = "FitDistr")
}

# cg_RunChart ----
cg_RunChart <- function (x, target, tolerance, ref.interval, facCg, facCgk,
                         n = 0.2, col = "black", pch = 19,
                         xlim = NULL, ylim = NULL, main = "Run Chart",
                         conf.level = 0.95, cgOut = TRUE){
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

# cg_HistChart ----
cg_HistChart <- function (x, target, tolerance, ref.interval, facCg, facCgk,
                          n = 0.2, col, xlim, ylim, main, conf.level = 0.95, cgOut = TRUE){
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

# cg_ToleranceChart ----
cg_ToleranceChart <- function (x, target, tolerance, ref.interval, facCg, facCgk,
                               n = 0.2, col, pch, xlim, ylim, main, conf.level = 0.95,
                               cgOut = TRUE){
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

# cg ----
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

# print.adtest ----
print.adtest <- function(x, digits = 4, quote = TRUE, prefix = "", ...) {
  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data: ", x$data.name, "\n")
  out <- character()
  if (!is.null(x$statistic))
    out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic[[1]], 4))))
  if (!is.null(x$parameter))
    out <- c(out, paste(names(x$parameter), "=", format(round(x$parameter, 3))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value[[1]], digits = digits)
    if (x$tableValue) {
      if (x$smaller)
        out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("<=", fp)))
      else out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "=") fp else paste(">", fp)))
    }
    else {
      out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
    }
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("alternative hypothesis: ")
  if (!is.null(x$null.value)) {
    if (length(x$null.value) == 1) {
      cat("true", names(x$null.value), "is not equal to", x$null.value, "\n")
    }
    else {
      cat(x$alternative, "\nnull values:\n")
      print(x$null.value, ...)
    }
  }
  if (!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")), "percent confidence interval:\n", format(c(x$conf.int[1L], x$conf.int[2L])), "\n")
  }
  if (!is.null(x$estimate)) {
    cat("sample estimates:\n")
    print(x$estimate, ...)
  }
  cat("\n")
  invisible(x)
}

# pcr ----
pcr <- function (x, distribution = "normal", lsl, usl, target, boxcox = FALSE,
                 lambda = c(-5, 5), main, xlim, ylim, grouping = NULL, std.dev = NULL,
                 conf.level = 0.9973002, lineWidth = 1, lineCol = "red",
                 lineType = "solid", specCol = "red3", specWidth = 1, cex.text = 2,
                 cex.val = 1.5, cex.col = "darkgray", plot = TRUE, ADtest = TRUE, bounds.lty = 3,
                 bounds.col = "red", ...){
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

  # PLOT ------------------
  {
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
  }

  if(ADtest){
    if(not3distr){
      print.adtest(adTestStats)
    }
  }

  if(plot==TRUE){
    show(main_plot)
    invisible(list(lambda = lambda, cp = cp, cpk = cpk,
                   cpl = cpl, cpu = cpu, ppt = ppt, ppl = ppl, ppu = ppu,
                   A = A, usl = usl, lsl = lsl, target = target,
                   adTest = adTestStats, plot = main_plot))
  }
  else{
    invisible(list(lambda = lambda, cp = cp, cpk = cpk,
                   cpl = cpl, cpu = cpu, ppt = ppt, ppl = ppl, ppu = ppu,
                   A = A, usl = usl, lsl = lsl, target = target,
                   adTest = adTestStats, plot = main_plot))
  }

}
set.seed(1234)
datos <- rnorm(20, mean = 20)
pcr(datos, "normal", lsl = 17, usl = 23)

set.seed(1234)
weib <- rweibull(20, shape = 2, scale = 8)
pcr(weib, "weibull", usl = 20)
