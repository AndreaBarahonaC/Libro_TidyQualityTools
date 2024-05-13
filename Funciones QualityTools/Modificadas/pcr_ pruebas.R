.sdSg = function(x, grouping = NULL, method = c("NOWEIGHT", "MVLUE", "RMSDF"), na.rm = TRUE, DB = TRUE) {
  DB = FALSE
  if (!is.data.frame(x) && !is.vector(x) && is.numeric(x))
    stop("x needs to be either a data.frame or a vector and numeric")
  if (is.null(grouping)) {
    if (is.data.frame(x))
      return(sd(x[, 1]))
    else return(sd(x))
  }
  else grouping = as.data.frame(grouping)
  group = unique(grouping)
  sdVec = numeric(length = length(group))
  for (i in 1:nrow(group)) {
    if (is.data.frame(x))
      temp = x[group[i, 1] == grouping[, 1], 1]
    if (is.vector(x))
      temp = x[group[i, 1] == grouping[, 1]]
    sdVec[i] = sd(temp, na.rm = T)/.c4(length(temp[!is.na(temp)]))
    if (DB) {
      print(group[i, 1])
      print(temp)
      print(length(temp[!is.na(temp)]))
    }
  }
  if (DB) {
    print(paste("std.dev: ", mean(sdVec)))
    print(sdVec)
  }
  return((mean(sdVec)))
}

.lfkp = function(wholeList, filterList) {
  if (!is.list(wholeList))
    stop(paste(deparse(substitute(wholeList)), "is not a list!"))
  if (length(wholeList) == 0)
    return(wholeList)
  if (!is.list(filterList))
    stop(paste(deparse(substitute(filterList)), "is not a list!"))
  if (length(filterList) == 0)
    return(filterList)
  logVec = lapply(names(wholeList), "%in%", names(filterList))
  filteredList = wholeList[unlist(logVec)]
  return(filteredList)
}



library(qualityTools)
library(ggplot2)

set.seed(1234)
x <- rnorm(20, mean = 20)
distribution <- "norm"
lsl = 17
usl = 23
boxcox = FALSE
lambda = c(-5, 5)
grouping = NULL
std.dev = NULL
conf.level = 0.9973002
lineWidth = 1
lineCol = "red"
lineType = "solid"
specCol = "red3"
specWidth = 1
cex.text = 2
cex.val = 1.5
cex.col = "darkgray"
plot = TRUE
bounds.lty = 3
bounds.col = "red"

#######################################33
DB = FALSE
data.name = deparse(substitute(x))[1]
par.orig <- par(c("mar", "oma", "mfrow"))
on.exit(par(par.orig))
parList = list()
parList$col = "lightblue"
parList$border = 1
parList$lwd = 1
parList$cex.axis = 1.5
target = NULL


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

# -------------
if (distribution == "weibull3" || distribution == "lognormal3" ||
    distribution == "gamma3")
  any3distr = TRUE
if (distribution != "weibull3" && distribution != "lognormal3" &&
    distribution != "gamma3")
  not3distr = TRUE
if (boxcox) {
  distribution = "norm"
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
if (!is.null(grouping))
  if (is.vector(grouping))
    grouping = as.data.frame(grouping)
center = colMeans(x)
if (!is.null(x) & !is.null(grouping)) {
  if (nrow(x) != nrow(grouping))
    stop(paste("length of ", deparse(substitute(grouping)),
               " differs from length of ", varName))
}
main = paste("Process Capability using", as.character(distribution),
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
if (DB) {
  print(paste("confHigh:", confHigh))
  print(paste("confLow:", confLow))
}
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

if (DB)
    print("TODO: Pass the estimated parameters correctly")
  fitList = vector(mode = "list", length = 0)
  fitList$x = x[, 1]
  fitList$densfun = dis
fitList$start = start
  if (DB)
    print(paste("parameter: ", paramsList))
#### Solo en este caso qFUN es eval(parse(text = paste("q", "norm", sep = "")))
qFun = eval(parse(text = paste("q", "norm", sep = "")))
paramsList = c(paramsList, .lfkp(parList, formals(qFun)))
if (distribution == "norm") {
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
### SOLO EN ESTE CASO pFUN = eval(parse(text = paste("p", "norm", sep = "")))
pFun = eval(parse(text = paste("p", "norm", sep = "")))
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
if (DB == TRUE) {
  print(cp)
  print(cpk)
  print(cpu)
  print(cpl)
  print(ppu)
  print(ppl)
  print(ppt)
}


if (missing(xlim)) {
  xlim <- range(x[, 1], usl, lsl)
  xlim <- xlim + diff(xlim) * c(-0.2, 0.2)
}
xVec <- seq(min(xlim), max(xlim), length = 200)
# SOLO EN ESTE CASO
dFun = eval(parse(text = paste("d", "norm", sep = "")))
dParamsList = .lfkp(paramsList, formals(dFun))
dParamsList$x = xVec
yVec = do.call(dFun, dParamsList)
histObj <- hist(x[, 1], plot = FALSE)
if (missing(ylim)) {
  ylim <- range(histObj$density, yVec)
  ylim <- ylim + diff(ylim) * c(0, 0.05)
}
## ----------------

# ----------------- plot
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
  geom_vline(aes(xintercept = usl, color = "Confidence interval"), linetype = "dashed", col = "red") + # USL
  geom_vline(aes(xintercept = lsl, color = "Confidence interval"), linetype = "dashed", col = "red") + # LSL
  theme(legend.position = "none")

#  etiquetas de los límites
if (!is.null(lsl) & !is.null(usl))
  p1 <- p1 +
  scale_x_continuous(limits = xlim, expand = c(0, 0),
                     sec.axis = sec_axis(~ ., breaks = c(lsl, usl),
                                         labels = c(paste("LSL =",format(lsl, digits = 3)), paste("USL =",format(usl, digits = 3)))
                     )) +
  theme(axis.text.y.right = element_text(size = 15))
else{
  if(!is.null(lsl))
    p1 <- p1 +
      scale_x_continuous(limits = xlim, expand = c(0, 0),
                         sec.axis = sec_axis(~ ., breaks = lsl,
                                             labels = paste("LSL =",format(lsl, digits = 3)) )
      )) +
  theme(axis.text.y.right = element_text(size = 15))
if(!is.null(usl))
  p1 <- p1 +
    scale_x_continuous(limits = xlim, expand = c(0, 0),
                       sec.axis = sec_axis(~ ., breaks = usl,
                                           labels = paste("USL =",format(usl, digits = 3)) )
    )) +
  theme(axis.text.y.right = element_text(size = 15))
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
  xlim(c(0.24, 0.27)) + ylim(c(0.24, 0.41))
{
  if(is.null(cpu))
    p2 <- p2 + annotate('text', x = 0.25, y = 0.40,
                        label = paste("C[pkU]==", "*"),
                        parse = TRUE, size = 4, hjust = 0)
  else p2 <- p2 + annotate('text', x = 0.25, y = 0.40,
                           label = paste("C[pkU]==", round(cpu, 2)),
                           parse = TRUE, size = 4, hjust = 0)
  if(is.null(cpl))
    p2 <- p2 + annotate('text', x = 0.25, y = 0.35,
                        label = paste("C[pkL]==", "*"),
                        parse = TRUE, size = 4, hjust = 0)
  else
    p2 <- p2 + annotate('text', x = 0.25, y = 0.35,
                        label = paste("C[pkL]==", round(cpl, 2)),
                        parse = TRUE, size = 4, hjust = 0)
  if(is.null(cpk))
    p2 <- p2 + annotate('text', x = 0.25, y = 0.30,
                        label = paste("C[pk]==", "*"),
                        parse = TRUE, size = 4, hjust = 0)
  else
    p2 <- p2 + annotate('text', x = 0.25, y = 0.30,
                        label = paste("C[pk]==", round(cpk, 2)),
                        parse = TRUE, size = 4, hjust = 0)
  if(is.null(cp))
    p2 <- p2 + annotate('text', x = 0.25, y = 0.25,
                        label = paste("C[p]==", "*"),
                        parse = TRUE, size = 4, hjust = 0)
  else
    p2 <- p2 + annotate('text',x = 0.25,y = 0.25,
                        label = paste("C[p]==", round(cp, 2)),
                        parse = TRUE,size = 4,hjust = 0)
  }


# 3. Cajita de info n, means, sd
index = 1:(length(estimates) + 3)
names(x) = data.name
if (not3distr) {
  names(x) = data.name
  adTestStats = .myADTest(x, distribution)
  A = numeric()
  p = numeric()
  if (class(adTestStats) == "adtest") {
    A = adTestStats$statistic
    p = adTestStats$p.value
  }
  else {
    A = NA
    p = NA
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
    xlim(c(0.245, 0.26)) + ylim(c(0.19, 0.41))
  {
    # n y A
    p3 <- p3 + annotate('text', x = 0.25, y = 0.40,
                        label = paste("n==", numObs),
                        parse = TRUE, size = 3, hjust = 0) +
      annotate('text', x = 0.25, y = 0.35,
               label = paste("A==", format(as.numeric(A), digits = 3)),
               parse = TRUE, size = 3, hjust = 0)

    # p
    if (!is.null(adTestStats$smaller) && adTestStats$smaller){
      p3 <- p3 + annotate(
        'text',
        x = 0.25,
        y = 0.30,
        label = paste("p<", format(as.numeric(p), digits =3)),
        parse = TRUE,
        size = 3,
        hjust = 0
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
        hjust = 0
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
        hjust = 0
      )
    }

    # mean y sd
    p3 <- p3 + annotate('text', x = 0.25, y = 0.25,
                        label = paste("mean==", format(estimates[[1]], digits = 3)),
                        parse = TRUE, size = 3, hjust = 0) +
      annotate('text', x = 0.25, y = 0.20,
               label = paste("sd==", format(estimates[[2]], digits = 3)),
               parse = TRUE, size = 3, hjust = 0)
    }
}
if (any3distr) {
  p3 <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    xlim(c(0.245, 0.26)) + ylim(c(0.19, 0.41))
  {
    # n y A
    p3 <- p3 + annotate('text', x = 0.25, y = 0.40,
                        label = paste("n==", numObs),
                        parse = TRUE, size = 3, hjust = 0) +
      annotate('text', x = 0.25, y = 0.35,
               label = paste("A = ", "*"),
               parse = FALSE, size = 3, hjust = 0) +
      annotate('text', x = 0.25, y = 0.30,
               label = paste("p = ", "*"),
               parse = FALSE, size = 3, hjust = 0) +
      annotate('text', x = 0.25, y = 0.25,
               label = paste("mean==", format(estimates[[1]], digits = 3)),
               parse = TRUE, size = 3, hjust = 0) +
      annotate('text', x = 0.25, y = 0.20,
               label = paste("sd==", format(estimates[[2]], digits = 3)),
               parse = TRUE, size = 3, hjust = 0)
    }
}

###################### Quantile Quantile PLot
# Para este necesitamos que funcione el archivo "Clase distribuciones y graficos"

qqPlot_o(x[, 1], y = distribution, ylab = "", main = "",
       axes = FALSE, bounds.lty = bounds.lty, bounds.col = bounds.col)


library(car)
car::qqPlot(x[,1], ylab = "", distribution = dis)

qqPlot(x[,1], distribution = distribution,
       ylab = "", xlab = "", )





nada <- FitDistr(x[,1],"normal")


nada$estimate

prueba <-

prueba
# Creamos un objeto DistrCollection
distr_coll <- DistrCollection$new()

# Creamos un objeto Distr para una distribución normal
distr_normal <- Distr$new(x = rnorm(100), name = "normal", parameters = c(mean = 0, sd = 1))

# Agregamos el objeto Distr a la colección
distr_coll$add(distr_normal)

# Creamos un objeto Distr para una distribución exponencial
distr_exp <- Distr$new(x = rexp(100), name = "exponential", parameters = c(rate = 1))

# Agregamos el objeto Distr a la colección
distr_coll$add(distr_exp)

# Llamamos a la función qqPlot_o con la colección de distribuciones
qqPlot_o(distr_coll, y = "normal", confbounds = TRUE, alpha = 0.05, main = "Q-Q Plot")



for (i in seq_along(distr_coll$distr)) {
  d <- distr_coll$distr[[i]]
  do.call(qqPlot, c(list(x = d$x, y = d$name), parList))
}

d <- distr_coll$distr[[1]]
do.call(qqPlot, c(list(x = d$x, y = d$name), parList))

