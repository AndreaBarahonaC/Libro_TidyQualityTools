library(R6)
library(ggplot2)
library(patchwork)
library(scales)
library(plotly)
library(gridExtra)
library(RColorBrewer)

############DISEÑOS FACTORIALES 2^k#############




###############USO DE facDesign#######################################################
dfac <- facDesign(k = 3, centerCube = 4)
#dfac$names()
dfac$names(c('Factor 1', 'Factor 2', 'Factor 3'))
#dfac$names()
dfac$lows(c(80,120,1))
#dfac$lows()
dfac$highs(c(120,140,2))
#dfac$highs()
dfac$summary()

###############################################################################

### lm #########################################################################
m1 <- dfac$lm(rend ~ A*B*C)
summary(m1)


#########DISEÑOS FACTORIALES FRACCIONARIOS 2^k-p#################

####Uso diseños factoriales fraccionarios####
dfacfrac <- fracDesign(k=3,gen='C=AB',centerCube = 4)
dfacfrac$summary()
aliasTable(dfacfrac)
confounds(dfacfrac)

#######DISEÑOS REPLICADOS Y PUNTOS CENTRALES###############
dfac2 <- facDesign(k = 3, centerCube = 2, replicates = 1)
dfac2$summary()

#######RESPUESTAS MÚLTIPLES###############################
set.seed(1234)
y2 <- rnorm(12,mean=120)
dfac$.response(data.frame(rend,y2))
dfac$.response()
#graficos
wirePlot(A, B, rend, data = dfac, form = "rend~A+B+C+A*B")
contourPlot(A, B, y2, data = dfac, form = "y2~A+B+C+A*B")

#AJUSTAR EL FACTOR
wirePlot(A, B, y2, data = dfac, factors = list(C=-1), form = "y2~A*B*C")
wirePlot(A, B, y2, data = dfac, factors = list(C=1), form = "y2~A*B*C")

#FITS
dfac$set.fits(dfac$lm(rend~A+B))
dfac$set.fits(dfac$lm(y2~A*B*C))
dfac$fits

########DISEÑOS DE SUPERFICIE DE RESPUESTA#################
set.seed(1234)
fdo2 <- facDesign(k = 2, centerCube = 3)
fdo2$names(c("Factor1","Factor2"))
fdo2$lows(c(134,155))
fdo2$highs(c(155,175))
rend <- c(simProc(134,175),simProc(144.5,165.5),simProc(155,155),simProc(144.5,165.5),simProc(155,175),simProc(144.5,165.5),simProc(134,155))
fdo2$.response(rend)


####Visualizacion####
wirePlot(A, B, rend2, form = "rend2 ~ A*B + I(A^2) + I(B^2)", data = rsdo)
contourPlot(A, B, rend2, form = "rend2 ~ A*B + I(A^2) + I(B^2)", data = rsdo)
####filled.contour####
A = seq(40, 210, length = 100)
B = seq(90, 190, length = 100)
C = seq(90, 190, length = 100)
# filled.contour(A, B,outer(A,B, simProc, noise = FALSE), xlab = "Factor 1", ylab = "Factor 2", color.palette = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))
contourPlot(x = A, y = B, z = outer(A,B, simProc, noise = FALSE), xlab = "Factor 1", ylab = "Factor 2")


####poner en orden estandar####
fdo <- randomize(fdo,so = TRUE)
fdo$summary()



###MONTAJE SECUENCIAL#######################
fdo3 <- facDesign(k = 6)
fdo3$summary()
rsdo <- starDesign(alpha = "orthogonal", data = fdo3)
rsdo$summary()
###ALEATORIZACIÓN##########################
randomize(fdo, random.seed = 123)
fdo$summary()
randomize(fdo, so = TRUE)
fdo$summary()
###BLOQUEO###################################
blocking(fdo,blocks = 1)
fdo$summary()
###DESEABILIDADES############################


##EJEMPLO:####


# d<-desirability$new()
# desirability


####Utilizacion de deseabilidades junto con experimentos desiñados#####
ddo <- rsmDesign(k = 3, alpha = 1.633, cc = 0, cs = 6)
ddo <- randomize(ddo,so=TRUE)
ddo$summary()
ddo$names(c("silica", "silan", "sulfur"))
ddo$highs(c(1.7, 60, 2.8))
ddo$lows(c(0.7, 40, 1.8))
ddo$summary()
#asignar response
ddo$.response(data.frame(y1, y2, y3, y4)[c(5, 2, 3, 8, 1, 6, 7, 4, 9:20), ])
d2 <- desirability(y2, 1000, 1300, target = "max")
d4 <- desirability(y4, 60, 75, target = 67.5)
#asignar desires
ddo$desires(d1)
ddo$desires(d2)
ddo$desires(d3)
ddo$desires(d4)
ddo$desires()
#asignar fits
ddo$set.fits(ddo$lm(y1 ~ A + B + C + A:B + A:C + B:C + I(A^2) + I(B^2) + I(C^2)))
ddo$set.fits(ddo$lm(y2 ~ A + B + C + A:B + A:C + B:C + I(A^2) + I(B^2) + I(C^2)))
ddo$set.fits(ddo$lm(y3 ~ A + B + C + A:B + A:C + B:C + I(A^2) + I(B^2) + I(C^2)))
ddo$set.fits(ddo$lm(y4 ~ A + B + C + A:B + A:C + B:C + I(A^2) + I(B^2) + I(C^2)))
ddo$fits



################################################################################################33



###DISEÑOS DE MEZCLAS####









