#1. Creación del Paquete ----
#New project, new package, crea:
#-DESCRIPTION: nombre, descripcion, autores, entre otros
#-NAMESPACE:
#las variables del paquete que se exportan (y son, por tanto, accesibles a los usuarios)
#las variables que se importan de otros paquetes.
#las clases y métodos que deben registrarse.
#-man: documentación
#-R: codigo

library(tidyqualitytools)
#2. Ejemplos distr

#3. Ejemplos facDesign
###Diseño factorial 2k####
set.seed(1234)
dfac <- facDesign(k = 3, centerCube = 4,random.seed=1234)
dfac$names(c('Factor 1', 'Factor 2', 'Factor 3'))
dfac$lows(c(80,120,1))
dfac$highs(c(120,140,2))
dfac$summary()

#Asignar rendimiento al diseño factorial
set.seed(1)
rend <- c(simProc(120,140,1),simProc(80,140,1),simProc(120,140,2),simProc(120,120,1),simProc(90,130,1.5),simProc(90,130,1.5),simProc(80,120,2),simProc(90,130,1.5),simProc(90,130,1.5),simProc(120,120,2),simProc(80,140,2),simProc(80,120,1))


dfac$.response(rend)
dfac$summary()
### effectPlot
dfac$effectPlot(classic=TRUE)

## interactionplot
interactionPlot(dfac)

### lm
m1 <- dfac$lm(rend ~ A*B*C)
summary(m1)

###paretoPlot y normalPlot
paretoPlot(dfac)
normalPlot(dfac)

#wirePlot y contourPlot
wirePlot(A,B,rend,data=dfac)
contourPlot(A,B,rend,data=dfac)

###Diseños factoriales fraccionarios 2k-p####
dfacfrac <- fracDesign(k=3,gen='C=AB',centerCube = 4)
dfacfrac$summary()
aliasTable(dfacfrac)
confounds(dfacfrac)

#fracChoose
dfac1<-fracChoose()

###Diseños replicados y puntos centrales####
dfac2 <- facDesign(k = 3, centerCube = 2, replicates = 2)
dfac2$summary()

###Respuestas múltiples####
set.seed(1234)
y2 <- rnorm(12,mean=120)
dfac$.response(data.frame(rend,y2))
dfac$summary()
#graficos
contourPlot(A, B, y2, data = dfac, form = "y2~A+B+C+A*B")
wirePlot(A, B, y2, data = dfac, factors = list(C=1), form = "y2~A*B*C")


###FITS####
dfac$set.fits(dfac$lm(rend~A+B))
dfac$set.fits(dfac$lm(y2~A*B*C))
dfac$fits
###Encontrar el mayor rendimiento esperado####
sao <- steepAscent(factors = c("A", "B"), response = "rend", data = dfac, steps = 20)
predicted <- simProc(sao$get(,5), sao$get(,6))
sao$.response(predicted)
sao
sao$plot()


###Diseños de superficie de respuesta(Modelizar relaciones no lineales y montaje secuencial)####
set.seed(1234)
fdo2 <- facDesign(k = 2, centerCube = 3)
fdo2$names(c("Factor1","Factor2"))
fdo2$lows(c(134,155))
fdo2$highs(c(155,175))
rend <- c(simProc(134,175),simProc(144.5,165.5),simProc(155,155),simProc(144.5,165.5),simProc(155,175),simProc(144.5,165.5),simProc(134,155))
fdo2$.response(rend)
fdo2$summary()

#Adicion de puntos estrella
rsdo <- starDesign(data=fdo2)
rsdo
rend2 <- c(rend,
           simProc(130, 165),
           simProc(155, 165),
           simProc(144, 155),
           simProc(144, 179),
           simProc(144, 165),
           simProc(144, 165),
           simProc(144, 165)
)

#Adicion secuencial
rsdo$.response(rend2)
rsdo

#Ajustar un modelo cuadrático
lm.3 <- rsdo$lm(rend2 ~ A*B + I(A^2) + I(B^2))
summary(lm.3)

#Graficos
wirePlot(A, B, rend2, form = "rend2 ~ A*B + I(A^2) + I(B^2)", data = rsdo)
contourPlot(A, B, rend2, form = "rend2 ~ A*B + I(A^2) + I(B^2)", data = rsdo)

#Otra forma
fdo <- rsmDesign(k=3, alpha=1.633, cc=0, cs=6)

#poner en orden estandar#
fdo <- randomize(fdo,so = TRUE)
fdo$summary()

#rsmChoose()
rsdo <- rsmChoose()



###Deseabilidades (optimizar más de una variable de respuesta)####
#desires: target: max, min, val
y1 <- c(102, 120, 117, 198, 103, 132, 132, 139, 102, 154, 96, 163, 116,
        153, 133, 133, 140, 142, 145, 142)
y2 <- c(900, 860, 800, 2294, 490, 1289, 1270, 1090, 770, 1690, 700, 1540,
        2184, 1784, 1300, 1300, 1145, 1090, 1260, 1344)
y3 <- c(470, 410, 570, 240, 640, 270, 410, 380, 590, 260, 520, 380, 520,
        290, 380, 380, 430, 430, 390, 390)
y4 <- c(67.5, 65, 77.5, 74.5, 62.5, 67, 78, 70, 76, 70, 63, 75, 65, 71,
        70, 68.5, 68, 68, 69, 70)

d1 <- desirability(y1, 120, 170, scale = c(1, 1), target = "max")
d3 <- desirability(y3, 400, 600, target = 500)
d1$print()
d1$plot()
d3$plot()
#Con experimentos diseñados
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
#Optimizar
optimum(ddo,type='optim')


