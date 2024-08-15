library(tidyqualitytools)
#as.data.frame.facDesign
dfrac <- fracDesign(k = 3, gen = "C = AB")
as.data.frame.facDesign(dfrac)

#aliastable
aliasTable(dfrac)
tdo <- taguchiDesign("L4_2")
aliasTable(tdo)
