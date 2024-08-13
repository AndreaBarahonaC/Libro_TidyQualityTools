#######################################################################################
####################### DISEÑO PLACKETT BURMAN - FUNCIONES ############################
#######################################################################################

# Función pbDesign ----
pbDesign <- function(n, k , randomize = TRUE, replicates = 1) {
  if(missing(n)&&missing(k))
    stop("Either n or k must be set!")
  if(missing(n)==FALSE && missing(k)==FALSE && k!=n-1 )
    stop("Wrong combination of n and k")
  if(missing(n))
    n=k+1
  if(missing(k))
    k=n-1
  DB = FALSE
  odo = NA
  if (DB)
    print(n)
  design = .pbDesign(n)
  repVec = rep(1, nrow(design))
  if (replicates > 1) {
    X = .pbDesign(n)
    for (i in 1:(replicates - 1)) {
      design = rbind(design, X)
      repVec = c(repVec, rep(i + 1, times = nrow(X)))
    }
  }
  Replicate = data.frame(Replicate = as.numeric(repVec))
  if (DB)
    print(Replicate)
  odo = pbDesign.c$new()
  odo$design = design
  names(odo$design) = .NAMES[1:ncol(design)]
  odo$replic = Replicate
  StandOrder = 1:nrow(odo$design)
  RunOrder = StandOrder
  if (randomize) {
    RunOrder = sample(1:nrow(odo$design), nrow(odo$design), replace = FALSE, prob = NULL)
  }
  odo$design = odo$design[order(RunOrder), ]
  odo$replic = data.frame(Replicate = odo$replic[order(RunOrder), 1])
  row.names(odo$design) = odo$design$RunOrder
  odo$runOrder = data.frame(RunOrder = data.frame(RunOrder = RunOrder)[order(RunOrder), ])
  odo$standardOrder = data.frame(StandOrder = data.frame(StandOrder = StandOrder)[order(RunOrder), ])
  odo$response = data.frame(y = rep(NA, nrow(odo$design)))
  tfList = vector("list", ncol(design))
  for (i in seq(along = tfList)) tfList[[i]] = pbFactor$new()
  names(tfList) = names(odo$design)
  odo$.factors(tfList)
  valList = list(length = length(odo$names()))
  for (i in names(odo$names())) valList[[i]] = sort(unique(odo$design[, i]))
  odo$values(valList)
  return(odo)
}
#Uso pbDesign####
pbdo<-pbDesign(n=5)
pbdo$summary()
