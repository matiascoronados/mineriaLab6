require(doParallel)
require(e1071)

datosA=read.csv("/Users/matiascoronado/Downloads/Taller_6_SVR/G2_002.csv")
datosB=read.csv("/Users/matiascoronado/Downloads/Taller_6_SVR/G4_002.csv")

# Funcion del profe
retardos_multi <- function(
  signalData,
  lags
)
{
  signal.uni <- signalData
  max.lag <- max(unlist(lags)) + 1
  indices <- 1:nrow(signal.uni)
  lag.mat <- embed(indices, max.lag)
  col.names <- list("PAMn","VFSCn")
  columns <- NULL
  lagged.columns.names <- c()
  for(colname in col.names){
    lag.order <- lags[[colname]]
    columns[[colname]] <- signal.uni[lag.mat[, 1], colname]
    if(!is.null(lag.order) && lag.order > 0)
      for(i in 1:lag.order){
        new.colname <- paste(colname, paste0("lag", i), sep = ".")
        lagged.columns.names <- c(lagged.columns.names, new.colname)
        columns[[new.colname]] <- signal.uni[lag.mat[, i+1], colname]
      }
  }
  folded.signal <- data.frame(columns)
  sorting <- order(lag.mat[, 1])
  folded.signal <- folded.signal[sorting, ]
  list(folded.signal = folded.signal, lagged.columns.names = lagged.columns.names)
}


#### Procesamiento DATOS A ####
## 70% of the sample size
datos <- datosA

attach(datos)
#Tiempo de muestreo
Ts=0.2
#Se genera la secuencia
Tiempo=seq(Ts,(length(VFSC))*Ts,Ts)
formula=VFSC ~ PAM
set.seed(123)

# Costos mas bajos 0.0625
# NU: 0 - 0.2
# Variar de 0.05
# Gamma: -3 - 2
# Variar de 1

#### ORIGINAL
registerDoParallel(cores = 6)
cost <- 2^seq(-4, 12, 2)
nu <- seq(0.1, 0.9, 0.4)
gamma<-2^seq(-4, 12, 2)
lagsList<-seq(1,5,1)

#######################################################
#######################################################
#######################################################
#######################################################
# SUJETO A
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 1. NU
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- seq(0.1, 0.2, 0.05)
gamma<-2^seq(-4, 12, 2)
lagsList<-seq(1,5,1)

# Se normalizan los datos en valores 0 - 1
PAMn<-(datos$PAM-min(datos$PAM))/(max(datos$PAM)-min(datos$PAM))
VFSCn<-(datos$VFSC-min(datos$VFSC))/(max(datos$VFSC)-min(datos$VFSC))
data <- data.frame(PAMn,VFSCn)

smp_size <- floor(0.50 * nrow(data))
## set the seed to make your partition reproducible
set.seed(123)
test <- data[1:smp_size,]
aux <- smp_size+1
train <- data[aux:nrow(data),]

#Se obtienen los mejores modelos
parms <- expand.grid(lagsList=lagsList, cost = cost, nu = nu, gamma=gamma)
salida <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
  c <- parms[i, ]$cost
  n <- parms[i, ]$nu
  g <- parms[i, ]$gamma
  l <- parms[i, ]$lagsList
  lag<-list(PAMn = l,VFSCn = 0)
  signal.train <- retardos_multi(train, lag)
  retDatos=signal.train$folded.signal
  x=subset(retDatos, select = -VFSCn)
  y=retDatos$VFSCn
  set.seed(123)
  modelo <- e1071::svm(x, y, type = "nu-regression", kernel = "radial", cost = c, nu = n, gamma=g)
  dataframe <- data.frame(PAMn = test$PAMn)
  colnames <- c("PAMn")
  aux <- ncol(x)-1
  for(i in 1:aux){
    colname <- paste('PAMn.lag',i,sep="")
    colnames <- append(colnames, colname)
    newcol <- data.frame(i = test$PAMn)
    dataframe <- cbind(dataframe, newcol)
  }
  colnames(dataframe) <-colnames
  pred <- predict(modelo, dataframe)
  corr_pred<-cor(pred,test$VFSCn,method = "pearson")
  dataframe <- NULL
  c(l, c, n, g, corr_pred)
}))

output0 <- matrix(unlist(salida), ncol = 5, byrow = TRUE)





























####### ####### #######  ####### ####### ####### 
# Para ver el tema de la capacidad de autoregulacion!.
####### ####### #######  ####### ####### ####### 

inverseStep=matrix(1,180/Ts,1)
inverseStep[(90/Ts):(180/Ts),1]=0
train <- train


# partes 2 son los mejores!!!!!!!!!!!
mejoresModelos.sujetoB.02 <- mejoresModelosDatosB.parte2
head(mejoresModelos.sujetoB.02)
#Mostrar cuales dan las mejores corr

mejoresModelos <- mejoresModelos.sujetoB.02

for (i in 1:6){
  PAMn<-(datos$PAM-min(datos$PAM))/(max(datos$PAM)-min(datos$PAM))
  VFSCn<-(datos$VFSC-min(datos$VFSC))/(max(datos$VFSC)-min(datos$VFSC))
  data <- data.frame(PAMn,VFSCn)
  lag<-list(PAMn = mejoresModelos[i,1],VFSCn = 0)
  signal.train <- retardos_multi(train , lag)
  retDatos=signal.train$folded.signal
  x=subset(retDatos, select = -VFSCn)
  y=retDatos$VFSCn
  mejorModelo <- svm(x, y, kernel = "radial",type = "nu-regression", cost = mejoresModelos[i,2], nu = mejoresModelos[i,3], gamma=mejoresModelos[i,4])
  PAMn=inverseStep
  VFSCn=inverseStep 
  data <- data.frame(PAMn,VFSCn)
  lag<-list(PAMn = mejoresModelos[i,1],VFSCn = 0)
  signal.train <- retardos_multi(data, lag)
  retDatos=signal.train$folded.signal
  x=subset(retDatos, select = -VFSCn)
  y=retDatos$VFSCn
  stepTime=seq(Ts,(length(retDatos$PAMn))*Ts,Ts)
  stepResponse <- predict(mejorModelo, x ) 
  plot(stepTime,retDatos$PAMn,type="l", col="red")
  lines(stepTime,stepResponse, col = "blue")
  legend("topright", c("Escalon de presión", "respuesta al escalon"), title = "autorregulacion", pch = 1, col=c("red","blue"),lty=c(1,1),inset = 0.01)
  print(paste("corr=",mejoresModelos[i,5]))
  readline(prompt="Press [enter] to continue")
}


#mejoresModelosA <- mejoresModelos
# SANO
#mejoresModelosB <- mejoresModelos
# ENFERMO





