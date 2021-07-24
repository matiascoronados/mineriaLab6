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

registerDoParallel(cores = 6)
cost <- 2^seq(-4, 12, 2)
nu <- seq(0.1, 0.9, 0.4)
gamma<-2^seq(-4, 12, 2)
lagsList<-seq(1,5,1)
# Se normalizan los datos en valores 0 - 1
PAMn<-(datos$PAM-min(datos$PAM))/(max(datos$PAM)-min(datos$PAM))
VFSCn<-(datos$VFSC-min(datos$VFSC))/(max(datos$VFSC)-min(datos$VFSC))
data <- data.frame(PAMn,VFSCn)

smp_size <- floor(0.50 * nrow(data))
## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)
train <- data[train_ind, ]
test <- data[-train_ind, ]



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
  modelo <- e1071::svm(x, y, type = "nu-regression", kernel = "radial", cost = c, nu = n, gamma=g)
  pred <- predict(modelo, test) 
  corr_pred<-cor(pred,test,method = "pearson")
   
  c(l, c, n, g, corr_pred)
}))

output <- matrix(unlist(salida), ncol = 5, byrow = TRUE)
# Se ordenan por correlacion.
mejoresModelos<-output[order(output[,5], decreasing = TRUE),]
print(mejoresModelos)









