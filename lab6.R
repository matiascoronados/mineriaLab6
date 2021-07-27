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
datos <- datosB

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
salida0 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output0 <- matrix(unlist(salida0), ncol = 5, byrow = TRUE)
mejoresModelos0<-output0[order(output0[,5], decreasing = TRUE),]
head(mejoresModelos0)




#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 2. NUv2
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- seq(0.05, 0.12, 0.025)
gamma<-2^seq(-4, 2,0.5)
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
salida1 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output1 <- matrix(unlist(salida1), ncol = 5, byrow = TRUE)
mejoresModelos1<-output1[order(output1[,5], decreasing = TRUE),]
head(mejoresModelos1)



#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 2. NUv3
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- seq(0.08, 0.11, 0.01)
gamma<-2^seq(-4, 2,0.5)
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
salida2 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output2 <- matrix(unlist(salida2), ncol = 5, byrow = TRUE)
mejoresModelos2<-output2[order(output2[,5], decreasing = TRUE),]
head(mejoresModelos2)




#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 2. GAMMA
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- c(0.10)
gamma<-2^seq(-7, 1,0.25)
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
salida3 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output3 <- matrix(unlist(salida3), ncol = 5, byrow = TRUE)
mejoresModelos3<-output3[order(output3[,5], decreasing = TRUE),]
head(mejoresModelos3)

#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 2. GAMMAv2
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- c(0.10)
gamma<-2^seq(-10, -5,0.15)
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
salida4 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output4 <- matrix(unlist(salida4), ncol = 5, byrow = TRUE)
mejoresModelos4<-output4[order(output4[,5], decreasing = TRUE),]
head(mejoresModelos4)


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 2. GAMMAv3
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- c(0.10)
gamma<-2^seq(-7, -3,0.1)
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
salida5 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output5 <- matrix(unlist(salida5), ncol = 5, byrow = TRUE)
mejoresModelos5<-output5[order(output5[,5], decreasing = TRUE),]
head(mejoresModelos5)

#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 3. COST
registerDoParallel(cores = 7)
cost <- 2^seq(-1, 6, 0.5)
nu <- c(0.10)
gamma<-2^seq(-7, -5, 0.10)
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
salida6 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output6 <- matrix(unlist(salida6), ncol = 5, byrow = TRUE)
mejoresModelos6<-output6[order(output6[,5], decreasing = TRUE),]
head(mejoresModelos6)


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 3. COSTv1
registerDoParallel(cores = 7)
cost <- 2^seq(3, 6, 0.1)
nu <- c(0.10)
gamma<-2^seq(-7, -5, 0.10)
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
salida7 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output7 <- matrix(unlist(salida7), ncol = 5, byrow = TRUE)
mejoresModelos7<-output7[order(output7[,5], decreasing = TRUE),]
head(mejoresModelos7)

#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 3. COSTv2
registerDoParallel(cores = 7)
cost <- 2^seq(3, 4, 0.05)
nu <- c(0.10)
gamma<-2^seq(-7, -5, 0.10)
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
salida8 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output8 <- matrix(unlist(salida8), ncol = 5, byrow = TRUE)
mejoresModelos8<-output8[order(output8[,5], decreasing = TRUE),]
head(mejoresModelos8)








#######################################################
#######################################################
#######################################################
#######################################################
# SUJETO B
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 1. NU
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- seq(0.4, 0.9, 0.1)
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
salida9 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output9 <- matrix(unlist(salida9), ncol = 5, byrow = TRUE)
mejoresModelos9<-output9[order(output9[,5], decreasing = TRUE),]
head(mejoresModelos9)

#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 1. NU v1
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- seq(0.85, 0.95, 0.05)
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
salida10 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output10 <- matrix(unlist(salida10), ncol = 5, byrow = TRUE)
mejoresModelos10<-output10[order(output10[,5], decreasing = TRUE),]
head(mejoresModelos10)


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 1. GAMMA
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- c(0.95)
gamma<-2^seq(-5, 3, 1)
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
salida11 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output11 <- matrix(unlist(salida11), ncol = 5, byrow = TRUE)
mejoresModelos11<-output11[order(output11[,5], decreasing = TRUE),]
head(mejoresModelos11)


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 1. GAMMAv1
registerDoParallel(cores = 7)
cost <- 2^seq(-4, 12, 2)
nu <- c(0.95)
gamma<-2^seq(-3, 1, 0.5)
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
salida12 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output12 <- matrix(unlist(salida12), ncol = 5, byrow = TRUE)
mejoresModelos12 <-output12[order(output12[,5], decreasing = TRUE),]
head(mejoresModelos12)

#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 2. COST
registerDoParallel(cores = 7)
cost <- 2^seq(2, 6, 0.5)
nu <- c(0.95)
gamma<-2^seq(-3, 1, 0.5)
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
salida13 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output13 <- matrix(unlist(salida13), ncol = 5, byrow = TRUE)
mejoresModelos13 <-output13[order(output13[,5], decreasing = TRUE),]
head(mejoresModelos13)


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
# TEST 2. COST v1
registerDoParallel(cores = 7)
cost <- 2^seq(11, 14, 0.5)
nu <- c(0.95)
gamma<-2^seq(-3, 1, 0.5)
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
salida14 <- (c( foreach(i = 1:nrow(parms),  combine = rbind, .inorder = FALSE) %dopar% {
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

output14 <- matrix(unlist(salida14), ncol = 5, byrow = TRUE)
mejoresModelos14 <-output14[order(output14[,5], decreasing = TRUE),]
head(mejoresModelos14)



####### ####### #######  ####### ####### ####### 
# Para ver el tema de la capacidad de autoregulacion!.
####### ####### #######  ####### ####### ####### 
 

inverseStep=matrix(1,180/Ts,1)
inverseStep[(90/Ts):(180/Ts),1]=0
train <- train
mejoresModelos <- mejoresModelos14

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





