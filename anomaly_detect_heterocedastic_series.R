library(hht)
library(daltoolbox)
library(harbinger)
library(e1071)
library(changepoint)
library(ggplot2)
library(dplyr)

#parâmetros do CEEMD
noise.amp=0.001
trials=5

obj <- harbinger()
obj$noise.amp <- noise.amp
obj$trials <- trials

#htt data

#serie volcano
serie = sig

#changepoint data
serie = head(wave.c44137,5000)
serie = head(HC1,5000)                       
serie = head(ftse100$V2,5000)  
     
id <- 1:length(serie)

#gerando IMFs
set.seed(123)#semente para rodar sempre o mesmo resultado (senão cada vez que roda vai ter um conjunto diferente de IMF)
ceemd.result <- try(CEEMD(serie, id, obj$noise.amp, obj$trials, verbose = TRUE))#roda CEEMD para gerar as IMFs

PlotIMFs(ceemd.result)#plota as IMFs com o sinal original e os ruídos

##Função que soma as IMFs dada IMF inicial e final
fc_somaIMF <- function(ceemd.result, inicio, fim){
  soma_imf <- rep(0, length(ceemd.result[["original.signal"]]))
  for (k in inicio:fim){
    soma_imf <- soma_imf + ceemd.result[["imf"]][,k]
  }
  return(soma_imf)
}

#armazena as IMFs
vec <- list()
for (n in 1:ceemd.result$nimf){
  vec[[n]] <- ceemd.result[["imf"]][,n]
}

#cria as IMFs acumuladas
cum.vec <- list()
for (i in 1:ceemd.result[["nimf"]]){
  cum.vec[[i]] <- fc_somaIMF(ceemd.result, 1, i)
}

##Função que divide IMFs (acha div) 
fc_div <- function(lista, mode="var"){
  if(mode=="var"){
    ##calculando var pra cada imf
    vec <- vector()
    for (n in 1:length(lista)){
      vec[n] <- var(lista[[n]])
    }
    #Curvatura mínima
    res <- daltoolbox::transform(daltoolbox::fit_curvature_min(), vec)
    div <- res$x
    #      print(vec)        
  }
  if(mode=="rug"){
    ##Função de rugosidade
    
    fc <- function(x){
      firstD = diff(x)
      normFirstD = (firstD - mean(firstD)) / sd(firstD)
      roughness = (diff(normFirstD) ** 2) / 4
      return(mean(roughness))
    }
    
    ##calculando rugosidade pra cada imf
    vec <- vector()
    for (n in 1:length(lista)){
      vec[n] <- fc(lista[[n]])
    }
    #Curvatura máxima
    res <- daltoolbox::transform(daltoolbox::fit_curvature_max(), vec)
    div <- res$x
    #      print(vec)
  }
  
  plot(vec, t="p", ylab="", xlab="", col=ifelse(vec==vec[div], 'red', 'blue'))
  return(div)
}

div <- fc_div(vec, "var")#normal com a variância (curvatura mínima)
div <- fc_div(cum.vec, "rug")#acumulado com a rugosidade (curvatura máxima)

#--- modificação curtose ---


##Função que soma as IMFs dada IMF inicial e final
fc_somaIMF <- function(ceemd.result, inicio, fim){
  soma_imf <- rep(0, length(ceemd.result[["original.signal"]]))
  for (k in inicio:fim){
    soma_imf <- soma_imf + ceemd.result[["imf"]][,k]
  }
  return(soma_imf)
}

#armazena as IMFs
vec <- list()
for (n in 1:ceemd.result$nimf){
  vec[[n]] <- ceemd.result[["imf"]][,n]
}

#cria as IMFs acumuladas
cum.vec <- list()
for (i in 1:ceemd.result[["nimf"]]){
  cum.vec[[i]] <- fc_somaIMF(ceemd.result, 1, i)
}

##Função que divide IMFs (acha div) 
fc_div2 <- function(lista, mode="var"){
  if(mode=="var"){
    ##calculando var pra cada imf
    vec <- vector()
    for (n in 1:length(lista)){
      vec[n] <- var(lista[[n]])
    }
    #Curvatura mínima
    res <- daltoolbox::transform(daltoolbox::fit_curvature_min(), vec)
    div <- res$x
    #      print(vec)        
  }
  if(mode=="curtose"){
    ##Função de rugosidade
    
    fc_curtose <- function(x){
      firstD = diff(x)
      normFirstD = (firstD - mean(firstD)) / sd(firstD)
      #kurtosis = abs((kurtosis(normFirstD)))
      kurtosis = (kurtosis(normFirstD))
      return(mean(kurtosis))
    }
    
    ##calculando rugosidade pra cada imf
    vec <- vector()
    for (n in 1:length(lista)){
      vec[n] <- fc_curtose(lista[[n]])
    }
    #Curvatura máxima
    res <- daltoolbox::transform(daltoolbox::fit_curvature_max(), vec)
    div <- res$x
    #      print(vec)
  }
  
  plot(vec, t="p", ylab="", xlab="", col=ifelse(vec==vec[div], 'red', 'blue'))
  return(div)
}


div <- fc_div2(vec, "var")#normal com a variância (curvatura mínima)
div <- fc_div2(cum.vec, "curtose")#acumulado com a rugosidade (curvatura máxima)

#--- Testando outras abordagens de normalização e discriminação das imf´s

##Função que soma as IMFs dada IMF inicial e final
#fc_somaIMF <- function(ceemd.result, inicio, fim){
#  soma_imf <- rep(0, length(ceemd.result[["original.signal"]]))
#  for (k in inicio:fim){
#    soma_imf <- soma_imf + ceemd.result[["imf"]][,k]
#  }
#  return(soma_imf)
#}
#
##armazena as IMFs
#vec <- list()
#for (n in 1:ceemd.result$nimf){
#  vec[[n]] <- ceemd.result[["imf"]][,n]
#}
#
##cria as IMFs acumuladas
#cum.vec <- list()
#for (i in 1:ceemd.result[["nimf"]]){
#  cum.vec[[i]] <- fc_somaIMF(ceemd.result, 1, i)
#}
#
###Função que divide IMFs (acha div) 
#fc_div2 <- function(lista, mode="var"){
#  if(mode=="var"){
#    ##calculando var pra cada imf
#    vec <- vector()
#    for (n in 1:length(lista)){
#      vec[n] <- var(lista[[n]])
#    }
#    #Curvatura mínima
#    res <- daltoolbox::transform(daltoolbox::fit_curvature_min(), vec)
#    div <- res$x
#    print(div)
#    #      print(vec)        
#  }
#  if(mode=="curtose"){
#    ##Função de rugosidade
#    
#    fc_curtose <- function(x){
#      firstD = diff(x)
#      #normFirstD = (firstD - mean(firstD)) / sd(firstD) #padronização
#      normFirstD = (firstD - min(firstD)) / ( max(firstD) - min(firstD) ) #normalização
#      kurtosis = abs((kurtosis(normFirstD)))
#      #kurtosis = (kurtosis(normFirstD))
#      return(mean(kurtosis))
#    }
#    
#    ##calculando rugosidade pra cada imf
#    vec <- vector()
#    for (n in 1:length(lista)){
#      vec[n] <- fc_curtose(lista[[n]])
#    }
#    #Curvatura máxima
#    res <- daltoolbox::transform(daltoolbox::fit_curvature_max(), vec)
#    div <- res$x
#    #      print(vec)
#  }
#  
#  plot(vec, t="p", ylab="", xlab="", col=ifelse(vec==vec[div], 'red', 'blue'))
#  return(div)
#}
#
#
#div <- fc_div2(vec, "var")#normal com a variância (curvatura mínima)
#div <- fc_div2(cum.vec, "curtose")#acumulado com a rugosidade (curvatura máxima)
#

#serie volcano
serie = sig

#changepoint data
serie = head(wave.c44137,5000)
serie = head(HC1,5000)                       
serie = head(ftse100$V2,5000)  

id <- 1:length(serie)

#gerando IMFs
set.seed(123)#semente para rodar sempre o mesmo resultado (senão cada vez que roda vai ter um conjunto diferente de IMF)
ceemd.result <- try(CEEMD(serie, id, obj$noise.amp, obj$trials, verbose = TRUE))#roda CEEMD para gerar as IMFs

PlotIMFs(ceemd.result)#plota as IMFs com o sinal original e os ruídos

source("teste_fun.R")

fc_div2(cum.vec, "curtose")

#--- estudo das imfs

##Função que soma as IMFs dada IMF inicial e final
fc_somaIMF <- function(ceemd.result, inicio, fim){
  soma_imf <- rep(0, length(ceemd.result[["original.signal"]]))
  for (k in inicio:fim){
    soma_imf <- soma_imf + ceemd.result[["imf"]][,k]
  }
  return(soma_imf)
}

#cria as IMFs acumuladas
cum.vec <- list()
for (i in 1:ceemd.result[["nimf"]]){
  cum.vec[[i]] <- fc_somaIMF(ceemd.result, 1, i)
}


# criando vetor das imfs
vec <- list()
for (n in 1:ceemd.result$nimf){
  vec[[n]] <- ceemd.result[["imf"]][,n]
}

par(mfrow=c(3,4))

metricas_dist_imfs <- function(dados) {
  # Cálculo das métricas
  #dados = ( dados - min(dados) )/ (max(dados) - min(dados))
  dados = ( dados - mean(dados) )/ sd(dados) 
  assimetria <- skewness(dados)
  curtose <- kurtosis(dados)
  rugosidade <- mean((diff(dados) ** 2) / 4)
  
  # Visualização dos dados
  plot(dados)
  hist(dados)
  
  # Armazenamento das métricas em um dataframe
  metricas <- data.frame(
    Assimetria = assimetria,
    Curtose = curtose,
    Rugosidade = rugosidade
  )
  
  
  # Retorno do dataframe com as métricas
  return(assign("metricas",metricas,envir = .GlobalEnv))
}

# métricas da distribuição das imfs
metricas_lista2 = sapply(vec,metricas_dist_imfs)

df_metricas_lista2 = data.frame(t(metricas_lista2))

# imfs acumuladas e métricas não acumuladas
for (m in 1:length(names(df_metricas_lista2)) ) {
  
  hist(as.numeric(df_metricas_lista2[m][1][[1]] ))
  
}  

# imfs acumuladas e métricas acumuladas
for (m in 1:length(names(df_metricas_lista2)) ) {
  
  hist(cumsum(as.numeric(df_metricas_lista2[m][1][[1]] )) )
  
}  

# métricas da distribuição acumulada das imfs não acumuladas

metricas_lista1 = sapply(cum.vec,metricas_dist_imfs)

df_metricas_lista1 = data.frame(t(metricas_lista1))

par(mfrow=c(3,4))

# sem estar ordenado, imfs sem acumular e sem acumular métricas

for (m in 1:length(names(df_metricas_lista1)) ) {
  
  hist(as.numeric(df_metricas_lista1[m][1][[1]] ))
  plot(as.numeric(df_metricas_lista1[m][1][[1]] ))
  
}  

# sem estar ordenado, imfs acumuladas e sem acumular métricas

for (m in 1:length(names(df_metricas_lista2)) ) {
  
  hist(as.numeric(df_metricas_lista2[m][1][[1]] )) 
  plot(as.numeric(df_metricas_lista2[m][1][[1]] ))
}  


# ordenado, imfs não acumuladas e sem acumular métricas

for (m in 1:length(names(df_metricas_lista1)) ) {
  
  hist(sort(as.numeric(df_metricas_lista1[m][1][[1]] )) )
  plot(sort(as.numeric(df_metricas_lista1[m][1][[1]] )) )
  
}  

# ordenado, imfs acumuladas e acumulando métricas

for (m in 1:length(names(df_metricas_lista2)) ) {
  
  hist(sort(cumsum(as.numeric(df_metricas_lista2[m][1][[1]] )) ))
  plot(sort(cumsum(as.numeric(df_metricas_lista2[m][1][[1]] ))) )
}  

