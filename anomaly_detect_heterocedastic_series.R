library(hht)
library(daltoolbox)
library(harbinger)
library(e1071)

#parâmetros do CEEMD
noise.amp=0.001
trials=5

obj <- harbinger()
obj$noise.amp <- noise.amp
obj$trials <- trials

#serie volcano
serie = sig
#serie sin_data
serie = sin_data$y
#serie sin_data
serie = tt



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
      kurtosis = abs((kurtosis(normFirstD)))
      #kurtosis = (kurtosis(normFirstD))
      return(mean(kurtosis))
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


div2 <- fc_div2(vec, "var")#normal com a variância (curvatura mínima)
div2 <- fc_div2(cum.vec, "curtose")#acumulado com a rugosidade (curvatura máxima)


