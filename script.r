### LIBRERÍAS ###

library(tidyverse)
library(broom)
library(rbokeh)
library(stringr)
library(rebus)
library(stringi)
library(Rcmdr)

### EDA ###

alz <- read.csv("meg_mci.csv", sep = ",")

View(alz[,1:20])

str(alz[,1:20])

tail(colnames(alz))

# Extraemos en un vector los diferentes tipos de estadísticos aplicados. Todos los nombres de las variables
# siguen la estructura "sync_X_X.Y" donde X representa el número de los sensores involucrados e Y el estadístico
# estudiado. Si separamos cada nombre por el ".", guardamos toda la información en un vector y se eliminan los
# valores repetidos se obtendrán los 5 estadísticos utilizados.
filtro <- sapply(colnames(alz)[3:length(colnames(alz))], 
                 function(x) strsplit(x, "[.]")[[1]][2],
                 USE.NAMES = F)

estadisticos <- unique(filtro)

### MANIPULACIÓN DE TEXTO ###

# Creamos una tabla para cada estadístico que relacione cada par de sensores con los estadísticos obtenidos
# ind_X es un vector que indica qué columnas hacen referencia a un estadístico determinado.
ind_media   <- str_detect(colnames(alz), "mean")
ind_std     <- str_detect(colnames(alz), "std")
ind_mediana <- str_detect(colnames(alz), "median")
ind_mad     <- str_detect(colnames(alz), "mad")
ind_cov     <- str_detect(colnames(alz), "cov")

# Extraemos todos los caracteres antes del número de los sensores
aux_colnames <- str_replace_all(colnames(alz), 
                                "sync_", 
                                replacement = "")

# En las columnas relativas a las medias de sincronización de sensores, eliminamos todo caracter que aparezca
# después del "."
aux_media <- str_replace(aux_colnames[ind_media],
                         DOT %R% "mean",
                         "")

# Aplicamos la prueba t a todas las medias en función del grupo (control o MCI) y extraemos el valor de p-value
# en el vector "resultados_ttest"
resultados_ttest <- c()

for(i in which(ind_media == T)){
  resultado <- testme(alz[,i])
  resultados_ttest[i-2] <- resultado$p.value
  if(i %% 100 == 0){ print(i) }
}

names(resultados_ttest) <- colnames(alz)[ind_media]
names(resultados_ttest) <- c()

# Comprobamos si las desviaciones estándar entre grupos difieren significativamente
leveneTest(alz$sync_1_2.std, alz$class, center = mean)$'Pr(>F)'[1]

resultados_levene <- c()
for (i in which(ind_std)){
  resultado <- leveneTest(alz[,i], alz$class, center = mean)
  resultados_levene[which(which(ind_std) == i)] <- resultado$'Pr(>F)'[1]
}

resultados_wilcox <- c()
for (i in which(ind_mediana)){
  resultado <- wilcox.test(alz[,i] ~ class, data = alz)
  resultados_wilcox[which(which(ind_mediana) == i)] <- resultado$'p.value'
}

### PREPROCESO DE DATOS ###

# Creamos una tabla con 2 columnas respectivas a los sensores estudiados
tabla <- str_split(aux_media, "_", simplify = TRUE)

# Añadimos el valor del p-value obtenido
tabla_p.value <- cbind(tabla, resultados_ttest)

# Convertimos la tabla de matriz a dataframe
tabla_p.value            <- as.data.frame(tabla_p.value)
tabla_p.value            <- as.data.frame(apply(tabla_p.value, 2, as.numeric))
tabla_p.value$sign       <- as.factor(ifelse(tabla_p.value$resultados_ttest < 0.05, 1, 0))
tabla_p.value$ttest_mod  <- ifelse(tabla_p.value$resultados_ttest > 0.05, 
                                   0.05, 
                                   tabla_p.value$resultados_ttest)
tabla_p.value            <- cbind(tabla_p.value, resultados_levene)
tabla_p.value$levene_mod <- ifelse(tabla_p.value$resultados_levene > 0.05, 
                                   0.05, 
                                   tabla_p.value$resultados_levene) 


# sync_mean <- matrix(data = nrow = 102, ncol = 102)
# 
# sync_p.value <- matrix(resultados_ttest, nrow = 102, ncol = 102)

# Esta tabla no la he utilizado aún
tabla_p.value_1 <- matrix(nrow = 102, ncol = 102)

for(i in 1:102){
  for (j in 1:102){
    if(i == j){
      next()
    }
    if(length(tabla_p.value[tabla_p.value$V1 == i & tabla_p.value$V2 == j, "resultados_ttest"]) == 0){
      tabla_p.value_1[i,j] = tabla_p.value[tabla_p.value$V1 == j & tabla_p.value$V2 == i, "resultados_ttest"] 
    } else {
      tabla_p.value_1[i,j] = tabla_p.value[tabla_p.value$V1 == i & tabla_p.value$V2 == j, "resultados_ttest"]      
    }
  }
}

tabla_p.value_1 <- as.data.frame(tabla_p.value_1)

### VISUALIZACIÓN ### 

ggplot(tabla_p.value, aes(x = V2, y = V1, fill = sign)) +  
  geom_tile(aes(fill = sign), colour = "grey")

ggplot(tabla_p.value[tabla_p.value$sign == 1,], aes(x = V2, y = V1, fill = resultados_ttest)) +  
  geom_tile(aes(fill = resultados_ttest), colour = "grey")

ggplot(tabla_p.value, aes(x = V2, y = V1, fill = ttest_mod)) +  
  geom_tile(aes(fill = ttest_mod), colour = "grey")

ggplot(tabla_p.value, aes(x = V2, y = V1, fill = levene_mod)) +  
  geom_tile(aes(fill = levene_mod), colour = "grey")

ggplot(alz, aes(x = sync_1_17.mean)) + geom_histogram() + facet_wrap(~class)




resultado_1 <- apply(alz[,3:ncol(alz)], 
                     2, 
                     function(x) t.test(x ~ class, data = alz) %>% 
                       tidy())
resultado_2 <- lapply(resultado_1, tidy)


### CREACIÓN DE FUNCIONES ### 
# "showme" te permite sacar el histograma de una variable por cada grupo de sujetos (class)
showme <- function(var) {
  df <- data.frame(sync = alz[,var], class = alz[, "class"]) 
  ggplot(df, aes(x = sync)) + geom_histogram() + facet_wrap(~class)
}

showme("sync_1_19.mean")

# "testme" te permite hacer un t.test entre grupo control y grupo afectado simplemente 
# indicando la variable
testme <- function(x){
  t.test(x ~ class, data = alz) %>% tidy()
}

for(i in 3:ncol(alz)){
  resultados[i-2] <- testme(alz[,i])
  if(i %% 100 == 0){ print(i) }
}
resultado$p.value



testme(alz$sync_1_2.mean)

colnames(alz)[3:length(colnames(alz))]


length(resultados[resultados< 0.05])
plot(head(resultados[resultados < 0.05]))
plot(resultados[resultados < 0.05])

head(resultados_ttest)

str_extract(colnames(alz), 
            one_or_more(DGT) %R% 
              one_or_more(DGT) %R% 
              DOT %R% 
              ANY_CHAR %R% 
              ANY_CHAR %R% 
              ANY_CHAR %R% 
              ANY_CHAR)

colnames(alz)[300:310]
