library(tidyverse)
library(janitor)
library(readxl)
library(distances)

url <- "https://github.com/jgbabativam/Curso_Multivariado/raw/main/Datos/"
url_data <- paste0(url, "ARWU_100_top.csv")
df <- read.csv2(url_data)

glimpse(df)
head(df)

top_20 <- df %>% arrange(df$World.Rank) %>% slice_head(n=20)
ggplot(data=top_20, aes(x= reorder(Institution, World.Rank), y=Award)) + geom_col()  


# Utilizar los datos del archivo r14_Sci_Qs_Webometrics.csv para elaborar gráficos de dispersión de

url2 <- "https://github.com/jgbabativam/Curso_Multivariado/raw/main/Datos/r14_Sci_Qs_Webometrics.csv"
dfweb <- read.csv2(url2, sep = ';')

ggplot(data=dfweb, aes(x=SC.Lac.Ranking, y=QS.Ranking)) + geom_point()
cor(dfweb$SC.Lac.Ranking, dfweb$QS.Ranking, use="complete.obs", method="kendall")

# Con los datos del ARWU_100_top.csv elaborar una matriz de dispersión con

table1 <- pairs(dfweb[, c("SC.Lac.Ranking", "SC.Ibe.Ranking","SC.Co.Ranking","SC.Productividad","SC.Colaboracion.Interciol","SC.Impacto.normalizado","SC.Publicaciones.de.alta.calidad","SC.Indice.de.especializacion","SC.Indice.de.excelencia","SC.Liderazgo.cientifico","SC.Excelencia.con.liderazgo")])

# Con los datos del archivo datos_ciudades.xlsx elaborar diagramas de cajas (Boxplots) para visualizar si hay datos atípicos en las variables

url <- "C:/Users/johan/Downloads/datos_ciudades.xlsx"
df3 <- read_excel(url)

boxplot(df3$CYT_21, df3$CYT_20, df3$CYT_19, df3$CYT_18, df3$CYT_17, names=c("CYT_21","CYT_20","CYT_19", "CYT_18", "CYT_17"))
    
# sí hay, en cyt20 y cyt18
quartiles <- quantile(df3$CYT_20)
quartiles
RIQ_CYT19_sup <-  quartiles["75%"] + (quartiles["75%"] - quartiles["25%"]) * 1.5
RIQ_CYT19_sup

# Calcular la matriz de correlación entre las variables del conjunto que le correspondió al grupo en el ejercicio 5 y escoger las dos variables que tienen mayor correlación.

cor(df3[, c("CYT_21","CYT_20","CYT_19", "CYT_18", "CYT_17")], use="complete.obs")

# las dos con más correlacion "CYT_20" y "CYT_21"

#Calcular la distancia euclidiana entre San Andrés y Riohacha con respecto a estas
#dos variables y luego calcular la distancia de Mahalanobis 
#entre las mismas ciudades respecto a las mismas dos variables.

datos_completos <- df3 %>% select(c("CIUDADES", "CYT_21", "CYT_20")) %>%
  mutate(across(c("CYT_21", "CYT_20"), as.numeric)) %>%
  column_to_rownames("CIUDADES")

dist_mahalanobis_todas <- distances(datos_completos, normalize = "mahalanobize", id_variable = row.names(datos_completos))
dist_euclidiana_todas  <- distances(datos_completos, normalize = "none",, id_variable = row.names(datos_completos))


dist_mahalanobis_todas
# la mahalanobis es 5.3143938

dist_euclidiana_todas
# la ecludiana es de 0.807051
# la difetencia sería que la correlacion entre estas no es fuerte


# matriz de correlacion a mano

rh_y_cyt <- merge(df3[,1:15], df3[, c("CIUDADES","CYT_21","CYT_20","CYT_19", "CYT_18", "CYT_17")], by="CIUDADES")
rh_y_cyt <- rh_y_cyt[2:15]
rh_y_cyt_m <- as.matrix(rh_y_cyt)
means <- colMeans(rh_y_cyt_m)

datos_centrados = sweep(rh_y_cyt_m, 2, means, FUN = "-")
n <- nrow(rh_y_cyt_m)


S <- ((t(datos_centrados) %*% datos_centrados) /n)
S

# matriz de covarianza con R
cov(rh_y_cyt_m)

# Difieren un poco en los resultados debido a que cov() usa un estimador 
# insesgado (1/n-1), mientras que manuealmente usamos uno asintoticamente insesgado


# ahora la de correacion
# la matriz de datos centrados y normalizados es estandarizar las matriz de datos centrados
# desviaciones

desviaciones <- sqrt(diag(S))
cen_es <- sweep(datos_centrados, 2, desviaciones, FUN = "/")

corr_matrix <- (t(cen_es) %*% cen_es) / n  
corr_matrix

# ahora corr en R
cor(rh_y_cyt_m)
# sí coinciden,

n <- 3
tmp <- S
diag(tmp) <- NA  
sort(tmp, decreasing = TRUE)
# los más son rh7 y rh1 y rh9 rh1
# los menos son rh10 y rh1, rh1 y rh14

tmp1 <- corr_matrix
diag(tmp1) <- NA
sort(tmp1, decreasing = TRUE)

# Convertir la matriz a tabla de pares
pares <- function(mat){
  diag(mat) <- NA          # ignorar diagonal
  mat[lower.tri(mat)] <- NA # evitar duplicados
  
  # Convertir a dataframe largo
  df <- as.data.frame(as.table(mat))
  df <- na.omit(df)
  colnames(df) <- c("Var1", "Var2", "Valor")
  df <- df[order(df$Valor, decreasing=TRUE), ]
  rownames(df) <- NULL
  df
}

tabla_cov  <- pares(S)
tabla_corr <- pares(corr_matrix)

# Mayor y menor covarianza
head(tabla_cov, 2)   # mayor
tail(tabla_cov, 2)   # menor

# Mayor y menor correlación
head(tabla_corr, 2)  # mayor
tail(tabla_corr, 2)  # menor

# RH_2 RH_5, RH_1 RH_7 la mayor
# RH_5  RH_9, RH_6 RH_16 la menor

# Construir una matriz de datos con los puntajes en los seis criterios del
# ranking ARWU para las universidades del grupo asignado y comprobar que con el TDVS se puede reconstruir.

X <- na.omit(df[4:length(df)])
descomposicion <- svd(X)

U <- descomposicion$u
D <- diag(descomposicion$d) # Convierte el vector en una matriz diagonal
V <- descomposicion$v

# Ahora las dimensiones encajarán perfectamente
X_reconstruida <- U %*% D %*% t(V)

round(X - X_reconstruida, 6)


#----------- no es parte del lab ---------------


# Datos de ejemplo
# genero <- c("M", "F", "M", "F", "M", "F", "M")
# respuesta <- c("Sí", "No", "Sí", "Sí", "No", "No", "Sí")

# Tabla de contingencia
# tabla <- table(genero, respuesta)
# tabla


#ggplot(datos, aes(x = peso, y = altura)) +
#geom_point()
