library(GEOquery)
library(knitr)


# Analizamos el GPL para obtener datos sobre el experimento:

gpl <- getGEO("GPL14932")

# exploramos el objeto gpl:
str(gpl)

# @dataTable contiene la informacion sobre la tabla de variantes:  
tabla_datos_gpl<-gpl@dataTable

# Accedemos a @table 
table_gpl<-(tabla_datos_gpl@table)

# Vemos las primeras filas
head(table_gpl)

# El numero de filas de este df es el numero de variantes del array:
nrow(table_gpl)

# Tambien podemos verlo en el objeto @header
gpl@header$data_row_count

# Comprobamos el número de individuos del experimento
length(gpl@header$sample_id)


gpl@header

# Exploramos ahora una muestra para ver las variables que contiene el experimento:
gsm <- getGEO("GSM895465")


# Vemos si la muestra tiene otras variables fenotipicas de interés:
gsm@header$characteristics_ch1
# Esta muestra era de sangre y el individuo tenía alzheimer

# Vemos además la tabla con los datos deque tiene cada muestra:
head(Table(gsm))

