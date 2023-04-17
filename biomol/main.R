# Crea un vector x que contenga 10 valores enteros aleatorios entre -100 y +50
set.seed(123);

x <- sample(-100:50, 10, replace = TRUE);

print(x)

# Calcula datos estadísticos simples de x: la media, la desviación estándar y la varianza.
media <- mean(x)
desviacionEstandar <- sd(x)
varianza <- var(x)

# Crear un objeto de tipo vector con el nombre est.x en el que guardes los 3 estadísticos calculados para x
est.x <- c(media,desviacionEstandar,varianza);
print(est.x)

# Calcula la suma de los números contenidos en el vector x.
sumX <- sum(x)
sumX

# Escribe una instrucción en R para encontrar el valor máximo y mínimo del vector x
valorMaximoX <- max(x)
valorMinimoX <- min(x)

print(valorMaximoX)
print(valorMinimoX)

# Crea dos vectores a y b de tipo entero y longitud 2, de la misma longitud ambos, multiplícalos y muestra el resultado de la multiplicación.
a = sample(-10:10, 2, replace = TRUE)
b = sample(-10:10, 2, replace = TRUE)

productoAB <- a*b

print(productoAB)

# Genera un vector con una secuencia de números del 1 al 10 de 0.1 en 0.1 y nómbralo z.
z = seq(1,10, by=0.1)

print(z)

# Genera una matriz llamada m de dimensiones 3, 4 que contenga los números del 1 al 12 empezando a llenar datos por las filas
m <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)

print(m)

"Genera un marco de datos llamado mdd de dimensiones 4, 3 que contenga: una
columna con datos 1, 2, 3, 4; una columna con datos “a”, “b”, “c”, “d” y una columna
con datos FALSE, FALSE, TRUE, TRUE. "
mdd <- data.frame(col1 = 1:4, col2 = letters[1:4] , col3 = c(FALSE, FALSE, TRUE, TRUE))

head(mdd)

# Genera una lista de nombre lista que contenga los objetos: x, est.x, a, b, z, m y mdd e imprime su contenido.
lista <- list(x = x, est.x = est.x, a = a, b = b, z = z, m = m, mdd = mdd)

print(lista)




m <- matrix(1:6, 2, 3)

m








