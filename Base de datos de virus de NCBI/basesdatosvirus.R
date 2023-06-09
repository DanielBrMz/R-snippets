install.packages("seqinr")
library(seqinr)
setwd("/Users/brmz/Downloads/BasicsExercises/Base de datos de virus de NCBI")

zika <- read.fasta("Zika.fasta")
dengue <- read.fasta("Dengue.fasta")
coronavirus <- read.fasta("SARS-Coronavirus.fasta")
sars_cov_2 <- read.fasta("SARS-COV-2.fasta")

virus <- c(zika, dengue, coronavirus, sars_cov_2)
names(virus) <- c("Zika", "Dengue", "Coronavirus", "Sars Cov 2")


# Mis funciones
calcular_porcentaje_bases <- function(secuencia){
  
  secuencia <- unlist(secuencia, use.names = FALSE)
  
  n <- length(secuencia)
  a <- 0
  c <- 0
  g <- 0
  t <- 0
  
  for (i in seq_along(secuencia)) {
    base <- secuencia[i]
    if (base == "a") {
      a <- a + 1
    } else if (base == "c") {
      c <- c + 1
    } else if (base == "g") {
      g <- g + 1
    } else if (base == "t") {
      t <- t + 1
    }
  }
  
  por_a <- a/n * 100
  por_c <- c/n * 100
  por_g <- g/n * 100
  por_t <- t/n * 100
  
  cat("Porcentaje de Adenina:", round(por_a, digits = 2), "%\n")
  cat("Porcentaje de Citosina:", round(por_c, digits = 2), "%\n")
  cat("Porcentaje de Guanina:", round(por_g, digits = 2), "%\n")
  cat("Porcentaje de Timina:", round(por_t, digits = 2), "%\n")
}

calcular_cantidad_bases <- function(secuencia){
  
  secuencia <- unlist(secuencia, use.names = FALSE)
  
  a <- 0
  c <- 0
  g <- 0
  t <- 0
  
  for (i in seq_along(secuencia)) {
    base <- secuencia[i]
    if (base == "a") {
      a <- a + 1
    } else if (base == "c") {
      c <- c + 1
    } else if (base == "g") {
      g <- g + 1
    } else if (base == "t") {
      t <- t + 1
    }
  }
  
  cat("Contenido de Adenina:", round(a, digits = 2), "%\n")
  cat("Contenido de Citocina:", round(c, digits = 2), "%\n")
  cat("Contenido de Guanina:", round(g, digits = 2), "%\n")
  cat("Contenido de Timina:", round(t, digits = 2), "%\n")
  
  return(list(c(a, c, g, t)))
}

contar_contenido_gc <- function(secuencia){
  
  secuencia <- unlist(secuencia, use.names = FALSE)
  
  gc <- 0
  
  for (i in seq_along(secuencia)){
    base <- secuencia[i]
    if(base == "g" || base == "c"){
      gc <- gc + 1
    } 
    
  }
  
}

calcular_hebra_complementaria <- function(secuencia){
  complementaria <- ""
  for (i in 1:nchar(secuencia)) {
    base <- substr(secuencia, i, i)
    if (base == "A") {
      complementaria <- paste0(complementaria, "T")
    } else if (base == "T") {
      complementaria <- paste0(complementaria, "A")
    } else if (base == "G") {
      complementaria <- paste0(complementaria, "C")
    } else if (base == "C") {
      complementaria <- paste0(complementaria, "G")
    }
  }
  return(complementaria)
}

hebra_complementaria <- calcular_hebra_complementaria(secuencia_aleatoria)

cat("Secuencia aleatoria: ", secuencia_aleatoria, "\n")
cat("Hebra complementaria: ", hebra_complementaria)


"1. ¿Cuál es el tamaño de cada genoma?"

calcular_tamaño_genoma <- function(genoma){
  
  length(genoma[[1]])
  
  calcular_porcentaje_bases(genoma[[1]])
}

calcular_tamaño_genoma(zika)
calcular_tamaño_genoma(dengue)
calcular_tamaño_genoma(coronavirus)
calcular_tamaño_genoma(sars_cov_2)
  
"2. ¿Cúal es la composición de nucleótidos (A, T, G y C) de cada genoma?"

calcular_composicion_genoma <- function(genoma) {
  nucleotidos <- c("Adenina: ", "Citosina: ", "Guanina: ", "Timina: ")
  composicion <- calcular_cantidad_bases(genoma[[1]], 1)
  
  for(i in seq_along(composicion)){
    cat("Contenido de ", nucleotidos[i], composicion[[i]], "\n")
  }
  
} 

calcular_composicion_genoma(zika)
calcular_composicion_genoma(dengue)
calcular_composicion_genoma(coronavirus)
names(calcular_composicion_genoma(sars_cov_2))

"3. ¿Cuál es el contenido de GC de cada virus?"

calcular_contenido_gc <- function(genoma){
  gc <- contar_contenido_gc(genoma)
  names(gc) <- ""
  cat("El contenido de gc del genoma es: ", gc)
}

calcular_contenido_gc(zika)
calcular_contenido_gc(dengue)
calcular_contenido_gc(coronavirus)
calcular_contenido_gc(sars_cov_2)

"4. Crear y aplica una función para obtener la secuencia complementaria e
imprimirla por cada secuencia, ejemplo:
Virus original: agttgttagt ctacgtggac cgacaagaac
Complementaria: tcaacaatca gatgcacctg gctgttcttg"

obtener_secuencia_complementaria <- function(genoma){
  
    cat("Secuencia complementaria: ", head(comp(genoma[[1]]), n =  15))
}

obtener_secuencia_complementaria(zika)
obtener_secuencia_complementaria(dengue)
obtener_secuencia_complementaria(coronavirus)
obtener_secuencia_complementaria(sars_cov_2)

"5. Crear una tabla comparativa a manera de resumen (marco de datos) en la que
se muestre la composición de nucleótidos de los 5 genomas virales."

mostrar_resumen_nucleotidos <- function(genomas){
  
  largo <- seq_along(genomas)
  tabla_resumen <- data.frame(adenina = largo,
                              citocina = largo,
                              guanina = largo,
                              timina = largo)
  
  rownames(tabla_resumen) <- names(genomas)
  
  for(i in seq_along(genomas)){
    tabla_resumen[i, ] <- count(genomas[i][[1]], 1)
  }
  
  return(tabla_resumen)
}

mostrar_resumen_nucleotidos(virus)


# Paquete seqinr

"1. ¿Cuál es el tamaño de cada genoma?"

calcular_tamaño_genoma <- function(genoma){
  
  length(genoma[[1]])
  
  calcular_porcentaje_bases(genoma[[1]])
}

calcular_tamaño_genoma(zika)
calcular_tamaño_genoma(dengue)
calcular_tamaño_genoma(coronavirus)
calcular_tamaño_genoma(sars_cov_2)

"2. ¿Cúal es la composición de nucleótidos (A, T, G y C) de cada genoma?"

calcular_composicion_genoma <- function(genoma) {
  nucleotidos <- c("Adenina: ", "Citosina: ", "Guanina: ", "Timina: ")
  composicion <- count(genoma[[1]], 1)
  seq_along(composicion)
  
  for(i in seq_along(composicion)){
    cat("Contenido de ", nucleotidos[i], composicion[[i]], "\n")
  }
}

"3. ¿Cuál es el contenido de GC de cada virus?"
  
calcular_contenido_gc <- function(genoma){
  gc <- count(genoma[[1]], 1)[2] + count(genoma[[1]], 1)[3]
  names(gc) <- ""
  cat("El contenido de gc del genoma es: ", gc)
}
  
calcular_contenido_gc(zika)
calcular_contenido_gc(dengue)
calcular_contenido_gc(coronavirus)
calcular_contenido_gc(sars_cov_2)

"4. Crear y aplica una función para obtener la secuencia complementaria e
imprimirla por cada secuencia, ejemplo:
Virus original: agttgttagt ctacgtggac cgacaagaac
Complementaria: tcaacaatca gatgcacctg gctgttcttg"

obtener_secuencia_complementaria <- function(genoma){
  
  cat("Secuencia complementaria: ", head(comp(genoma[[1]]), n =  15))
}

obtener_secuencia_complementaria(zika)
obtener_secuencia_complementaria(dengue)
obtener_secuencia_complementaria(coronavirus)
obtener_secuencia_complementaria(sars_cov_2)

"5. Crear una tabla comparativa a manera de resumen (marco de datos) en la que
se muestre la composición de nucleótidos de los 5 genomas virales."

mostrar_resumen_nucleotidos <- function(genomas){
  
  largo <- seq_along(genomas)
  tabla_resumen <- data.frame(adenina = largo,
                              citocina = largo,
                              guanina = largo,
                              timina = largo)
  
  rownames(tabla_resumen) <- names(genomas)
  
  for(i in seq_along(genomas)){
    tabla_resumen[i, ] <- count(genomas[i][[1]], 1)
  }
  
  return(tabla_resumen)
}

mostrar_resumen_nucleotidos(virus)
  
