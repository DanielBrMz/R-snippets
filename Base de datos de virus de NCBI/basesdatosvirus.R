install.packages("seqinr")
library(seqinr)
setwd("/Users/danie/Downloads/BasicsExercises/Base de datos de virus de NCBI")

zika <- read.fasta("Zika.fasta")
dengue <- read.fasta("Dengue.fasta")
coronavirus <- read.fasta("SARS-Coronavirus.fasta")
sars_cov_2 <- read.fasta("SARS-COV-2.fasta")

virus <- c(zika, dengue, coronavirus, sars_cov_2)
names(virus) <- c("Zika", "Dengue", "Coronavirus", "Sars Cov 2")

"1. ¿Cuál es el tamaño de cada genoma?"

calcular_tamaño_genoma <- function(genoma){
  return(length(genoma[[1]]))
}

calcular_tamaño_genoma(zika)
calcular_tamaño_genoma(dengue)
calcular_tamaño_genoma(coronavirus)
calcular_tamaño_genoma(sars_cov_2)

"2. ¿Cúal es la composición de nucleótidos (A, T, G y C) de cada genoma?"

calcular_composicion_genoma <- function(genoma) {
  return(count(genoma[[1]], 1))
} 

calcular_composicion_genoma(zika)
calcular_composicion_genoma(dengue)
calcular_composicion_genoma(coronavirus)
names(calcular_composicion_genoma(sars_cov_2))

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


