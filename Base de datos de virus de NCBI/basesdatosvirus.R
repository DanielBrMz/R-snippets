library(seqinr)

setwd("C:\Users\danie\Downloads\BasicsExercises\Base de datos de virus de NCBI")

calcular_tamaño_genoma <- function(file_path) {

  secuencia <- read.fasta(file_path)
  
  tamaño <- sum(nchar(as.character(seq[[1]])))
  

  return(tamaño)
}

obtener_composicion_nucleotidos <- function(archivo_fasta) {
  seq <- read.fasta(archivo_fasta)[[1]]$seq
  n <- length(seq)
  A <- count(seq, "A") / n * 100
  T <- count(seq, "T") / n * 100
  G <- count(seq, "G") / n * 100
  C <- count(seq, "C") / n * 100
  
  composicion <- data.frame(A, T, G, C)
  
  return(composicion)
}

