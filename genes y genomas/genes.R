secuencia_establecida <- "ATGCTTGACGCTCAAACCATCGC"


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

calcular_hebra_reversa <- function(secuencia){
  inversa <- ""
  for (i in nchar(secuencia):1) {
    inversa <- paste0(inversa, substr(secuencia, i, i))
  }
  return(reversa)
}

hebra_complementaria <- calcular_hebra_complementaria(secuencia_establecida)
hebra_reversa <- calcular_hebra_inversa(secuencia_establecida)
hebra_reversa_complementaria <- calculaHebraComplementaria(hebra_reversa)

cat("Secuencia original: ", secuencia_establecida, "\n")
cat("Hebra complementaria: ", hebra_complementaria)
cat("Hebra reversa: ", hebra_inversa)
cat("Hebra reversa complementaria: ", hebra_reversa_complementaria)


