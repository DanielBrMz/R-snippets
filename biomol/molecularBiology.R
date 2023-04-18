# PARTE 1

"Crea una función que genere una secuencia aleatoria de nucleótidos de ADN (A, T,
G y C) de tamaño “n”.
Ejecútala solicitando una secuencia aleatoria de ADN de 30 nucleótidos y muestra
el resultado impreso en consola."

crear_secuencia_adn <- function(n = 30) {
  
  nucleotidos <- c("A", "T", "G", "C")
  secuencia <- sample(nucleotidos, n, replace=TRUE)
  
  
  return(paste(secuencia, collapse=""))
}


secuencia_aleatoria <- crear_secuencia_adn(30);
print(secuencia_aleatoria)

"Crea una función que calcule el tamaño de una secuencia de ADN.
Utilízala para calcular el tamaño de la secuencia que generaste en el punto 1 y
muestra el resultado impreso en consola."

calcular_tamaño_secuencia <- function(secuencia) {
  return(nchar(secuencia))
}

tamaño_secuencia <- calcular_tamaño_secuencia(secuencia_aleatoria)

cat("El tamaño de la secuencia de ADN es:", tamaño_secuencia)

"Crea una función que recibe una secuencia de DNA e imprime el porcentaje de
cada base (A, C, G y T) en la secuencia. Ejecútala sobre la secuencia que generaste
en el punto 1 y muestra el resultado impreso en consola."

calcular_porcentaje_bases <- function(secuencia){
  
  n <- nchar(secuencia)
  a <- 0
  c <- 0
  g <- 0
  t <- 0
  
  for (i in 1:n) {
    base <- substr(secuencia, i, i)
    if (base == "A") {
      a <- a + 1
    } else if (base == "C") {
      c <- c + 1
    } else if (base == "G") {
      g <- g + 1
    } else if (base == "T") {
      t <- t + 1
    }
  }

  por_a <- a/n * 100
  por_c <- c/n * 100
  por_g <- g/n * 100
  por_t <- t/n * 100
  
  cat("Porcentaje de bases A:", round(por_a, digits = 2), "%\n")
  cat("Porcentaje de bases C:", round(por_c, digits = 2), "%\n")
  cat("Porcentaje de bases G:", round(por_g, digits = 2), "%\n")
  cat("Porcentaje de bases T:", round(por_t, digits = 2), "%\n")
}

calcular_porcentaje_bases(secuencia_aleatoria)

"Crea una función que recibe una hebra directa y regresa la hebra inversa.
Ejecútala sobre la secuencia que generaste en el punto 1 y muestra el resultado
impreso en consola"

calcular_hebra_inversa <- function(secuencia){
  inversa <- ""
  for (i in nchar(secuencia):1) {
    inversa <- paste0(inversa, substr(secuencia, i, i))
  }
  return(inversa)
}

hebra_inversa <- calcular_hebra_inversa(secuencia_aleatoria)
cat("Secuencia aleatoria: ", secuencia_aleatoria, "\n")
cat("Hebra inversa: ", hebra_inversa)


# Parte 2

"Crea una función qué recibe una hebra directa y obtiene la hebra complementaria.
Ejecútala sobre la secuencia que generaste en el punto 1 y muestra el resultado
impreso en consola"

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
cat("Hebra complementaria:", hebra_complementaria)


"Crea una función que transcribe ADN a ARN. Es decir, que recibe una secuencia
molde de nucleótidos de ADN y devuelve la secuencia del transcrito
correspondiente que se transcribiría a partir de dicha secuencia.
Ejecútala sobre la secuencia que generaste en el punto 1 y muestra el resultado
impreso en consola."

transcribir_adn_arn <- function(secuencia) {
  arn <- ""
  for(i in 1:nchar(secuencia)){
    base <- substr(secuencia, i, i)
    if(base == "T"){
      arn <- paste0(arn, "U")
    }else{
      arn <- paste0(arn, base)
    }
  }
  return(arn)
}

arn_transcrito <- transcribir_adn_arn(secuencia_aleatoria)
cat("Secuencia aleatoria: ", secuencia_aleatoria, "\n")
cat("Secuencia transcrita a ARN: ", arn_transcrito)

"Crea una función que te sirva para encontrar codones de inicio y de terminación
en una secuencia dada de ARN. A la secuencia presente entre un codón de
inicio y uno de terminación en un gen se le conoce como marco de lectura.
Ejecútala sobre la secuencia que generaste en el punto anterior (6) y muestra el
resultado impreso en consola."


encontrar_codones <- function(secuencia){
  codon_inicio = "AUG"
  codon_final = c("UAA", "UAG", "UGA")
  
  secuencia <- paste0(codon_inicio, secuencia, "UAAUAGUGA")

  codones <- c();
  
  
  for(i in  seq(1, nchar(secuencia), by=3)){
    
    j <- (i+2)/3;
    
    codones[j] <- paste0(substr(secuencia, i, i), substr(secuencia, i+1, i+1), substr(secuencia, i+2, i+2));
  }
  print(codones)
  print(length(codones))
  
  for(i in 1:length(codones)){
    if(codones[i] == codon_inicio){
      print("START")
    }else if(codones[i] %in% codon_final){
      print("STOP")
      break
    } else {
      print(codones[i])
    }
  }
  
}

encontrar_codones(arn_transcrito)