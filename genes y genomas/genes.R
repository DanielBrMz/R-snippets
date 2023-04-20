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


hmp <- "ATGCTTGACGCTCAAACCATCGCTACAGTAAAAGCCACCATCCCTTTACTGGTGGAAACGGGGCCAAAGT
TAACCGCCCATTTCTACGACCGTATGTTTACTCATAACCCAGAACTCAAAGAAATTTTTAACATGAGTAA
CCAGCGTAATGGCGATCAACGTGAAGCCCTGTTTAACGCTATTGCCGCCTACGCCAGTAATATTGAAAAC
CTGCCTGCGCTGCTGCCAGCGGTAGAAAAAATCGCGCAGAAGCACACCAGCTTCCAGATCAAACCGGAAC
AGTACAACATCGTCGGTGAACACCTGTTGGCAACGCTGGACGAAATGTTCAGCCCGGGGCAGGAAGTGCT
GGACGCGTGGGGTAAAGCCTATGGTGTACTGGCTAATGTATTTATCAATCGCGAGGCGGAAATCTATAAC
GAAAACGCCAGCAAAGCCGGTGGTTGGGAAGGTACTCGCGATTTCCGCATTGTGGCTAAAACACCGCGCA
GCGCGCTTATCACCAGCTTCGAACTGGAGCCGGTCGACGGTGGCGCAGTGGCAGAATACCGTCCGGGGCA
ATATCTCGGCGTCTGGCTGAAGCCGGAAGGTTTCCCACATCAGGAAATTCGTCAGTACTCTTTGACTCGC
AAACCGGATGGCAAAGGCTATCGTATTGCGGTGAAACGCGAAGAGGGTGGGCAGGTATCCAACTGGTTGC
ACAATCACGCCAATGTTGGCGATGTCGTGAAACTGGTCGCTCCGGCAGGTGATTTCTTTATGGCTGTCGC
AGATGACACACCAGTGACGTTAATCTCTGCCGGTGTTGGTCAAACGCCAATGCTGGCAATGCTCGACACG
CTGGCAAAAGCAGGCCACACAGCACAAGTGAACTGGTTCCATGCGGCAGAAAATGGCGATGTTCACGCCT
TTGCCGATGAAGTTAAGGAACTGGGGCAGTCACTGCCGCGCTTTACCGCGCACACCTGGTATCGTCAGCC
GAGCGAAGCCGATCGCGCTAAAGGTCAGTTTGATAGCGAAGGTCTGATGGATTTGAGCAAACTGGAAGGT
GCGTTCAGCGATCCGACAATGCAGTTCTATCTCTGCGGCCCGGTTGGCTTCATGCAGTTTACCGCGAAAC
AGTTAGTGGATCTGGGCGTGAAGCAGGAAAACATTCATTACGAATGCTTTGGCCCGCATAAGGTGCTGTA
A"

gene_hmp <- gsub(pattern = "\n", replacement = "", hmp)


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

encontrar_codones <- function(secuencia){
  codon_inicio = "AUG"
  codon_final = c("UAA", "UAG", "UGA")
  
  secuencia <- paste0(codon_inicio, secuencia, "UAAUAGUGA")
  
  codones <- c();
  
  
  for(i in  seq(1, nchar(secuencia), by=3)){
    
    j <- (i+2)/3;
    
    codones[j] <- paste0(substr(secuencia, i, i), substr(secuencia, i+1, i+1), substr(secuencia, i+2, i+2));
  }
  
  return(codones)
}


traducir_arn <- function(codones){
  tabla_aminoacidos <- data.frame(
    codon = c("UUU", "UUC", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "AUU", "AUC", "AUA", "AUG",
              "GUU", "GUC", "GUA", "GUG", "UCU", "UCC", "UCA", "UCG", "CCU", "CCC", "CCA", "CCG",
              "ACU", "ACC", "ACA", "ACG", "GCU", "GCC", "GCA", "GCG", "UAU", "UAC", "UAA", "UAG",
              "CAU", "CAC", "CAA", "CAG", "AAU", "AAC", "AAA", "AAG", "GAU", "GAC", "GAA", "GAG",
              "UGU", "UGC", "UGA", "UGG", "CGU", "CGC", "CGA", "CGG", "AGU", "AGC", "AGA", "AGG",
              "GGU", "GGC", "GGA", "GGG"),
    amino_acido = c("F", "F", "L", "L", "L", "L", "L", "L", "I", "I", "I", "M",
                    "V", "V", "V", "V", "S", "S", "S", "S", "P", "P", "P", "P",
                    "T", "T", "T", "T", "A", "A", "A", "A", "Y", "Y", "STOP", "STOP",
                    "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E",
                    "C", "C", "STOP", "W", "R", "R", "R", "R", "S", "S", "R", "R",
                    "G", "G", "G", "G")
  )
  
  cat("Aminoacidos del gen hmp: \nSTART ")
  for(i in 2:length(codones)){
    
    codon_encontrado <- codones[i] == tabla_aminoacidos$codon
    posicion <- which(codon_encontrado == TRUE)
    

    cat(tabla_aminoacidos$amino_acido[posicion], " ")
  } 
  
}

arn_hmp <- transcribir_adn_arn(gene_hmp)
codones_hmp <- encontrar_codones(arn_hmp)
traducir_arn(codones_hmp)
traducir_arn(codones_hmp)

