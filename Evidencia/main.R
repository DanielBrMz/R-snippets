setwd("/Users/danie/Downloads/BasicsExercises/Evidencia")


library(viridis)
library(Biostrings)
library(DECIPHER)
library(ade4)
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(viridis)
library(ggplot2)


names_sars_cov_2 <- c("SARS-CoV-2", "RaTG13", "RmYN02",
                      "Pangolin-CoV", "HKU2", "MERS-CoV",
                      "HCoV-NL63", "HCoV-OC43", "HCoV-HKU1", "SARS-CoV")


variants_coronavirus <- c("MN908947", "MT121215", "NC_045512", 
                           "MT121216", "NC_004718", "NC_019843",
                           "AY567487", "NC_006213", "NC_006577", "NC_004718")

names(variants_coronavirus) <- names_sars_cov_2


"Calcula la longitud de las secuencias de cada variante. Muestra el resultado 
impreso en consola y el código que utilizaste para realizar este punto.`"

coronavirus_sequences <- read.GenBank(variants_coronavirus)
sequence_lengths <- sapply(coronavirus_sequences, length)
names(sequence_lengths) <- variantes_coronavirus
sequence_lengths

"Genera una gráfica en la que se visualice la longitud de los 10 genomas
virales que analizaste. Muestra el código que utilizaste para realizar este punto."

barplot(sequence_lengths, main = "Sequence Lengths of Viral Genomes", xlab = "Viral Genome", ylab = "Sequence Length", col = rainbow(length(sequence_lengths)))

"Calcula la composición de nucleótidos a cada una de las variantes del virus. 
Muestra el resultado impreso en consola y el código que utilizaste para realizar este punto."

nucleotide_composition <- sapply(coronavirus_sequences, function(x) {
  table(unlist(strsplit(as.character(x), "")))
})

nucleotide_composition

 
