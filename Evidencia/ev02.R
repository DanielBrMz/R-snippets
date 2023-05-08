setwd("/Users/USER/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/Documents/Proyectos/RProjects/R-snippets/Evidencia")


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

"Obtén las secuencias de los genomas de los virus elegidos según la investigación que hayas decidido realizar. 
Para ello utiliza la función read.GenBank para obtener los genomas directamente del NCBI desde R Studio. 
Muestra el código empleado para obtenerlos."

coronavirus_sequences <- read.GenBank(variants_coronavirus)


"Realiza el alineamiento de los genomas virales y visualiza el resultado de tu alineamiento en tu navegador web. 
Muestra el código empleado para realizar lo anterior e incluye dos imágenes con el resultado del alineamiento, 
una de los primeros 150 nucleótidos y otra de los nucleótidos 500 al 650."

# Escribir las secuencias en formato fasta
write.dna(coronavirus_sequences, file = "coronavirus_sequences.fasta", format = "fasta")

# Leer las secuencias en formato fasta
coronavirus_string_set <- readDNAStringSet("coronavirus_sequences.fasta", format = "fasta")

# Orientar las secuencias
oriented_coronavirus_string_set <- OrientNucleotides(coronavirus_string_set)

# Alinear las secuencias 
aligned_coronavirus_string_set <- AlignSeqs(oriented_coronavirus_string_set)



first_150_nucleotides <- subseq(coronavirus_string_set, start = 1, end = 150)

first_150_nucleotides

nucleotides_500_650 <-  subseq(coronavirus_string_set, start = 500, end = 650)

nucleotides_500_650


BrowseSeqs(first_150_nucleotides)
BrowseSeqs(nucleotides_500_650)

"Agrega una interpretación escrita, desde el punto de vista biológico, para esta gráfica (de 6 a 12 renglones)."


"Las funciones BrowseSeqs permiten visualizar secuencias de virus relacionados con SARS-CoV-2 encontrados en animales. 
Al observar los primeros 150 nucleótidos, podemos identificar regiones conservadas y estructuras importantes. 
La visualización de nucleótidos 500-650 destaca regiones codificantes y regulatorias. Estos análisis ayudan a comprender 
la relación genética y funcional entre los virus y su posible implicación en la transmisión animal-humano. 
La investigación se centra en examinar secuencias de coronavirus en animales para explicar la posibilidad de transmisión. 
Interpretar las secuencias en el navegador proporciona información clave para identificar similitudes y diferencias entre 
los genomas, lo que puede ser relevante para comprender la transmisión y la adaptación de virus entre especies. Herramientas 
adicionales y análisis más detallados son necesarios para una comprensión completa."

"Genera una matriz de distancia a partir de los genomas alineados. Crea una tabla en escala de grises en la que observes de manera 
visual el resultado de la matriz de distancia e inclúyela en tu reporte. Muestra el código empleado para obtener lo anterior e incluye la matriz de distancia y la tabla que obtuviste."
