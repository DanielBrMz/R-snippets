setwd("/Users/danie/Downloads/BasicsExercises/Árboles Filogenéticos en R")


"Instala los paquetes y carga las librerías necesarias para hacer un análisis filogenético de genomas virales en R.
Por ejemplo, ejecuta los siguientes comandos una vez que hayas instalado los paquetes necesarios:"

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

"Obtén las secuencias de 10 genomas virales a partir del código de acceso al genoma en GenBank usando la
función read.GenBank del paquete ape y muestra el código empleado para hacerlo."

coronavirus <- c("MT810758", "AY278489", "MW030193", "AY508724", "MT594401",
                  
                  "AY485277", "MT292571", "AY390556", "MN908947","JX869059",

                  "MW133981", "MT577009", "MT835383", "MW056032", "MT470219")


coronavirus_sequences <- read.GenBank(coronavirus)


"Imprime en consola la estructura del objeto de clase DNABin que obtuviste con la función read.GenBank.
Muestra el código empleado para hacerlo y la estructura del objeto DNAbin que obtuviste."

str(coronavirus_sequences)

"Concentra en un archivo todas las secuencias de los genomas con la función write.dna del paquete ape y
muestra el código empleado para hacerlo."

write.dna(secuencias_coronavirus, file = "coronavirus_sequences.fasta", format = "fasta")

"Carga las secuencias concentradas en el archivo del punto anterior con la función readDNAStringSet de
Biostrings. Muestra el código empleado para hacerlo, imprime en consola el contenido del objeto tipo
DNAStringSet e inclúyelo en tu entregable."

my_dna_string_set <- readDNAStringSet("coronavirus_sequences.fasta", format = "fasta")

my_dna_string_set

"Orienta los nucleótidos de los genomas con la función OrientNucleotides del paquete DECIPHER y muestra
el código empleado para hacerlo."

oriented_dna_string_set <- OrientNucleotides(my_dna_string_set)

oriented_dna_string_set

"Realiza el alineamiento de las secuencias de los genomas virales con la función AlignSeqs y visualiza el resultado
del alineamiento en tu navegador con la función BrowseSeqs de DECIPHER"

aligned_dna_string_set <- AlignSeqs(oriented_dna_string_set)

BrowseSeqs(aligned_dna_string_set)

"Guardar el resultado del alineamiento en un archivo en formato .fasta con la función writeXStringSet de
Biostrings"

writeXStringSet(aligned_dna_string_set, "coronavirus_sequences_aligned.fasta", format="fasta")

"Carga el archivo .fasta que construiste en el paso anterior con la función read.alignment de seqinr."

alignment <- read.alignment(file = "coronavirus_sequences_aligned.fasta", format = "fasta")

alignment

"Crear una matriz de distancia con la función dist.alignment de seqinr y obtén una tabla en escala de grises con la
función table.paint del paquete ape4 en donde sombras más oscuras de gris significan una mayor distancia.
Muestra el contenido de la matriz de distancia que creaste, la imagen de la tabla en escala de grises que
construiste, y el código empleado para conseguir todo lo anterior"

dist_matrix <- dist.alignment(alignment, matrix = "similarity")

dist_matrix

temp <- as.data.frame(as.matrix(dist_matrix))

class(temp)

typeof(temp)

dim(temp)

grayscale_table <- table.paint(temp, cleg  = 0, clabel.row = .5, clabel.col = .5) + scale_color_viridis()

"Construye un objeto de tipo phylo con la función nj del paquete ape a partir de la matriz de distancia que
obtuviste en el paso anterior. Muestra el código que empleaste para hacerlo"

virus_tree <- nj(dist_matrix)
class(virus_tree)

