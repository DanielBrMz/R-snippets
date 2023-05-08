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


# Se eligen nucleotidos determinados
first_150_nucleotides <- subseq(coronavirus_string_set, start = 1, end = 150)

names(first_150_nucleotides) <- names_sars_cov_2

nucleotides_500_650 <-  subseq(coronavirus_string_set, start = 500, end = 650)

names(nucleotides_500_650) <- names_sars_cov_2

#Se muestran en el navegador
BrowseSeqs(first_150_nucleotides)
BrowseSeqs(nucleotides_500_650)

"Agrega una interpretación escrita, desde el punto de vista biológico, para esta gráfica (de 6 a 12 renglones)."


"Las funciones BrowseSeqs permiten visualizar secuencias de virus relacionados con
SARS-CoV-2 encontrados en animales. Al observar los primeros 150 nucleótidos, podemos
identificar regiones conservadas y estructuras importantes. La visualización de nucleótidos
500-650 destaca regiones codificantes y regulatorias. Estos análisis ayudan a comprender 
la relación genética y funcional entre los virus y su posible implicación en la transmisión
animal-humano. La investigación se centra en examinar secuencias de coronavirus en
animales para explicar la posibilidad de transmisión. Interpretar las secuencias en el
navegador proporciona información clave para identificar similitudes y diferencias entre los
genomas, lo que puede ser relevante para comprender la transmisión y la adaptación de
virus entre especies."

"Genera una matriz de distancia a partir de los genomas alineados. Crea una tabla en escala de grises en la que observes de manera 
visual el resultado de la matriz de distancia e inclúyela en tu reporte. Muestra el código empleado para obtener lo anterior e 
incluye la matriz de distancia y la tabla que obtuviste."

writeXStringSet(aligned_coronavirus_string_set, "coronavirus_sequences_aligned.fasta", format="fasta")

alignment_coronavirus <- read.alignment(file = "coronavirus_sequences_aligned.fasta", format = "fasta")

distance_coronavirus <- dist.alignment(alignment_coronavirus, matrix = "similarity")

temp <- as.data.frame(as.matrix(distance_coronavirus))

rownames(temp) <- names_sars_cov_2
colnames(temp) <- names_sars_cov_2

grayscale_table <- table.paint(temp, cleg  = 0, clabel.row = .5, clabel.col = .5) + scale_color_viridis()


"Agrega una interpretación escrita, desde el punto de vista biológico, para esta gráfica (de 6 a 12 renglones)."

"La gráfica generada muestra una representación visual de la matriz de distancias entre los
genomas virales alineados. Cada casilla de la tabla representa la similitud entre dos 
genomas, donde los colores más oscuros indican una mayor similitud y los colores más
claros indican una menor similitud.
Al observar la gráfica, podemos identificar patrones de similitud y diferencias entre los
genomas. Las regiones con colores oscuros indican una alta similitud, lo que sugiere que
esos genomas comparten secuencias conservadas. Por otro lado, las regiones con colores
claros indican diferencias entre los genomas, lo que puede corresponder a regiones únicas 
o variantes genéticas.
El análisis de la matriz de distancias nos permite identificar relaciones evolutivas entre los
genomas virales y comprender su diversidad genética. Esto es relevante para estudiar la
transmisión y adaptación de los virus entre especies, así como para comprender la posible
emergencia de nuevas variantes. Una comprensión más detallada de la matriz de distancias
requeriría análisis adicionales y comparaciones con genomas de referencia."


"Construye un árbol filogenético a partir de la matriz de distancia obtenida e incluye en el árbol los números de accesión de 
los genomas utilizados, sus nombres comunes o cualquier otra leyenda que te permita indicar la ubicación de cada virus en el árbol.
Muestra el código empleado para realizar lo anterior e incluye la imagen del árbol filogenético."


phylogenetic_tree_coronavirus <- nj(distance_coronavirus)

phylogenetic_tree_coronavirus <- ladderize(phylogenetic_tree_coronavirus)

original_tip_labels <- rownames(phylogenetic_tree_coronavirus$tip.label)

tip_labels <- setNames(names_sars_cov_2, original_tip_labels)


modified_phylogenetic_tree <- phylogenetic_tree_coronavirus
modified_phylogenetic_tree$tip.label <- tip_labels


plot.phylo(modified_phylogenetic_tree, show.tip.label = TRUE)


plot_virus_phylogeny <- ggtree(modified_phylogenetic_tree) + geom_tiplab() + ggtitle("Análisis filogenético de genomas de variantes del SARS-CoV-2 encontrados en animales")

plot_virus_phylogeny


"Agrega una interpretación escrita, desde el punto de vista biológico, para esta gráfica"

"En el árbol filogenético, cada rama representa una conexión evolutiva entre las variantes del 
virus. Las ramas más cercanas entre sí indican una mayor similitud genética, lo que sugiere 
una relación más estrecha en términos de descendencia común. Por otro lado, las ramas 
más separadas representan una mayor divergencia genética y pueden indicar eventos de 
evolución independiente.
Las etiquetas en los extremos de las ramas representan las variantes del virus y están 
asociadas con los nombres comunes asignados a cada una de ellas. Estos nombres 
proporcionan una forma más fácil de identificar y referirse a las diferentes variantes del 
SARS-CoV-2 encontradas en animales.
Al observar el árbol filogenético, podemos identificar agrupaciones o clados de variantes del
virus que comparten una mayor similitud genética. Estos clados pueden indicar la existencia 
de linajes virales específicos que se han propagado en poblaciones animales y pueden tener 
implicaciones en la transmisión y la adaptación del virus."
