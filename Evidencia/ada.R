## ---- Settings ----
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.align = "center",
  fig.width = 7,
  fig.height = 5,
  fig.path = "images/"
)


#Cargar librerias
library(viridis)
library(Biostrings)
library(DECIPHER)
library(ade4)
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(ggplot2)
#added packages
# library(ggthemes)
# library(reshape2)


"Obtén las secuencias de los genomas de las variantes de SARS-CoV-2 y otros virus que tu decidas (8-10 genomas en total) 
desde el NCBI Links to an external site. o el buscador de virus del NCBI Links to an external site."

variantes_coronavirus <- c("MN908947", "MT121215", "NC_019843", 
                           "MT121216", "NC_004718", "NC_019843",
                           "AY567487", "NC_006213", "NC_006577",
                           "NC_004718")


#raw GenBank
variant_secuences <- read.GenBank(variantes_coronavirus, as.character = FALSE)

#Transform to Biostrings
write.FASTA(variant_secuences, "variant_sequences.fasta")
fasta_variant_secuences <- readDNAStringSet("variant_sequences.fasta", format = "fasta")


"Calcula la longitud de las secuencias de cada variante.
Muestra el resultado impreso en consola y el código que utilizaste para 
realizar este punto."

variant_sizes <- lengths(variant_secuences) #Get lenghts

#Display in table form
my_dataframe <- data.frame("Variante" = variantes_coronavirus,
                        "Longitud" = variant_sizes)
table_data <- data.frame(variantes_coronavirus, variant_sizes)
print(table_data)

"Genera una gráfica en la que se visualice la longitud de los 10 genomas 
virales que analizaste. Muestra el código que utilizaste para realizar 
este punto."

names(variant_sizes)
barplot(variant_sizes, col = "blue", xlab = "Longitud", ylab = "Variantes", width = 0.5, horiz = TRUE)

"Agrega una interpretación escrita, desde el 
punto de vista biológico, para esta gráfica (de 6 a 12 renglones)."

#Explicacion
"Se puede observar principalmente el hecho de que todas las secuencias tienen una longitud
bastante similar, esto se debe a que las variantes comparten una cepa comun, pero aun asi varian de manera
sigmificativa devido las mutaciones generadas en la replicacion del virus.
A pesar de que en la mayoria de casos las mutaciones de este tipo de virus son aleatorias
gracias a la escala en la que este tipo de organismos tienden a replicarse ademas de 
factores como la selccion natural tienden a fomentar mutaciones que provocan que el virus se propage
de manera mas eficiente.

En estos casos incluso diferencias muy pequenas pueden tener un impacto sigmificativo
en la capacidad del virus para propagarse y imapctar al ser humano y otras especies."

#Paso 5
"Calcula la composición de nucleótidos a cada una de las variantes del virus.
Muestra el resultado impreso en consola y el código que utilizaste para realizar 
este punto."
 
# Get nucleotide composition for each sequence using lapply and alphabetFrequency
nt_comp_list <- lapply(fasta_variant_secuences, function(x) alphabetFrequency(x)[intersect(c("A", "C", "G", "T"), names(alphabetFrequency(x)))])
# Convert the list to a matrix
nt_comp <- do.call(rbind, nt_comp_list)

# Get matrix for percentages
nt_pct <- nt_comp / variant_sizes * 100

# Round the matrix percentages
nt_pct <- round(nt_pct, 2)

# Print the matrix
cat("\nNucleotids by number:\n")
print(nt_comp)
cat("\nNucleotids by percentage:\n")
print(nt_pct)

"Genera una gráfica en la que se visualice la composición de nucleótidos de los 
10 genomas virales que analizaste. Muestra el código que utilizaste para realizar 
este punto."


#Convert to dataframe
nt_pct_df <- as.data.frame(nt_pct)
nt_pct_df$seqid <- rownames(nt_pct_df)

#Melt for display
nt_pct_long <- melt(nt_pct_df, id.vars = "seqid", variable.name = "neucleotide", value.name = "percentage")

#Create bidimensional heatmap
X11()
p <- ggplot(nt_pct_long, aes(x = neucleotide, y = seqid, fill = percentage)) +
  geom_tile() +
  scale_fill_gradientn(colors = viridis::viridis(10)) +
  theme_minimal()


print(p)


"Agrega una interpretación escrita, desde el punto de vista biológico, 
para esta gráfica (de 6 a 12 renglones)."

"Esta gráfica es una representación visual de la variación en la composición de
 nucleótidos de diferentes variantes del virus del COVID-19. Los nucleótidos son 
 los bloques de construcción que forman el material genético del virus, el ARN.

Cada variante del virus se muestra en la gráfica como una columna, mientras que
las filas representan posiciones específicas en la secuencia de ARN del virus. 
El calor en cada celda indica la frecuencia relativa de un nucleótido específico 
en esa posición de la secuencia.

Las variaciones en la composición de nucleótidos pueden tener importantes 
consecuencias biológicas, como la capacidad del virus para infectar células, 
evadir el sistema inmunológico y replicarse. Por lo tanto, esta gráfica puede 
ayudar a los investigadores a identificar diferencias significativas entre 
variantes del virus y comprender mejor cómo se propagan y afectan a los pacientes.

En resumen, esta gráfica proporciona una representación visual de la variación 
en la secuencia de ARN del COVID-19 y puede ser útil para los investigadores 
que estudian la biología y la evolución del virus."

#Siguiente paso
"Calcula el %GC de cada variante. Muestra el resultado impreso en consola 
y el código que utilizaste para realizar este punto."


gc_frecuences <- c()
for (i in seq_along(variantes_coronavirus)){
  gc_frecuences <- append(gc_frecuences, GC.content(c(variant_secuences[i])))
}

barplot(gc_frecuences,
        names.arg = variantes_coronavirus,
        xlab = "Variant",
        ylab = "Percentage",
        main = "%GC",
        col = "steelblue"
)

"Agrega una interpretación escrita, desde el punto de vista biológico, 
para esta gráfica (de 6 a 12 renglones)."

"Esta gráfica muestra la comparación del contenido de guanina y citosina (%GC) 
en las secuencias de nucleótidos de diferentes variantes del virus COVID-19. 
El %GC es una medida de la proporción de guanina y citosina en relación a adenina 
y timina en la secuencia de nucleótidos de un organismo.

En general, el %GC es un factor importante en la estructura y función del material 
genético. En esta gráfica, se puede observar que el %GC varía ligeramente entre cada 
variante del virus, lo que sugiere que cada variante tiene una secuencia única que 
puede influir en sus características biológicas.

Aunque los cambios en el %GC pueden parecer pequeños, pueden tener un gran impacto 
en la estructura y función de las proteínas codificadas por los genes del virus. 
Por ejemplo, los cambios en el %GC pueden afectar la estabilidad de la estructura 
secundaria del ARN viral, lo que a su vez podría influir en la capacidad del virus 
para replicarse y propagarse.

En resumen, esta gráfica proporciona información sobre la variación en el contenido 
de GC en diferentes variantes del virus COVID-19. El análisis de esta variación puede 
ser útil para los investigadores que estudian la biología y la evolución del virus, 
y para el desarrollo de tratamientos y vacunas específicos para cada variante.

Crea una tabla (marco de datos) en la que se inluya el nombre común de cada virus, 
el ID (identificador en base de datos) de cada genoma, la longitud y el contenido GC 
de cada genoma. Muestra el resultado impreso en consola y el código que utilizaste 
para realizar este punto."

nombres_virus <- c("SARS-CoV-2", "RaTG13", "RmYN02", "Pangolin-CoV", "HKU2", "MERS-CoV", "HCoV-NL63", "HCoV-OC43", "HCoV-HKU1", "SARS-CoV")


punto_final <- data.frame(Nombres = nombres_virus, ID = variantes_coronavirus, Longitud = variant_sizes, Contenido_GC = round(gc_frecuences * 100, 2))
print("\n")
print(punto_final)