###############################################################################
#                                    TFM 
###############################################################################
#           Script 01 B: Reconstrucción del objeto phyloseq
###############################################################################

# Objetivo:
# Reconstrucción del objeto phyloseq utilizando las secuencias reales de los 
# ASVs obtenidas en Galaxy
#
# Entrada:
#   - phyloseqObject.phyloseq
#   - Galaxy-removeBimera.dada2_sequencetable
#
# Salida:
#   - phyloseqObject_sequences.rds

###############################################################################

rm(list = ls())

###############################################################################
# Configuración de rutas
###############################################################################

# El script funciona tanto si R se inicia en la raíz del repositorio como si
# se inicia en la carpeta inmediatamente superior a "VaginalMicrobiome".
project_dir <- if (dir.exists("cleaned")) {
  "."
} else if (dir.exists(file.path("VaginalMicrobiome", "cleaned"))) {
  "VaginalMicrobiome"
} else {
  stop(
    "ERROR: no se encontró la carpeta 'cleaned'. ",
    "Ejecute el script desde la raíz del repositorio VaginalMicrobiome ",
    "o desde su carpeta padre."
  )
}

project_path <- function(...) file.path(project_dir, ...)

###############################################################################
# 1. Librerías
###############################################################################

library(phyloseq)
library(readr)
library(dplyr)

###############################################################
# 2. Cargar phyloseq exportado desde Galaxy
###############################################################

ps <- readRDS(
  project_path("cleaned", "phyloseqObject.phyloseq"))

###############################################################
# 3. Explorar tabla removeBimeraDenovo
###############################################################

seqtab <-
  read.delim(
    project_path("cleaned", "Galaxy-removeBimera.dada2_sequencetable"),
    check.names = FALSE,
    stringsAsFactors = FALSE)

# Explorar la tabla
dim(seqtab) # [1] 8298 filas (ASVs = sec de nucleotidos) x 2368 columnas (muestras = SRRs)
head(seqtab[, 1:5])

###############################################################
# 4. Extraer las secuencias ASV
###############################################################

asv_sequences <- seqtab[, 1]
head(asv_sequences)

# Comprobar longitud de las secuencias
table(nchar(asv_sequences))
summary(nchar(asv_sequences))
#    Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
#   240.0   253.0   253.0   253.4   253.0   425.0 

# Comparar número de ASVs con phyloseq
cat("ASVs en phyloseq:", ntaxa(ps), "\n") # ASVs en phyloseq: 8298 
cat("Secuencias en tabla Galaxy:", length(asv_sequences), "\n") 
# Secuencias en tabla Galaxy: 8298 

###############################################################
# 5. Extraer únicamente la matriz de conteos
###############################################################

counts_galaxy <- as.matrix(seqtab[, -1])

cat("Dimensiones de counts_galaxy:", dim(counts_galaxy), "\n")
# Dimensiones de counts_galaxy: 8298 2367

###############################################################
# 6. Extraer tabla de abundancias desde phyloseq
###############################################################

otu_ps <- as(otu_table(ps), "matrix")

# Asegurar que las filas sean taxa/ASVs
if (!taxa_are_rows(ps)) {otu_ps <- t(otu_ps)}

cat("Dimensiones de otu_ps:", dim(otu_ps), "\n") # Dimensiones de otu_ps: 8298 2367 
cat("Dimensiones de counts_galaxy:", dim(counts_galaxy), "\n") # Dimensiones de counts_galaxy: 8298 2367 

###############################################################
# 7. Verificar dimensiones
###############################################################

if (!identical(dim(otu_ps), dim(counts_galaxy))) {
  stop("ERROR: otu_ps y counts_galaxy tienen dimensiones diferentes.")}

cat("OK: las dimensiones coinciden.\n") # OK: las dimensiones coinciden.

###############################################################
# 8. Verificar muestras
###############################################################

if (!identical(colnames(otu_ps), colnames(counts_galaxy))) {
  stop("ERROR: el orden/nombre de las muestras no coincide entre ",
    "phyloseq y la tabla Galaxy.")}

cat("OK: las muestras coinciden y están en el mismo orden.\n") # OK: las muestras coinciden y están en el mismo orden.

###############################################################
# 9. Verificar que las abundancias sean idénticas
###############################################################

if (all(otu_ps == counts_galaxy)) {
  cat("OK: las abundancias son idénticas.\n\n")
} else {cat("ERROR: las abundancias NO son idénticas.\n\n")
  # Número de diferencias
    n_differences <- sum(otu_ps != counts_galaxy)
    cat("Número de diferencias:", n_differences, "\n\n")
  # Ubicación de las primeras diferencias
    differences <- which(otu_ps != counts_galaxy, arr.ind = TRUE)
    print(differences[1:min(20, nrow(differences)), , drop = FALSE])
  # Detener el script
    stop("Las matrices de abundancias no son idénticas. ",
      "Revisar antes de continuar.")}

###############################################################
# 10. Verificar número de secuencias
###############################################################

if (length(asv_sequences) != nrow(otu_ps)) {
  stop(
    "ERROR: el número de secuencias en asv_sequences (",length(asv_sequences),")
    no coincide con el número de ASVs (",nrow(otu_ps),").")}

cat("OK: hay una secuencia para cada ASV.\n") # OK: hay una secuencia para cada ASV.

###############################################################
# 11. Reemplazar nombres de taxa por las secuencias ASV
###############################################################

taxa_names(ps) <- asv_sequences

cat("\nPrimeros nombres de taxa después del cambio:\n")

print(head(taxa_names(ps)))

###############################################################
# 12. Verificar que no haya secuencias duplicadas
###############################################################

n_unique <- length(unique(taxa_names(ps)))
n_total  <- length(taxa_names(ps))

cat("\nNúmero total de ASVs:", n_total, "\n") # Número total de ASVs: 8298 
cat("Número de secuencias únicas:", n_unique, "\n") # Número de secuencias únicas: 8298

if (n_unique != n_total) {
  stop("ERROR: existen secuencias ASV duplicadas.")
} else {cat("OK: no hay secuencias duplicadas.\n")} # OK: no hay secuencias duplicadas.

###############################################################################
# 13. Reconstruir tabla de abundancias con las secuencias como nombres de fila
###############################################################################

counts_galaxy <- as(otu_table(ps), "matrix")

if (!taxa_are_rows(ps)) {counts_galaxy <- t(counts_galaxy)}

###############################################################
# 14. Guardar tabla de secuencias + abundancias
###############################################################

write.table(
  data.frame(
    Sequence = rownames(counts_galaxy),
    counts_galaxy,
    check.names = FALSE),
  file = project_path("results", "Script01B_ASVtable_Galaxy.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

cat("\nTabla ASV exportada correctamente.\n")

###############################################################
# 15. Guardar nuevo objeto phyloseq
###############################################################

saveRDS(ps,
  project_path("cleaned", "phyloseqObject_sequences.rds"))

cat("Objeto phyloseq guardado correctamente.\n")

###############################################################
# 16. Resumen
###############################################################

cat("\n")
cat("=========================================\n")
cat("SCRIPT 01B FINALIZADO\n")
cat("=========================================\n\n")

cat("ASVs:", ntaxa(ps), "\n") # ASVs: 8298
cat("Muestras:", nsamples(ps), "\n") # Muestras: 2367 
cat("Identificadores reemplazados:", ntaxa(ps), "\n") # Identificadores reemplazados: 8298 

# Ahora todos los scripts utilizarán las secuencias reales de los ASVs.
