###############################################################################
#                                    TFM
###############################################################################
#                       Script 04: Filtrado taxonómico
###############################################################################

# Objetivo:
#   Eliminar ASVs identificados explícitamente como no bacterianos antes de
#   realizar los análisis ecológicos.

# Se eliminan:
#   - Eukaryota
#   - Archaea
#   - secuencias asignadas a mitocondrias
#   - secuencias asignadas a cloroplastos

# Los ASVs cuya clasificación taxonómica es incompleta o NA se conservan
# de forma conservadora, salvo que exista evidencia explícita de que pertenecen
# a alguno de los grupos anteriores.

# Entrada:
#   - ps_decontam.rds

# Salidas:
#   - ps_taxFiltered.rds
#   - Taxonomic_Filtering_Report.txt
#   - Taxonomic_Filtering_Removed_ASVs.csv
#   - Taxonomic_Filtering_Summary.csv

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
# 1. Cargar librerías
###############################################################################

library(phyloseq)
library(dplyr)
library(readr)

dir.create(project_path("results"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("cleaned"), recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2. Cargar objeto phyloseq
###############################################################################

ps <- readRDS(project_path("cleaned", "ps_decontam.rds"))

cat("\nObjeto cargado\n")
cat("Número de muestras:", nsamples(ps), "\n") # 2179
cat("Número de ASVs:", ntaxa(ps), "\n") # 6539

###############################################################################
# 3. Extraer clasificación taxonómica antes del filtrado
###############################################################################

tax_before <- as.data.frame(
  tax_table(ps),
  stringsAsFactors = FALSE)

cat("\nDistribución inicial por Reino:\n")
print(
  table(
    tax_before$Kingdom,
    useNA = "ifany"))
#   Archaea  Bacteria Eukaryota      <NA> 
#        6      6391       123        19 

cat("\nDistribución inicial por Filo:\n")
print(
  table(
    tax_before$Phylum,
    useNA = "ifany"))
# ...              <NA> 201 

###############################################################################
# 4. Identificar taxones a eliminar
###############################################################################

# Se buscan ASVs asignados explícitamente a:
#   - mitocondrias
#   - cloroplastos
#   - Eukaryota
#   - Archaea

mitochondria <- grepl(
  "mitochond",
  tax_before$Family,
  ignore.case = TRUE)

chloroplast <- grepl(
  "chloroplast",
  tax_before$Order,
  ignore.case = TRUE)

eukaryota <- grepl(
  "eukary",
  tax_before$Kingdom,
  ignore.case = TRUE)

archaea <- grepl(
  "archaea",
  tax_before$Kingdom,
  ignore.case = TRUE)

###############################################################################
# 5. Reemplazar valores NA de los vectores lógicos por FALSE
###############################################################################

# Si un ASV tiene NA en alguna categoría taxonómica, no se lo elimina
# automáticamente. Solo se elimina cuando existe una asignación explícita
# a alguno de los grupos no bacterianos definidos arriba.

mitochondria[is.na(mitochondria)] <- FALSE
chloroplast[is.na(chloroplast)] <- FALSE
eukaryota[is.na(eukaryota)] <- FALSE
archaea[is.na(archaea)] <- FALSE

###############################################################################
# 6. Identificar ASVs a eliminar y ASVs a conservar
###############################################################################

remove <- (
  mitochondria |
    chloroplast |
    eukaryota |
    archaea)

keep <- !remove

cat("\n========================================\n")
cat("TAXONES IDENTIFICADOS PARA ELIMINAR\n")
cat("========================================\n")

cat("\nMitocondrias:", sum(mitochondria), "\n") # Mitocondrias: 51 
cat("Cloroplastos:", sum(chloroplast), "\n") # Cloroplastos: 22 
cat("Eucariotas:", sum(eukaryota), "\n") # Eucariotas: 123 
cat("Arqueas:", sum(archaea), "\n") # Arqueas: 6 
cat("ASVs únicos a eliminar:", sum(remove), "\n")
# ASVs únicos a eliminar: 202 

###############################################################################
# 7. Guardar tabla de ASVs eliminados
###############################################################################

# Esta tabla permite documentar qué ASVs fueron eliminados y por qué motivo.

removed_taxa <- tax_before[remove, , drop = FALSE]

removed_table <- data.frame(
  ASV = rownames(removed_taxa),
  removed_taxa,
  Mitochondria = mitochondria[remove],
  Chloroplast = chloroplast[remove],
  Eukaryota = eukaryota[remove],
  Archaea = archaea[remove],
  stringsAsFactors = FALSE)

write_csv(
  removed_table,
  project_path("results", "Taxonomic_Filtering_Removed_ASVs.csv"))

###############################################################################
# 8. Aplicar filtrado taxonómico
###############################################################################

ps_filtered <- prune_taxa(keep, ps)

removed_asvs <- ntaxa(ps) - ntaxa(ps_filtered)

cat("\nASVs iniciales:", ntaxa(ps), "\n")
# ASVs iniciales: 6539 
cat("ASVs finales:", ntaxa(ps_filtered), "\n")
# ASVs finales: 6337
cat("ASVs eliminados:", removed_asvs, "\n")
# ASVs eliminados: 202

###############################################################################
# 9. Comprobar si alguna muestra quedó con 0 reads
###############################################################################

# Después de eliminar ASVs no bacterianos, se comprueba si alguna muestra
# quedó completamente vacía.

# Si una muestra tiene 0 reads, se elimina del objeto porque no puede
# utilizarse en análisis posteriores.

zero_samples <- sample_names(ps_filtered)[
  sample_sums(ps_filtered) == 0]

cat("\nMuestras con 0 reads después del filtrado:",
  length(zero_samples),
  "\n")
# Muestras con 0 reads después del filtrado: 0

if (length(zero_samples) > 0) {
  cat("Estas muestras serán eliminadas del objeto filtrado.\n")
  ps_filtered <- prune_samples(
    sample_sums(ps_filtered) > 0,
    ps_filtered)}

###############################################################################
# 10. Extraer clasificación taxonómica después del filtrado
###############################################################################

tax_after <- as.data.frame(
  tax_table(ps_filtered),
  stringsAsFactors = FALSE)

cat("\nDistribución final por Reino:\n")
print(
  table(
    tax_after$Kingdom,
    useNA = "ifany"))
# Bacteria     <NA> 
# 6318       19 

cat("\nDistribución final por Filo:\n")
print(
  table(
    tax_after$Phylum,
    useNA = "ifany"))
# ... <NA> 87 

###############################################################################
# 11. Comprobar que no hayan quedado taxones explícitamente no bacterianos
###############################################################################

remaining_non_bacterial <- (
  grepl(
    "eukary",
    tax_after$Kingdom,
    ignore.case = TRUE) |
    grepl(
      "archaea",
      tax_after$Kingdom,
      ignore.case = TRUE) |
    grepl(
      "mitochond",
      tax_after$Family,
      ignore.case = TRUE) |
    grepl(
      "chloroplast",
      tax_after$Order,
      ignore.case = TRUE))

remaining_non_bacterial[is.na(remaining_non_bacterial)] <- FALSE

cat("\nASVs explícitamente no bacterianos restantes:",
  sum(remaining_non_bacterial),
  "\n") # 0

if (sum(remaining_non_bacterial) > 0) {
  stop("ERROR: quedaron ASVs explícitamente no bacterianos después del filtrado.")}

###############################################################################
# 12. Crear tabla resumen del filtrado
###############################################################################

filter_summary <- data.frame(
  Metric = c(
    "Muestras antes del filtrado",
    "Muestras después del filtrado",
    "ASVs antes del filtrado",
    "ASVs después del filtrado",
    "ASVs eliminados",
    "Mitocondrias detectadas",
    "Cloroplastos detectados",
    "Eucariotas detectados",
    "Arqueas detectadas",
    "Muestras con 0 reads eliminadas"),
  Value = c(
    nsamples(ps),
    nsamples(ps_filtered),
    ntaxa(ps),
    ntaxa(ps_filtered),
    removed_asvs,
    sum(mitochondria),
    sum(chloroplast),
    sum(eukaryota),
    sum(archaea),
    length(zero_samples)))

print(filter_summary)
#                             Metric Value
# 1      Muestras antes del filtrado  2179
# 2    Muestras después del filtrado  2179
# 3          ASVs antes del filtrado  6539
# 4        ASVs después del filtrado  6337
# 5                  ASVs eliminados   202
# 6          Mitocondrias detectadas    51
# 7          Cloroplastos detectados    22
# 8            Eucariotas detectados   123
# 9               Arqueas detectadas     6
# 10 Muestras con 0 reads eliminadas     0

write_csv(
  filter_summary,
  project_path("results", "Taxonomic_Filtering_Summary.csv"))

###############################################################################
# 13. Guardar objeto phyloseq filtrado
###############################################################################

saveRDS(
  ps_filtered,
  project_path("cleaned", "ps_taxFiltered.rds"))

###############################################################################
# 14. Guardar informe del filtrado taxonómico
###############################################################################

sink(project_path("results", "Taxonomic_Filtering_Report.txt"))

cat("========================================\n")
cat("INFORME DE FILTRADO TAXONOMICO\n")
cat("========================================\n\n")

cat("Numero de muestras antes del filtrado:\n")
cat(nsamples(ps), "\n\n")

cat("Numero de muestras despues del filtrado:\n")
cat(nsamples(ps_filtered), "\n\n")

cat("Numero de ASVs antes del filtrado:\n")
cat(ntaxa(ps), "\n\n")

cat("Numero de ASVs despues del filtrado:\n")
cat(ntaxa(ps_filtered), "\n\n")

cat("Numero de ASVs eliminados:\n")
cat(removed_asvs, "\n\n")

cat("Mitocondrias detectadas:\n")
cat(sum(mitochondria), "\n\n")

cat("Cloroplastos detectados:\n")
cat(sum(chloroplast), "\n\n")

cat("Eucariotas detectados:\n")
cat(sum(eukaryota), "\n\n")

cat("Arqueas detectadas:\n")
cat(sum(archaea), "\n\n")

cat("Muestras con 0 reads eliminadas:\n")
cat(length(zero_samples), "\n\n")

cat("Distribucion por Reino ANTES del filtrado:\n")
print(
  table(
    tax_before$Kingdom,
    useNA = "ifany"))

cat("\nDistribucion por Reino DESPUES del filtrado:\n")
print(
  table(
    tax_after$Kingdom,
    useNA = "ifany"))

cat("\nDistribucion por Filo ANTES del filtrado:\n")
print(
  table(
    tax_before$Phylum,
    useNA = "ifany"))

cat("\nDistribucion por Filo DESPUES del filtrado:\n")
print(
  table(
    tax_after$Phylum,
    useNA = "ifany"))

cat("\nInterpretacion:\n")
cat(
  paste(
    "El filtrado elimino los ASVs identificados explicitamente como",
    "no bacterianos, incluyendo Eukaryota, Archaea, mitocondrias y",
    "cloroplastos. Los ASVs con clasificacion taxonomica incompleta o",
    "ausente se conservaron, salvo que existiera evidencia explicita de",
    "pertenencia a alguno de los grupos eliminados."),
  "\n")

sink()

###############################################################################
# 15. Fin del script
###############################################################################

cat("\n========================================\n")
cat("SCRIPT 04 FINALIZADO\n")
cat("========================================\n")

cat("\nObjeto guardado en:\n",
  project_path("cleaned", "ps_taxFiltered.rds"), "\n")
# Objeto guardado en:
# VaginalMicrobiome/cleaned/ps_taxFiltered.rds

cat("\nScript 04 finalizado correctamente.\n")

