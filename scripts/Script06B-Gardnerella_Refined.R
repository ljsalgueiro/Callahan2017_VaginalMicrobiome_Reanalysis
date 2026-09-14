###############################################################################
#                                    TFM
###############################################################################
#         Script 06: Refinamiento taxonómico de Gardnerella
###############################################################################

# Objetivo:
#   1. Recuperar las secuencias G1, G2 y G3 del objeto processed.rda utilizado
#      en el análisis original de Callahan et al. (2017).
#   2. Identificar los ASVs de Gardnerella de nuestra réplica compatibles
#      exactamente con dichas variantes.
#   3. Clasificar los ASVs de Gardnerella como G1, G2, G3 o G_other.
#   4. Incorporar GardnerellaVariant a una copia del objeto phyloseq.
#
# Fundamento:
#   En el código original de Callahan:
#     gg <- which(tax[, "Genus"] %in% "Gardnerella")
#     G1 <- gg[[1]]
#     G2 <- gg[[2]]
#     G3 <- gg[[3]]

#   En processed.rda se recuperaron 9 variantes Gardnerella:
#     7, 8, 12, 39, 40, 80, 664, 1017, 1236

#   Por tanto:
#     G1 = rownames(tax)[7]
#     G2 = rownames(tax)[8]
#     G3 = rownames(tax)[12]

#   Las referencias de Callahan tienen 235 nt y nuestros ASVs 253 nt.
#   Se busca una coincidencia exacta de los 235 nt dentro de cada ASV.

#   IMPORTANTE:
#   Más de un ASV de 253 nt puede corresponder a la misma variante de Callahan
#   si difieren únicamente en posiciones externas a sus 235 nt de referencia.

# Entradas:
#   - ps_species_byBlast.rds
#   - processed.rda

# Salidas:
#   - Callahan_Gardnerella_G1_G2_G3_reference.csv
#   - Gardnerella_G1_G2_G3_matches.csv
#   - Script06_Gardnerella_ByHand.csv
#   - Gardnerella_ByHand.rds
#   - Script06_DominanciaGard.png
#   - Script06_Gardnerella_variant_summary.csv
#   - ps_gard_refined.rds

###############################################################################

rm(list = ls())

# Permite ejecutar el script desde la raíz del repositorio o desde su carpeta
# inmediatamente superior.
project_dir <- if (dir.exists("VaginalMicrobiome")) {
  "VaginalMicrobiome"
} else {
  "."
}


###############################################################################
# 1. Librerías
###############################################################################

library(phyloseq)
library(dplyr)
library(ggplot2)
library(readr)

###############################################################################
# 2. Directorios
###############################################################################
results_dir <- file.path(project_dir, "results", "Gardnerella")
images_dir <- file.path(project_dir, "images", "Gardnerella")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 3. Cargar objeto phyloseq
###############################################################################

ps <- readRDS(file.path(project_dir, "cleaned", "ps_lacto_refined.rds"))

cat("\n========================================\n")
cat("OBJETO PHYLOSEQ\n")
cat("========================================\n")
print(ps)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6337 taxa and 2179 samples ]
# sample_data() Sample Data:       [ 2179 samples by 54 sample variables ]
# tax_table()   Taxonomy Table:    [ 6337 taxa by 10 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 6337 reference sequences ]

###############################################################################
# 4. Recuperar secuencias G1, G2 y G3 originales de Callahan
###############################################################################

if (!file.exists(
  file.path(project_dir, "Callahan2017", "RepRefine_Scripts", "input", "processed.rda"))) {
  stop("ERROR: no se encontró processed.rda.")}

# Se carga en un entorno independiente para no sobrescribir objetos propios.
callahan_env <- new.env()
load(file.path(project_dir, "Callahan2017", "RepRefine_Scripts", "input", "processed.rda"), 
     envir = callahan_env)

if (!exists("tax", envir = callahan_env, inherits = FALSE)) {
  stop("ERROR: processed.rda no contiene el objeto 'tax'.")}

tax_callahan <- get("tax", envir = callahan_env)

# Mismo procedimiento utilizado en el código original de Callahan.
gg_callahan <- which(tax_callahan[, "Genus"] %in% "Gardnerella")

cat("\nVariantes Gardnerella en processed.rda:", length(gg_callahan), "\n")
cat("Índices:\n")
print(gg_callahan)
# [1]    7    8   12   39   40   80  664 1017 1236

if (length(gg_callahan) < 3) {
  stop("ERROR: processed.rda contiene menos de tres variantes Gardnerella.")}

# Callahan identifica como G1, G2 y G3 las tres primeras variantes Gardnerella.
G_reference <- c(
  G1 = rownames(tax_callahan)[gg_callahan[[1]]],
  G2 = rownames(tax_callahan)[gg_callahan[[2]]],
  G3 = rownames(tax_callahan)[gg_callahan[[3]]])

cat("\nSecuencias G1/G2/G3 recuperadas:\n")
print(G_reference)

reference_lengths <- nchar(G_reference)

cat("\nLongitud de las referencias:\n")
print(reference_lengths)
#  G1  G2  G3 
# 235 235 235 

if (!all(reference_lengths == 235)) {
  warning("Las referencias G1/G2/G3 no tienen todas 235 nt.")}

###############################################################################
# 5. Guardar referencias originales
###############################################################################

callahan_reference_table <- data.frame(
  GardnerellaVariant = names(G_reference),
  Callahan_index = gg_callahan[1:3],
  Sequence = unname(G_reference),
  Length = unname(reference_lengths),
  stringsAsFactors = FALSE)

write_csv(
  callahan_reference_table,
  file.path(project_dir, "results", "Gardnerella", "Callahan_Gardnerella_G1_G2_G3_reference.csv"))

###############################################################################
# 6. Seleccionar ASVs Gardnerella de nuestra réplica
###############################################################################

tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

gard_idx <- tax %>% filter(Genus == "Gardnerella")

n_gard <- nrow(gard_idx)

cat("\nASVs Gardnerella en nuestra réplica:", n_gard, "\n")

if (n_gard == 0) stop("ERROR: no se encontraron ASVs de Gardnerella.")

cat("\nAsignación SILVA:\n")
print(table(gard_idx$Species, useNA = "always"))
# vaginalis      <NA> 
#        7        30 

# Los nombres de taxa son las secuencias ASV.
gard_sequences <- rownames(gard_idx)

stopifnot(!anyNA(gard_sequences))
stopifnot(all(nchar(gard_sequences) > 0))

cat("\nLongitud de los ASVs Gardnerella:\n")
print(table(nchar(gard_sequences)))
# 253

###############################################################################
# 7. Buscar coincidencias exactas con G1/G2/G3
###############################################################################

# Las referencias Callahan tienen 235 nt y nuestros ASVs 253 nt.
# Por ello se busca la secuencia completa de 235 nt dentro del ASV.

find_exact_reference <- function(reference_sequence, sequences) {
  sequences[
    grepl(reference_sequence, sequences, fixed = TRUE)]}

G1_match <- find_exact_reference(G_reference["G1"], gard_sequences)
G2_match <- find_exact_reference(G_reference["G2"], gard_sequences)
G3_match <- find_exact_reference(G_reference["G3"], gard_sequences)

cat("\n========================================\n")
cat("COINCIDENCIAS EXACTAS\n")
cat("========================================\n")

cat("\nG1:", length(G1_match), "ASV(s)\n")
print(G1_match) # G1: 2 ASV(s)

cat("\nG2:", length(G2_match), "ASV(s)\n")
print(G2_match) # G2: 2 ASV(s)

cat("\nG3:", length(G3_match), "ASV(s)\n")
print(G3_match) # G3: 1 ASV(s)

###############################################################################
# 8. Verificaciones
###############################################################################

# No se exige una única coincidencia por variante.
# Diferentes ASVs de 253 nt pueden compartir exactamente la región de
# referencia de 235 nt de Callahan.

if (length(G1_match) == 0) {
  stop("ERROR: no se encontró ningún ASV compatible con G1.")}

if (length(G2_match) == 0) {
  stop("ERROR: no se encontró ningún ASV compatible con G2.")}

if (length(G3_match) == 0) {
  stop("ERROR: no se encontró ningún ASV compatible con G3.")}

# Ningún ASV debe corresponder simultáneamente a dos variantes.
if (length(intersect(G1_match, G2_match)) > 0) {
  stop("ERROR: existen ASVs asignables simultáneamente a G1 y G2.")}

if (length(intersect(G1_match, G3_match)) > 0) {
  stop("ERROR: existen ASVs asignables simultáneamente a G1 y G3.")}

if (length(intersect(G2_match, G3_match)) > 0) {
  stop("ERROR: existen ASVs asignables simultáneamente a G2 y G3.")}

###############################################################################
###############################################################################

G_match_table <- bind_rows(
  data.frame(
    GardnerellaVariant = "G1",
    ASV = G1_match,
    stringsAsFactors = FALSE),
  data.frame(
    GardnerellaVariant = "G2",
    ASV = G2_match,
    stringsAsFactors = FALSE),
  data.frame(
    GardnerellaVariant = "G3",
    ASV = G3_match,
    stringsAsFactors = FALSE))

cat("\nCorrespondencias G1/G2/G3:\n")
print(G_match_table)

###############################################################################
# 10. Calcular abundancia total de los ASVs Gardnerella
###############################################################################

ps.gard <- prune_taxa(gard_sequences, ps)

gard.otu <- as(otu_table(ps.gard), "matrix")

# Dejar muestras en filas y ASVs en columnas.
if (taxa_are_rows(ps.gard)) {
  gard.otu <- t(gard.otu)}

gard.rank <- sort(
  colSums(gard.otu),
  decreasing = TRUE)

top.gard <- data.frame(
  ASV = names(gard.rank),
  Reads = as.numeric(gard.rank),
  stringsAsFactors = FALSE)

top.gard$Species <- tax[top.gard$ASV, "Species"]

top.gard <- top.gard %>%
  mutate(
    Rank = row_number(),
    Prop = Reads / sum(Reads),
    CumProp = cumsum(Prop))

cat("\nASVs Gardnerella más abundantes:\n")
print(head(top.gard, 20))
#       Reads   Species Rank         Prop   CumProp
# 1  10883817 vaginalis    1 4.455160e-01 0.4455160
# 2   7267043 vaginalis    2 2.974677e-01 0.7429837
# 3   3687349 vaginalis    3 1.509372e-01 0.8939209

###############################################################################
# 11. Añadir abundancia a la tabla G1/G2/G3
###############################################################################

G_match_table <- G_match_table %>%
  left_join(top.gard %>% select(ASV, Rank, Reads, Prop), by = "ASV") %>%
    arrange(GardnerellaVariant, Rank)

cat("\nG1/G2/G3 y abundancias:\n")
print(G_match_table)
#   GardnerellaVariant
# 1                 G1
# 2                 G1
# 3                 G2
# 4                 G2
# 5                 G3

#   Rank    Reads         Prop
# 1    1 10883817 4.455160e-01
# 2   10     2950 1.207547e-04
# 3    2  7267043 2.974677e-01
# 4   13      677 2.771218e-05
# 5    3  3687349 1.509372e-01

write_csv(
  G_match_table,
  file.path(project_dir, "results", "Gardnerella", "Gardnerella_G1_G2_G3_matches.csv"))

###############################################################################
# 12. Crear clasificación de todos los ASVs Gardnerella
###############################################################################

Gardnerella_ByHand <- data.frame(
  ASV = gard_sequences,
  GardnerellaVariant = "G_other",
  stringsAsFactors = FALSE)

Gardnerella_ByHand$GardnerellaVariant[
  Gardnerella_ByHand$ASV %in% G1_match] <- "G1"

Gardnerella_ByHand$GardnerellaVariant[
  Gardnerella_ByHand$ASV %in% G2_match] <- "G2"

Gardnerella_ByHand$GardnerellaVariant[
  Gardnerella_ByHand$ASV %in% G3_match] <- "G3"

cat("\nNúmero de ASVs por variante:\n")
print(table(Gardnerella_ByHand$GardnerellaVariant))
# G_other      G1      G2      G3 
#     32       2       2       1 

# Controles.
stopifnot(
  sum(Gardnerella_ByHand$GardnerellaVariant == "G1") == length(G1_match))

stopifnot(
  sum(Gardnerella_ByHand$GardnerellaVariant == "G2") == length(G2_match))

stopifnot(
  sum(Gardnerella_ByHand$GardnerellaVariant == "G3") == length(G3_match))

stopifnot(
  nrow(Gardnerella_ByHand) == nrow(gard_idx))

###############################################################################
# 13. Calcular abundancia total por variante
###############################################################################

gard.summary <- top.gard %>%
  left_join(Gardnerella_ByHand, by = "ASV") %>%
  group_by(GardnerellaVariant) %>%
  summarise(N_ASVs = n(), Reads = sum(Reads), .groups = "drop") %>%
  mutate(Prop = Reads / sum(Reads)) %>%
  arrange(desc(Reads))

cat("\n========================================\n")
cat("RESUMEN POR VARIANTE\n")
cat("========================================\n")
print(gard.summary)
# GardnerellaVariant  N_ASVs   Reads  Prop
#  <chr>               <int>    <dbl> <dbl>
# 1 G1                    2 10886767 0.446
# 2 G2                    2  7267720 0.297
# 3 G3                    1  3687349 0.151
# 4 G_other              32  2587851 0.106

###############################################################################
# 14. Proporción total G1 + G2 + G3
###############################################################################

dominant_prop <- gard.summary %>%
  filter(GardnerellaVariant %in% c("G1", "G2", "G3")) %>%
  summarise(Prop = sum(Prop)) %>%
  pull(Prop)

cat("\nProporción de Gardnerella representada por G1 + G2 + G3: ",
  round(dominant_prop * 100, 2),
  "%\n",
  sep = "")
# Proporción de Gardnerella representada por G1 + G2 + G3: 89.41%

###############################################################################
# 15. Curva de dominancia
###############################################################################

# La curva se utiliza como descripción secundaria de la estructura de
# abundancias. NO determina la asignación G1/G2/G3.

p.gard <- ggplot(
  top.gard,
  aes(x = Rank, y = CumProp)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.90, linetype = 2) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(
    x = "ASV rank",
    y = "Abundancia acumulativa",
    title = "Curva de dominancia - Gardnerella")

p.gard

ggsave(
  file.path(project_dir, "images", "Gardnerella", "DominanciaGard.png"),
  p.gard,
  width = 8,
  height = 5)

###############################################################################
# 16. Guardar clasificación y resumen
###############################################################################

saveRDS(
  Gardnerella_ByHand,
  file.path(project_dir, "cleaned", "Gardnerella_refined.rds"))

write_csv(
  Gardnerella_ByHand,
  file.path(project_dir, "results", "Gardnerella", "Gardnerella_refined.csv"))

write_csv(
  gard.summary,
  file.path(project_dir, "results", "Gardnerella", "Gardnerella_refinement_summary.csv"))

###############################################################################
# 17. Crear copia del objeto phyloseq
###############################################################################

ps_gard_corrected <- ps

###############################################################################
# 18. Incorporar GardnerellaVariant a tax_table
###############################################################################

tax_new <- as.data.frame(
  tax_table(ps_gard_corrected),
  stringsAsFactors = FALSE)

if (!"gardnerella_refined" %in% colnames(tax_new)) {
  tax_new$gardnerella_refined <- NA_character_}

idx <- match(rownames(tax_new), Gardnerella_ByHand$ASV)

tax_new$gardnerella_refined[!is.na(idx)] <-
  Gardnerella_ByHand$GardnerellaVariant[idx[!is.na(idx)]]

tax_table(ps_gard_corrected) <-
  tax_table(as.matrix(tax_new))

###############################################################################
# 19. Comprobaciones finales
###############################################################################

cat("\n========================================\n")
cat("CONTROL DEL NUEVO PHYLOSEQ\n")
cat("========================================\n")

cat("\nTaxones totales:", ntaxa(ps_gard_corrected), "\n")

cat("\nColumnas taxonómicas:\n")
print(colnames(tax_table(ps_gard_corrected)))

tax_check <- as.data.frame(
  tax_table(ps_gard_corrected),
  stringsAsFactors = FALSE)

gard_check <- tax_check %>% filter(Genus == "Gardnerella")

stopifnot(nrow(gard_check) == n_gard)
stopifnot(!anyNA(gard_check$gardnerella_refined))

cat("\nClasificación dentro de Gardnerella:\n")
print(table(gard_check$gardnerella_refined))

###############################################################################
# 20. Guardar nuevo objeto phyloseq
###############################################################################

saveRDS(ps_gard_corrected, file.path(project_dir, "cleaned", "ps_gard_refined.rds"))

###############################################################################
# 21. Resumen final
###############################################################################

cat("\n========================================\n")
cat("PROCESO COMPLETADO\n")
cat("========================================\n")

cat("\nG1/G2/G3 fueron definidos mediante coincidencia exacta con las ",
  "secuencias recuperadas de processed.rda de Callahan.\n",
  sep = "")

cat("Número de ASVs G1: ", length(G1_match),
  "\nNúmero de ASVs G2: ", length(G2_match),
  "\nNúmero de ASVs G3: ", length(G3_match),
  "\n",
  sep = "")

cat("\nNuevo objeto phyloseq:\n",
  file.path(project_dir, "cleaned", "ps_gard_refined.rds\n"),
  sep = "")

cat("========================================\n")


###############################################################################
# Reproducibilidad
###############################################################################

software_versions <- data.frame(
  Package = c("R", "phyloseq", "dplyr", "ggplot2", "readr"),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("phyloseq")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("readr"))),
  stringsAsFactors = FALSE)

write_csv(
  software_versions,
  file.path(dirname(results_dir), "software_versions.csv"))

cat("\n=== VERSIONES ===\n")
print(software_versions)
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

