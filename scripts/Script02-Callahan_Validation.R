###############################################################################
#                                    TFM
###############################################################################
#         Script 02: Comparación y validación frente a Callahan 2017
###############################################################################

# Objetivo:
# Comparar el objeto phyloseq reconstruido a partir del reprocesamiento en
# Galaxy con los objetos publicados por Callahan 2017.

# El script evalúa:
#   1. Número de muestras y ASVs.
#   2. Profundidad de secuenciación.
#   3. Profundidad y riqueza en las muestras presentes en ambos datasets.
#   4. Composición taxonómica general.
#   5. Coincidencia de secuencias ASV tras armonizar la región comparada.
#   6. Concordancia de abundancias de los ASVs compartidos.
#   7. Concordancia taxonómica de los ASVs compartidos.

# IMPORTANTE:
# Este script es exclusivamente de validación y NO modifica el phyloseq.

# Entradas:
#   - phyloseqObject_sequences.rds
#   - processed.rda
#   - ps_allsams.rds

# Salidas:
#   - tablas CSV de comparación y validación
#   - gráficos de profundidad y abundancia

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
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)

dir.create(project_path("results"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("images"), recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2. Cargar objetos
###############################################################################

ps <- readRDS(project_path("cleaned", "phyloseqObject_sequences.rds"))

load(project_path("Callahan2017", "RepRefine_Scripts", "input", "processed.rda"))
# processed.rda contiene, entre otros:
# st  = tabla de conteos de Callahan
# ft  = tabla utilizada en los análisis publicados
# tax = tabla taxonómica

ps_allsams <- readRDS(
  project_path("Callahan2017", "RepRefine_Scripts", "input", "ps_allsams.rds"))

###############################################################################
# 3. Comparación general
###############################################################################

comparison <- data.frame(
  Dataset = c("Galaxy", "processed.rda", "ps_allsams"),
  Samples = c(nsamples(ps), nrow(st), nsamples(ps_allsams)),
  ASVs = c(ntaxa(ps), ncol(st), ntaxa(ps_allsams)))

print(comparison)
#         Dataset Samples ASVs
# 1        Galaxy    2367 8298
# 2 processed.rda    2179 2699
# 3    ps_allsams    6088 5481

write_csv(
  comparison,
  project_path("results", "Script02_Comparison_General.csv"))

###############################################################################
# 4. Profundidad de secuenciación
###############################################################################

# Profundidad = número total de reads de una muestra sumando todos sus ASVs.

reads_galaxy <- sample_sums(ps)
reads_processed <- rowSums(st)
reads_ps_allsams <- sample_sums(ps_allsams)

depth_table <- data.frame(
  Dataset = c("Galaxy", "processed.rda", "ps_allsams"),
  N = c(length(reads_galaxy), length(reads_processed), length(reads_ps_allsams)),
  Mean = c(mean(reads_galaxy), mean(reads_processed), mean(reads_ps_allsams)),
  Median = c(median(reads_galaxy), median(reads_processed), median(reads_ps_allsams)),
  Min = c(min(reads_galaxy), min(reads_processed), min(reads_ps_allsams)),
  Max = c(max(reads_galaxy), max(reads_processed), max(reads_ps_allsams)))

print(depth_table)
#         Dataset    N     Mean   Median   Min    Max
# 1        Galaxy 2367 147013.7 153696.0     4 592866
# 2 processed.rda 2179 159697.6 157543.0 28580 592382
# 3    ps_allsams 6088 142660.9 145741.5     0 592382

# Esto sugiere fuertemente que el conjunto processed.rda de Callahan ya 
# contiene una selección/filtrado de muestras, mientras que el objeto Galaxy 
# conserva muestras adicionales, incluidas algunas extremadamente poco profundas.

write_csv(
  depth_table,
  project_path("results", "Script02_SequencingDepth_General.csv"))

depth_df <- bind_rows(
  data.frame(Reads = reads_galaxy, Dataset = "Galaxy"),
  data.frame(Reads = reads_processed, Dataset = "processed.rda"),
  data.frame(Reads = reads_ps_allsams, Dataset = "ps_allsams"))

zero_depth <- depth_df %>%
  group_by(Dataset) %>%
  summarise(
    N_zero_reads = sum(Reads == 0),
    .groups = "drop")

cat("\nMuestras con 0 reads:\n")
print(zero_depth)
#   Dataset       N_zero_reads
#    <chr>                <int>
# 1 Galaxy                   0
# 2 processed.rda            0
# 3 ps_allsams               8

write_csv(
  zero_depth,
  project_path("results", "Script02_ZeroDepthSamples.csv"))

p_depth <- depth_df %>%
  filter(Reads > 0) %>%
  ggplot(aes(x = Dataset, y = Reads)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "Sequencing depth",
    subtitle = "Samples with 0 reads excluded from the logarithmic visualization",
    x = NULL,
    y = "Reads per sample (log10 scale)")
p_depth

ggsave(
  project_path("images", "Script02_SequencingDepth.png"),
  p_depth,
  width = 7,
  height = 5,
  dpi = 300)

###############################################################################
# 5. Identificar muestras compartidas entre Galaxy y processed.rda
###############################################################################

# Esta comparación permite analizar exactamente las mismas muestras en ambos
# pipelines y evita que las diferencias se deban a muestras adicionales.

# Buscar correspondencias entre rownames(st) y metadata Galaxy
#-------------------------------------------------------------
head(sample_names(ps), 20)
head(rownames(st), 20)
head(rownames(df), 20)

meta_galaxy <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

meta_galaxy$phyloseq_sample_name <- sample_names(ps)

overlap_st_vs_galaxy <- sapply(
  meta_galaxy,
  function(x) {length(intersect(as.character(rownames(st)), as.character(x)))})

overlap_st_vs_galaxy <- sort(overlap_st_vs_galaxy, decreasing = TRUE)

cat("\nCoincidencias entre rownames(st) y columnas Galaxy:\n")
print(overlap_st_vs_galaxy)

nrow(st) # muestras st: 2179
nrow(ft) # muestras ft: 2179
nrow(df) # filas df: 2179

identical(rownames(st), rownames(df)) # rownames(st) == rownames(df): TRUE
identical(rownames(ft), rownames(df)) # rownames(ft) == rownames(df): TRUE
identical(rownames(st), rownames(ft)) # rownames(st) == rownames(ft): TRUE

meta_galaxy <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

sum(meta_galaxy$Library.Name == meta_galaxy$Sample.Name, na.rm = TRUE) # 2367

sum(!is.na(meta_galaxy$Library.Name)) # 2367

sum(!is.na(meta_galaxy$Sample.Name)) # 2367

head(meta_galaxy[, c("Library.Name", "Sample.Name")], 20) 
# Library.Name y Sample.Name son equivalentes

# sample_names(ps) contiene accesiones SRR.
# Los objetos de Callahan st, ft y df utilizan como rownames identificadores
# que corresponden a Sample.Name/Library.Name de la metadata Galaxy.
# Por tanto, el emparejamiento debe realizarse utilizando Sample.Name.

meta_galaxy$SRR <- sample_names(ps)

shared_sample_names <- intersect(
  as.character(meta_galaxy$Sample.Name),
  rownames(st))

length(shared_sample_names) # Muestras compartidas Galaxy / Callahan: 2179

# Se esperan las 2179 muestras utilizadas en processed.rda
if (length(shared_sample_names) != 2179) {
  warning(
    paste0(
      "Se esperaban 2179 muestras compartidas y se encontraron ",
      length(shared_sample_names), "."))}

# Construir tabla de correspondencia
#-----------------------------------
sample_mapping <- meta_galaxy %>%
  filter(Sample.Name %in% shared_sample_names) %>%
  select(SRR, Sample.Name, Library.Name)

# Comprobar que no existan identificadores Sample.Name duplicados
if (anyDuplicated(sample_mapping$Sample.Name)) {
  stop("ERROR: existen Sample.Name duplicados en la metadata Galaxy.")}

cat("Filas de la tabla de correspondencia:", nrow(sample_mapping), "\n") # 2179

write_csv(
  sample_mapping,
  project_path("results", "Script02_Sample_Mapping_Galaxy_Callahan.csv"))

###############################################################################
# 6. Comparación pareada de profundidad de secuenciación
###############################################################################

# Para cada muestra biológica comparamos:
#   total de reads en Galaxy vs total de reads en processed.rda

galaxy_depth <- data.frame(
  SRR = names(sample_sums(ps)),
  Galaxy_Depth = as.numeric(sample_sums(ps)),
  stringsAsFactors = FALSE)

depth_shared <- sample_mapping %>%
  left_join(galaxy_depth, by = "SRR") %>%
  mutate(Callahan_Depth = as.numeric(rowSums(st)[Sample.Name]))

stopifnot(!anyNA(depth_shared$Galaxy_Depth))

stopifnot(!anyNA(depth_shared$Callahan_Depth))

depth_cor <- cor.test(
  depth_shared$Galaxy_Depth,
  depth_shared$Callahan_Depth,
  method = "spearman",
  exact = FALSE)

cat("\nSpearman profundidad muestras compartidas:",
  unname(depth_cor$estimate), "\n") 
# Spearman profundidad muestras compartidas: 0.9956144

write_csv(
  depth_shared,
  project_path("results", "Script02_SharedSamples_Depth.csv"))

###############################################################################
# 7. Comparación pareada de riqueza observada
###############################################################################

# Riqueza observada = número de ASVs con conteo > 0 en una muestra.

ps_shared <- prune_samples(sample_mapping$SRR, ps)

rich_galaxy_all <- estimate_richness(ps_shared, measures = "Observed")

rich_galaxy_table <- data.frame(
  SRR = rownames(rich_galaxy_all),
  Galaxy_Richness = rich_galaxy_all$Observed,
  stringsAsFactors = FALSE)

# Riqueza Callahan
rich_callahan_all <- rowSums(ft > 0)

# Unir mediante Sample.Name
richness_shared <- sample_mapping %>%
  left_join(rich_galaxy_table, by = "SRR") %>%
  mutate(
    Callahan_Richness = as.numeric(
      rich_callahan_all[Sample.Name]))

stopifnot(!anyNA(richness_shared$Galaxy_Richness))

stopifnot(!anyNA(richness_shared$Callahan_Richness))

richness_cor <- cor.test(
  richness_shared$Galaxy_Richness,
  richness_shared$Callahan_Richness,
  method = "spearman",
  exact = FALSE)

cat("Spearman riqueza muestras compartidas:",
  unname(richness_cor$estimate), "\n")
# Spearman riqueza muestras compartidas: 0.997853

write_csv(
  richness_shared,
  project_path("results", "Script02_SharedSamples_Richness.csv"))

###############################################################################
# 8. Muestras Galaxy no incluidas en processed.rda
###############################################################################

galaxy_not_callahan <- meta_galaxy %>%
  filter(!Sample.Name %in% rownames(st)) %>%
  mutate(Sequencing_Depth = as.numeric(sample_sums(ps)[SRR]))

cat("\nMuestras Galaxy no presentes en processed.rda:",
  nrow(galaxy_not_callahan), "\n") # 188
 
cat("\nProfundidad de las muestras no incluidas:\n")
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#     4.0    182.8    423.5   7360.5   2080.8 229960.0 

print(summary(galaxy_not_callahan$Sequencing_Depth))

write_csv(
  galaxy_not_callahan,
  project_path("results", "Script02_Galaxy_Samples_NotIn_Callahan.csv"))

###############################################################################
# 9. Comparación taxonómica general
###############################################################################

tax_galaxy <- as.data.frame(tax_table(ps))

gen_galaxy <- unique(na.omit(tax_galaxy$Genus))
gen_callahan <- unique(na.omit(tax[, "Genus"]))

species_galaxy <- unique(na.omit(tax_galaxy$Species))
species_callahan <- unique(na.omit(tax[, "Species"]))

shared_genus <- intersect(gen_galaxy, gen_callahan)
shared_species <- intersect(species_galaxy, species_callahan)

# Índice de Jaccard:
# número de taxones compartidos / número total de taxones presentes entre
# ambos pipelines.
jaccard_genus <- length(shared_genus) /
  length(union(gen_galaxy, gen_callahan))

jaccard_species <- length(shared_species) /
  length(union(species_galaxy, species_callahan))

taxonomy_general <- data.frame(
  Rank = c("Genus", "Species"),
  Galaxy_unique = c(length(gen_galaxy), length(species_galaxy)),
  Callahan_unique = c(length(gen_callahan), length(species_callahan)),
  Shared = c(length(shared_genus), length(shared_species)),
  Jaccard = c(jaccard_genus, jaccard_species))

print(taxonomy_general)
#      Rank Galaxy_unique Callahan_unique Shared    Jaccard
# 1   Genus           636             345    284 0.40746055
# 2 Species           449              13      7 0.01538462

write_csv(
  taxonomy_general,
  project_path("results", "Script02_Taxonomy_General.csv"))

###############################################################################
# 10. Comparar longitud de las secuencias ASV
###############################################################################

# Antes de comparar directamente las secuencias de ambos pipelines se evalúa
# su longitud, ya que las ASVs publicadas por Callahan y las obtenidas en
# Galaxy no necesariamente representan exactamente la misma extensión del
# amplicón.

asv_galaxy <- taxa_names(ps)
asv_callahan <- rownames(tax)

galaxy_lengths <- nchar(asv_galaxy)
callahan_lengths <- nchar(asv_callahan)

cat("\n========================================\n")
cat("LONGITUD DE SECUENCIAS ASV\n")
cat("========================================\n")

cat("\nGalaxy:\n")
print(summary(galaxy_lengths))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   240.0   253.0   253.0   253.4   253.0   425.0 

print(table(galaxy_lengths))

cat("\nCallahan:\n")
print(summary(callahan_lengths))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   235.0   235.0   235.0   235.1   235.0   295.0 

print(table(callahan_lengths))

length_summary <- data.frame(
  Dataset = c("Galaxy", "Callahan"),
  N_ASVs = c(length(asv_galaxy), length(asv_callahan)),
  Min = c(min(galaxy_lengths), min(callahan_lengths)),
  Q1 = c(
    unname(quantile(galaxy_lengths, 0.25)),
    unname(quantile(callahan_lengths, 0.25))),
  Median = c(median(galaxy_lengths), median(callahan_lengths)),
  Mean = c(mean(galaxy_lengths), mean(callahan_lengths)),
  Q3 = c(
    unname(quantile(galaxy_lengths, 0.75)),
    unname(quantile(callahan_lengths, 0.75))),
  Max = c(max(galaxy_lengths), max(callahan_lengths)))

print(length_summary)
#    Dataset N_ASVs Min  Q1 Median     Mean  Q3 Max
# 1   Galaxy   8298 240 253    253 253.3606 253 425
# 2 Callahan   2699 235 235    235 235.0889 235 295

write_csv(
  length_summary,
  project_path("results", "Script02_ASV_Length_Summary.csv"))

# Longitudes más frecuentes
#---------------------------
galaxy_length_mode <- as.numeric(
  names(which.max(table(galaxy_lengths))))

callahan_length_mode <- as.numeric(
  names(which.max(table(callahan_lengths))))

cat("\nLongitud más frecuente Galaxy:",
  galaxy_length_mode, "nt\n") # 253 nt

cat("Longitud más frecuente Callahan:",
  callahan_length_mode, "nt\n") # 235 nt

# Las ASVs del procesamiento Galaxy son generalmente de 253 nt, mientras que
# las secuencias publicadas por Callahan tienen mayormente 235 nt.

###############################################################################
# 11. Comprobar posición de la región Callahan dentro de Galaxy
###############################################################################

gal_example <- asv_galaxy[1] 
cal_example <- asv_callahan[1]

cat("\nEjemplo Galaxy:\n")
cat(gal_example, "\n")

cat("\nEjemplo Callahan:\n")
cat(cal_example, "\n")

cat("\nLongitud Galaxy:", nchar(gal_example), "nt\n")
# Longitud Galaxy: 253 nt

cat("Longitud Callahan:", nchar(cal_example), "nt\n")
# Longitud Callahan: 235 nt

match_position <- regexpr(cal_example, gal_example, fixed = TRUE)

cat("\nPosición inicial de Callahan dentro de Galaxy:", match_position, "\n")
# Posición inicial de Callahan dentro de Galaxy: 11 

if (match_position > 0) {
  extra_5prime <- match_position - 1
  extra_3prime <- nchar(gal_example) - extra_5prime - nchar(cal_example)
  cat("Bases adicionales Galaxy en extremo 5':", extra_5prime, "\n")
  cat("Bases adicionales Galaxy en extremo 3':", extra_3prime, "\n")
} else {
  warning(
    "La secuencia Callahan seleccionada no está contenida exactamente ",
    "en la secuencia Galaxy seleccionada.")}
# Bases adicionales Galaxy en extremo 5': 10 
# Bases adicionales Galaxy en extremo 3': 8 

# Esto demuestra que la región Callahan está contenida dentro
# de la secuencia Galaxy, dejando 10 nt adicionales en el extremo 5' y 8 nt
# en el extremo 3'.

###############################################################################
# 12. Evaluar si el recorte 10 nt + 8 nt reproduce la región Callahan
###############################################################################

# La comparación inicial sugiere que las secuencias Galaxy contienen una
# región interna equivalente a la publicada por Callahan, con 10 nt
# adicionales en el extremo 5' y 8 nt adicionales en el extremo 3'.
# Para comparar exactamente la misma región eliminamos esos extremos.
# Se prueba esta armonización sobre todas las secuencias Galaxy compatibles
# con esa longitud.

galaxy_trim <- substring(
  asv_galaxy,
  first = 11,
  last = nchar(asv_galaxy) - 8)

cat("\nLongitud después del trimming:\n")
print(summary(nchar(galaxy_trim)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 222.0   235.0   235.0   235.4   235.0   407.0 

length(asv_galaxy) # ASVs Galaxy originales: 8298

cat("Secuencias únicas después del trimming:",
  length(unique(galaxy_trim)), "\n") 
# Secuencias únicas después del trimming: 8072 

cat("ASVs que colapsan después del trimming:",
  length(galaxy_trim) - length(unique(galaxy_trim)), "\n")
# ASVs que colapsan después del trimming: 226 

###############################################################################
# 13. ASVs compartidos
###############################################################################

shared_asvs <- intersect(
  unique(galaxy_trim),
  asv_callahan)

n_shared_asvs <- length(shared_asvs) 
n_shared_asvs
# Secuencias armonizadas Galaxy presentes exactamente en Callahan: 1875

percent_callahan_recovered <- round( 
  100 * n_shared_asvs / length(unique(asv_callahan)), 2)
percent_callahan_recovered
# Porcentaje de ASVs Callahan recuperados: 69.47

percent_galaxy_shared <- 100 *
  n_shared_asvs /
  length(unique(galaxy_trim))
percent_galaxy_shared
# ASVs Galaxy armonizados compartidos: 23.22844

###############################################################################
# 14. Abundancia total de ASVs compartidos
###############################################################################

# taxa_sums() = número total de reads asignados a cada ASV sumando todas las
# muestras del dataset.

# Si varios ASVs Galaxy se vuelven idénticos después del trimming, sus
# conteos deben sumarse antes de compararlos con Callahan.

galaxy_counts <- taxa_sums(ps)

galaxy_df <- data.frame(
  ASV = galaxy_trim,
  Galaxy = as.numeric(galaxy_counts), stringsAsFactors = FALSE) %>%
  group_by(ASV) %>%
  summarise(
    Galaxy = sum(Galaxy),
    N_original_ASVs = n(),
    .groups = "drop")

callahan_df <- data.frame(
  ASV = colnames(st),
  Callahan = colSums(st),
  stringsAsFactors = FALSE)

shared_counts <- inner_join(
  galaxy_df,
  callahan_df,
  by = "ASV")

stopifnot(nrow(shared_counts) == n_shared_asvs)

###############################################################################
# 12. Concordancia de abundancias
###############################################################################

# Spearman evalúa si los ASVs que son más abundantes en Galaxy también tienden
# a ser los más abundantes en Callahan.

# rho cercano a 1 = orden de abundancias muy similar.

abundance_cor <- cor.test(
  shared_counts$Galaxy,
  shared_counts$Callahan,
  method = "spearman",
  exact = FALSE)

cat("\nSpearman abundancia ASVs compartidos:",
    unname(abundance_cor$estimate), "\n")
# Spearman abundancia ASVs compartidos: 0.9706608

cat("p-value:", abundance_cor$p.value, "\n") # p-value: 0 

write_csv(
  shared_counts,
  project_path("results", "Script02_SharedASV_Abundances.csv"))

p_shared <- ggplot(
  shared_counts,
  aes(x = log10(Galaxy + 1), y = log10(Callahan + 1))) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  labs(
    title = "Abundance of shared ASVs",
    x = "Galaxy log10(total reads + 1)",
    y = "Callahan log10(total reads + 1)")
p_shared

ggsave(
  project_path("images", "Script02_Shared_ASV_Abundances.png"),
  p_shared,
  width = 8,
  height = 5,
  dpi = 300)

###############################################################################
# 13. Evaluar el colapso taxonómico provocado por el trimming
###############################################################################

# Varios ASVs Galaxy pueden convertirse en la misma secuencia de 235 nt.
# No elegimos arbitrariamente la taxonomía del primer ASV.

# Para cada secuencia recortada conservamos el género/especie únicamente si
# todas las asignaciones no-NA son concordantes.

gal_tax_raw <- data.frame(
  ASV = galaxy_trim,
  Genus = tax_galaxy$Genus,
  Species = tax_galaxy$Species,
  stringsAsFactors = FALSE)

collapse_taxonomy <- function(x) {
  x <- unique(na.omit(x))
  if (length(x) == 0) return(NA_character_)
  if (length(x) == 1) return(x)
  return("Discordant")}

gal_tax <- gal_tax_raw %>%
  group_by(ASV) %>%
  summarise(
    Genus_Galaxy = collapse_taxonomy(Genus),
    Species_Galaxy = collapse_taxonomy(Species),
    N_Galaxy_ASVs = n(),
    .groups = "drop")

call_tax <- data.frame(
  ASV = rownames(tax),
  Genus_Callahan = tax[, "Genus"],
  Species_Callahan = tax[, "Species"],
  stringsAsFactors = FALSE)

shared_tax <- inner_join(
  gal_tax,
  call_tax,
  by = "ASV")

###############################################################################
# 14. Concordancia taxonómica de ASVs compartidos
###############################################################################

genus_comparable <- shared_tax %>%
  filter(
    !is.na(Genus_Galaxy),
    !is.na(Genus_Callahan),
    Genus_Galaxy != "Discordant")

species_comparable <- shared_tax %>%
  filter(
    !is.na(Species_Galaxy),
    !is.na(Species_Callahan),
    Species_Galaxy != "Discordant")

genus_concordance <- mean(
  genus_comparable$Genus_Galaxy ==
    genus_comparable$Genus_Callahan)

species_concordance <- mean(
  species_comparable$Species_Galaxy ==
    species_comparable$Species_Callahan)

taxonomy_concordance <- data.frame(
  Rank = c("Genus", "Species"),
  Comparable_ASVs = c(
    nrow(genus_comparable),
    nrow(species_comparable)),
  Concordant_ASVs = c(
    sum(genus_comparable$Genus_Galaxy ==
          genus_comparable$Genus_Callahan),
    sum(species_comparable$Species_Galaxy ==
          species_comparable$Species_Callahan)),
  Concordance = c(
    genus_concordance,
    species_concordance))

print(taxonomy_concordance)
#      Rank Comparable_ASVs Concordant_ASVs Concordance
# 1   Genus            1455            1119   0.7690722
# 2 Species               8               6   0.7500000

write_csv(
  shared_tax,
  project_path("results", "Script02_SharedASV_Taxonomy.csv"))

write_csv(
  taxonomy_concordance,
  project_path("results", "Script02_Taxonomic_Concordance.csv"))

###############################################################################
# 15. Tabla final de validación
###############################################################################

validation_table <- data.frame(
  Metric = c(
    "Galaxy samples",
    "Callahan processed samples",
    "Shared samples",
    "Galaxy ASVs",
    "Callahan ASVs",
    "Shared ASVs after sequence harmonization",
    "Callahan ASVs recovered (%)",
    "Galaxy harmonized ASVs shared (%)",
    "Spearman shared-ASV abundance",
    "Spearman shared-sample depth",
    "Spearman shared-sample richness",
    "Genus Jaccard",
    "Species Jaccard",
    "Genus concordance shared ASVs",
    "Species concordance shared ASVs"),
  Value = c(
    nsamples(ps),
    nrow(st),
    length(shared_sample_names),
    ntaxa(ps),
    ncol(st),
    n_shared_asvs,
    round(percent_callahan_recovered, 2),
    round(percent_galaxy_shared, 2),
    round(unname(abundance_cor$estimate), 4),
    round(unname(depth_cor$estimate), 4),
    round(unname(richness_cor$estimate), 4),
    round(jaccard_genus, 4),
    round(jaccard_species, 4),
    round(genus_concordance, 4),
    round(species_concordance, 4)))

print(validation_table)

write_csv(
  validation_table,
  project_path("results", "Script02_Pipeline_Validation.csv"))

###############################################################################
# 16. Fin
###############################################################################

cat("\n========================================\n")
cat("SCRIPT 02 FINALIZADO\n")
cat("========================================\n")
cat("Muestras Galaxy:", nsamples(ps), "\n")
# Muestras Galaxy: 2367 
cat("Muestras Callahan processed:", nrow(st), "\n")
# Muestras Callahan processed: 2179 
cat("Muestras compartidas:", length(shared_sample_names), "\n")
# Muestras compartidas: 2179 
cat("ASVs Galaxy:", ntaxa(ps), "\n") # ASVs Galaxy: 8298 
cat("ASVs Callahan:", ncol(st), "\n") # ASVs Callahan: 2699 
cat("ASVs compartidos:", n_shared_asvs, "\n") # ASVs compartidos: 1875 
cat("Spearman abundancias:", round(unname(abundance_cor$estimate), 4), "\n")
# Spearman abundancias: 0.9707 
cat("Concordancia género:", round(genus_concordance * 100, 2), "%\n")
# Concordancia género: 76.91 %
cat("Concordancia especie:", round(species_concordance * 100, 2), "%\n")
# Concordancia especie: 75 %
