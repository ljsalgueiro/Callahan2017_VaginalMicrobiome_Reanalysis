###############################################################################
#                                    TFM
###############################################################################
# Script 07B: Prevotella - asignación BLAST y preparación para vSpeciateDB
###############################################################################

# Objetivo:
#   Filtrar hits BLAST, construir el Best Hit Set, generar prevotella_byBlast,
#   incorporar la clasificación al phyloseq y corregir orientación antes de
#   ejecutar SpeciateIT/vSpeciateDB.

# Entrada:
#   ps_gard_refined.rds
#   Prevotella_ASV_mapping.csv
#   Prevotella_BLAST_results_all.csv

# Salidas:
#   Prevotella_BestHitSet.csv
#   Prevotella_species_byBlast.csv
#   Prevotella_BLAST_true_ties.csv
#   Prevotella_BLAST_tie_patterns.csv
#   ps_prevotella_byBlast.rds
#   Prevotella_vSpeciate_orientation.csv
#   Prevotella_ASVs_vSpeciateDB.fasta
#   vSpeciateDB_metadata.csv
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
# Librerías
###############################################################################

library(phyloseq)
library(Biostrings)
library(dplyr)
library(readr)
library(stringr)

###############################################################################
# 1. Rutas y parámetros
###############################################################################

prev_dir <- file.path(project_dir, "results", "Prevotella")
blast_dir <- file.path(prev_dir, "Blast")
vs_dir <- file.path(project_dir, "results", "vSpeciateDB", "Prevotella")
dir.create(vs_dir, recursive = TRUE, showWarnings = FALSE)

input_ps <- file.path(project_dir, "cleaned", "ps_gard_refined.rds")
mapping_file <- file.path(prev_dir, "Prevotella_ASV_mapping.csv")
blast_file <- file.path(blast_dir, "Prevotella_BLAST_results_all.csv")

MIN_QUERY_COVER <- 95
MIN_PER_IDENT <- 99


required_files <- c(
  input_ps,
  mapping_file,
  blast_file)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Faltan archivos de entrada:",
      paste(missing_files, collapse = ", ")))
}

###############################################################################
# 2. Leer datos
###############################################################################

ps <- readRDS(input_ps)
mapping <- read_csv(mapping_file, show_col_types = FALSE)
blast <- read_csv(blast_file, show_col_types = FALSE)

required <- c("ASV_ID", "Accession", "Description", "Scientific_Name",
  "Query_Cover", "Per_Ident", "E_value", "Bit_Score")

if (!all(required %in% names(blast))) {
  stop("ERROR: faltan columnas necesarias en Prevotella_BLAST_results_all.csv.")}

###############################################################################
# 3. Filtro de calidad y Best Hit Set
###############################################################################

blast_qc <- blast %>%
  filter(!is.na(Query_Cover), !is.na(Per_Ident), !is.na(E_value),
    Query_Cover >= MIN_QUERY_COVER, Per_Ident >= MIN_PER_IDENT)

best_hits <- blast_qc %>%
  group_by(ASV_ID) %>%
  filter(
    round(Query_Cover, 6) ==
      round(max(Query_Cover, na.rm = TRUE), 6)) %>%
  filter(
    round(Per_Ident, 6) ==
      round(max(Per_Ident, na.rm = TRUE), 6)) %>%
  filter(
    E_value == min(E_value, na.rm = TRUE)) %>%
  ungroup()

normalize_prevotella <- function(x) {
  x <- str_squish(x)
  unresolved <- is.na(x) | x == "" |
    str_detect(str_to_lower(x), "uncultured|unidentified|unclassified|metagenome|bacterium") |
    str_detect(str_to_lower(x), "^prevotella sp\\.")
  result <- rep(NA_character_, length(x))
  valid <- !unresolved
  result[valid] <- str_extract(x[valid], "^Prevotella\\s+[A-Za-z0-9_-]+")
  result}

best_hits <- best_hits %>%
  mutate(Species_normalized = normalize_prevotella(Scientific_Name),
    Taxon_type = if_else(is.na(Species_normalized), "unresolved", "species"))

write_csv(best_hits, file.path(blast_dir, "Prevotella_BestHitSet.csv"))

###############################################################################
# 4. Resumir y asignar BLAST
###############################################################################

format_species_tie <- function(x) {
  if (is.na(x)) return("Prevotella sp.")
  species <- str_split(x, "\\s*\\|\\s*")[[1]]
  species <- str_squish(species)
  if (length(species) == 1) return(species)
  epithet <- str_remove(species, "^Prevotella\\s+")
  paste0("Prevotella ", paste(epithet, collapse = "/"))
}

best_summary <- best_hits %>%
  group_by(ASV_ID) %>%
  summarise(N_best_hits = n(),
    N_valid_species = n_distinct(Species_normalized[Taxon_type == "species"]),
    N_unresolved_hits = sum(Taxon_type == "unresolved"),
    Species_clean = {
      sp <- sort(unique(Species_normalized[Taxon_type == "species"]))
      sp <- sp[!is.na(sp)]
      if (length(sp) == 0) NA_character_ else paste(sp, collapse = " | ")
    }, .groups = "drop")

all_asvs <- mapping %>% dplyr::select(ASV_ID) %>% distinct()

best_summary <- all_asvs %>%
  left_join(best_summary, by = "ASV_ID") %>%
  mutate(Has_high_quality_hit = !is.na(N_best_hits),
    N_best_hits = coalesce(N_best_hits, 0L),
    N_valid_species = coalesce(N_valid_species, 0L),
    N_unresolved_hits = coalesce(N_unresolved_hits, 0L),
    species_byBlast = vapply(seq_len(n()), function(i) {
      if (!Has_high_quality_hit[i] || N_valid_species[i] == 0) return("Prevotella sp.")
      format_species_tie(Species_clean[i])
    }, character(1)),
    Assignment_type = case_when(!Has_high_quality_hit ~ "No_high_quality_hit",
      N_valid_species == 0 ~ "Unresolved_best_hit",
      N_valid_species == 1 ~ "Single_species",
      N_valid_species > 1 ~ "True_species_tie"),
    ASV_number = as.numeric(str_remove(ASV_ID, "^Prev_ASV_"))) %>%
  arrange(ASV_number) %>%
  dplyr::select(-ASV_number)

write_csv(best_summary, file.path(blast_dir, "Prevotella_species_byBlast.csv"))

true_ties <- best_summary %>% filter(Assignment_type == "True_species_tie")
tie_patterns <- true_ties %>% count(species_byBlast, name = "N_ASVs", sort = TRUE)

write_csv(true_ties, file.path(blast_dir, "Prevotella_BLAST_true_ties.csv"))
write_csv(tie_patterns, file.path(blast_dir, "Prevotella_BLAST_tie_patterns.csv"))

cat("\nTIPOS DE ASIGNACIÓN BLAST\n")
print(table(best_summary$Assignment_type))
# No_high_quality_hit      Single_species    True_species_tie Unresolved_best_hit 
#                 67                 101                  11                 186 

cat("\nPatrones de empate:\n")
print(tie_patterns, n = Inf)
# # A tibble: 6 × 2
#   species_byBlast                                N_ASVs
#   <chr>                                           <int>
#1 Prevotella histicola/melaninogenica                 5
#2 Prevotella melaninogenica/scopos                    2
#3 Prevotella fusca/melaninogenica                     1
#4 Prevotella jejuni/melaninogenica                    1
#5 Prevotella melaninogenica/veroralis                 1
#6 Prevotella melaninogenica/veroralis/vespertina      1

###############################################################################
# 5. Incorporar prevotella_byBlast al phyloseq
###############################################################################

blast_mapping <- mapping %>%
  left_join(best_summary %>% dplyr::select(ASV_ID, species_byBlast), by = "ASV_ID")

if (anyNA(blast_mapping$species_byBlast)) stop("ERROR: faltan asignaciones BLAST.")

tax_matrix <- as(tax_table(ps), "matrix")
prevotella_byBlast_vector <- rep(NA_character_, nrow(tax_matrix))
names(prevotella_byBlast_vector) <- rownames(tax_matrix)

lookup <- setNames(blast_mapping$species_byBlast, blast_mapping$Sequence)
matching <- intersect(rownames(tax_matrix), names(lookup))
prevotella_byBlast_vector[matching] <- lookup[matching]

tax_table(ps) <- tax_table(cbind(tax_matrix, prevotella_byBlast = prevotella_byBlast_vector))
saveRDS(ps, file.path(project_dir, "cleaned", "ps_prevotella_byBlast.rds"))

###############################################################################
# 6. Corregir orientación antes de vSpeciateDB
###############################################################################

orientation <- mapping %>%
  mutate(Reverse_Complement = as.character(reverseComplement(DNAStringSet(Sequence))),
    Forward_start = substr(Sequence, 1, 10),
    RevComp_start = substr(Reverse_Complement, 1, 10),
    Orientation = case_when(str_detect(Sequence, "^TACGG") ~ "Forward",
      str_detect(Reverse_Complement, "^TACGG") ~ "Reverse",
      TRUE ~ "Review"),
    Sequence_vSpeciate = if_else(Orientation == "Reverse", Reverse_Complement, Sequence))

write_csv(orientation, file.path(vs_dir, "Prevotella_vSpeciate_orientation.csv"))

vs_seqs <- DNAStringSet(orientation$Sequence_vSpeciate)
names(vs_seqs) <- orientation$ASV_ID
writeXStringSet(vs_seqs, filepath = file.path(vs_dir, "Prevotella_ASVs_vSpeciateDB.fasta"), format = "fasta")

vs_metadata <- data.frame(
  Tool = "SpeciateIT",
  Database = "vSpeciateDB V4 models",
  Model_directory = "vSpeciateIT_V4V4",
  Reference_database_source = "GTDB-SSU release 214.1",
  Publication = "Holm JB, Gajer P, Ravel J. BMC Bioinformatics. 2024;25:313.",
  DOI = "10.1186/s12859-024-05930-3",
  Models_DOI = "10.6084/m9.figshare.25254229",
  stringsAsFactors = FALSE)

write_csv(vs_metadata, file.path(vs_dir, "vSpeciateDB_metadata.csv"))

cat("\nOrientación para vSpeciateDB:\n")
print(table(orientation$Orientation))
# Forward Reverse  Review 
#     315      46       4 

cat("\nEjecutar SpeciateIT con:\n")
cat(file.path(vs_dir, "Prevotella_ASVs_vSpeciateDB.fasta"), "\n")
# VaginalMicrobiome/results/vSpeciateDB/Prevotella/Prevotella_ASVs_vSpeciateDB.fasta


###############################################################################
# Reproducibilidad
###############################################################################

software_versions <- data.frame(
  Package = c("R", "phyloseq", "Biostrings", "dplyr", "readr", "stringr"),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("phyloseq")),
    as.character(packageVersion("Biostrings")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("readr")),
    as.character(packageVersion("stringr"))),
  stringsAsFactors = FALSE)

cat("\n=== VERSIONES ===\n")
print(software_versions)
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

###############################################################################
# FIN
###############################################################################

