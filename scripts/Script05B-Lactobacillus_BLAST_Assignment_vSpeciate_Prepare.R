###############################################################################
#                                    TFM
###############################################################################
# Script 05B: Lactobacillus - asignación por BLASTn y preparación vSpeciateDB
###############################################################################

# Objetivo:
#   1. Aplicar filtros Query Cover >=95% y Percent Identity >=99%.
#   2. Seleccionar el Best Hit Set de forma lexicográfica: cobertura > identidad > E-value.
#   3. Aplicar reglas conservadoras de asignación taxonómica.
#   4. Incorporar species_byBlast al phyloseq preservando SILVA.
#   5. Revisar orientación y generar un FASTA ya orientado para vSpeciateDB.
#
# Entrada:
#   - ps_taxFiltered.rds
#   - Lactobacillus_ASV_mapping.csv
#   - BLAST_results_all.csv
#
# Salidas:
#   - BLAST_species_byBlast_final.csv
#   - species_byBlast_for_phyloseq.csv
#   - ps_lacto_byBlast.rds
#   - Lactobacillus_vSpeciate_orientation.csv
#   - Lactobacillus_ASVs_vSpeciateDB.fasta
#   - vSpeciateDB_metadata.csv

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

input_ps <- file.path(project_dir, "cleaned", "ps_taxFiltered.rds")
base_dir <- file.path(project_dir, "results", "Lactobacillus")
blast_dir <- file.path(base_dir, "Blast")
vs_dir <- file.path(project_dir, "results", "vSpeciateDB", "Lactobacillus")
dir.create(vs_dir, recursive = TRUE, showWarnings = FALSE)

MIN_QUERY_COVER <- 95
MIN_PER_IDENT <- 99
SPECIAL_DOMINANCE_THRESHOLD <- 90


required_files <- c(
  input_ps,
  file.path(base_dir, "Lactobacillus_ASV_mapping.csv"),
  file.path(blast_dir, "BLAST_results_all.csv"))

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

ps_original <- readRDS(input_ps)
mapping <- read_csv(file.path(base_dir, "Lactobacillus_ASV_mapping.csv"), show_col_types = FALSE)
blast <- read_csv(file.path(blast_dir, "BLAST_results_all.csv"), show_col_types = FALSE)
all_asvs <- mapping %>% distinct(ASV_ID)

required <- c("ASV_ID", "Hit_Rank", "Accession", "Description", "Scientific_Name",
  "Query_Cover", "Per_Ident", "E_value", "Bit_Score")
if (!all(required %in% names(blast))) stop("ERROR: faltan columnas necesarias en BLAST_results_all.csv.")

###############################################################################
# 3. Normalización y filtros BLAST
###############################################################################

blast <- blast %>%
  mutate(Scientific_Name_original = Scientific_Name, Scientific_Name = str_squish(Scientific_Name),
    Scientific_Name = case_when(Scientific_Name == "Lactobacillus leichmannii" ~ "Lactobacillus delbrueckii",
      Scientific_Name == "Lactobacillus mulieri" ~ "Lactobacillus mulieris", TRUE ~ Scientific_Name),
    Scientific_Name = if_else(str_detect(Scientific_Name, "^Lactobacillus iners\\s"), "Lactobacillus iners", Scientific_Name),
    Scientific_Name = str_replace(Scientific_Name, "^(Lactobacillus\\s+[^\\s]+)\\s+subsp\\..*$", "\\1"),
    Taxon_type = case_when(is.na(Scientific_Name) | Scientific_Name == "" ~ "unresolved",
      str_detect(str_to_lower(Scientific_Name), "^uncultured|^unidentified") ~ "unresolved",
      str_detect(Scientific_Name, "^Lactobacillus sp\\.") ~ "unresolved",
      Scientific_Name == "Lactobacillus casei group sp." ~ "unresolved",
      str_detect(Scientific_Name, "^Lactobacillus\\s+[a-z][a-z-]+$") ~ "species", TRUE ~ "unresolved"))

eligible_hits <- blast %>%
  filter(!is.na(Query_Cover), !is.na(Per_Ident), !is.na(E_value),
    Query_Cover >= MIN_QUERY_COVER, Per_Ident >= MIN_PER_IDENT)

best_hits <- eligible_hits %>%
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

best_species_hits <- best_hits %>%
  filter(Taxon_type == "species") %>%
  count(ASV_ID, Scientific_Name, name = "N_hits_species") %>%
  group_by(ASV_ID) %>%
  mutate(N_species_hits_total = sum(N_hits_species), Percent_hits_species = 100 * N_hits_species / N_species_hits_total) %>%
  ungroup()

best_species_groups <- best_species_hits %>%
  mutate(Species_group_test = if_else(Scientific_Name %in% c("Lactobacillus jensenii", "Lactobacillus mulieris"),
    "Lactobacillus jensenii/mulieris", Scientific_Name)) %>%
  group_by(ASV_ID, Species_group_test) %>% summarise(N_hits_group = sum(N_hits_species), .groups = "drop") %>%
  group_by(ASV_ID) %>% mutate(N_species_hits_total = sum(N_hits_group),
    Percent_hits_group = 100 * N_hits_group / N_species_hits_total) %>% ungroup()

jm_stats <- best_species_groups %>%
  filter(Species_group_test == "Lactobacillus jensenii/mulieris") %>%
  select(ASV_ID, Jensenii_Mulieris_N_hits = N_hits_group, Jensenii_Mulieris_Percent = Percent_hits_group)

best_summary <- best_hits %>%
  group_by(ASV_ID) %>%
  summarise(Best_Query_Cover = dplyr::first(Query_Cover), Best_Per_Ident = dplyr::first(Per_Ident),
    Best_E_value = dplyr::first(E_value), N_best_hits = n(),
    N_valid_species = n_distinct(Scientific_Name[Taxon_type == "species"], na.rm = TRUE),
    Species_candidates = {x <- sort(unique(Scientific_Name[Taxon_type == "species"])); x <- x[!is.na(x)];
      if (length(x) == 0) NA_character_ else paste(x, collapse = " | ")},
    Has_unresolved_best_hits = any(Taxon_type == "unresolved"), .groups = "drop") %>%
  left_join(jm_stats, by = "ASV_ID") %>% mutate(Pattern = Species_candidates)

###############################################################################
# 4. Reglas de asignación BLAST
###############################################################################

pattern_fgJM <- paste(c("Lactobacillus fornicalis", "Lactobacillus gasseri",
  "Lactobacillus jensenii", "Lactobacillus mulieris"), collapse = " | ")
pattern_gasseri <- paste(c("Lactobacillus gasseri", "Lactobacillus hominis",
  "Lactobacillus johnsonii", "Lactobacillus paragasseri"), collapse = " | ")
pattern_acidophilus <- paste(c("Lactobacillus acidophilus", "Lactobacillus delbrueckii",
  "Lactobacillus gasseri", "Lactobacillus helveticus"), collapse = " | ")
pattern_conflict <- paste(c("Lactobacillus crispatus", "Lactobacillus gallinarum",
  "Lactobacillus jensenii"), collapse = " | ")

best_summary <- best_summary %>%
  mutate(species_byBlast = case_when(
    N_valid_species == 0 ~ "Lactobacillus sp.",
    N_valid_species == 1 ~ Species_candidates,
    Pattern == pattern_fgJM & !is.na(Jensenii_Mulieris_Percent) &
      Jensenii_Mulieris_Percent >= SPECIAL_DOMINANCE_THRESHOLD ~ "Lactobacillus jensenii/mulieris",
    Pattern == "Lactobacillus acidophilus | Lactobacillus crispatus" ~ "Lactobacillus acidophilus/crispatus",
    Pattern == pattern_gasseri ~ "Lactobacillus gasseri group",
    Pattern == "Lactobacillus delbrueckii | Lactobacillus helveticus" ~ "Lactobacillus delbrueckii/helveticus",
    Pattern == "Lactobacillus acidophilus | Lactobacillus delbrueckii" ~ "Lactobacillus acidophilus/delbrueckii",
    Pattern == pattern_acidophilus ~ "Lactobacillus acidophilus group",
    Pattern == "Lactobacillus acidophilus | Lactobacillus amylovorus" ~ "Lactobacillus acidophilus/amylovorus",
    Pattern == pattern_conflict ~ "Lactobacillus sp.",
    Pattern == "Lactobacillus gasseri | Lactobacillus johnsonii" ~ "Lactobacillus gasseri/johnsonii group",
    Pattern == "Lactobacillus jensenii | Lactobacillus mulieris" ~ "Lactobacillus jensenii/mulieris",
    TRUE ~ "Lactobacillus sp."),
    Assignment_Level = case_when(species_byBlast == "Lactobacillus sp." ~ "Genus_only",
      str_detect(species_byBlast, " group$") ~ "Species_group",
      str_detect(species_byBlast, "/") ~ "Species_ambiguous", TRUE ~ "Species"),
    Assignment_Rule = case_when(N_valid_species == 0 ~ "no_species_resolved_among_best_hits",
      N_valid_species == 1 ~ "unique_species_among_best_hits",
      Pattern == pattern_fgJM & species_byBlast == "Lactobacillus jensenii/mulieris" ~ "special_dominant_jensenii_mulieris",
      species_byBlast == "Lactobacillus gasseri group" ~ "A_species_group_gasseri",
      species_byBlast == "Lactobacillus acidophilus group" ~ "A_species_group_acidophilus",
      species_byBlast == "Lactobacillus gasseri/johnsonii group" ~ "A_species_group_gasseri_johnsonii",
      Assignment_Level == "Species_ambiguous" ~ "B_V4_species_ambiguity",
      species_byBlast == "Lactobacillus sp." & N_valid_species > 1 ~ "C_unresolved_species_conflict",
      TRUE ~ "other"))

final_assignment <- all_asvs %>%
  left_join(best_summary, by = "ASV_ID") %>%
  mutate(species_byBlast = if_else(is.na(species_byBlast), "Lactobacillus sp.", species_byBlast),
    Assignment_Level = if_else(is.na(Assignment_Level), "Genus_only", Assignment_Level),
    Assignment_Rule = if_else(is.na(Assignment_Rule), "no_hit_passed_quality_thresholds", Assignment_Rule),
    ASV_number = as.numeric(str_remove(ASV_ID, "^ASV_"))) %>%
  arrange(ASV_number) %>% select(-ASV_number)

write_csv(final_assignment, file.path(blast_dir, "BLAST_species_byBlast_final.csv"))
write_csv(final_assignment %>% select(ASV_ID, species_byBlast),
  file.path(blast_dir, "species_byBlast_for_phyloseq.csv"))

###############################################################################
# 5. Incorporar species_byBlast al phyloseq
###############################################################################

blast_mapping <- mapping %>% left_join(final_assignment %>% select(ASV_ID, species_byBlast), by = "ASV_ID")
if (anyNA(blast_mapping$species_byBlast)) stop("ERROR: faltan asignaciones species_byBlast.")

ps_blast <- ps_original
tax_matrix <- as(tax_table(ps_blast), "matrix")
lactobacillus_byBlast_vector <- rep(NA_character_, nrow(tax_matrix))
names(lactobacillus_byBlast_vector) <- rownames(tax_matrix)
lookup <- setNames(blast_mapping$species_byBlast, blast_mapping$Sequence)
matching <- intersect(rownames(tax_matrix), names(lookup))
lactobacillus_byBlast_vector[matching] <- lookup[matching]
tax_table(ps_blast) <- tax_table(
  cbind(tax_matrix, lactobacillus_byBlast = lactobacillus_byBlast_vector))
saveRDS(ps_blast, file.path(project_dir, "cleaned", "ps_lacto_byBlast.rds"))

###############################################################################
# 6. Orientación previa a vSpeciateDB
###############################################################################

orientation <- mapping %>%
  mutate(Reverse_Complement = as.character(reverseComplement(DNAStringSet(Sequence))),
    Forward_start = substr(Sequence, 1, 10), RevComp_start = substr(Reverse_Complement, 1, 10),
    Orientation = case_when(str_detect(Sequence, "^TACGT") ~ "Forward",
      str_detect(Reverse_Complement, "^TACGT") ~ "Reverse", TRUE ~ "Review"),
    Sequence_vSpeciate = if_else(Orientation == "Reverse", Reverse_Complement, Sequence))

write_csv(orientation, file.path(vs_dir, "Lactobacillus_vSpeciate_orientation.csv"))

vs_seqs <- DNAStringSet(orientation$Sequence_vSpeciate)
names(vs_seqs) <- orientation$ASV_ID
writeXStringSet(vs_seqs, filepath = file.path(vs_dir, "Lactobacillus_ASVs_vSpeciateDB.fasta"), format = "fasta")

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

cat("\nASIGNACIÓN BLAST COMPLETADA\n")
cat("ASVs Lactobacillus:", nrow(mapping), "\n") # 202
cat("Asignaciones:\n")
print(sort(table(final_assignment$species_byBlast), decreasing = TRUE))
cat("\nOrientación para vSpeciateDB:\n")
print(table(orientation$Orientation))
# Forward Reverse  Review 
#     172      23       7 
cat("\nEjecutar vSpeciateDB una sola vez con:\n")
cat(file.path(vs_dir, "Lactobacillus_ASVs_vSpeciateDB.fasta"), "\n")
cat("Guardar el resultado como:\n")
cat(file.path(vs_dir, "MC_order7_results.txt"), "\n")


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

