###############################################################################
#                                    TFM
###############################################################################
# Script 07C: Prevotella - integración vSpeciateDB y taxonomía final
###############################################################################

# Objetivo:
#   Integrar SpeciateIT/vSpeciateDB V4, compararlo con BLASTn, aplicar reglas
#   conservadoras, incorporar prevotella_byVSpeciate y prevotella_refined al
#   phyloseq con una estrategia conservadora.

# Referencia:
#   Holm JB, Gajer P, Ravel J. BMC Bioinformatics. 2024;25:313.
#   DOI: 10.1186/s12859-024-05930-3

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
library(dplyr)
library(readr)
library(stringr)

###############################################################################
# 1. Rutas
###############################################################################

prev_dir <- file.path(project_dir, "results", "Prevotella")
blast_dir <- file.path(prev_dir, "Blast")
vs_dir <- file.path(project_dir, "results", "vSpeciateDB", "Prevotella")

mapping <- read_csv(file.path(prev_dir, "Prevotella_ASV_mapping.csv"), show_col_types = FALSE)
blast <- read_csv(file.path(blast_dir, "Prevotella_species_byBlast.csv"), show_col_types = FALSE)
orientation <- read_csv(file.path(vs_dir, "Prevotella_vSpeciate_orientation.csv"), show_col_types = FALSE)
ps <- readRDS(file.path(project_dir, "cleaned", "ps_prevotella_byBlast.rds"))

vs_file <- file.path(vs_dir, "MC_order7_results.txt")
if (!file.exists(vs_file)) stop("ERROR: falta MC_order7_results.txt. Ejecute SpeciateIT en Ubuntu.")


required_files <- c(
  file.path(prev_dir, "Prevotella_ASV_mapping.csv"),
  file.path(blast_dir, "Prevotella_species_byBlast.csv"),
  file.path(vs_dir, "Prevotella_vSpeciate_orientation.csv"),
  file.path(project_dir, "cleaned", "ps_prevotella_byBlast.rds"),
  vs_file)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Faltan archivos de entrada:",
      paste(missing_files, collapse = ", ")))
}

###############################################################################
# 2. Leer vSpeciateDB
###############################################################################

vspeciate <- read_tsv(vs_file,
  col_names = c("ASV_ID", "vSpeciate_raw", "vSpeciate_Posterior", "vSpeciate_N_Decisions"),
  show_col_types = FALSE, trim_ws = TRUE)

if (anyDuplicated(vspeciate$ASV_ID)) stop("ERROR: vSpeciateDB contiene ASV_ID duplicados.")
if (!setequal(mapping$ASV_ID, vspeciate$ASV_ID)) stop("ERROR: los ASV_ID de vSpeciateDB no coinciden con el mapping.")

vspeciate <- vspeciate %>%
  mutate(vSpeciate_label = str_replace_all(vSpeciate_raw, "_", " "),
    vSpeciate_Level = case_when(str_detect(vSpeciate_raw, "^d_") ~ "Domain",
      str_detect(vSpeciate_raw, "^p_") ~ "Phylum",
      str_detect(vSpeciate_raw, "^c_") ~ "Class",
      str_detect(vSpeciate_raw, "^o_") ~ "Order",
      str_detect(vSpeciate_raw, "^f_") ~ "Family",
      str_detect(vSpeciate_raw, "^g_") ~ "Genus",
      TRUE ~ "Species"),
    species_byVSpeciate = if_else(vSpeciate_Level == "Species", vSpeciate_label, NA_character_))

###############################################################################
# 3. Tabla maestra
###############################################################################

master <- mapping %>%
  left_join(blast %>% dplyr::select(ASV_ID, species_byBlast, Assignment_type), by = "ASV_ID") %>%
  left_join(orientation %>% dplyr::select(ASV_ID, Orientation, Sequence_vSpeciate), by = "ASV_ID") %>%
  left_join(vspeciate, by = "ASV_ID")

vspeciate_in_blast_tie <- function(blast_assignment, vspeciate_assignment) {
  if (is.na(vspeciate_assignment) || !str_detect(blast_assignment, "/")) return(FALSE)
  candidates <- str_split(str_remove(blast_assignment, "^Prevotella\\s+"), "/")[[1]]
  vspecies <- str_remove(vspeciate_assignment, "^Prevotella\\s+")
  vspecies %in% candidates}

master$vSpeciate_in_BLAST_tie <- mapply(vspeciate_in_blast_tie,
  master$species_byBlast, master$species_byVSpeciate)

###############################################################################
# 4. Refinamiento conservador
###############################################################################

master <- master %>%
  mutate(BLAST_is_tie = str_detect(species_byBlast, "/"),
    prevotella_refined = case_when(
      !BLAST_is_tie & species_byBlast != "Prevotella sp." ~ species_byBlast,
      species_byBlast == "Prevotella sp." ~ "Prevotella sp.",
      species_byBlast == "Prevotella jejuni/melaninogenica" &
        species_byVSpeciate == "Prevotella jejuni" &
        vSpeciate_Posterior >= 0.90 ~ "Prevotella jejuni",
      species_byBlast == "Prevotella fusca/melaninogenica" &
        species_byVSpeciate == "Prevotella fusca" &
        vSpeciate_Posterior >= 0.90 ~ "Prevotella fusca",
      species_byBlast == "Prevotella melaninogenica/scopos" ~ "Prevotella melaninogenica/scopos",
      BLAST_is_tie & is.na(species_byVSpeciate) ~ species_byBlast,
      BLAST_is_tie & !vSpeciate_in_BLAST_tie ~ "Prevotella sp.",
      BLAST_is_tie ~ species_byBlast,
      TRUE ~ species_byBlast),
    Refinement_rule = case_when(
      species_byBlast == "Prevotella jejuni/melaninogenica" &
        prevotella_refined == "Prevotella jejuni" ~ "vSpeciate_valid_resolution_jejuni",
      species_byBlast == "Prevotella fusca/melaninogenica" &
        prevotella_refined == "Prevotella fusca" ~ "vSpeciate_valid_resolution_fusca",
      species_byBlast == "Prevotella melaninogenica/scopos" ~ "16S_V4_ambiguity_retained",
      BLAST_is_tie & is.na(species_byVSpeciate) ~ "vSpeciate_unresolved_BLAST_tie_retained",
      BLAST_is_tie & prevotella_refined == "Prevotella sp." ~ "BLAST_vSpeciate_strong_discordance",
      BLAST_is_tie ~ "BLAST_tie_retained",
      species_byBlast == "Prevotella sp." ~ "BLAST_unresolved_retained",
      TRUE ~ "BLAST_assignment_retained"),
    Refined_by_vSpeciate = prevotella_refined != species_byBlast)

###############################################################################
# 5. Incorporar al phyloseq
###############################################################################

tax_final <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
idx <- match(master$Sequence, rownames(tax_final))
if (anyNA(idx)) stop("ERROR: no se pudieron mapear todos los ASVs de Prevotella.")

tax_final$prevotella_byVSpeciate <- NA_character_
tax_final$prevotella_refined <- NA_character_
tax_final$prevotella_byVSpeciate[idx] <- master$species_byVSpeciate
tax_final$prevotella_refined[idx] <- master$prevotella_refined

# Crear taxonomía final integrada
tax_final$species_refined <- dplyr::case_when(
  tax_final$Genus == "Lactobacillus" &
    !is.na(tax_final$lactobacillus_refined) ~ tax_final$lactobacillus_refined,
  tax_final$Genus == "Gardnerella" &
    !is.na(tax_final$gardnerella_refined) ~ tax_final$gardnerella_refined,
  tax_final$Genus == "Prevotella" &
    !is.na(tax_final$prevotella_refined) ~ tax_final$prevotella_refined,
  !is.na(tax_final$Species) ~
    paste(tax_final$Genus, tax_final$Species),
  !is.na(tax_final$Genus) ~
    paste0(tax_final$Genus, " sp."),
  TRUE ~ NA_character_)

tax_table(ps) <- tax_table(as.matrix(tax_final))
saveRDS(ps, file.path(project_dir, "cleaned", "ps_prevotella_refined.rds"))

# Exploración de species_refined
ps <- readRDS(file.path(project_dir, "cleaned", "ps_prevotella_refined.rds"))
tax <- as.data.frame(tax_table(ps))

# ASVs sin resolución al menos a género
unresolved <- is.na(tax_final$species_refined)

cat("ASVs sin species_refined:", sum(unresolved), "\n") # 1343
cat("ASVs con species_refined:", sum(!unresolved), "\n")# 4994
cat("Porcentaje ASVs sin resolver:",
    round(100 * sum(unresolved) / nrow(tax_final), 2), "%\n") # 21.19 %

# Comprobar qué proporción de las lecturas representan
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)

reads_unresolved <- sum(otu[unresolved, , drop = FALSE])
reads_total <- sum(otu)

cat("Lecturas de ASVs sin resolver:", reads_unresolved, "\n") #  2529956
cat("Lecturas totales:", reads_total, "\n") # 344812265
cat("Porcentaje de lecturas sin resolver:",
    round(100 * reads_unresolved / reads_total, 4), "%\n") # 0.7337 %

# Prevalencia: en cuántas muestras aparecen
prev_unresolved <- rowSums(otu[unresolved, , drop = FALSE] > 0)

summary(prev_unresolved)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.000   1.000   7.509   2.000 854.000 

unresolved_info <- tax_final[unresolved, , drop = FALSE]
unresolved_info$Prevalence <- prev_unresolved
unresolved_info$Reads <- rowSums(otu[unresolved, , drop = FALSE])

unresolved_info %>%
  arrange(desc(Prevalence)) %>%
  dplyr::select(Kingdom, Phylum, Class, Order, Family,
                Genus, Species, Prevalence, Reads) %>%
  head(10)

###############################################################################
# 6. Guardar tablas y resumen
###############################################################################

master_final <- master %>%
  dplyr::select(ASV_ID, Sequence, species_byBlast, species_byVSpeciate,
    vSpeciate_Posterior, vSpeciate_N_Decisions, vSpeciate_Level, Orientation,
    prevotella_refined, Refined_by_vSpeciate, Refinement_rule)

ties_final <- master %>%
  filter(BLAST_is_tie) %>%
  count(species_byBlast, species_byVSpeciate, prevotella_refined,
    Refinement_rule, name = "N_ASVs", sort = TRUE)

summary_table <- data.frame(
  Metric = c("N_Prevotella_ASVs", "BLAST_species_or_tie_assigned",
    "vSpeciate_species_assigned", "Final_species_or_tie_assigned", "Changed_by_vSpeciate"),
  Value = c(nrow(master), sum(master$species_byBlast != "Prevotella sp."),
    sum(!is.na(master$species_byVSpeciate)),
    sum(master$prevotella_refined != "Prevotella sp."),
    sum(master$Refined_by_vSpeciate, na.rm = TRUE)),
  stringsAsFactors = FALSE)

write_csv(master_final, file.path(vs_dir, "Prevotella_taxonomy_master_final.csv"))
write_csv(ties_final, file.path(vs_dir, "Prevotella_ties_final.csv"))
write_csv(summary_table, file.path(vs_dir, "Prevotella_refinement_summary.csv"))

cat("\nRESUMEN TAXONÓMICO FINAL\n")
print(summary_table)
# RESUMEN TAXONÓMICO FINAL

#                         Metric Value
#1             N_Prevotella_ASVs   365
#2 BLAST_species_or_tie_assigned   112
#3    vSpeciate_species_assigned   251
#4 Final_species_or_tie_assigned   106
#5          Changed_by_vSpeciate     8

cat("\nReglas aplicadas:\n")
print(master %>% count(Refinement_rule, sort = TRUE), n = Inf)
# Reglas aplicadas:
# A tibble: 7 × 2
#  Refinement_rule                        n
#   <chr>                              <int>
#1 BLAST_unresolved_retained            253
#2 BLAST_assignment_retained            101
#3 BLAST_vSpeciate_strong_discordance     6
#4 16S_V4_ambiguity_retained              2
#5 BLAST_tie_retained                     1
#6 vSpeciate_valid_resolution_fusca       1
#7 vSpeciate_valid_resolution_jejuni      1

cat("\nClasificación final:\n")
print(sort(table(master$prevotella_refined), decreasing = TRUE))
# Clasificación final:
#Prevotella sp.                    Prevotella bivia 
#259                                  23 
#Prevotella melaninogenica         Prevotella corporis 
#16                                  15 
#Prevotella amnii                  Prevotella disiens 
#7                                   6 
#Prevotella denticola               Prevotella intermedia 
#5                                   5 
#Prevotella massiliensis            Prevotella pallens 
#5                                   5 
#Prevotella brunnea              Prevotella multiformis 
#3                                   3 
#Prevotella melaninogenica/scopos   Prevotella dentalis 
#2                                   1 
#Prevotella fusca                   Prevotella ihumii 
#1                                   1 
#Prevotella illustrans              Prevotella jejuni 
#1                                   1 
#Prevotella koreensis Prevotella melaninogenica/veroralis 
#1                                   1 
#Prevotella merdae                   Prevotella micans 
#1                                   1 
#Prevotella nigrescens              Prevotella phocaeensis 
#1                                   1 


###############################################################################
# Reproducibilidad
###############################################################################

software_versions <- data.frame(
  Package = c("R", "phyloseq", "dplyr", "readr", "stringr"),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("phyloseq")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("readr")),
    as.character(packageVersion("stringr"))),
  stringsAsFactors = FALSE)

cat("\n=== VERSIONES ===\n")
print(software_versions)
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

