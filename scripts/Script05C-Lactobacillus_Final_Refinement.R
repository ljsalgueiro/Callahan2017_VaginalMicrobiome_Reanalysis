###############################################################################
#                                    TFM
###############################################################################
# Script 05C: Lactobacillus - integración vSpeciateDB y taxonomía final
###############################################################################

# Objetivo:
#   1. Integrar los resultados de SpeciateIT/vSpeciateDB V4.
#   2. Compararlos con BLASTn.
#   3. Aplicar reglas conservadoras de refinamiento.
#   4. Preservar SILVA, lactobacillus_byBlast y lactobacillus_byVSpeciate.
#   5. Incorporar lactobacillus_refined al objeto phyloseq.

# Referencia:
#   Holm JB, Gajer P, Ravel J. SpeciateIT and vSpeciateDB: novel, fast, and
#   accurate per sequence 16S rRNA gene taxonomic classification of vaginal
#   microbiota. BMC Bioinformatics. 2024;25:313.
#   DOI: 10.1186/s12859-024-05930-3

# Entrada:
#   - ps_lacto_byBlast.rds
#   - Lactobacillus_ASV_mapping.csv
#   - species_byBlast_for_phyloseq.csv
#   - Lactobacillus_vSpeciate_orientation.csv
#   - MC_order7_results.txt

# Salidas:
#   - Lactobacillus_taxonomy_master_final.csv
#   - Lactobacillus_refinement_summary.csv
#   - ps_lacto_refined.rds

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

base_dir <- file.path(project_dir, "results", "Lactobacillus")
blast_dir <- file.path(base_dir, "Blast")
vs_dir <- file.path(project_dir, "results", "vSpeciateDB", "Lactobacillus")

ps <- readRDS(file.path(project_dir, "cleaned", "ps_lacto_byBlast.rds"))
mapping <- read_csv(file.path(base_dir, "Lactobacillus_ASV_mapping.csv"), show_col_types = FALSE)
blast <- read_csv(file.path(blast_dir, "species_byBlast_for_phyloseq.csv"), show_col_types = FALSE)
orientation <- read_csv(file.path(vs_dir, "Lactobacillus_vSpeciate_orientation.csv"), show_col_types = FALSE)


required_files <- c(
  file.path(project_dir, "cleaned", "ps_lacto_byBlast.rds"),
  file.path(base_dir, "Lactobacillus_ASV_mapping.csv"),
  file.path(blast_dir, "species_byBlast_for_phyloseq.csv"),
  file.path(vs_dir, "Lactobacillus_vSpeciate_orientation.csv"),
  file.path(vs_dir, "MC_order7_results.txt"))

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Faltan archivos de entrada:",
      paste(missing_files, collapse = ", ")))
}

###############################################################################
# 2. Leer y normalizar vSpeciateDB
###############################################################################

vspeciate <- read_tsv(file.path(vs_dir, "MC_order7_results.txt"),
  col_names = c("ASV_ID", "vSpeciate_raw", "vSpeciate_Posterior", "vSpeciate_N_Decisions"),
  show_col_types = FALSE, trim_ws = TRUE)

if (anyDuplicated(vspeciate$ASV_ID)) stop("ERROR: vSpeciateDB contiene ASV_ID duplicados.")
if (!setequal(mapping$ASV_ID, vspeciate$ASV_ID)) stop("ERROR: los ASV_ID de vSpeciateDB no coinciden con el mapping.")

vspeciate <- vspeciate %>%
  mutate(vSpeciate_label = str_replace_all(vSpeciate_raw, "_", " "),
    vSpeciate_Level = case_when(str_detect(vSpeciate_raw, "^d_") ~ "Domain",
      str_detect(vSpeciate_raw, "^p_") ~ "Phylum", str_detect(vSpeciate_raw, "^c_") ~ "Class",
      str_detect(vSpeciate_raw, "^o_") ~ "Order", str_detect(vSpeciate_raw, "^f_") ~ "Family",
      str_detect(vSpeciate_raw, "^g_") ~ "Genus", TRUE ~ "Species"),
    species_byVSpeciate = if_else(vSpeciate_Level == "Species", vSpeciate_label, NA_character_))

###############################################################################
# 3. Tabla maestra SILVA + BLAST + vSpeciateDB
###############################################################################

tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
tax$Sequence <- rownames(tax)
silva <- tax %>%
  dplyr::select(Sequence, Genus, Species) %>%
  dplyr::rename(Species_SILVA = Species)

master <- mapping %>%
  left_join(silva, by = "Sequence") %>%
  left_join(blast, by = "ASV_ID") %>%
  left_join(orientation %>% select(ASV_ID, Orientation, Sequence_vSpeciate), by = "ASV_ID") %>%
  left_join(vspeciate, by = "ASV_ID") %>%
  mutate(BLAST_vs_vSpeciate = case_when(is.na(species_byVSpeciate) ~ "vSpeciate_unresolved",
    species_byBlast == species_byVSpeciate ~ "Exact_concordance",
    species_byBlast == "Lactobacillus sp." ~ "vSpeciate_more_resolved",
    TRUE ~ "Different_assignment"))

###############################################################################
# 4. Reglas conservadoras de refinamiento
###############################################################################

master <- master %>%
  mutate(species_refined = case_when(
    !str_detect(species_byBlast, "/") & !str_detect(species_byBlast, " group$") &
      species_byBlast != "Lactobacillus sp." ~ species_byBlast,
    species_byBlast == "Lactobacillus jensenii/mulieris" &
      species_byVSpeciate == "Lactobacillus jensenii" ~ "Lactobacillus jensenii",
    species_byBlast == "Lactobacillus jensenii/mulieris" &
      species_byVSpeciate == "Lactobacillus mulieris" ~ "Lactobacillus mulieris",
    species_byBlast == "Lactobacillus jensenii/mulieris" ~ "Lactobacillus jensenii/mulieris",
    species_byBlast == "Lactobacillus acidophilus/crispatus" ~ "Lactobacillus acidophilus/crispatus",
    species_byBlast == "Lactobacillus gasseri group" ~ "Lactobacillus gasseri group",
    species_byBlast == "Lactobacillus gasseri/johnsonii group" ~ "Lactobacillus gasseri/johnsonii group",
    species_byBlast == "Lactobacillus acidophilus/amylovorus" ~ "Lactobacillus acidophilus/amylovorus",
    species_byBlast == "Lactobacillus acidophilus group" &
      species_byVSpeciate == "Lacticaseibacillus paracasei" ~ "Lactobacillus sp.",
    species_byBlast == "Lactobacillus acidophilus/delbrueckii" &
      species_byVSpeciate == "Lacticaseibacillus paracasei" ~ "Lactobacillus sp.",
    species_byBlast == "Lactobacillus acidophilus/delbrueckii" ~ "Lactobacillus acidophilus/delbrueckii",
    species_byBlast == "Lactobacillus delbrueckii/helveticus" &
      species_byVSpeciate == "Limosilactobacillus fermentum" ~ "Lactobacillus sp.",
    species_byBlast == "Lactobacillus delbrueckii/helveticus" ~ "Lactobacillus delbrueckii/helveticus",
    species_byBlast == "Lactobacillus sp." ~ "Lactobacillus sp.",
    TRUE ~ species_byBlast),
    Refinement_rule = case_when(
      species_byBlast == "Lactobacillus jensenii/mulieris" &
        species_refined == "Lactobacillus jensenii" ~ "vSpeciate_valid_resolution_jensenii",
      species_byBlast == "Lactobacillus jensenii/mulieris" &
        species_refined == "Lactobacillus mulieris" ~ "vSpeciate_valid_resolution_mulieris",
      species_byBlast == "Lactobacillus acidophilus/crispatus" ~ "V4_indistinguishable_keep_BLAST_ambiguity",
      species_byBlast == "Lactobacillus gasseri group" ~ "V4_gasseri_paragasseri_not_resolved",
      species_byBlast == "Lactobacillus gasseri/johnsonii group" ~ "insufficient_resolution_keep_group",
      species_byBlast == "Lactobacillus acidophilus/amylovorus" ~ "vSpeciate_not_valid_for_this_tie_keep_BLAST",
      species_byBlast == "Lactobacillus acidophilus group" &
        species_refined == "Lactobacillus sp." ~ "BLAST_vSpeciate_strong_discordance",
      species_byBlast == "Lactobacillus acidophilus/delbrueckii" &
        species_refined == "Lactobacillus sp." ~ "BLAST_vSpeciate_strong_discordance",
      species_byBlast == "Lactobacillus delbrueckii/helveticus" &
        species_refined == "Lactobacillus sp." ~ "BLAST_vSpeciate_strong_discordance",
      species_byBlast == "Lactobacillus sp." ~ "BLAST_unresolved_retained",
      species_refined == species_byBlast ~ "BLAST_assignment_retained",
      TRUE ~ "other"),
    Refined_by_vSpeciate = species_refined != species_byBlast)

if (anyNA(master$species_refined)) stop("ERROR: existen ASVs sin species_refined.")

###############################################################################
# 5. Incorporar taxonomía final al phyloseq
###############################################################################

ps_final <- ps
tax_final <- as.data.frame(tax_table(ps_final), stringsAsFactors = FALSE)
idx <- match(master$Sequence, rownames(tax_final))
if (anyNA(idx)) stop("ERROR: no se pudieron mapear todos los ASVs de Lactobacillus al phyloseq.")

tax_final$lactobacillus_byVSpeciate <- NA_character_
tax_final$lactobacillus_refined <- NA_character_
tax_final$lactobacillus_byVSpeciate[idx] <- master$species_byVSpeciate
tax_final$lactobacillus_refined[idx] <- master$species_refined
tax_table(ps_final) <- tax_table(as.matrix(tax_final))

###############################################################################
# 6. Resumen y guardado
###############################################################################

summary_table <- data.frame(
  Metric = c("N_Lactobacillus_ASVs", "SILVA_species_assigned", "BLAST_species_or_group_assigned",
    "vSpeciate_species_assigned", "Final_species_or_group_assigned", "Changed_by_vSpeciate"),
  Value = c(nrow(master), sum(!is.na(master$Species_SILVA)),
    sum(master$species_byBlast != "Lactobacillus sp."),
    sum(!is.na(master$species_byVSpeciate)),
    sum(master$species_refined != "Lactobacillus sp."),
    sum(master$Refined_by_vSpeciate, na.rm = TRUE)),
  stringsAsFactors = FALSE)

write_csv(master, file.path(vs_dir, "Lactobacillus_taxonomy_master_final.csv"))
write_csv(summary_table, file.path(vs_dir, "Lactobacillus_refinement_summary.csv"))
saveRDS(ps_final, file.path(project_dir, "cleaned", "ps_lacto_refined.rds"))

cat("\nREFINAMIENTO FINAL COMPLETADO\n")
print(summary_table)
#                            Metric Value
# 1            N_Lactobacillus_ASVs   202
# 2          SILVA_species_assigned    19
# 3 BLAST_species_or_group_assigned    88
# 4      vSpeciate_species_assigned   188
# 5 Final_species_or_group_assigned    83
# 6            Changed_by_vSpeciate    17

cat("\nClasificación final:\n")
print(sort(table(master$species_refined), decreasing = TRUE))
#   Lactobacillus sp.                 Lactobacillus iners 
#                 119                                  41 
#Lactobacillus gasseri              Lactobacillus mulieris 
#                 17                                   9 
#Lactobacillus crispatus              Lactobacillus jensenii 
#                     8                                   3 
#Lactobacillus acidophilus Lactobacillus acidophilus/crispatus 
#                       2                                   1 
#Lactobacillus delbrueckii            Lactobacillus kalixensis 
#                       1                                   1 

master %>%
  filter(Refined_by_vSpeciate) %>%
  dplyr::select(
    ASV_ID,
    Species_SILVA,
    species_byBlast,
    species_byVSpeciate,
    vSpeciate_Posterior,
    species_refined,
    Refinement_rule) %>%
  arrange(Refinement_rule) %>%
  print(n = Inf)

master %>%
  count(Refinement_rule, sort = TRUE) %>%
  print(n = Inf)


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

