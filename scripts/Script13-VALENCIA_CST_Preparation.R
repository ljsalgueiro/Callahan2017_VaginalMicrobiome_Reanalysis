###############################################################################
#                                    TFM
###############################################################################
#             Script 13: Preparación e integración de VALENCIA
###############################################################################

# Objetivo:
# Preparar desde el objeto phyloseq definitivo la tabla de conteos requerida
# por VALENCIA, auditar el mapeo contra los centroides oficiales e integrar
# CST, sub-CST y score cuando VALENCIA_output.csv esté disponible.
#
# Armonización taxonómica:
# - conservar species_refined cuando coincide exactamente con la referencia;
# - usar SILVA a especie solo si species_refined está realmente ausente;
# - Gardnerella G1/G2/G3/G_other se agrega a Gardnerella_vaginalis solo para
#   compatibilidad con VALENCIA, sin modificar species_refined;
# - Lactobacillus, Prevotella, Atopobium y Sneathia sin especie compatible se
#   agregan al género si dicho género existe en la referencia;
# - el resto se mapea jerárquicamente a género, familia, orden, clase o filo;
# - si no existe correspondencia, se asigna a d_NA.
#
# El script no calcula los CST en R. Genera VALENCIA_input.csv para ejecutar
# Valencia.py en Python. Si el output ya existe, crea ps_valencia.rds.

# ============================================================
# VALENCIA reproducibility information
# ============================================================
#
# VALENCIA repository:
# https://github.com/ravel-lab/VALENCIA
#
# VALENCIA commit:
# 8559d454387479f7155333693d854961463c3b15
#
# Valencia.py MD5:
# 46f695e12331be15ad303602ab9fbf7e
#
# Reference centroids:
# CST_centroids_012920.csv
#
# Centroids MD5:
# 40029ff8c590f919a8506ad3bff4269f
#
# Execution environment:
# Ubuntu 20.04.6 LTS
# Python 3.8.10
# pandas 2.0.3
# numpy 1.24.4
# ============================================================

###############################################################################

rm(list = ls())

library(phyloseq)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

project_dir <- if (
  dir.exists("VaginalMicrobiome")) {
  "VaginalMicrobiome"
} else {
  "."
}

###############################################################################
# 1. Rutas
###############################################################################

ps_file <- file.path(
  project_dir,
  "cleaned",
  "ps_prevotella_refined.rds")

valencia_dir <- file.path(
  project_dir,
  "results",
  "VALENCIA")

ref_file <- file.path(
  project_dir,
  "databases",
  "VALENCIA",
  "CST_centroids_012920.csv")

input_file <- file.path(
  valencia_dir,
  "VALENCIA_input.csv")

mapping_file <- file.path(
  valencia_dir,
  "VALENCIA_ASV_mapping.csv")

metadata_file <- file.path(
  valencia_dir,
  "VALENCIA_metadata.csv")

output_file <- file.path(
  valencia_dir,
  "VALENCIA_output.csv")

ps_output_file <- file.path(
  project_dir,
  "cleaned",
  "ps_valencia.rds")

dir.create(
  valencia_dir,
  recursive = TRUE,
  showWarnings = FALSE)

###############################################################################
# 2. Comprobar archivos y cargar datos
###############################################################################

if (!file.exists(ps_file))
  stop("No se encontró ps_prevotella_refined.rds.")

if (!file.exists(ref_file))
  stop("No se encontró CST_centroids_012920.csv.")

ps <- readRDS(ps_file)

centroids <- read_csv(
  ref_file,
  show_col_types = FALSE)

if (!"sub_CST" %in% names(centroids))
  stop("El archivo de centroides no contiene sub_CST.")

ref_taxa <- setdiff(
  names(centroids),
  "sub_CST")

if (!"d_NA" %in% ref_taxa)
  stop("La referencia de VALENCIA no contiene d_NA.")

###############################################################################
# 3. Conservar muestras clínicas
###############################################################################

meta <- data.frame(
  sample_data(ps),
  stringsAsFactors = FALSE)

meta$sampleID <- rownames(meta)

if ("Groups" %in% names(meta)) {

  keep_samples <- rownames(meta)[
    meta$Groups %in% c(
      "ST","SP","UT","UP")]

  ps <- prune_samples(
    keep_samples,
    ps)

  meta <- data.frame(
    sample_data(ps),
    stringsAsFactors = FALSE)

  meta$sampleID <- rownames(meta)
}

cat("\n=== OBJETO DE ENTRADA ===\n")
cat("Muestras clínicas:",nsamples(ps),"\n") # 2179
cat("ASVs:",ntaxa(ps),"\n") # 6337
cat("Filotipos VALENCIA:",length(ref_taxa),"\n")

###############################################################################
# 4. Extraer conteos y taxonomía
###############################################################################

otu <- as(
  otu_table(ps),
  "matrix")

if (!taxa_are_rows(ps))
  otu <- t(otu)

tax <- as.data.frame(
  tax_table(ps),
  stringsAsFactors = FALSE)

tax <- tax[
  match(
    rownames(otu),
    rownames(tax)),
  ,
  drop = FALSE]

stopifnot(
  identical(
    rownames(otu),
    rownames(tax)))

required_tax <- c(
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species",
  "species_refined")

missing_tax <- setdiff(
  required_tax,
  names(tax))

if (length(missing_tax) > 0)
  stop(
    paste(
      "Faltan columnas taxonómicas:",
      paste(
        missing_tax,
        collapse = ", ")))

###############################################################################
# 5. Funciones de normalización y matching
###############################################################################

normalize_name <- function(x) {

  x <- ifelse(
    is.na(x),
    "",
    x)

  tolower(
    gsub(
      "[^A-Za-z0-9]",
      "",
      x))
}

ref_norm <- normalize_name(
  ref_taxa)

if (anyDuplicated(ref_norm)) {

  duplicated_names <- ref_taxa[
    duplicated(ref_norm) |
      duplicated(
        ref_norm,
        fromLast = TRUE)]

  stop(
    paste(
      "Nombres VALENCIA ambiguos tras normalización:",
      paste(
        duplicated_names,
        collapse = ", ")))
}

match_reference <- function(candidate) {

  if (
    length(candidate) == 0 ||
      is.na(candidate) ||
      candidate == "")
    return(NA_character_)

  idx <- match(
    normalize_name(candidate),
    ref_norm)

  if (is.na(idx))
    return(NA_character_)

  ref_taxa[idx]
}

is_ambiguous_species <- function(x) {

  if (is.na(x) || x == "")
    return(TRUE)

  str_detect(
    x,
    "/|\\bsp\\.$|\\bsp\\b|group$|uncultured|unclassified")
}

###############################################################################
# 6. Mapear cada ASV a un filotipo VALENCIA
###############################################################################

map_asv <- function(i) {

  genus <- tax$Genus[i]
  species_refined <- tax$species_refined[i]
  species_silva <- tax$Species[i]

  # Gardnerella: armonización exclusiva para VALENCIA.
  if (
    !is.na(genus) &&
      genus == "Gardnerella") {

    return(
      c(
        Taxon = "Gardnerella_vaginalis",
        Level = "VALENCIA_compatibility",
        Rule = "Gardnerella_all_to_vaginalis"))
  }

  # Prioridad 1: species_refined con coincidencia exacta.
  if (!is_ambiguous_species(
    species_refined)) {

    hit <- match_reference(
      gsub(
        " ",
        "_",
        species_refined))

    if (!is.na(hit)) {

      return(
        c(
          Taxon = hit,
          Level = "Species_refined",
          Rule = "Exact_species_refined"))
    }
  }

  # Prioridad 2: SILVA solo cuando species_refined está ausente.
  refined_missing <- is.na(
    species_refined) ||
    species_refined == "" ||
    species_refined == "NA"

  if (
    refined_missing &&
      !is.na(genus) &&
      genus != "" &&
      !is.na(species_silva) &&
      species_silva != "") {

    silva_full <- paste(
      genus,
      species_silva)

    if (!is_ambiguous_species(
      silva_full)) {

      hit <- match_reference(
        gsub(
          " ",
          "_",
          silva_full))

      if (!is.na(hit)) {

        return(
          c(
            Taxon = hit,
            Level = "Species_SILVA",
            Rule = "Exact_species_SILVA"))
      }
    }
  }

  # Géneros vaginales principales sin especie compatible.
  genus_fallback <- c(
    "Lactobacillus",
    "Prevotella",
    "Atopobium",
    "Sneathia")

  if (
    !is.na(genus) &&
      genus %in% genus_fallback) {

    hit <- match_reference(
      paste0(
        "g_",
        genus))

    if (!is.na(hit)) {

      return(
        c(
          Taxon = hit,
          Level = "Genus",
          Rule = paste0(
            "Fallback_g_",
            genus)))
    }
  }

  # Fallback jerárquico general.
  candidates <- list(
    Genus = ifelse(
      is.na(tax$Genus[i]),
      NA_character_,
      paste0(
        "g_",
        tax$Genus[i])),
    Family = ifelse(
      is.na(tax$Family[i]),
      NA_character_,
      paste0(
        "f_",
        tax$Family[i])),
    Order = ifelse(
      is.na(tax$Order[i]),
      NA_character_,
      paste0(
        "o_",
        tax$Order[i])),
    Class = ifelse(
      is.na(tax$Class[i]),
      NA_character_,
      paste0(
        "c_",
        tax$Class[i])),
    Phylum = ifelse(
      is.na(tax$Phylum[i]),
      NA_character_,
      paste0(
        "p_",
        tax$Phylum[i])))

  for (level in names(
    candidates)) {

    hit <- match_reference(
      candidates[[level]])

    if (!is.na(hit)) {

      return(
        c(
          Taxon = hit,
          Level = level,
          Rule = paste0(
            "Fallback_",
            level)))
    }
  }

  c(
    Taxon = "d_NA",
    Level = "Unmapped",
    Rule = "Fallback_d_NA")
}

mapped <- t(
  vapply(
    seq_len(nrow(tax)),
    map_asv,
    FUN.VALUE = c(
      Taxon = "",
      Level = "",
      Rule = "")))

mapped <- as.data.frame(
  mapped,
  stringsAsFactors = FALSE)

asv_mapping <- data.frame(
  ASV_ID = rownames(tax),
  Genus = tax$Genus,
  Species = tax$Species,
  species_refined = tax$species_refined,
  VALENCIA_taxon = mapped$Taxon,
  Mapping_level = mapped$Level,
  Mapping_rule = mapped$Rule,
  stringsAsFactors = FALSE)

###############################################################################
# 7. Auditoría del mapeo
###############################################################################

target_genera <- c(
  "Lactobacillus",
  "Gardnerella",
  "Prevotella",
  "Atopobium",
  "Sneathia")

target_mapping <- asv_mapping %>%
  filter(
    Genus %in% target_genera) %>%
  count(
    Genus,
    species_refined,
    VALENCIA_taxon,
    Mapping_level,
    Mapping_rule,
    name = "N_ASVs") %>%
  arrange(
    Genus,
    VALENCIA_taxon)

mapping_audit <- asv_mapping %>%
  count(
    Mapping_level,
    Mapping_rule,
    name = "N_ASVs",
    sort = TRUE)

cat(
  "\n=== MAPEO DE TAXONES VAGINALES PRINCIPALES ===\n")

print(
  as.data.frame(
    target_mapping),
  row.names = FALSE)

cat("\n=== AUDITORÍA GENERAL DEL MAPEO ===\n")
print(mapping_audit)

cat(
  "\nASVs asignados a d_NA:",
  sum(
    asv_mapping$VALENCIA_taxon ==
      "d_NA"),
  "\n")

write_csv(
  asv_mapping,
  mapping_file)

write_csv(
  target_mapping,
  file.path(
    valencia_dir,
    "VALENCIA_TargetTaxa_FinalMapping.csv"))

write_csv(
  mapping_audit,
  file.path(
    valencia_dir,
    "VALENCIA_mapping_audit.csv"))

###############################################################################
# 8. Agregar conteos por filotipo VALENCIA
###############################################################################

counts_valencia <- rowsum(
  otu,
  group =
    asv_mapping$VALENCIA_taxon,
  reorder = FALSE)

missing_ref_taxa <- setdiff(
  ref_taxa,
  rownames(counts_valencia))

if (length(
  missing_ref_taxa) > 0) {

  zeros <- matrix(
    0,
    nrow =
      length(missing_ref_taxa),
    ncol =
      ncol(counts_valencia),
    dimnames = list(
      missing_ref_taxa,
      colnames(counts_valencia)))

  counts_valencia <- rbind(
    counts_valencia,
    zeros)
}

counts_valencia <- counts_valencia[
  ref_taxa,
  ,
  drop = FALSE]

###############################################################################
# 9. Construir input para VALENCIA
###############################################################################

sample_counts <- t(
  counts_valencia)

original_read_count <- colSums(
  otu)

mapped_read_count <- rowSums(
  sample_counts)

if (!all(
  mapped_read_count ==
    original_read_count[
      rownames(sample_counts)])) {

  stop(
    "Los conteos mapeados no coinciden con los conteos originales.")
}

valencia_input <- data.frame(
  sampleID =
    rownames(sample_counts),
  read_count = as.numeric(
    original_read_count[
      rownames(sample_counts)]),
  sample_counts,
  check.names = FALSE,
  stringsAsFactors = FALSE)

if (
  names(valencia_input)[1] !=
    "sampleID" ||
    names(valencia_input)[2] !=
      "read_count") {

  stop(
    "Las dos primeras columnas no cumplen el formato VALENCIA.")
}

if (anyDuplicated(
  valencia_input$sampleID))
  stop("Existen sampleID duplicados.")

if (any(
  valencia_input$read_count <= 0))
  stop("Existen muestras con read_count <= 0.")

write_csv(
  valencia_input,
  input_file)

###############################################################################
# 10. Auditoría de lecturas no resueltas
###############################################################################

dna_reads <- sample_counts[
  ,
  "d_NA"]

dna_audit <- data.frame(
  sampleID =
    rownames(sample_counts),
  read_count =
    original_read_count[
      rownames(sample_counts)],
  d_NA_reads = dna_reads,
  d_NA_fraction =
    dna_reads /
    original_read_count[
      rownames(sample_counts)],
  stringsAsFactors = FALSE)

cat(
  "\n=== LECTURAS NO RESUELTAS PARA VALENCIA ===\n")

print(
  summary(
    dna_audit$d_NA_fraction))

write_csv(
  dna_audit,
  file.path(
    valencia_dir,
    "VALENCIA_dNA_by_sample.csv"))

###############################################################################
# 11. Guardar metadata
###############################################################################

meta_out <- meta %>%
  as.data.frame() %>%
  rownames_to_column(
    var = "rowname_internal") %>%
  mutate(
    sampleID =
      rowname_internal) %>%
  select(
    -rowname_internal)

write_csv(
  meta_out,
  metadata_file)

###############################################################################
# 12. Resumen del input
###############################################################################

cat("\n========================================\n")
cat("INPUT PARA VALENCIA CREADO\n")
cat("========================================\n")
cat("Archivo:",input_file,"\n")

cat("Muestras:",nrow(valencia_input),"\n") # 2179
cat(
  "Filotipos:",
  ncol(valencia_input) - 2,
  "\n")
cat(
  "Lecturas totales:",
  sum(valencia_input$read_count),
  "\n")

###############################################################################
# 13. Integrar resultados de VALENCIA
###############################################################################

if (file.exists(output_file)) {

  cat(
    "\n=== INTEGRANDO RESULTADOS DE VALENCIA ===\n")

  valencia_out <- read_csv(
    output_file,
    show_col_types = FALSE)

  required_out <- c(
    "sampleID",
    "subCST",
    "score",
    "CST")

  missing_out <- setdiff(
    required_out,
    names(valencia_out))

  if (length(missing_out) > 0)
    stop(
      paste(
        "Faltan columnas en VALENCIA_output.csv:",
        paste(
          missing_out,
          collapse = ", ")))

  if (anyDuplicated(
    valencia_out$sampleID))
    stop(
      "VALENCIA_output.csv contiene sampleID duplicados.")

  if (!setequal(
    sample_names(ps),
    valencia_out$sampleID)) {

    stop(
      "Los sampleID del output no coinciden con phyloseq.")
  }

  if (anyNA(
    suppressWarnings(
      as.numeric(
        valencia_out$score)))) {

    stop(
      "VALENCIA_output.csv contiene score ausente o no numérico.")
  }

  valencia_out <- valencia_out[
    match(
      sample_names(ps),
      valencia_out$sampleID),
    ]

  sd <- data.frame(
    sample_data(ps),
    stringsAsFactors = FALSE)

  sd$VALENCIA_subCST <-
    valencia_out$subCST

  sd$VALENCIA_CST <-
    valencia_out$CST

  sd$VALENCIA_score <-
    as.numeric(
      valencia_out$score)

  rownames(sd) <-
    sample_names(ps)

  sample_data(ps) <-
    sample_data(sd)

  saveRDS(
    ps,
    ps_output_file)

  cst_summary <- sd %>%
    rownames_to_column(
      "sampleID") %>%
    count(
      cohort,
      term_vs_preterm_delivery,
      VALENCIA_CST,
      name = "N_samples") %>%
    group_by(
      cohort,
      term_vs_preterm_delivery) %>%
    mutate(
      Proportion =
        N_samples /
        sum(N_samples)) %>%
    ungroup()

  write_csv(
    cst_summary,
    file.path(
      valencia_dir,
      "VALENCIA_CST_summary.csv"))

  write_csv(
    valencia_out,
    file.path(
      valencia_dir,
      "VALENCIA_assignments_full.csv"))

  cat(
    "Objeto con CST guardado en:",
    ps_output_file,
    "\n")

} else {

  cat(
    "\nVALENCIA_output.csv todavía no existe.\n")

  cat(
    "Ejecute Valencia.py en Ubuntu y vuelva a ejecutar este script.\n")
}

###############################################################################
# 14. Reproducibilidad
###############################################################################

versions <- data.frame(
  Package = c(
    "R",
    "phyloseq",
    "dplyr",
    "readr",
    "stringr",
    "tibble"),
  Version = c(
    paste(
      R.version$major,
      R.version$minor,
      sep = "."),
    as.character(
      packageVersion("phyloseq")),
    as.character(
      packageVersion("dplyr")),
    as.character(
      packageVersion("readr")),
    as.character(
      packageVersion("stringr")),
    as.character(
      packageVersion("tibble"))),
  stringsAsFactors = FALSE)

write_csv(
  versions,
  file.path(
    valencia_dir,
    "VALENCIA_Preparation_SoftwareVersions.csv"))

reproducibility <- data.frame(
  Item = c(
    "Input_phyloseq",
    "Input_phyloseq_MD5",
    "VALENCIA_centroids",
    "VALENCIA_centroids_MD5",
    "VALENCIA_repository",
    "VALENCIA_input_MD5",
    "VALENCIA_output_MD5"),
  Value = c(
    ps_file,
    unname(
      tools::md5sum(
        ps_file)),
    ref_file,
    unname(
      tools::md5sum(
        ref_file)),
    "https://github.com/ravel-lab/VALENCIA",
    unname(
      tools::md5sum(
        input_file)),
    if (
      file.exists(
        output_file)) {
      unname(
        tools::md5sum(
          output_file))
    } else {
      NA_character_
    }),
  stringsAsFactors = FALSE)

write_csv(
  reproducibility,
  file.path(
    valencia_dir,
    "VALENCIA_Preparation_Reproducibility.csv"))

writeLines(
  capture.output(
    sessionInfo()),
  file.path(
    valencia_dir,
    "VALENCIA_Preparation_sessionInfo.txt"))

cat("\n=== VERSIONES ===\n")
print(versions)

cat("\n=== REPRODUCIBILIDAD ===\n")
print(reproducibility)

cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

###############################################################################
# FIN
###############################################################################
