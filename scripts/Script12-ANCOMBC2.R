###############################################################################
#                                    TFM
###############################################################################
#              Script 12: Análisis de abundancia diferencial
###############################################################################

# Objetivo:
# Identificar taxones diferencialmente abundantes entre parto a término y
# pretérmino mediante ANCOM-BC2, a nivel de género y especie.
#
# Diseño principal:
# - ventana gestacional 15-33 semanas;
# - Global:
#   Delivery + cohort + GestWeek + (1 | host_subject_id);
# - Stanford/UAB:
#   Delivery + GestWeek + (1 | host_subject_id);
# - sin interacción Delivery x cohort;
# - sensibilidad adicional por edad materna.
#
# ANCOM-BC2 se utiliza porque permite ajustar covariables y medidas repetidas
# mediante efectos aleatorios y aplica bias correction para datos
# composicionales. No se añade profundidad de secuenciación como covariable.
#
# Prevalencia:
# prv_cut = 0.05 se aplica dentro de cada análisis sobre las muestras 15-33.
#
# Precisión:
# LFC mínimo detectable optimista = 1.96 * SE. Esta cota no incorpora el FDR.
#
# Referencia:
# Lin H, Peddada SD. Nat Methods. 2024;21:83-91.
#
# Ejecutar desde la raíz del repositorio VaginalMicrobiome o desde su carpeta
# inmediatamente superior.

###############################################################################

rm(list = ls())

library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(readr)
library(ggplot2)

###############################################################################
# 1. Directorios y parámetros
###############################################################################

project_dir <- if (
  dir.exists("VaginalMicrobiome")) {
  "VaginalMicrobiome"
} else {
  "."
}

results_dir <- file.path(
  project_dir,
  "results",
  "Differential_abundance")

images_dir <- file.path(
  project_dir,
  "images",
  "Differential_abundance")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  images_dir,
  recursive = TRUE,
  showWarnings = FALSE)

PREVALENCE_CUTOFF <- 0.05
ALPHA <- 0.05
GW_MIN <- 15
GW_MAX <- 33
SEED <- 20260913

analyses <- c(
  "Global",
  "Stanford",
  "UAB")

set.seed(SEED)

###############################################################################
# 2. Cargar objeto y preparar metadata
###############################################################################

ps_file <- file.path(
  project_dir,
  "cleaned",
  "ps_prevotella_refined.rds")

if (!file.exists(ps_file))
  stop(
    paste(
      "No se encontró:",
      ps_file))

ps <- readRDS(ps_file)

cat("\n=== OBJETO PHYLOSEQ ===\n")
cat("Muestras:",nsamples(ps),"\n")
cat("ASVs:",ntaxa(ps),"\n")

meta <- data.frame(
  sample_data(ps),
  stringsAsFactors = FALSE)

meta$SampleID <- rownames(meta)

required_meta <- c(
  "host_subject_id",
  "Groups",
  "cohort",
  "term_vs_preterm_delivery",
  "gest_week_collection",
  "gest_wk_delivery",
  "AGE")

missing_meta <- setdiff(
  required_meta,
  names(meta))

if (length(missing_meta) > 0)
  stop(
    paste(
      "Faltan variables de metadata:",
      paste(
        missing_meta,
        collapse = ", ")))

meta <- meta %>%
  filter(
    Groups %in% c(
      "ST","SP","UT","UP")) %>%
  mutate(
    cohort = case_when(
      cohort %in% c(
        "S","Stanford") ~
        "Stanford",
      cohort %in% c(
        "U","UAB") ~
        "UAB",
      TRUE ~
        as.character(cohort)),
    Delivery = case_when(
      term_vs_preterm_delivery ==
        "Term" ~
        "Term",
      term_vs_preterm_delivery ==
        "Preterm" ~
        "Preterm",
      Groups %in% c(
        "ST","UT") ~
        "Term",
      Groups %in% c(
        "SP","UP") ~
        "Preterm",
      TRUE ~
        NA_character_),
    GestWeek = suppressWarnings(
      parse_number(
        as.character(
          gest_week_collection))),
    GestWeek_Delivery =
      suppressWarnings(
        parse_number(
          as.character(
            gest_wk_delivery))),
    AGE = suppressWarnings(
      as.numeric(
        as.character(AGE))),
    Temporal_inconsistency =
      !is.na(GestWeek) &
      !is.na(GestWeek_Delivery) &
      GestWeek >
        GestWeek_Delivery,
    cohort = factor(
      cohort,
      levels = c(
        "Stanford","UAB")),
    Delivery = factor(
      Delivery,
      levels = c(
        "Term","Preterm")),
    host_subject_id =
      factor(host_subject_id))

if (
  anyNA(meta$cohort) ||
    anyNA(meta$Delivery))
  stop(
    "Existen muestras sin cohorte o desenlace correctamente definido.")

###############################################################################
# 3. Ventana gestacional
###############################################################################

meta <- meta %>%
  filter(
    !Temporal_inconsistency,
    !is.na(GestWeek),
    GestWeek >= GW_MIN,
    GestWeek <= GW_MAX)

cat("\n=== VENTANA GESTACIONAL ===\n")
cat("Semanas:",GW_MIN,"-",GW_MAX,"\n")
cat("Muestras:",nrow(meta),"\n")
cat(
  "Mujeres:",
  n_distinct(
    meta$host_subject_id),
  "\n")

print(
  meta %>%
    distinct(
      host_subject_id,
      cohort,
      Delivery) %>%
    count(
      cohort,
      Delivery))

rownames(meta) <-
  meta$SampleID

ps <- prune_samples(
  meta$SampleID,
  ps)

ps <- prune_taxa(
  taxa_sums(ps) > 0,
  ps)

###############################################################################
# 4. Comprobar estructura longitudinal
###############################################################################

subject_check <- meta %>%
  group_by(
    host_subject_id) %>%
  summarise(
    N_samples = n(),
    N_cohort =
      n_distinct(cohort),
    N_delivery =
      n_distinct(Delivery),
    N_age =
      n_distinct(
        AGE[
          !is.na(AGE)]),
    .groups = "drop")

if (
  any(
    subject_check$N_cohort != 1) ||
    any(
      subject_check$N_delivery != 1) ||
    any(
      subject_check$N_age > 1))
  stop(
    "Alguna mujer presenta metadata clínica inconsistente.")

cat("\nMuestras por mujer:\n")
print(
  summary(
    subject_check$N_samples))

###############################################################################
# 5. Extraer taxonomía y conteos
###############################################################################

tax <- as.data.frame(
  tax_table(ps),
  stringsAsFactors = FALSE)

otu <- as(
  otu_table(ps),
  "matrix")

if (!taxa_are_rows(ps))
  otu <- t(otu)

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
# 6. Crear taxones de género y especie
###############################################################################

tax$GenusTaxon <-
  tax$Genus

tax$GenusTaxon[
  is.na(
    tax$GenusTaxon) |
    tax$GenusTaxon == "" |
    tax$GenusTaxon == "NA"] <-
  NA_character_

tax$SpeciesTaxon <-
  NA_character_

idx_refined <-
  !is.na(tax$Genus) &
  tax$Genus != "" &
  !is.na(
    tax$species_refined) &
  tax$species_refined != "" &
  tax$species_refined != "NA"

idx_full <- idx_refined &
  startsWith(
    tax$species_refined,
    paste0(
      tax$Genus,
      " "))

tax$SpeciesTaxon[
  idx_full] <-
  tax$species_refined[
    idx_full]

idx_short <- idx_refined &
  !idx_full

tax$SpeciesTaxon[
  idx_short] <- paste(
    tax$Genus[
      idx_short],
    tax$species_refined[
      idx_short])

idx_sp <-
  is.na(
    tax$SpeciesTaxon) &
  !is.na(tax$Genus) &
  tax$Genus != ""

tax$SpeciesTaxon[
  idx_sp] <- paste(
    tax$Genus[
      idx_sp],
    "sp.")

tax$SpeciesTaxon[
  is.na(tax$Genus) |
    tax$Genus == ""] <-
  NA_character_

cat("\n=== TAXONOMÍA PARA ANCOM-BC2 ===\n")
cat(
  "ASVs con género:",
  sum(
    !is.na(
      tax$GenusTaxon)),
  "\n")
cat(
  "ASVs sin género excluidos:",
  sum(
    is.na(
      tax$GenusTaxon)),
  "\n")
cat(
  "Géneros únicos:",
  n_distinct(
    tax$GenusTaxon,
    na.rm = TRUE),
  "\n")
cat(
  "Taxones de especie únicos:",
  n_distinct(
    tax$SpeciesTaxon,
    na.rm = TRUE),
  "\n")

###############################################################################
# 7. Agregar ASVs por nivel taxonómico
###############################################################################

aggregate_taxa <- function(
  otu,
  grouping) {

  valid <-
    !is.na(grouping) &
    grouping != ""

  counts <- rowsum(
    otu[
      valid,
      ,
      drop = FALSE],
    group =
      grouping[
        valid],
    reorder = FALSE)

  counts[
    rowSums(counts) > 0,
    ,
    drop = FALSE]
}

counts_list <- list(
  Genus =
    aggregate_taxa(
      otu,
      tax$GenusTaxon),
  Species =
    aggregate_taxa(
      otu,
      tax$SpeciesTaxon))

###############################################################################
# 8. Auditoría de prevalencia
###############################################################################

audit_prevalence <- function(
  counts,
  metadata,
  analysis,
  tax_level,
  adjustment) {

  metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE]

  n_samples <- ncol(counts)

  n_subjects <- n_distinct(
    metadata$host_subject_id)

  sample_presence <-
    counts > 0

  sample_prev <- rowSums(
    sample_presence) /
    n_samples

  subjects <- unique(
    as.character(
      metadata$host_subject_id))

  subject_presence <- sapply(
    subjects,
    function(subject) {

      samples <- rownames(
        metadata)[
          as.character(
            metadata$host_subject_id) ==
            subject]

      rowSums(
        counts[
          ,
          samples,
          drop = FALSE]) >
        0
    })

  if (is.vector(
    subject_presence)) {

    subject_presence <- matrix(
      subject_presence,
      ncol = 1)
  }

  subject_prev <- rowSums(
    subject_presence) /
    n_subjects

  pass <- sample_prev >=
    PREVALENCE_CUTOFF

  detail <- data.frame(
    Adjustment = adjustment,
    Taxonomic_Level =
      tax_level,
    Analysis = analysis,
    Taxon =
      rownames(counts),
    N_samples =
      n_samples,
    Samples_present =
      rowSums(
        sample_presence),
    Sample_Prevalence =
      sample_prev,
    N_subjects =
      n_subjects,
    Subjects_present =
      rowSums(
        subject_presence),
    Subject_Prevalence =
      subject_prev,
    Pass_prv_cut =
      pass,
    stringsAsFactors = FALSE)

  summary <- data.frame(
    Adjustment = adjustment,
    Taxonomic_Level =
      tax_level,
    Analysis = analysis,
    Samples = n_samples,
    Subjects = n_subjects,
    Taxa_before_prv_cut =
      nrow(counts),
    Taxa_after_prv_cut =
      sum(pass),
    Taxa_removed_prv_cut =
      sum(!pass),
    prv_cut =
      PREVALENCE_CUTOFF,
    stringsAsFactors = FALSE)

  list(
    table = detail,
    summary = summary)
}

###############################################################################
# 9. Preparar subconjuntos
###############################################################################

prepare_subset <- function(
  counts,
  metadata,
  analysis,
  age_adjusted = FALSE) {

  d <- metadata

  if (analysis == "Stanford")
    d <- d %>%
      filter(
        cohort == "Stanford")

  if (analysis == "UAB")
    d <- d %>%
      filter(
        cohort == "UAB")

  if (age_adjusted)
    d <- d %>%
      filter(
        !is.na(AGE))

  d <- droplevels(d)

  counts_sub <- counts[
    ,
    rownames(d),
    drop = FALSE]

  list(
    counts = counts_sub,
    metadata = d)
}

###############################################################################
# 10. Ejecutar ANCOM-BC2
###############################################################################

run_ancombc2 <- function(
  counts,
  metadata,
  fixed_formula,
  analysis,
  tax_level,
  adjustment) {

  metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE]

  ps_temp <- phyloseq(
    otu_table(
      counts,
      taxa_are_rows = TRUE),
    sample_data(
      metadata))

  cat("\n========================================\n")
  cat(
    "ANCOM-BC2 |",
    adjustment,"|",
    tax_level,"|",
    analysis,"\n")
  cat("========================================\n")
  cat("Muestras:",ncol(counts),"\n")
  cat(
    "Mujeres:",
    n_distinct(
      metadata$host_subject_id),
    "\n")
  cat(
    "Fórmula fija:",
    fixed_formula,
    "\n")

  ancombc2(
    data = ps_temp,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula =
      fixed_formula,
    rand_formula =
      "(1|host_subject_id)",
    p_adj_method = "BH",
    pseudo_sens = TRUE,
    prv_cut =
      PREVALENCE_CUTOFF,
    lib_cut = 0,
    group = NULL,
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = ALPHA,
    n_cl = 1,
    verbose = TRUE)
}

###############################################################################
# 11. Extraer efecto Delivery
###############################################################################

extract_delivery <- function(
  result,
  analysis,
  tax_level,
  adjustment) {

  res <- result$res

  if (!is.data.frame(res))
    stop(
      "resultado$res no es un data.frame.")

  lfc_col <- grep(
    "^lfc_.*DeliveryPreterm",
    colnames(res),
    value = TRUE)

  if (length(lfc_col) != 1) {

    cat("\nColumnas disponibles:\n")
    print(
      colnames(res))

    stop(
      paste(
        "No se pudo identificar DeliveryPreterm en",
        tax_level,
        analysis,
        adjustment))
  }

  coefficient <- sub(
    "^lfc_",
    "",
    lfc_col)

  se_col <- paste0(
    "se_",
    coefficient)

  w_col <- paste0(
    "W_",
    coefficient)

  p_col <- paste0(
    "p_",
    coefficient)

  q_col <- paste0(
    "q_",
    coefficient)

  diff_col <- paste0(
    "diff_",
    coefficient)

  passed_col <- paste0(
    "passed_ss_",
    coefficient)

  robust_col <- paste0(
    "diff_robust_",
    coefficient)

  required <- c(
    lfc_col,
    se_col,
    w_col,
    p_col,
    q_col,
    diff_col)

  missing <- setdiff(
    required,
    colnames(res))

  if (length(missing) > 0)
    stop(
      paste(
        "Faltan columnas ANCOM-BC2:",
        paste(
          missing,
          collapse = ", ")))

  taxon <- if (
    "taxon" %in%
      colnames(res)) {
    res$taxon
  } else {
    rownames(res)
  }

  passed_sensitivity <- if (
    passed_col %in%
      colnames(res)) {
    res[[passed_col]]
  } else {
    rep(
      NA,
      nrow(res))
  }

  diff_robust <- if (
    robust_col %in%
      colnames(res)) {
    res[[robust_col]]
  } else {
    res[[diff_col]] &
      !is.na(
        passed_sensitivity) &
      passed_sensitivity
  }

  data.frame(
    Adjustment = adjustment,
    Taxonomic_Level =
      tax_level,
    Analysis = analysis,
    Taxon = taxon,
    LFC =
      res[[lfc_col]],
    SE =
      res[[se_col]],
    W =
      res[[w_col]],
    P_value =
      res[[p_col]],
    FDR =
      res[[q_col]],
    Differential_abundance =
      res[[diff_col]],
    Passed_sensitivity =
      passed_sensitivity,
    Robust_significant =
      diff_robust,
    stringsAsFactors = FALSE) %>%
    mutate(
      IC95_low =
        LFC -
        1.96 * SE,
      IC95_high =
        LFC +
        1.96 * SE,
      LFC_min_detectable =
        1.96 * SE,
      Below_optimistic_MDL =
        abs(LFC) <
        LFC_min_detectable,
      Significant =
        !is.na(FDR) &
        FDR < ALPHA,
      Direction = case_when(
        is.na(LFC) ~
          NA_character_,
        LFC > 0 ~
          "Higher_in_Preterm",
        LFC < 0 ~
          "Higher_in_Term",
        TRUE ~
          "No_direction"))
}

###############################################################################
# 12. Ejecutar bloque principal o ajustado por edad
###############################################################################

run_block <- function(
  counts_list,
  metadata,
  age_adjusted = FALSE) {

  adjustment <- if (
    age_adjusted) {
    "Age_adjusted"
  } else {
    "Primary"
  }

  results <- list()
  prev_detail <- list()
  prev_summary <- list()
  models <- list()

  k <- 1

  for (tax_level in
    names(counts_list)) {

    for (analysis in
      analyses) {

      sub <- prepare_subset(
        counts_list[[
          tax_level]],
        metadata,
        analysis,
        age_adjusted)

      fixed_formula <- if (
        analysis == "Global") {

        if (age_adjusted) {
          paste(
            "Delivery + cohort +",
            "GestWeek + AGE")
        } else {
          paste(
            "Delivery + cohort +",
            "GestWeek")
        }

      } else {

        if (age_adjusted) {
          paste(
            "Delivery + GestWeek +",
            "AGE")
        } else {
          "Delivery + GestWeek"
        }
      }

      prev <- audit_prevalence(
        sub$counts,
        sub$metadata,
        analysis,
        tax_level,
        adjustment)

      set.seed(SEED)

      model <- run_ancombc2(
        sub$counts,
        sub$metadata,
        fixed_formula,
        analysis,
        tax_level,
        adjustment)

      results[[k]] <-
        extract_delivery(
          model,
          analysis,
          tax_level,
          adjustment)

      prev_detail[[k]] <-
        prev$table

      prev_summary[[k]] <-
        prev$summary

      models[[
        paste(
          adjustment,
          tax_level,
          analysis,
          sep = "_")
      ]] <- model

      k <- k + 1
    }
  }

  list(
    results =
      bind_rows(results),
    prevalence =
      bind_rows(
        prev_detail),
    prevalence_summary =
      bind_rows(
        prev_summary),
    models = models)
}

###############################################################################
# 13. Análisis principal
###############################################################################

cat("\n\n========================================\n")
cat("ANÁLISIS PRINCIPAL\n")
cat("========================================\n")

primary <- run_block(
  counts_list,
  meta,
  age_adjusted = FALSE)

Primary_Results <-
  primary$results

Primary_Prevalence <-
  primary$prevalence

Primary_Prevalence_Summary <-
  primary$prevalence_summary

cat("\n=== PREVALENCIA: PRINCIPAL ===\n")
print(
  Primary_Prevalence_Summary)

###############################################################################
# 14. Sensibilidad ajustada por edad materna
###############################################################################

cat("\n\n========================================\n")
cat("SENSIBILIDAD AJUSTADA POR EDAD\n")
cat("========================================\n")

age_adjusted <- run_block(
  counts_list,
  meta,
  age_adjusted = TRUE)

AgeAdjusted_Results <-
  age_adjusted$results

AgeAdjusted_Prevalence <-
  age_adjusted$prevalence

AgeAdjusted_Prevalence_Summary <-
  age_adjusted$prevalence_summary

cat(
  "\n=== PREVALENCIA: AJUSTADO POR EDAD ===\n")

print(
  AgeAdjusted_Prevalence_Summary)

###############################################################################
# 15. Resumen de resultados
###############################################################################

create_summary <- function(
  results,
  prevalence_summary) {

  results %>%
    group_by(
      Adjustment,
      Taxonomic_Level,
      Analysis) %>%
    summarise(
      Taxa_tested = n(),
      Significant_taxa =
        sum(
          Significant,
          na.rm = TRUE),
      Robust_significant_taxa =
        sum(
          Robust_significant,
          na.rm = TRUE),
      Higher_in_Preterm =
        sum(
          Significant &
          Direction ==
            "Higher_in_Preterm",
          na.rm = TRUE),
      Higher_in_Term =
        sum(
          Significant &
          Direction ==
            "Higher_in_Term",
          na.rm = TRUE),
      .groups = "drop") %>%
    left_join(
      prevalence_summary,
      by = c(
        "Adjustment",
        "Taxonomic_Level",
        "Analysis"))
}

Primary_Summary <- create_summary(
  Primary_Results,
  Primary_Prevalence_Summary)

AgeAdjusted_Summary <- create_summary(
  AgeAdjusted_Results,
  AgeAdjusted_Prevalence_Summary)

DA_Summary <- bind_rows(
  Primary_Summary,
  AgeAdjusted_Summary)

cat("\n=== RESUMEN ANCOM-BC2 ===\n")
print(
  DA_Summary)

###############################################################################
# 16. Taxones significativos y dirigidos
###############################################################################

All_Results <- bind_rows(
  Primary_Results,
  AgeAdjusted_Results)

Significant_Results <- All_Results %>%
  filter(
    Significant) %>%
  arrange(
    Adjustment,
    Taxonomic_Level,
    Analysis,
    FDR)

Robust_Results <- All_Results %>%
  filter(
    Robust_significant) %>%
  arrange(
    Adjustment,
    Taxonomic_Level,
    Analysis,
    FDR)

genera_interest <- c(
  "Lactobacillus",
  "Gardnerella",
  "Prevotella")

species_interest <- c(
  "Lactobacillus crispatus",
  "Lactobacillus iners",
  "Lactobacillus jensenii",
  "Lactobacillus mulieris",
  "Lactobacillus gasseri",
  "Gardnerella G1",
  "Gardnerella G2",
  "Gardnerella G3",
  "Gardnerella G_other")

Taxa_Interest <- All_Results %>%
  filter(
    Taxon %in%
      genera_interest |
    Taxon %in%
      species_interest |
    (
      grepl(
        "^Prevotella ",
        Taxon) &
      Taxon !=
        "Prevotella sp."
    )) %>%
  arrange(
    Adjustment,
    Taxonomic_Level,
    Analysis,
    Taxon)

Target_Precision <- Taxa_Interest %>%
  select(
    Adjustment,
    Taxonomic_Level,
    Analysis,
    Taxon,
    LFC,
    SE,
    IC95_low,
    IC95_high,
    LFC_min_detectable,
    Below_optimistic_MDL,
    P_value,
    FDR,
    Significant,
    Passed_sensitivity,
    Robust_significant,
    Direction)

cat("\n=== PRECISIÓN DE TAXONES DIRIGIDOS ===\n")
print(
  Target_Precision)

###############################################################################
# 17. Comparar principal vs ajustado por edad
###############################################################################

primary_cmp <- Primary_Results %>%
  select(
    Taxonomic_Level,
    Analysis,
    Taxon,
    LFC_primary = LFC,
    SE_primary = SE,
    P_primary =
      P_value,
    FDR_primary = FDR,
    Significant_primary =
      Significant,
    Robust_primary =
      Robust_significant,
    Direction_primary =
      Direction)

age_cmp <- AgeAdjusted_Results %>%
  select(
    Taxonomic_Level,
    Analysis,
    Taxon,
    LFC_age = LFC,
    SE_age = SE,
    P_age =
      P_value,
    FDR_age = FDR,
    Significant_age =
      Significant,
    Robust_age =
      Robust_significant,
    Direction_age =
      Direction)

Age_Comparison <- full_join(
  primary_cmp,
  age_cmp,
  by = c(
    "Taxonomic_Level",
    "Analysis",
    "Taxon")) %>%
  mutate(
    Comparable =
      !is.na(
        LFC_primary) &
      !is.na(
        LFC_age),
    Delta_LFC =
      LFC_age -
      LFC_primary,
    Direction_changed =
      Comparable &
      !is.na(
        Direction_primary) &
      !is.na(
        Direction_age) &
      Direction_primary !=
        Direction_age,
    Significance_changed =
      Comparable &
      !is.na(
        Significant_primary) &
      !is.na(
        Significant_age) &
      Significant_primary !=
        Significant_age)

Age_Comparison_Summary <-
  Age_Comparison %>%
  group_by(
    Taxonomic_Level,
    Analysis) %>%
  summarise(
    Taxa_compared =
      sum(
        Comparable),
    Spearman_LFC = ifelse(
      sum(
        Comparable) >= 3,
      cor(
        LFC_primary[
          Comparable],
        LFC_age[
          Comparable],
        method = "spearman",
        use = "complete.obs"),
      NA_real_),
    Median_abs_Delta_LFC =
      median(
        abs(
          Delta_LFC[
            Comparable]),
        na.rm = TRUE),
    N_direction_changes =
      sum(
        Direction_changed,
        na.rm = TRUE),
    N_significance_changes =
      sum(
        Significance_changed,
        na.rm = TRUE),
    .groups = "drop")

cat(
  "\n=== PRINCIPAL VS AJUSTADO POR EDAD ===\n")

print(
  Age_Comparison_Summary)

###############################################################################
# 18. Gráficos principal vs ajustado por edad
###############################################################################

plot_data <- Age_Comparison %>%
  filter(
    Comparable)

for (tax_level in c(
  "Genus","Species")) {

  for (analysis in analyses) {

    dat <- plot_data %>%
      filter(
        Taxonomic_Level ==
          tax_level,
        Analysis ==
          analysis)

    if (nrow(dat) < 3)
      next

    rho <- cor(
      dat$LFC_primary,
      dat$LFC_age,
      method = "spearman",
      use = "complete.obs")

    p <- ggplot(
      dat,
      aes(
        x =
          LFC_primary,
        y =
          LFC_age)) +
      geom_point(
        alpha = .65) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = 2) +
      theme_bw(
        base_size = 12) +
      labs(
        title = paste(
          tax_level,
          "-",
          analysis),
        subtitle = paste0(
          "Spearman rho = ",
          round(
            rho,
            3)),
        x =
          "LFC principal",
        y =
          "LFC ajustado por edad")

    ggsave(
      file.path(
        images_dir,
        paste0(
          "ANCOMBC2_LFC_",
          tax_level,
          "_",
          analysis,
          ".png")),
      p,
      width = 6,
      height = 6,
      dpi = 300)
  }
}

###############################################################################
# 19. Guardar resultados
###############################################################################

write_csv(
  Primary_Results,
  file.path(
    results_dir,
    "ANCOMBC2_Primary_All_Results.csv"))

write_csv(
  AgeAdjusted_Results,
  file.path(
    results_dir,
    "ANCOMBC2_AgeAdjusted_All_Results.csv"))

write_csv(
  All_Results,
  file.path(
    results_dir,
    "ANCOMBC2_All_Results.csv"))

write_csv(
  DA_Summary,
  file.path(
    results_dir,
    "ANCOMBC2_Summary.csv"))

write_csv(
  Significant_Results,
  file.path(
    results_dir,
    "ANCOMBC2_Significant.csv"))

write_csv(
  Robust_Results,
  file.path(
    results_dir,
    "ANCOMBC2_Robust_Significant.csv"))

write_csv(
  Taxa_Interest,
  file.path(
    results_dir,
    "ANCOMBC2_Taxa_Interest.csv"))

write_csv(
  Target_Precision,
  file.path(
    results_dir,
    "ANCOMBC2_TargetTaxa_Precision.csv"))

write_csv(
  Primary_Prevalence,
  file.path(
    results_dir,
    "Prevalence_Primary_All_Taxa.csv"))

write_csv(
  AgeAdjusted_Prevalence,
  file.path(
    results_dir,
    "Prevalence_AgeAdjusted_All_Taxa.csv"))

write_csv(
  bind_rows(
    Primary_Prevalence_Summary,
    AgeAdjusted_Prevalence_Summary),
  file.path(
    results_dir,
    "Prevalence_Summary.csv"))

write_csv(
  Age_Comparison,
  file.path(
    results_dir,
    "ANCOMBC2_Primary_vs_AgeAdjusted.csv"))

write_csv(
  Age_Comparison_Summary,
  file.path(
    results_dir,
    "ANCOMBC2_Primary_vs_AgeAdjusted_Summary.csv"))

if (
  requireNamespace(
    "openxlsx",
    quietly = TRUE)) {

  openxlsx::write.xlsx(
    list(
      Summary =
        DA_Summary,
      Primary =
        Primary_Results,
      Age_adjusted =
        AgeAdjusted_Results,
      Significant =
        Significant_Results,
      Robust =
        Robust_Results,
      Taxa_interest =
        Taxa_Interest,
      Target_precision =
        Target_Precision,
      Prevalence =
        bind_rows(
          Primary_Prevalence_Summary,
          AgeAdjusted_Prevalence_Summary),
      Age_comparison =
        Age_Comparison,
      Age_comparison_summary =
        Age_Comparison_Summary),
    file = file.path(
      results_dir,
      "ANCOMBC2_Final_Results.xlsx"),
    rowNames = FALSE,
    overwrite = TRUE)

} else {

  warning(
    paste(
      "openxlsx no está instalado;",
      "se omite la exportación XLSX."))
}

saveRDS(
  list(
    Primary_models =
      primary$models,
    AgeAdjusted_models =
      age_adjusted$models,
    Primary_results =
      Primary_Results,
    AgeAdjusted_results =
      AgeAdjusted_Results,
    Target_precision =
      Target_Precision,
    Prevalence_primary =
      Primary_Prevalence_Summary,
    Prevalence_age =
      AgeAdjusted_Prevalence_Summary),
  file.path(
    results_dir,
    "ANCOMBC2_Final_Objects.rds"))

###############################################################################
# 20. Reproducibilidad
###############################################################################

versions <- data.frame(
  Package = c(
    "R",
    "phyloseq",
    "ANCOMBC",
    "dplyr",
    "readr",
    "ggplot2"),
  Version = c(
    paste(
      R.version$major,
      R.version$minor,
      sep = "."),
    as.character(
      packageVersion(
        "phyloseq")),
    as.character(
      packageVersion(
        "ANCOMBC")),
    as.character(
      packageVersion(
        "dplyr")),
    as.character(
      packageVersion(
        "readr")),
    as.character(
      packageVersion(
        "ggplot2"))),
  stringsAsFactors = FALSE)

write_csv(
  versions,
  file.path(
    results_dir,
    "ANCOMBC2_SoftwareVersions.csv"))

cat("\n=== VERSIONES ===\n")
print(versions)

cat("\n========================================\n")
cat("ANÁLISIS DE ABUNDANCIA DIFERENCIAL FINALIZADO\n")
cat("========================================\n")
cat("Ventana:",GW_MIN,"-",GW_MAX,"semanas\n")
cat("Muestras:",nrow(meta),"\n")
cat(
  "Mujeres:",
  n_distinct(
    meta$host_subject_id),
  "\n")
cat(
  "No se ajustó por profundidad de secuenciación.\n")

cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

###############################################################################
# FIN
###############################################################################

