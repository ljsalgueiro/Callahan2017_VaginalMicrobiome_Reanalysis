###############################################################################
#                                    TFM
###############################################################################
#        Script 11: Estabilidad temporal intraindividual de la microbiota
###############################################################################

# Objetivo:
# Evaluar si la microbiota vaginal cambia más entre muestras consecutivas en
# mujeres con parto pretérmino que en mujeres con parto a término.
#
# Principal:
# - ventana gestacional 15-33 semanas;
# - modelo mixto a nivel de par;
# - logit(Bray-Curtis) como respuesta;
# - estructura varExp para la heterocedasticidad residual;
# - intercepto aleatorio por participante.
#
# Sensibilidades:
# - todo el embarazo con el mismo modelo;
# - número mínimo de pares por participante.
#
# Global:
# logit(BC) ~ Delivery + cohort + Interval + GestWeek_mid +
#   (1 | host_subject_id)
#
# Stanford/UAB:
# logit(BC) ~ Delivery + Interval + GestWeek_mid +
#   (1 | host_subject_id)
#
# Los pares con Interval = 0 se conservan. Los pares consecutivos pueden
# compartir una muestra; el intercepto aleatorio por participante no captura
# por completo esa dependencia adicional.
#
# Ejecutar desde la raíz del repositorio VaginalMicrobiome o desde su carpeta
# inmediatamente superior.

###############################################################################

rm(list = ls())

library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(nlme)

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
  "Temporal_stability",
  "FINAL")

images_dir <- file.path(
  project_dir,
  "images",
  "Temporal_stability",
  "FINAL")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  images_dir,
  recursive = TRUE,
  showWarnings = FALSE)

GW_MIN <- 15
GW_MAX <- 33
EPS <- 1e-3
MIN_PAIRS_VALUES <- c(
  1,2,3,5,10)

analyses <- c(
  "Global",
  "Stanford",
  "UAB")

###############################################################################
# 2. Cargar objeto phyloseq definitivo
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

###############################################################################
# 3. Preparar metadata
###############################################################################

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
  "gest_wk_delivery")

missing_meta <- setdiff(
  required_meta,
  names(meta))

if (length(missing_meta) > 0)
  stop(
    paste(
      "Faltan variables:",
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
# 4. Excluir inconsistencias y preparar abundancias relativas
###############################################################################

meta_valid <- meta %>%
  filter(
    !Temporal_inconsistency,
    !is.na(GestWeek))

cat("\n=== DISPONIBILIDAD TEMPORAL ===\n")
cat("Muestras clínicas:",nrow(meta),"\n")
cat(
  "Muestras con semana válida:",
  nrow(meta_valid),
  "\n")
cat(
  "Inconsistencias temporales excluidas:",
  sum(
    meta$Temporal_inconsistency,
    na.rm = TRUE),
  "\n")

ps_valid <- prune_samples(
  meta_valid$SampleID,
  ps)

ps_valid <- prune_taxa(
  taxa_sums(ps_valid) > 0,
  ps_valid)

ps_rel <- transform_sample_counts(
  ps_valid,
  function(x) {
    if (sum(x) == 0) {
      x
    } else {
      x / sum(x)
    }
  })

otu_rel <- as(
  otu_table(ps_rel),
  "matrix")

if (taxa_are_rows(ps_rel))
  otu_rel <- t(otu_rel)

if (!all(
  meta_valid$SampleID %in%
    rownames(otu_rel)))
  stop(
    "No todas las muestras de metadata están en la matriz OTU.")

otu_rel <- otu_rel[
  meta_valid$SampleID,
  ,
  drop = FALSE]

###############################################################################
# 5. Funciones auxiliares
###############################################################################

logit_bc <- function(
  p,
  eps = EPS) {

  p <- p *
    (1 - 2 * eps) +
    eps

  log(
    p /
      (1 - p))
}

skewness_simple <- function(x) {

  x <- x[
    is.finite(x)]

  if (length(x) < 3)
    return(NA_real_)

  s <- sd(x)

  if (
    !is.finite(s) ||
      s == 0)
    return(NA_real_)

  mean(
    (x - mean(x))^3) /
    s^3
}

make_pairs <- function(
  metadata,
  otu_matrix,
  window_only = TRUE) {

  d <- metadata

  if (window_only) {

    d <- d %>%
      filter(
        GestWeek >= GW_MIN,
        GestWeek <= GW_MAX)
  }

  # SampleID rompe empates de forma determinista cuando hay dos muestras
  # recogidas en la misma semana gestacional.
  d <- d %>%
    arrange(
      host_subject_id,
      GestWeek,
      SampleID)

  pair_index <- d %>%
    group_by(
      host_subject_id) %>%
    mutate(
      SampleID_next =
        lead(SampleID),
      GestWeek_next =
        lead(GestWeek)) %>%
    filter(
      !is.na(
        SampleID_next)) %>%
    ungroup() %>%
    mutate(
      Interval =
        GestWeek_next -
        GestWeek,
      GestWeek_mid =
        (
          GestWeek +
            GestWeek_next
        ) / 2)

  if (nrow(pair_index) == 0)
    return(pair_index)

  pair_index$BC <- vapply(
    seq_len(
      nrow(pair_index)),
    function(i) {

      a <- pair_index$SampleID[i]
      b <- pair_index$SampleID_next[i]

      as.numeric(
        vegdist(
          otu_matrix[
            c(a,b),
            ,
            drop = FALSE],
          method = "bray"))
    },
    numeric(1))

  pair_index %>%
    mutate(
      BC_logit =
        logit_bc(BC))
}

fit_pair_model <- function(
  pairs,
  analysis = "Global",
  transformed = TRUE) {

  d <- pairs

  if (analysis == "Stanford")
    d <- d %>%
      filter(
        cohort == "Stanford")

  if (analysis == "UAB")
    d <- d %>%
      filter(
        cohort == "UAB")

  d <- droplevels(d)

  response <- if (
    transformed) {
    "BC_logit"
  } else {
    "BC"
  }

  rhs <- if (
    analysis == "Global") {
    paste(
      "Delivery + cohort +",
      "Interval + GestWeek_mid")
  } else {
    paste(
      "Delivery + Interval +",
      "GestWeek_mid")
  }

  formula_model <- as.formula(
    paste(
      response,
      "~",
      rhs))

  ctrl <- lmeControl(
    opt = "optim",
    maxIter = 200,
    msMaxIter = 500,
    niterEM = 100,
    tolerance = 1e-6)

  if (transformed) {

    lme(
      formula_model,
      random =
        ~1 | host_subject_id,
      weights =
        varExp(
          form =
            ~fitted(.)),
      data = d,
      method = "REML",
      na.action = na.omit,
      control = ctrl)

  } else {

    lme(
      formula_model,
      random =
        ~1 | host_subject_id,
      data = d,
      method = "REML",
      na.action = na.omit,
      control = ctrl)
  }
}

model_diagnostics <- function(
  model,
  dataset,
  analysis,
  scale_name) {

  r <- residuals(
    model,
    type = "normalized")

  f <- fitted(model)
  d <- getData(model)

  data.frame(
    Dataset = dataset,
    Analysis = analysis,
    Response = scale_name,
    N_pairs = nobs(model),
    N_subjects =
      n_distinct(
        d$host_subject_id),
    Residual_skew =
      skewness_simple(r),
    Cor_absResidual_Fitted =
      suppressWarnings(
        cor(
          abs(r),
          f,
          method = "spearman",
          use = "complete.obs")),
    Cor_Residual_Interval =
      suppressWarnings(
        cor(
          r,
          d$Interval,
          method = "spearman",
          use = "complete.obs")),
    Cor_Residual_GestWeek =
      suppressWarnings(
        cor(
          r,
          d$GestWeek_mid,
          method = "spearman",
          use = "complete.obs")),
    stringsAsFactors = FALSE)
}

save_diagnostic_plots <- function(
  model,
  dataset,
  analysis,
  scale_name) {

  df <- data.frame(
    Fitted =
      fitted(model),
    Residual =
      residuals(
        model,
        type = "normalized"))

  p1 <- ggplot(
    df,
    aes(
      Fitted,
      Residual)) +
    geom_point(
      alpha = .45) +
    geom_smooth(
      method = "loess",
      se = FALSE) +
    geom_hline(
      yintercept = 0,
      linetype = 2) +
    theme_bw() +
    labs(
      title =
        paste(
          dataset,
          analysis,
          scale_name),
      x =
        "Valores ajustados",
      y =
        "Residuos normalizados")

  p2 <- ggplot(
    df,
    aes(
      sample =
        Residual)) +
    stat_qq() +
    stat_qq_line() +
    theme_bw() +
    labs(
      title =
        paste(
          "QQ",
          dataset,
          analysis,
          scale_name))

  prefix <- paste(
    dataset,
    analysis,
    scale_name,
    sep = "_")

  ggsave(
    file.path(
      images_dir,
      paste0(
        "Residuals_",
        prefix,
        ".png")),
    p1,
    width = 6,
    height = 5,
    dpi = 300)

  ggsave(
    file.path(
      images_dir,
      paste0(
        "QQ_",
        prefix,
        ".png")),
    p2,
    width = 6,
    height = 5,
    dpi = 300)
}

extract_delivery <- function(
  model,
  dataset,
  analysis) {

  tt <- summary(
    model)$tTable

  term <- "DeliveryPreterm"

  if (!term %in%
    rownames(tt))
    stop(
      paste(
        "No se encontró",
        term,
        "en",
        dataset,
        analysis))

  ci <- intervals(
    model,
    which = "fixed")$fixed

  data.frame(
    Dataset = dataset,
    Analysis = analysis,
    Contrast =
      "Preterm_vs_Term",
    Beta =
      tt[
        term,
        "Value"],
    SE =
      tt[
        term,
        "Std.Error"],
    DF =
      tt[
        term,
        "DF"],
    T_value =
      tt[
        term,
        "t-value"],
    Pvalue =
      tt[
        term,
        "p-value"],
    CI_low =
      ci[
        term,
        "lower"],
    CI_high =
      ci[
        term,
        "upper"],
    N_pairs =
      nobs(model),
    N_subjects =
      n_distinct(
        getData(model)$host_subject_id),
    Direction = ifelse(
      tt[
        term,
        "Value"] > 0,
      "Higher_instability_in_Preterm",
      "Higher_instability_in_Term"),
    stringsAsFactors = FALSE)
}

###############################################################################
# 6. Construir pares
###############################################################################

pairs_window <- make_pairs(
  meta_valid,
  otu_rel,
  TRUE)

pairs_full <- make_pairs(
  meta_valid,
  otu_rel,
  FALSE)

cat("\n=== PARES: VENTANA 15-33 ===\n")
cat("Pares:",nrow(pairs_window),"\n")
cat(
  "Mujeres con pares:",
  n_distinct(
    pairs_window$host_subject_id),
  "\n")

cat("\n=== PARES: TODO EL EMBARAZO ===\n")
cat("Pares:",nrow(pairs_full),"\n")
cat(
  "Mujeres con pares:",
  n_distinct(
    pairs_full$host_subject_id),
  "\n")

###############################################################################
# 7. Disponibilidad longitudinal
###############################################################################

base_subjects <- meta %>%
  distinct(
    host_subject_id,
    cohort,
    Delivery)

sample_counts_window <- meta_valid %>%
  filter(
    GestWeek >= GW_MIN,
    GestWeek <= GW_MAX) %>%
  count(
    host_subject_id,
    name =
      "N_samples_window")

sample_counts_full <- meta_valid %>%
  count(
    host_subject_id,
    name =
      "N_samples_full")

pair_counts_window <- pairs_window %>%
  count(
    host_subject_id,
    name =
      "N_pairs_window")

pair_counts_full <- pairs_full %>%
  count(
    host_subject_id,
    name =
      "N_pairs_full")

Pair_Availability <- base_subjects %>%
  left_join(
    sample_counts_window,
    by =
      "host_subject_id") %>%
  left_join(
    sample_counts_full,
    by =
      "host_subject_id") %>%
  left_join(
    pair_counts_window,
    by =
      "host_subject_id") %>%
  left_join(
    pair_counts_full,
    by =
      "host_subject_id") %>%
  mutate(
    across(
      c(
        N_samples_window,
        N_samples_full,
        N_pairs_window,
        N_pairs_full),
      ~replace_na(
        .x,
        0L)),
    Eligible_window =
      N_samples_window >= 2,
    Eligible_full =
      N_samples_full >= 2)

Pair_Availability_Summary <- Pair_Availability %>%
  group_by(
    cohort,
    Delivery) %>%
  summarise(
    N_subjects = n(),
    N_window_lt2 =
      sum(
        N_samples_window < 2),
    N_window_with_pairs =
      sum(
        Eligible_window),
    N_full_lt2 =
      sum(
        N_samples_full < 2),
    N_full_with_pairs =
      sum(
        Eligible_full),
    Median_pairs_window =
      median(
        N_pairs_window),
    Median_pairs_full =
      median(
        N_pairs_full),
    .groups = "drop")

Zero_Interval_Summary <- bind_rows(
  data.frame(
    Dataset =
      "Window_15_33",
    N_pairs =
      nrow(pairs_window),
    N_interval_zero =
      sum(
        pairs_window$Interval == 0),
    Percent_interval_zero =
      100 *
      mean(
        pairs_window$Interval == 0)),
  data.frame(
    Dataset =
      "Full_pregnancy",
    N_pairs =
      nrow(pairs_full),
    N_interval_zero =
      sum(
        pairs_full$Interval == 0),
    Percent_interval_zero =
      100 *
      mean(
        pairs_full$Interval == 0)))

cat("\n=== DISPONIBILIDAD DE PARES ===\n")
print(
  Pair_Availability_Summary)

cat("\n=== INTERVALOS DE 0 SEMANAS ===\n")
print(
  Zero_Interval_Summary)

###############################################################################
# 8. Distribución de Bray-Curtis
###############################################################################

diagnose_response <- function(
  pairs,
  label) {

  data.frame(
    Dataset = label,
    N_pairs =
      nrow(pairs),
    N_zero =
      sum(
        pairs$BC == 0),
    N_one =
      sum(
        pairs$BC == 1),
    Min =
      min(
        pairs$BC),
    Q1 =
      quantile(
        pairs$BC,
        .25),
    Median =
      median(
        pairs$BC),
    Mean =
      mean(
        pairs$BC),
    Q3 =
      quantile(
        pairs$BC,
        .75),
    Max =
      max(
        pairs$BC),
    Skew_BC =
      skewness_simple(
        pairs$BC),
    Skew_logit =
      skewness_simple(
        pairs$BC_logit),
    stringsAsFactors = FALSE)
}

BC_Distribution <- bind_rows(
  diagnose_response(
    pairs_window,
    "Window_15_33"),
  diagnose_response(
    pairs_full,
    "Full_pregnancy"))

cat("\n=== DISTRIBUCIÓN BRAY-CURTIS ===\n")
print(
  BC_Distribution)

###############################################################################
# 9. Diagnósticos: BC crudo vs logit(BC) + varExp
###############################################################################

diagnostic_models <- list()
diagnostic_results <- list()

k <- 1

for (dataset in c(
  "Window_15_33",
  "Full_pregnancy")) {

  pairs <- if (
    dataset ==
      "Window_15_33") {
    pairs_window
  } else {
    pairs_full
  }

  for (analysis in analyses) {

    for (transformed in c(
      FALSE,TRUE)) {

      scale_name <- if (
        transformed) {
        "Logit_varExp"
      } else {
        "Raw"
      }

      model <- fit_pair_model(
        pairs,
        analysis,
        transformed)

      key <- paste(
        dataset,
        analysis,
        scale_name,
        sep = "_")

      diagnostic_models[[key]] <-
        model

      diagnostic_results[[k]] <-
        model_diagnostics(
          model,
          dataset,
          analysis,
          scale_name)

      save_diagnostic_plots(
        model,
        dataset,
        analysis,
        scale_name)

      k <- k + 1
    }
  }
}

Diagnostic_Results <- bind_rows(
  diagnostic_results)

cat("\n=== DIAGNÓSTICOS DE MODELOS ===\n")
print(
  Diagnostic_Results)

###############################################################################
# 10. Análisis principal
###############################################################################

Primary_Models <- list(
  Global =
    diagnostic_models[[
      "Window_15_33_Global_Logit_varExp"]],
  Stanford =
    diagnostic_models[[
      "Window_15_33_Stanford_Logit_varExp"]],
  UAB =
    diagnostic_models[[
      "Window_15_33_UAB_Logit_varExp"]])

Primary_Results <- bind_rows(
  lapply(
    names(Primary_Models),
    function(a)
      extract_delivery(
        Primary_Models[[a]],
        "Window_15_33",
        a))) %>%
  mutate(
    Significant =
      Pvalue < 0.05)

cat(
  "\n=== MODELOS PRINCIPALES: VENTANA 15-33 ===\n")

print(
  Primary_Results)

###############################################################################
# 11. Sensibilidad: todo el embarazo
###############################################################################

Full_Models <- list(
  Global =
    diagnostic_models[[
      "Full_pregnancy_Global_Logit_varExp"]],
  Stanford =
    diagnostic_models[[
      "Full_pregnancy_Stanford_Logit_varExp"]],
  UAB =
    diagnostic_models[[
      "Full_pregnancy_UAB_Logit_varExp"]])

FullPregnancy_Results <- bind_rows(
  lapply(
    names(Full_Models),
    function(a)
      extract_delivery(
        Full_Models[[a]],
        "Full_pregnancy",
        a))) %>%
  mutate(
    Significant =
      Pvalue < 0.05)

cat(
  "\n=== SENSIBILIDAD: TODO EL EMBARAZO ===\n")

print(
  FullPregnancy_Results)

###############################################################################
# 12. Descriptivos de pares
###############################################################################

summarise_pairs <- function(
  pairs,
  dataset) {

  pairs %>%
    group_by(
      cohort,
      Delivery) %>%
    summarise(
      Dataset = dataset,
      N_pairs = n(),
      N_subjects =
        n_distinct(
          host_subject_id),
      Median_BC =
        median(BC),
      IQR_BC =
        IQR(BC),
      Mean_BC =
        mean(BC),
      Median_interval =
        median(Interval),
      Mean_interval =
        mean(Interval),
      Median_GestWeek_mid =
        median(
          GestWeek_mid),
      .groups = "drop") %>%
    select(
      Dataset,
      everything())
}

Pair_Summary <- bind_rows(
  summarise_pairs(
    pairs_window,
    "Window_15_33"),
  summarise_pairs(
    pairs_full,
    "Full_pregnancy"))

###############################################################################
# 13. Sensibilidad por número mínimo de pares
###############################################################################

min_pairs_sensitivity <- function(
  pairs,
  dataset) {

  subject_pairs <- pairs %>%
    count(
      host_subject_id,
      name = "N_pairs")

  out <- list()
  k <- 1

  for (min_pairs in
    MIN_PAIRS_VALUES) {

    eligible <- subject_pairs %>%
      filter(
        N_pairs >= min_pairs) %>%
      pull(
        host_subject_id)

    d <- pairs %>%
      filter(
        host_subject_id %in%
          eligible)

    for (analysis in analyses) {

      da <- d

      if (analysis == "Stanford")
        da <- da %>%
          filter(
            cohort == "Stanford")

      if (analysis == "UAB")
        da <- da %>%
          filter(
            cohort == "UAB")

      n_term <- n_distinct(
        da$host_subject_id[
          da$Delivery ==
            "Term"])

      n_preterm <- n_distinct(
        da$host_subject_id[
          da$Delivery ==
            "Preterm"])

      if (
        n_term >= 2 &&
          n_preterm >= 2) {

        model <- fit_pair_model(
          d,
          analysis,
          TRUE)

        res <- extract_delivery(
          model,
          dataset,
          analysis)

        res$Min_pairs <-
          min_pairs

        res$N_Term_subjects <-
          n_term

        res$N_Preterm_subjects <-
          n_preterm

        out[[k]] <- res

      } else {

        out[[k]] <- data.frame(
          Dataset = dataset,
          Analysis = analysis,
          Contrast =
            "Preterm_vs_Term",
          Beta = NA_real_,
          SE = NA_real_,
          DF = NA_real_,
          T_value = NA_real_,
          Pvalue = NA_real_,
          CI_low = NA_real_,
          CI_high = NA_real_,
          N_pairs = nrow(da),
          N_subjects =
            n_distinct(
              da$host_subject_id),
          Direction =
            "Insufficient_sample_size",
          Min_pairs = min_pairs,
          N_Term_subjects =
            n_term,
          N_Preterm_subjects =
            n_preterm,
          stringsAsFactors = FALSE)
      }

      k <- k + 1
    }
  }

  bind_rows(out)
}

MinPairs_Results <- bind_rows(
  min_pairs_sensitivity(
    pairs_window,
    "Window_15_33"),
  min_pairs_sensitivity(
    pairs_full,
    "Full_pregnancy"))

cat(
  "\n=== SENSIBILIDAD: NÚMERO MÍNIMO DE PARES ===\n")

print(
  MinPairs_Results)

###############################################################################
# 14. Gráficos descriptivos
###############################################################################

plot_window <- pairs_window %>%
  mutate(
    Group =
      interaction(
        cohort,
        Delivery,
        sep = "-"))

p_bc <- ggplot(
  plot_window,
  aes(
    Group,
    BC)) +
  geom_boxplot(
    outlier.shape = NA) +
  geom_jitter(
    width = .15,
    alpha = .35) +
  theme_bw() +
  labs(
    title =
      "Cambio composicional entre muestras consecutivas",
    subtitle =
      "Ventana gestacional 15-33 semanas",
    x = NULL,
    y = "Bray-Curtis")

p_interval <- ggplot(
  plot_window,
  aes(
    Interval,
    BC,
    colour =
      Delivery)) +
  geom_point(
    alpha = .35) +
  geom_smooth(
    method = "loess",
    se = TRUE) +
  facet_wrap(
    ~cohort) +
  theme_bw() +
  labs(
    title =
      "Bray-Curtis según intervalo entre muestras",
    subtitle =
      "Ventana gestacional 15-33 semanas",
    x =
      "Intervalo entre muestras, semanas",
    y =
      "Bray-Curtis")

p_position <- ggplot(
  plot_window,
  aes(
    GestWeek_mid,
    BC,
    colour =
      Delivery)) +
  geom_point(
    alpha = .35) +
  geom_smooth(
    method = "loess",
    se = TRUE) +
  facet_wrap(
    ~cohort) +
  theme_bw() +
  labs(
    title =
      "Bray-Curtis según posición gestacional del par",
    subtitle =
      "Ventana gestacional 15-33 semanas",
    x =
      "Semana gestacional media del par",
    y =
      "Bray-Curtis")

ggsave(
  file.path(
    images_dir,
    "TemporalStability_Bray_Window.png"),
  p_bc,
  width = 8,
  height = 5,
  dpi = 300)

ggsave(
  file.path(
    images_dir,
    "TemporalStability_Interval_vs_Bray.png"),
  p_interval,
  width = 8,
  height = 5,
  dpi = 300)

ggsave(
  file.path(
    images_dir,
    "TemporalStability_GestWeekMid_vs_Bray.png"),
  p_position,
  width = 8,
  height = 5,
  dpi = 300)

###############################################################################
# 15. Guardar resultados
###############################################################################

write_csv(
  pairs_window,
  file.path(
    results_dir,
    "TemporalStability_Pairs_Window15_33.csv"))

write_csv(
  pairs_full,
  file.path(
    results_dir,
    "TemporalStability_Pairs_FullPregnancy.csv"))

write_csv(
  Pair_Availability,
  file.path(
    results_dir,
    "TemporalStability_PairAvailability.csv"))

write_csv(
  Pair_Availability_Summary,
  file.path(
    results_dir,
    "TemporalStability_PairAvailability_Summary.csv"))

write_csv(
  Zero_Interval_Summary,
  file.path(
    results_dir,
    "TemporalStability_ZeroIntervals.csv"))

write_csv(
  BC_Distribution,
  file.path(
    results_dir,
    "TemporalStability_BC_Distribution.csv"))

write_csv(
  Diagnostic_Results,
  file.path(
    results_dir,
    "TemporalStability_ModelDiagnostics.csv"))

write_csv(
  Primary_Results,
  file.path(
    results_dir,
    "TemporalStability_Primary_Window15_33.csv"))

write_csv(
  FullPregnancy_Results,
  file.path(
    results_dir,
    "TemporalStability_Sensitivity_FullPregnancy.csv"))

write_csv(
  Pair_Summary,
  file.path(
    results_dir,
    "TemporalStability_PairSummary.csv"))

write_csv(
  MinPairs_Results,
  file.path(
    results_dir,
    "TemporalStability_Sensitivity_MinPairs.csv"))

saveRDS(
  list(
    pairs_window =
      pairs_window,
    pairs_full =
      pairs_full,
    Primary_Models =
      Primary_Models,
    Full_Models =
      Full_Models,
    Primary_Results =
      Primary_Results,
    FullPregnancy_Results =
      FullPregnancy_Results,
    Diagnostic_Results =
      Diagnostic_Results,
    MinPairs_Results =
      MinPairs_Results),
  file.path(
    results_dir,
    "TemporalStability_Final_Objects.rds"))

###############################################################################
# 16. Reproducibilidad
###############################################################################

versions <- data.frame(
  Package = c(
    "R",
    "phyloseq",
    "vegan",
    "dplyr",
    "tidyr",
    "readr",
    "ggplot2",
    "nlme"),
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
        "vegan")),
    as.character(
      packageVersion(
        "dplyr")),
    as.character(
      packageVersion(
        "tidyr")),
    as.character(
      packageVersion(
        "readr")),
    as.character(
      packageVersion(
        "ggplot2")),
    as.character(
      packageVersion(
        "nlme"))),
  stringsAsFactors = FALSE)

write_csv(
  versions,
  file.path(
    results_dir,
    "TemporalStability_SoftwareVersions.csv"))

cat("\n=== VERSIONES ===\n")
print(versions)

cat("\n========================================\n")
cat("ANÁLISIS DE ESTABILIDAD TEMPORAL FINALIZADO\n")
cat("========================================\n")

cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

###############################################################################
# FIN
###############################################################################
