###############################################################################
#                                    TFM
###############################################################################
#              Script 14: CST y dinámica temporal con VALENCIA
###############################################################################

# Objetivo:
# Evaluar la distribución de CST y sus transiciones durante el embarazo.
#
# Análisis principal:
# - ventana gestacional común de 15-33 semanas;
# - modelos logísticos mixtos one-vs-rest para cada CST con prevalencia >=5 %;
# - modelos globales aditivo y con interacción Delivery * cohort;
# - análisis estratificado por Stanford/UAB;
# - ajuste por semana gestacional;
# - FDR BH dentro de cada bloque Global / Stanford / UAB.
#
# Dinámica temporal:
# - pares de muestras consecutivas de una misma participante;
# - cambio de CST (sí/no) como respuesta;
# - modelos globales aditivo y con interacción Delivery * cohort;
# - ajuste por intervalo entre muestras y semana gestacional media del par.
#
# Sensibilidades:
# - todo el embarazo con los mismos modelos;
# - score VALENCIA >=0.1, usado solo como sensibilidad y no como criterio
#   primario de exclusión;
# - en transiciones, la sensibilidad por score conserva únicamente pares
#   originales en los que ambas muestras tienen score >=0.1.
#
# Las covariables continuas se estandarizan para mejorar la estabilidad
# numérica de los modelos. Los modelos no convergentes, singulares o con
# indicios claros de separación se conservan en la auditoría, pero no se usan
# para inferencia ni para el cálculo de FDR.

###############################################################################

rm(list = ls())

library(phyloseq)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(lme4)

###############################################################################
# 1. Directorios y parámetros
###############################################################################

results_dir <- "VaginalMicrobiome/results/VALENCIA/Analysis"
images_dir <- "VaginalMicrobiome/images/VALENCIA"

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
PREV_CUT <- 0.05
SCORE_CUT <- 0.1
ALPHA <- 0.05

analyses <- c(
  "Global_additive",
  "Global_interaction",
  "Stanford",
  "UAB")

###############################################################################
# 2. Cargar objeto phyloseq
###############################################################################

ps_file <- "VaginalMicrobiome/cleaned/ps_valencia.rds"

if (!file.exists(ps_file))
  stop("No se encontró ps_valencia.rds.")

ps <- readRDS(ps_file)

cat("\n=== OBJETO PHYLOSEQ ===\n")
cat("Muestras:",nsamples(ps),"\n")# 2179
cat("ASVs:",ntaxa(ps),"\n")# 6337

###############################################################################
# 3. Preparar metadata
###############################################################################

meta <- data.frame(
  sample_data(ps),
  stringsAsFactors = FALSE)

meta$SampleID <- rownames(meta)

required <- c(
  "host_subject_id",
  "Groups",
  "cohort",
  "term_vs_preterm_delivery",
  "gest_week_collection",
  "gest_wk_delivery",
  "VALENCIA_CST",
  "VALENCIA_subCST",
  "VALENCIA_score")

missing <- setdiff(
  required,
  names(meta))

if (length(missing) > 0)
  stop(
    paste(
      "Faltan variables:",
      paste(
        missing,
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
      factor(host_subject_id),
    CST = factor(
      VALENCIA_CST),
    subCST = factor(
      VALENCIA_subCST),
    VALENCIA_score =
      as.numeric(
        VALENCIA_score))

if (
  anyNA(meta$cohort) ||
    anyNA(meta$Delivery))
  stop(
    "Hay muestras sin cohorte o desenlace correctamente definido.")

if (anyNA(meta$CST))
  stop(
    "Hay muestras sin CST asignado por VALENCIA.")

###############################################################################
# 4. Control descriptivo del score VALENCIA
###############################################################################

Score_Overall <- meta %>%
  summarise(
    N_samples = n(),
    Min = min(
      VALENCIA_score,
      na.rm = TRUE),
    Q1 = quantile(
      VALENCIA_score,
      .25,
      na.rm = TRUE),
    Median = median(
      VALENCIA_score,
      na.rm = TRUE),
    Mean = mean(
      VALENCIA_score,
      na.rm = TRUE),
    Q3 = quantile(
      VALENCIA_score,
      .75,
      na.rm = TRUE),
    Max = max(
      VALENCIA_score,
      na.rm = TRUE),
    N_score_lt_0.1 =
      sum(
        VALENCIA_score <
          SCORE_CUT,
        na.rm = TRUE),
    Pct_score_lt_0.1 =
      100 * mean(
        VALENCIA_score <
          SCORE_CUT,
        na.rm = TRUE),
    N_score_lt_0.5 =
      sum(
        VALENCIA_score < .5,
        na.rm = TRUE),
    Pct_score_lt_0.5 =
      100 * mean(
        VALENCIA_score < .5,
        na.rm = TRUE))

Score_By_CST <- meta %>%
  group_by(CST) %>%
  summarise(
    N_samples = n(),
    Median_score =
      median(
        VALENCIA_score,
        na.rm = TRUE),
    IQR_score =
      IQR(
        VALENCIA_score,
        na.rm = TRUE),
    Mean_score =
      mean(
        VALENCIA_score,
        na.rm = TRUE),
    Pct_score_lt_0.1 =
      100 * mean(
        VALENCIA_score <
          SCORE_CUT,
        na.rm = TRUE),
    .groups = "drop")

Score_By_Group <- meta %>%
  group_by(
    cohort,
    Delivery) %>%
  summarise(
    N_samples = n(),
    Median_score =
      median(
        VALENCIA_score,
        na.rm = TRUE),
    IQR_score =
      IQR(
        VALENCIA_score,
        na.rm = TRUE),
    Pct_score_lt_0.1 =
      100 * mean(
        VALENCIA_score <
          SCORE_CUT,
        na.rm = TRUE),
    .groups = "drop")

cat("\n=== SCORE VALENCIA: GLOBAL ===\n")
print(Score_Overall)

cat("\n=== SCORE VALENCIA POR CST ===\n")
print(Score_By_CST)

cat("\n=== SCORE VALENCIA POR GRUPO ===\n")
print(Score_By_Group)

p_score <- ggplot(
  meta,
  aes(
    x = CST,
    y = VALENCIA_score)) +
  geom_boxplot(
    outlier.shape = NA) +
  geom_jitter(
    width = .15,
    alpha = .2,
    size = .7) +
  facet_wrap(~cohort) +
  theme_bw() +
  labs(
    title =
      "Score de similitud de VALENCIA por CST",
    x = "CST",
    y = "VALENCIA score")

ggsave(
  file.path(
    images_dir,
    "VALENCIA_Score_ByCST.png"),
  p_score,
  width = 9,
  height = 5.5,
  dpi = 300)

###############################################################################
# 5. Crear conjuntos analíticos
###############################################################################

meta_valid <- meta %>%
  filter(
    !Temporal_inconsistency,
    !is.na(GestWeek))

meta_window <- meta_valid %>%
  filter(
    GestWeek >= GW_MIN,
    GestWeek <= GW_MAX)

meta_full <- meta_valid

meta_window_score <- meta_window %>%
  filter(
    VALENCIA_score >=
      SCORE_CUT)

meta_full_score <- meta_full %>%
  filter(
    VALENCIA_score >=
      SCORE_CUT)

cat("\n=== CONJUNTOS ANALÍTICOS ===\n")

cat(
  "Ventana 15-33: ",
  nrow(meta_window),
  " muestras | ",
  n_distinct(
    meta_window$host_subject_id),
  " mujeres\n",
  sep = "")
# Ventana 15-33: 1658 muestras | 135 mujeres

cat(
  "Todo embarazo: ",
  nrow(meta_full),
  " muestras | ",
  n_distinct(
    meta_full$host_subject_id),
  " mujeres\n",
  sep = "")
# Todo embarazo: 2176 muestras | 135 mujeres

cat(
  "Ventana 15-33 con score >= ",
  SCORE_CUT,": ",
  nrow(meta_window_score),
  " muestras | ",
  n_distinct(
    meta_window_score$host_subject_id),
  " mujeres\n",
  sep = "")
# Ventana 15-33 con score >= 0.1: 1210 muestras | 119 mujeres

cat(
  "Todo embarazo con score >= ",
  SCORE_CUT,": ",
  nrow(meta_full_score),
  " muestras | ",
  n_distinct(
    meta_full_score$host_subject_id),
  " mujeres\n",
  sep = "")
# Todo embarazo con score >= 0.1: 1555 muestras | 120 mujeres

###############################################################################
# 6. Funciones auxiliares de estandarización y diagnóstico
###############################################################################

safe_scale <- function(x) {

  sx <- sd(
    x,
    na.rm = TRUE)

  if (
    is.na(sx) ||
      sx == 0)
    return(
      rep(
        0,
        length(x)))

  as.numeric(
    scale(x))
}

prepare_cst_data <- function(dat) {

  dat %>%
    mutate(
      GestWeek_z =
        safe_scale(
          GestWeek))
}

prepare_pair_data <- function(dat) {

  dat %>%
    mutate(
      Interval_z =
        safe_scale(
          Interval),
      GestWeek_mid_z =
        safe_scale(
          GestWeek_mid))
}

diagnose_glmer <- function(
  model,
  dataset,
  analysis,
  outcome) {

  singular <- isSingular(
    model,
    tol = 1e-4)

  conv_messages <-
    model@optinfo$conv$lme4$messages

  converged <- is.null(
    conv_messages)

  fitted_values <- fitted(
    model)

  beta_max <- max(
    abs(
      fixef(model)),
    na.rm = TRUE)

  pct_extreme <- 100 * mean(
    fitted_values < 1e-6 |
      fitted_values >
        1 - 1e-6)

  separation_warning <-
    beta_max > 10 ||
    pct_extreme > 10

  valid_inference <-
    converged &&
    !singular &&
    !separation_warning

  data.frame(
    Dataset = dataset,
    Analysis = analysis,
    Outcome = outcome,
    N_obs = nobs(model),
    Singular = singular,
    Converged = converged,
    Min_fitted =
      min(fitted_values),
    Max_fitted =
      max(fitted_values),
    Pct_fitted_extreme =
      pct_extreme,
    Max_abs_beta =
      beta_max,
    Separation_warning =
      separation_warning,
    Valid_inference =
      valid_inference,
    Convergence_message =
      ifelse(
        converged,
        NA_character_,
        paste(
          conv_messages,
          collapse = " | ")),
    stringsAsFactors = FALSE)
}

###############################################################################
# 7. Distribución descriptiva de CST y sub-CST
###############################################################################

cst_distribution <- function(
  dat,
  dataset) {

  dat %>%
    count(
      cohort,
      Delivery,
      CST,
      name = "N_samples") %>%
    group_by(
      cohort,
      Delivery) %>%
    mutate(
      Proportion =
        N_samples /
        sum(N_samples),
      Dataset = dataset) %>%
    ungroup() %>%
    select(
      Dataset,
      everything())
}

subcst_distribution <- function(
  dat,
  dataset) {

  dat %>%
    count(
      cohort,
      Delivery,
      subCST,
      name = "N_samples") %>%
    group_by(
      cohort,
      Delivery) %>%
    mutate(
      Proportion =
        N_samples /
        sum(N_samples),
      Dataset = dataset) %>%
    ungroup() %>%
    select(
      Dataset,
      everything())
}

CST_Distribution <- bind_rows(
  cst_distribution(
    meta_window,
    "Window_15_33"),
  cst_distribution(
    meta_full,
    "Full_pregnancy"),
  cst_distribution(
    meta_window_score,
    "Window_15_33_score_ge_0.1"))

subCST_Distribution <- bind_rows(
  subcst_distribution(
    meta_window,
    "Window_15_33"),
  subcst_distribution(
    meta_full,
    "Full_pregnancy"),
  subcst_distribution(
    meta_window_score,
    "Window_15_33_score_ge_0.1"))

p_cst <- meta_window %>%
  count(
    cohort,
    Delivery,
    CST,
    name = "N") %>%
  group_by(
    cohort,
    Delivery) %>%
  mutate(
    Proportion =
      N / sum(N)) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = Delivery,
      y = Proportion,
      fill = CST)) +
  geom_col() +
  facet_wrap(~cohort) +
  scale_y_continuous(
    labels =
      scales::percent) +
  theme_bw() +
  labs(
    title =
      "Distribución de CST por cohorte y desenlace",
    subtitle =
      "Ventana gestacional 15-33 semanas",
    x = NULL,
    y = "Proporción de muestras",
    fill = "CST")

ggsave(
  file.path(
    images_dir,
    "VALENCIA_CST_Distribution_Window15_33.png"),
  p_cst,
  width = 8,
  height = 5.5,
  dpi = 300)

###############################################################################
# 8. Modelos logísticos mixtos por CST
###############################################################################

fit_cst_model <- function(
  dat,
  cst,
  analysis) {

  d <- dat

  if (analysis == "Stanford")
    d <- d %>%
      filter(
        cohort == "Stanford")

  if (analysis == "UAB")
    d <- d %>%
      filter(
        cohort == "UAB")

  d <- droplevels(d)

  d$CST_present <- as.integer(
    as.character(d$CST) ==
      cst)

  prevalence <- mean(
    d$CST_present)

  if (
    prevalence < PREV_CUT ||
      prevalence >
        1 - PREV_CUT)
    return(NULL)

  formula_model <- if (
    analysis == "Global_additive") {

    CST_present ~
      Delivery +
      cohort +
      GestWeek_z +
      (1|host_subject_id)

  } else if (
    analysis == "Global_interaction") {

    CST_present ~
      Delivery * cohort +
      GestWeek_z +
      (1|host_subject_id)

  } else {

    CST_present ~
      Delivery +
      GestWeek_z +
      (1|host_subject_id)
  }

  model <- tryCatch(
    glmer(
      formula_model,
      data = d,
      family = binomial,
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(
          maxfun = 2e5))),
    error = function(e)
      NULL)

  if (is.null(model))
    return(NULL)

  list(
    model = model,
    prevalence = prevalence,
    n_samples = nrow(d),
    n_subjects =
      n_distinct(
        d$host_subject_id))
}

extract_cst_effect <- function(
  fit,
  cst,
  analysis,
  dataset) {

  if (is.null(fit))
    return(NULL)

  tt <- summary(
    fit$model)$coefficients

  term <- if (
    analysis == "Global_interaction") {
    "DeliveryPreterm:cohortUAB"
  } else {
    "DeliveryPreterm"
  }

  if (!term %in%
    rownames(tt))
    return(NULL)

  beta <- tt[
    term,
    "Estimate"]

  se <- tt[
    term,
    "Std. Error"]

  p <- tt[
    term,
    "Pr(>|z|)"]

  data.frame(
    Dataset = dataset,
    Analysis = analysis,
    CST = cst,
    Prevalence =
      fit$prevalence,
    N_samples =
      fit$n_samples,
    N_subjects =
      fit$n_subjects,
    Beta = beta,
    SE = se,
    OR = exp(beta),
    CI_low =
      exp(
        beta -
          1.96 * se),
    CI_high =
      exp(
        beta +
          1.96 * se),
    Pvalue_raw = p,
    Direction = if (
      analysis == "Global_interaction") {
      ifelse(
        beta > 0,
        "Stronger_Preterm_association_in_UAB",
        "Weaker_Preterm_association_in_UAB")
    } else {
      ifelse(
        beta > 0,
        "Higher_odds_in_Preterm",
        "Lower_odds_in_Preterm")
    },
    stringsAsFactors = FALSE)
}

run_cst_models <- function(
  dat,
  dataset) {

  d <- prepare_cst_data(
    dat)

  csts <- levels(
    droplevels(
      d$CST))

  results <- list()
  diagnostics <- list()
  models <- list()

  k <- 1
  j <- 1

  for (analysis in analyses) {

    for (cst in csts) {

      fit <- fit_cst_model(
        d,
        cst,
        analysis)

      if (is.null(fit))
        next

      res <- extract_cst_effect(
        fit,
        cst,
        analysis,
        dataset)

      diag <- diagnose_glmer(
        fit$model,
        dataset,
        analysis,
        paste0(
          "CST_",cst))

      if (!is.null(res)) {

        res$Valid_inference <-
          diag$Valid_inference

        res$Pvalue <- ifelse(
          diag$Valid_inference,
          res$Pvalue_raw,
          NA_real_)

        results[[k]] <- res
        k <- k + 1
      }

      diagnostics[[j]] <- diag

      models[[
        paste(
          dataset,
          analysis,
          cst,
          sep = "_")
      ]] <- fit$model

      j <- j + 1
    }
  }

  results_df <- bind_rows(
    results) %>%
    group_by(Analysis) %>%
    mutate(
      FDR = p.adjust(
        Pvalue,
        method = "BH"),
      Significant_FDR =
        !is.na(FDR) &
        FDR < ALPHA) %>%
    ungroup()

  list(
    results = results_df,
    diagnostics =
      bind_rows(
        diagnostics),
    models = models)
}

CST_Primary <- run_cst_models(
  meta_window,
  "Window_15_33")

CST_Full <- run_cst_models(
  meta_full,
  "Full_pregnancy")

CST_Score <- run_cst_models(
  meta_window_score,
  "Window_15_33_score_ge_0.1")

CST_Models_Primary <-
  CST_Primary$results

CST_Models_Full <-
  CST_Full$results

CST_Models_Score <-
  CST_Score$results

CST_Diagnostics <- bind_rows(
  CST_Primary$diagnostics,
  CST_Full$diagnostics,
  CST_Score$diagnostics)

cat("\n=== CST: MODELOS PRINCIPALES ===\n")
print(CST_Models_Primary)

cat("\n=== CST: SENSIBILIDAD SCORE ===\n")
print(CST_Models_Score)

cat("\n=== CST: DIAGNÓSTICOS ===\n")
print(
  as.data.frame(
    CST_Diagnostics),
  row.names = FALSE)

###############################################################################
# 9. Construir pares consecutivos para transiciones CST
###############################################################################

make_cst_pairs <- function(dat) {

  dat %>%
    arrange(
      host_subject_id,
      GestWeek,
      SampleID) %>%
    group_by(
      host_subject_id) %>%
    mutate(
      SampleID_next =
        lead(SampleID),
      CST_next =
        lead(CST),
      subCST_next =
        lead(subCST),
      Score_next =
        lead(VALENCIA_score),
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
        (GestWeek +
          GestWeek_next) / 2,
      CST_change =
        as.integer(
          CST != CST_next),
      subCST_change =
        as.integer(
          subCST !=
            subCST_next),
      Transition =
        paste(
          CST,
          CST_next,
          sep = " -> "))
}

pairs_window <- make_cst_pairs(
  meta_window)

pairs_full <- make_cst_pairs(
  meta_full)

# No se reconstruyen pares tras excluir muestras de score bajo.
pairs_window_score <- pairs_window %>%
  filter(
    VALENCIA_score >=
      SCORE_CUT,
    Score_next >=
      SCORE_CUT)

pairs_full_score <- pairs_full %>%
  filter(
    VALENCIA_score >=
      SCORE_CUT,
    Score_next >=
      SCORE_CUT)

cat(
  "\n=== PARES CON SCORE >= ",
  SCORE_CUT,
  " EN AMBAS MUESTRAS ===\n",
  sep = "")

cat(
  "Ventana 15-33:",
  nrow(pairs_window_score),
  "pares |",
  n_distinct(
    pairs_window_score$host_subject_id),
  "mujeres\n")
# === PARES CON SCORE >= 0.1 EN AMBAS MUESTRAS ===
# Ventana 15-33: 1047 pares | 109 mujeres

cat(
  "Todo embarazo:",
  nrow(pairs_full_score),
  "pares |",
  n_distinct(
    pairs_full_score$host_subject_id),
  "mujeres\n")
# Todo embarazo: 1371 pares | 113 mujeres

###############################################################################
# 10. Disponibilidad de pares
###############################################################################

pair_availability <- function(
  dat,
  pairs,
  dataset) {

  subjects <- dat %>%
    distinct(
      host_subject_id,
      cohort,
      Delivery)

  pair_counts <- pairs %>%
    count(
      host_subject_id,
      name = "N_pairs")

  subjects %>%
    left_join(
      pair_counts,
      by =
        "host_subject_id") %>%
    mutate(
      N_pairs =
        replace_na(
          N_pairs,
          0L),
      Has_pair =
        N_pairs >= 1,
      Dataset = dataset)
}

Pair_Availability <- bind_rows(
  pair_availability(
    meta_window,
    pairs_window,
    "Window_15_33"),
  pair_availability(
    meta_full,
    pairs_full,
    "Full_pregnancy"),
  pair_availability(
    meta_window,
    pairs_window_score,
    "Window_15_33_score_ge_0.1"),
  pair_availability(
    meta_full,
    pairs_full_score,
    "Full_pregnancy_score_ge_0.1"))

Pair_Availability_Summary <-
  Pair_Availability %>%
  group_by(
    Dataset,
    cohort,
    Delivery) %>%
  summarise(
    N_subjects = n(),
    N_with_pairs =
      sum(Has_pair),
    N_without_pairs =
      sum(!Has_pair),
    Median_pairs =
      median(N_pairs),
    .groups = "drop")

###############################################################################
# 11. Modelos mixtos de cambio de CST
###############################################################################

fit_transition_model <- function(
  pairs,
  analysis) {

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

  formula_model <- if (
    analysis == "Global_additive") {

    CST_change ~
      Delivery +
      cohort +
      Interval_z +
      GestWeek_mid_z +
      (1|host_subject_id)

  } else if (
    analysis == "Global_interaction") {

    CST_change ~
      Delivery * cohort +
      Interval_z +
      GestWeek_mid_z +
      (1|host_subject_id)

  } else {

    CST_change ~
      Delivery +
      Interval_z +
      GestWeek_mid_z +
      (1|host_subject_id)
  }

  tryCatch(
    glmer(
      formula_model,
      data = d,
      family = binomial,
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(
          maxfun = 2e5))),
    error = function(e)
      NULL)
}

extract_transition_effect <- function(
  model,
  dataset,
  analysis) {

  if (is.null(model))
    return(NULL)

  tt <- summary(
    model)$coefficients

  term <- if (
    analysis == "Global_interaction") {
    "DeliveryPreterm:cohortUAB"
  } else {
    "DeliveryPreterm"
  }

  if (!term %in%
    rownames(tt))
    return(NULL)

  beta <- tt[
    term,
    "Estimate"]

  se <- tt[
    term,
    "Std. Error"]

  p <- tt[
    term,
    "Pr(>|z|)"]

  data.frame(
    Dataset = dataset,
    Analysis = analysis,
    Beta = beta,
    SE = se,
    OR = exp(beta),
    CI_low =
      exp(
        beta -
          1.96 * se),
    CI_high =
      exp(
        beta +
          1.96 * se),
    Pvalue_raw = p,
    N_pairs = nobs(model),
    N_subjects =
      n_distinct(
        getME(
          model,
          "flist")[[1]]),
    Direction = if (
      analysis == "Global_interaction") {
      ifelse(
        beta > 0,
        "Stronger_Preterm_transition_association_in_UAB",
        "Weaker_Preterm_transition_association_in_UAB")
    } else {
      ifelse(
        beta > 0,
        "More_transitions_in_Preterm",
        "Fewer_transitions_in_Preterm")
    },
    stringsAsFactors = FALSE)
}

run_transition_models <- function(
  pairs,
  dataset) {

  d <- prepare_pair_data(
    pairs)

  models <- list()
  results <- list()
  diagnostics <- list()

  k <- 1

  for (analysis in analyses) {

    model <- fit_transition_model(
      d,
      analysis)

    if (is.null(model))
      next

    res <- extract_transition_effect(
      model,
      dataset,
      analysis)

    diag <- diagnose_glmer(
      model,
      dataset,
      analysis,
      "CST_change")

    if (!is.null(res)) {

      res$Valid_inference <-
        diag$Valid_inference

      res$Pvalue <- ifelse(
        diag$Valid_inference,
        res$Pvalue_raw,
        NA_real_)

      results[[k]] <- res
    }

    diagnostics[[k]] <- diag
    models[[analysis]] <- model

    k <- k + 1
  }

  list(
    models = models,
    results =
      bind_rows(results),
    diagnostics =
      bind_rows(diagnostics))
}

Transitions_Primary <- run_transition_models(
  pairs_window,
  "Window_15_33")

Transitions_Full <- run_transition_models(
  pairs_full,
  "Full_pregnancy")

Transitions_Score <- run_transition_models(
  pairs_window_score,
  "Window_15_33_score_ge_0.1")

Transitions_Full_Score <- run_transition_models(
  pairs_full_score,
  "Full_pregnancy_score_ge_0.1")

Transition_Diagnostics <- bind_rows(
  Transitions_Primary$diagnostics,
  Transitions_Full$diagnostics,
  Transitions_Score$diagnostics,
  Transitions_Full_Score$diagnostics)

cat("\n=== TRANSICIONES CST: PRINCIPAL ===\n")
print(
  Transitions_Primary$results)

cat("\n=== TRANSICIONES CST: SENSIBILIDAD SCORE ===\n")
print(
  Transitions_Score$results)

cat("\n=== DIAGNÓSTICOS DE TRANSICIONES ===\n")
print(
  as.data.frame(
    Transition_Diagnostics),
  row.names = FALSE)

###############################################################################
# 12. Matrices descriptivas de transición
###############################################################################

transition_table <- function(
  pairs,
  dataset) {

  pairs %>%
    count(
      cohort,
      Delivery,
      CST,
      CST_next,
      name = "N") %>%
    group_by(
      cohort,
      Delivery,
      CST) %>%
    mutate(
      Transition_probability =
        N / sum(N),
      Dataset = dataset) %>%
    ungroup() %>%
    select(
      Dataset,
      everything())
}

Transition_Table <- bind_rows(
  transition_table(
    pairs_window,
    "Window_15_33"),
  transition_table(
    pairs_full,
    "Full_pregnancy"),
  transition_table(
    pairs_window_score,
    "Window_15_33_score_ge_0.1"))

###############################################################################
# 13. Estabilidad individual descriptiva
###############################################################################

individual_stability <- function(
  pairs,
  dataset) {

  pairs %>%
    group_by(
      host_subject_id,
      cohort,
      Delivery) %>%
    summarise(
      Dataset = dataset,
      N_pairs = n(),
      N_changes =
        sum(CST_change),
      Proportion_changes =
        mean(CST_change),
      .groups = "drop")
}

Individual_Stability <- bind_rows(
  individual_stability(
    pairs_window,
    "Window_15_33"),
  individual_stability(
    pairs_full,
    "Full_pregnancy"),
  individual_stability(
    pairs_window_score,
    "Window_15_33_score_ge_0.1"))

p_change <- Individual_Stability %>%
  filter(
    Dataset ==
      "Window_15_33") %>%
  ggplot(
    aes(
      x = Delivery,
      y = Proportion_changes)) +
  geom_boxplot(
    outlier.shape = NA) +
  geom_jitter(
    width = .15,
    alpha = .45) +
  facet_wrap(~cohort) +
  theme_bw() +
  labs(
    title =
      "Transiciones de CST por participante",
    subtitle =
      "Ventana gestacional 15-33 semanas",
    x = NULL,
    y =
      "Proporción de pares con cambio de CST")

ggsave(
  file.path(
    images_dir,
    "VALENCIA_CST_Transitions_Window15_33.png"),
  p_change,
  width = 8,
  height = 5.5,
  dpi = 300)

###############################################################################
# 14. Guardar resultados
###############################################################################

write_csv(
  Score_Overall,
  file.path(
    results_dir,
    "VALENCIA_Score_Overall.csv"))

write_csv(
  Score_By_CST,
  file.path(
    results_dir,
    "VALENCIA_Score_ByCST.csv"))

write_csv(
  Score_By_Group,
  file.path(
    results_dir,
    "VALENCIA_Score_ByGroup.csv"))

write_csv(
  CST_Distribution,
  file.path(
    results_dir,
    "VALENCIA_CST_Distribution.csv"))

write_csv(
  subCST_Distribution,
  file.path(
    results_dir,
    "VALENCIA_subCST_Distribution.csv"))

write_csv(
  CST_Models_Primary,
  file.path(
    results_dir,
    "VALENCIA_CST_Models_Primary.csv"))

write_csv(
  CST_Models_Full,
  file.path(
    results_dir,
    "VALENCIA_CST_Models_FullPregnancy.csv"))

write_csv(
  CST_Models_Score,
  file.path(
    results_dir,
    "VALENCIA_CST_Models_ScoreSensitivity.csv"))

write_csv(
  CST_Diagnostics,
  file.path(
    results_dir,
    "VALENCIA_CST_ModelDiagnostics.csv"))

write_csv(
  Pair_Availability_Summary,
  file.path(
    results_dir,
    "VALENCIA_PairAvailability.csv"))

write_csv(
  Transitions_Primary$results,
  file.path(
    results_dir,
    "VALENCIA_TransitionModels_Primary.csv"))

write_csv(
  Transitions_Full$results,
  file.path(
    results_dir,
    "VALENCIA_TransitionModels_FullPregnancy.csv"))

write_csv(
  Transitions_Score$results,
  file.path(
    results_dir,
    "VALENCIA_TransitionModels_ScoreSensitivity.csv"))

write_csv(
  Transitions_Full_Score$results,
  file.path(
    results_dir,
    "VALENCIA_TransitionModels_FullScoreSensitivity.csv"))

write_csv(
  Transition_Diagnostics,
  file.path(
    results_dir,
    "VALENCIA_TransitionModelDiagnostics.csv"))

write_csv(
  Transition_Table,
  file.path(
    results_dir,
    "VALENCIA_TransitionTable.csv"))

write_csv(
  Individual_Stability,
  file.path(
    results_dir,
    "VALENCIA_IndividualStability.csv"))

###############################################################################
# 15. Guardar objetos
###############################################################################

VALENCIA_results <- list(
  meta_window = meta_window,
  meta_full = meta_full,
  meta_window_score = meta_window_score,
  pairs_window = pairs_window,
  pairs_full = pairs_full,
  pairs_window_score = pairs_window_score,
  CST_Diagnostics = CST_Diagnostics,
  Transition_Diagnostics = Transition_Diagnostics
)

saveRDS(
  VALENCIA_results,
  file.path(
    results_dir,
    "VALENCIA_Analysis_Objects.rds"
  )
)

saveRDS(
  CST_Models_Primary,
  file.path(
    results_dir,
    "CST_Models_Primary.rds"
  )
)

saveRDS(
  CST_Models_Full,
  file.path(
    results_dir,
    "CST_Models_Full.rds"
  )
)

saveRDS(
  CST_Models_Score,
  file.path(
    results_dir,
    "CST_Models_Score.rds"
  )
)

###############################################################################
# 16. Reproducibilidad
###############################################################################

versions <- data.frame(
  Package = c(
    "R",
    "phyloseq",
    "lme4"),
  Version = c(
    paste(
      R.version$major,
      R.version$minor,
      sep = "."),
    as.character(
      packageVersion("phyloseq")),
    as.character(
      packageVersion("lme4"))),
  stringsAsFactors = FALSE)

write_csv(
  versions,
  file.path(
    results_dir,
    "VALENCIA_SoftwareVersions.csv"))

cat("\n=== VERSIONES ===\n")
print(versions)
#    Package Version
# 1        R   4.4.2
# 2 phyloseq  1.50.0
# 3     lme4   2.0.1

cat("\n========================================\n")
cat("ANÁLISIS CST / VALENCIA FINALIZADO\n")
cat("========================================\n")

cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

###############################################################################
# FIN
###############################################################################
