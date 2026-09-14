###############################################################################
#                                    TFM
###############################################################################
#                       Script 09: Diversidad alfa
###############################################################################

# Objetivo:
# Evaluar Shannon y Gini-Simpson (1-D) entre parto a término y pretérmino
# dentro de una ventana gestacional común, conservando las medidas repetidas.
#
# Análisis principal:
# - Unidad de análisis: muestra.
# - Ventana gestacional: 15-33 semanas.
# - Shannon: log1p(Shannon).
# - Gini-Simpson: 1-D sin transformar, como análisis complementario.
# - Global:
#   índice ~ Delivery * cohort + GestWeek + log(profundidad).
# - Stanford/UAB:
#   índice ~ Delivery + GestWeek + log(profundidad).
# - Intercepto aleatorio por participante mediante nlme::lme().
#
# Sensibilidades:
# - ajuste adicional por edad materna;
# - promedio por participante dentro de la ventana 15-33 semanas.
#
# Potencia/MDE:
# - Shannon únicamente;
# - Global: MDE de la interacción Delivery x cohort;
# - Stanford/UAB: MDE del efecto Preterm vs Term;
# - simulación paramétrica basada en el modelo mixto ajustado;
# - comprobación obligatoria bajo H0: tasa de rechazo ~5 %.
#
# No se realiza rarefacción. La profundidad de secuenciación se incorpora
# como covariable log-transformada.
#
# Referencia para no rarefacción:
# McMurdie PJ, Holmes S. PLoS Comput Biol. 2014;10(4):e1003531.

###############################################################################

rm(list = ls())

library(phyloseq)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(nlme)

###############################################################################
# 1. Directorios y parámetros
###############################################################################

results_dir <- "VaginalMicrobiome/results/Alpha_diversity/FINAL"
images_dir <- "VaginalMicrobiome/images/Alpha_diversity/FINAL"

diag_dir <- file.path(
  images_dir,
  "Model_diagnostics_final")

dist_dir <- file.path(
  images_dir,
  "Distribution_diagnostics_final")

power_dir <- file.path(
  results_dir,
  "Power")

power_images <- file.path(
  images_dir,
  "Power")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  images_dir,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  diag_dir,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  dist_dir,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  power_dir,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  power_images,
  recursive = TRUE,
  showWarnings = FALSE)

GW_MIN <- 15
GW_MAX <- 33

POWER_TARGET <- 0.80
N_SIM_POWER <- 1000
SEED_POWER <- 20260913

MDE_GRID <- seq(
  0,
  0.80,
  by = 0.05)

analyses <- c(
  "Global",
  "Stanford",
  "UAB")

indices <- c(
  "Shannon",
  "Gini-Simpson")

###############################################################################
# 2. Cargar objeto y preparar metadata
###############################################################################

ps_file <- "VaginalMicrobiome/cleaned/ps_prevotella_refined.rds"

if (!file.exists(ps_file))
  stop("No se encontró el objeto phyloseq refinado.")

ps <- readRDS(ps_file)

cat("\n=== OBJETO DE ENTRADA ===\n")
cat("Muestras:",nsamples(ps),"\n")
cat("ASVs:",ntaxa(ps),"\n")

meta <- data.frame(
  sample_data(ps),
  stringsAsFactors = FALSE)

meta$SampleID <- rownames(meta)

required <- c(
  "host_subject_id",
  "cohort",
  "term_vs_preterm_delivery",
  "Groups",
  "gest_week_collection",
  "gest_wk_delivery",
  "AGE")

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
  anyNA(meta$Delivery) ||
    anyNA(meta$cohort))
  stop(
    "Existen muestras clínicas sin Delivery o cohorte.")

cat(
  "Registros con muestra posterior al parto:",
  sum(
    meta$Temporal_inconsistency,
    na.rm = TRUE),
  "\n")

###############################################################################
# 3. Calcular diversidad alfa y profundidad
###############################################################################

alpha <- estimate_richness(
  ps,
  measures = c(
    "Shannon","Simpson"))

alpha <- data.frame(
  SampleID = rownames(alpha),
  alpha,
  row.names = NULL)

depth <- data.frame(
  SampleID = sample_names(ps),
  Sequencing_Depth =
    sample_sums(ps),
  stringsAsFactors = FALSE)

alpha_sample <- meta %>%
  left_join(
    alpha,
    by = "SampleID") %>%
  left_join(
    depth,
    by = "SampleID") %>%
  mutate(
    Log_Sequencing_Depth =
      log(Sequencing_Depth),
    Shannon_log1p =
      log1p(Shannon))

if (any(
  alpha_sample$Sequencing_Depth <= 0,
  na.rm = TRUE))
  stop(
    "Existen muestras con profundidad <= 0.")

###############################################################################
# 4. Distribución temporal y ventana gestacional
###############################################################################

gw_samples <- alpha_sample %>%
  filter(
    !Temporal_inconsistency,
    !is.na(GestWeek)) %>%
  mutate(
    GestWeek_plot =
      floor(GestWeek),
    StudyGroup =
      interaction(
        cohort,
        Delivery,
        sep = "-")) %>%
  count(
    StudyGroup,
    GestWeek_plot,
    name = "N_samples")

gw_participants <- alpha_sample %>%
  filter(
    !Temporal_inconsistency,
    !is.na(GestWeek)) %>%
  mutate(
    GestWeek_plot =
      floor(GestWeek),
    StudyGroup =
      interaction(
        cohort,
        Delivery,
        sep = "-")) %>%
  distinct(
    host_subject_id,
    StudyGroup,
    GestWeek_plot) %>%
  count(
    StudyGroup,
    GestWeek_plot,
    name = "N_participants")

group_levels <- c(
  "Stanford-Term",
  "Stanford-Preterm",
  "UAB-Term",
  "UAB-Preterm")

gw_samples$StudyGroup <- factor(
  gw_samples$StudyGroup,
  levels = group_levels)

gw_participants$StudyGroup <- factor(
  gw_participants$StudyGroup,
  levels = group_levels)

p_gw_samples <- ggplot(
  gw_samples,
  aes(
    x = GestWeek_plot,
    y = N_samples,
    colour = StudyGroup,
    group = StudyGroup)) +
  annotate(
    "rect",
    xmin = GW_MIN,
    xmax = GW_MAX,
    ymin = -Inf,
    ymax = Inf,
    alpha = .08) +
  geom_line(
    linewidth = .8) +
  geom_point(
    size = 2) +
  labs(
    title =
      "Distribución de muestras por semana gestacional",
    subtitle = paste0(
      "Área sombreada: ventana analítica ",
      GW_MIN,"-",GW_MAX,
      " semanas"),
    x = "Semana gestacional",
    y = "Número de muestras",
    colour = "Grupo") +
  theme_bw()

p_gw_participants <- ggplot(
  gw_participants,
  aes(
    x = GestWeek_plot,
    y = N_participants,
    colour = StudyGroup,
    group = StudyGroup)) +
  annotate(
    "rect",
    xmin = GW_MIN,
    xmax = GW_MAX,
    ymin = -Inf,
    ymax = Inf,
    alpha = .08) +
  geom_hline(
    yintercept = 4,
    linetype = 2) +
  geom_line(
    linewidth = .8) +
  geom_point(
    size = 2) +
  labs(
    title =
      "Distribución de participantes por semana gestacional",
    subtitle = paste0(
      "Área sombreada: ",
      GW_MIN,"-",GW_MAX,
      " semanas; línea discontinua: n = 4"),
    x = "Semana gestacional",
    y = "Número de participantes",
    colour = "Grupo") +
  theme_bw()

ggsave(
  file.path(
    images_dir,
    "GestationalWeek_Samples.png"),
  p_gw_samples,
  width = 9,
  height = 6,
  dpi = 300)

ggsave(
  file.path(
    images_dir,
    "GestationalWeek_Participants.png"),
  p_gw_participants,
  width = 9,
  height = 6,
  dpi = 300)

write_csv(
  gw_samples,
  file.path(
    results_dir,
    "GestationalWeek_Samples.csv"))

write_csv(
  gw_participants,
  file.path(
    results_dir,
    "GestationalWeek_Participants.csv"))

alpha_window <- alpha_sample %>%
  filter(
    !Temporal_inconsistency,
    !is.na(GestWeek),
    GestWeek >= GW_MIN,
    GestWeek <= GW_MAX)

if (nrow(alpha_window) == 0)
  stop(
    "La ventana 15-33 no contiene muestras.")

cat("\n=== VENTANA GESTACIONAL ===\n")
cat("Semanas:",GW_MIN,"-",GW_MAX,"\n")
cat("Muestras:",nrow(alpha_window),"\n")
cat(
  "Mujeres:",
  n_distinct(
    alpha_window$host_subject_id),
  "\n")

print(
  alpha_window %>%
    distinct(
      host_subject_id,
      cohort,
      Delivery) %>%
    count(
      cohort,
      Delivery))

write_csv(
  alpha_window,
  file.path(
    results_dir,
    "AlphaDiversity_Final_Samples.csv"))

###############################################################################
# 5. Descriptivos
###############################################################################

depth_summary <- alpha_window %>%
  group_by(
    cohort,
    Delivery) %>%
  summarise(
    N_samples = n(),
    N_subjects =
      n_distinct(
        host_subject_id),
    Median_depth =
      median(
        Sequencing_Depth),
    IQR_depth =
      IQR(
        Sequencing_Depth),
    Mean_depth =
      mean(
        Sequencing_Depth),
    SD_depth =
      sd(
        Sequencing_Depth),
    Min_depth =
      min(
        Sequencing_Depth),
    Max_depth =
      max(
        Sequencing_Depth),
    .groups = "drop")

alpha_long <- alpha_window %>%
  pivot_longer(
    cols = c(
      Shannon,
      Simpson),
    names_to = "Index",
    values_to = "Value") %>%
  mutate(
    Index = recode(
      Index,
      Simpson =
        "Gini-Simpson"))

skewness_manual <- function(x) {

  x <- x[
    is.finite(x)]

  n <- length(x)

  if (n < 3)
    return(NA_real_)

  s <- sd(x)

  if (s == 0)
    return(0)

  n / (
    (n - 1) *
      (n - 2)) *
    sum(
      (
        (x - mean(x)) /
          s
      )^3)
}

desc_df <- alpha_long %>%
  group_by(
    cohort,
    Delivery,
    Index) %>%
  summarise(
    N_samples = n(),
    N_subjects =
      n_distinct(
        host_subject_id),
    Mean =
      mean(
        Value,
        na.rm = TRUE),
    Median =
      median(
        Value,
        na.rm = TRUE),
    SD =
      sd(
        Value,
        na.rm = TRUE),
    IQR =
      IQR(
        Value,
        na.rm = TRUE),
    Min =
      min(
        Value,
        na.rm = TRUE),
    Max =
      max(
        Value,
        na.rm = TRUE),
    Skewness =
      skewness_manual(
        Value),
    .groups = "drop")

write_csv(
  depth_summary,
  file.path(
    results_dir,
    "AlphaDiversity_Final_DepthSummary.csv"))

write_csv(
  desc_df,
  file.path(
    results_dir,
    "AlphaDiversity_Final_DescriptiveStatistics.csv"))

###############################################################################
# 6. Funciones de modelado y diagnóstico
###############################################################################

fit_lmm <- function(
  dat,
  index,
  analysis,
  age_adjusted = FALSE) {

  response <- if (
    index == "Shannon") {
    "Shannon_log1p"
  } else {
    "Simpson"
  }

  d <- dat %>%
    filter(
      !is.na(
        .data[[response]]),
      !is.na(GestWeek),
      !is.na(
        Log_Sequencing_Depth))

  if (age_adjusted)
    d <- d %>%
      filter(
        !is.na(AGE))

  age_term <- if (
    age_adjusted) {
    " + AGE"
  } else {
    ""
  }

  if (analysis == "Global") {

    rhs <- paste0(
      "Delivery * cohort + GestWeek + ",
      "Log_Sequencing_Depth",
      age_term)

  } else {

    d <- d %>%
      filter(
        cohort == analysis) %>%
      droplevels()

    rhs <- paste0(
      "Delivery + GestWeek + ",
      "Log_Sequencing_Depth",
      age_term)
  }

  model <- lme(
    fixed = as.formula(
      paste(
        response,
        "~",
        rhs)),
    random =
      ~1 | host_subject_id,
    data = d,
    method = "REML",
    na.action = na.omit,
    control = lmeControl(
      opt = "optim"))

  list(
    model = model,
    data = d)
}

extract_main_effect <- function(
  model,
  index,
  analysis,
  adjustment) {

  tt <- summary(
    model)$tTable

  term <- if (
    analysis == "Global") {
    "DeliveryPreterm:cohortUAB"
  } else {
    "DeliveryPreterm"
  }

  if (!term %in%
    rownames(tt))
    stop(
      paste(
        "No se encontró",
        term,
        "en",
        analysis,
        index))

  ci <- intervals(
    model,
    which = "fixed")$fixed[
      term,
      ]

  data.frame(
    Analysis = analysis,
    Index = index,
    Adjustment = adjustment,
    Contrast = ifelse(
      analysis == "Global",
      "Delivery_x_Cohort",
      "Preterm_vs_Term"),
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
      ci["lower"],
    CI_high =
      ci["upper"],
    Direction = if (
      analysis == "Global") {
      NA_character_
    } else if (
      tt[
        term,
        "Value"] > 0) {
      "Higher_in_Preterm"
    } else {
      "Higher_in_Term"
    },
    N_samples =
      nrow(model$data),
    N_subjects =
      n_distinct(
        model$data$host_subject_id),
    stringsAsFactors = FALSE)
}

save_diagnostics <- function(
  model,
  index,
  analysis,
  adjustment) {

  d <- data.frame(
    Fitted =
      fitted(model),
    Residual =
      residuals(
        model,
        type = "normalized"),
    GestWeek =
      model$data$GestWeek)

  re <- data.frame(
    Random_intercept =
      ranef(model)[,1])

  tag <- paste(
    analysis,
    index,
    adjustment,
    sep = "_")

  p1 <- ggplot(
    d,
    aes(
      Fitted,
      Residual)) +
    geom_point(
      alpha = .3) +
    geom_hline(
      yintercept = 0,
      linetype = 2) +
    geom_smooth(
      method = "loess",
      se = FALSE) +
    theme_bw()

  p2 <- ggplot(
    d,
    aes(
      sample = Residual)) +
    stat_qq() +
    stat_qq_line() +
    theme_bw()

  p3 <- ggplot(
    d,
    aes(
      GestWeek,
      Residual)) +
    geom_point(
      alpha = .3) +
    geom_hline(
      yintercept = 0,
      linetype = 2) +
    geom_smooth(
      method = "loess",
      se = FALSE) +
    theme_bw()

  p4 <- ggplot(
    re,
    aes(
      sample =
        Random_intercept)) +
    stat_qq() +
    stat_qq_line() +
    theme_bw()

  ggsave(
    file.path(
      diag_dir,
      paste0(
        tag,
        "_Residuals_vs_Fitted.pdf")),
    p1,
    width = 5.5,
    height = 4.5)

  ggsave(
    file.path(
      diag_dir,
      paste0(
        tag,
        "_QQ_Residuals.pdf")),
    p2,
    width = 5.5,
    height = 4.5)

  ggsave(
    file.path(
      diag_dir,
      paste0(
        tag,
        "_Residuals_vs_GestWeek.pdf")),
    p3,
    width = 5.5,
    height = 4.5)

  ggsave(
    file.path(
      diag_dir,
      paste0(
        tag,
        "_QQ_RandomIntercepts.pdf")),
    p4,
    width = 5.5,
    height = 4.5)

  data.frame(
    Analysis = analysis,
    Index = index,
    Adjustment = adjustment,
    Cor_absResidual_Fitted =
      cor(
        abs(d$Residual),
        d$Fitted,
        method = "spearman",
        use = "complete.obs"),
    Cor_Residual_GestWeek =
      cor(
        d$Residual,
        d$GestWeek,
        method = "spearman",
        use = "complete.obs"),
    Residual_SD =
      sd(
        d$Residual,
        na.rm = TRUE),
    Residual_abs_max =
      max(
        abs(d$Residual),
        na.rm = TRUE),
    N_random_intercepts =
      nrow(re),
    stringsAsFactors = FALSE)
}

###############################################################################
# 7. Análisis principal
###############################################################################

primary_results <- list()
primary_diagnostics <- list()

k <- 1
j <- 1

for (idx in indices) {

  for (a in analyses) {

    fit <- fit_lmm(
      alpha_window,
      idx,
      a,
      age_adjusted = FALSE)

    primary_results[[k]] <-
      extract_main_effect(
        fit$model,
        idx,
        a,
        "Primary")

    primary_diagnostics[[j]] <-
      save_diagnostics(
        fit$model,
        idx,
        a,
        "Primary")

    k <- k + 1
    j <- j + 1
  }
}

Primary_Results <- bind_rows(
  primary_results) %>%
  mutate(
    Significant =
      Pvalue < 0.05)

Primary_Diagnostics <- bind_rows(
  primary_diagnostics)

cat("\n=== ANÁLISIS PRINCIPAL ===\n")
print(Primary_Results)

###############################################################################
# 8. Potencia y MDE para Shannon
###############################################################################

# Global corresponde al término de interacción del modelo global.
# Stanford/UAB corresponden a DeliveryPreterm en los modelos estratificados.

fit_mde_base <- function(
  dat,
  analysis) {

  if (analysis == "Global") {

    d <- dat %>%
      filter(
        !is.na(
          Shannon_log1p),
        !is.na(GestWeek),
        !is.na(
          Log_Sequencing_Depth)) %>%
      droplevels()

    formula_model <-
      Shannon_log1p ~
      Delivery * cohort +
      GestWeek +
      Log_Sequencing_Depth

    target_term <-
      "DeliveryPreterm:cohortUAB"

  } else {

    d <- dat %>%
      filter(
        cohort == analysis,
        !is.na(
          Shannon_log1p),
        !is.na(GestWeek),
        !is.na(
          Log_Sequencing_Depth)) %>%
      droplevels()

    formula_model <-
      Shannon_log1p ~
      Delivery +
      GestWeek +
      Log_Sequencing_Depth

    target_term <-
      "DeliveryPreterm"
  }

  model <- lme(
    fixed = formula_model,
    random =
      ~1 | host_subject_id,
    data = d,
    method = "REML",
    na.action = na.omit,
    control = lmeControl(
      opt = "optim"))

  list(
    data = d,
    model = model,
    formula = formula_model,
    target_term = target_term)
}

simulate_alpha_power <- function(
  base,
  effect,
  n_sim = N_SIM_POWER) {

  d <- base$data
  model <- base$model

  X <- model.matrix(
    formula(
      model),
    data = d)

  beta <- fixed.effects(
    model)

  if (!base$target_term %in%
    names(beta))
    stop(
      "No se encontró el término objetivo del MDE.")

  beta[
    base$target_term] <-
    effect

  mu_fixed <- as.numeric(
    X %*% beta)

  var_random <- as.numeric(
    VarCorr(model)[1,"Variance"])

  sd_random <- sqrt(
    var_random)

  sigma_res <- model$sigma

  subjects <- levels(
    droplevels(
      d$host_subject_id))

  detected <- rep(
    NA,
    n_sim)

  for (i in seq_len(
    n_sim)) {

    u <- rnorm(
      length(subjects),
      0,
      sd_random)

    names(u) <- subjects

    sim <- d

    sim$Shannon_sim <-
      mu_fixed +
      u[
        as.character(
          d$host_subject_id)] +
      rnorm(
        nrow(d),
        0,
        sigma_res)

    fit <- try(
      lme(
        fixed =
          update(
            base$formula,
            Shannon_sim ~ .),
        random =
          ~1 | host_subject_id,
        data = sim,
        method = "REML",
        na.action = na.omit,
        control =
          lmeControl(
            opt = "optim")),
      silent = TRUE)

    if (!inherits(
      fit,
      "try-error")) {

      tt <- summary(
        fit)$tTable

      if (base$target_term %in%
        rownames(tt)) {

        detected[i] <-
          tt[
            base$target_term,
            "p-value"] <
          0.05
      }
    }
  }

  mean(
    detected,
    na.rm = TRUE)
}

estimate_mde <- function(
  dat,
  analysis,
  effects = MDE_GRID) {

  base <- fit_mde_base(
    dat,
    analysis)

  set.seed(
    SEED_POWER)

  power <- vapply(
    effects,
    function(x)
      simulate_alpha_power(
        base,
        x),
    numeric(1))

  out <- data.frame(
    Analysis = analysis,
    Beta = effects,
    Power = power,
    stringsAsFactors = FALSE)

  i <- which(
    out$Power >=
      POWER_TARGET)[1]

  if (is.na(i)) {

    mde <- NA_real_

  } else if (i == 1) {

    mde <- out$Beta[i]

  } else {

    mde <- approx(
      x = out$Power[
        c(i - 1,i)],
      y = out$Beta[
        c(i - 1,i)],
      xout =
        POWER_TARGET)$y
  }

  observed <- abs(
    fixed.effects(
      base$model)[
        base$target_term])

  list(
    Results = out,
    MDE80 = mde,
    Observed = observed,
    Target_term =
      base$target_term)
}

MDE_Global <- estimate_mde(
  alpha_window,
  "Global")

MDE_Stanford <- estimate_mde(
  alpha_window,
  "Stanford")

MDE_UAB <- estimate_mde(
  alpha_window,
  "UAB")

MDE_Results <- bind_rows(
  MDE_Global$Results,
  MDE_Stanford$Results,
  MDE_UAB$Results)

MDE_Summary <- data.frame(
  Analysis = c(
    "Global",
    "Stanford",
    "UAB"),
  Target_term = c(
    MDE_Global$Target_term,
    MDE_Stanford$Target_term,
    MDE_UAB$Target_term),
  Beta_observed_abs = c(
    MDE_Global$Observed,
    MDE_Stanford$Observed,
    MDE_UAB$Observed),
  Beta_MDE80 = c(
    MDE_Global$MDE80,
    MDE_Stanford$MDE80,
    MDE_UAB$MDE80),
  stringsAsFactors = FALSE) %>%
  mutate(
    Fraction_of_MDE =
      Beta_observed_abs /
      Beta_MDE80)

MDE_TypeI_Check <- MDE_Results %>%
  filter(
    Beta == 0) %>%
  transmute(
    Analysis,
    Rejection_rate = Power,
    Expected = 0.05,
    Difference =
      Rejection_rate -
      Expected)

cat("\n=== MDE SHANNON: GLOBAL, STANFORD Y UAB ===\n")
print(MDE_Summary)

cat("\n=== COMPROBACIÓN H0 DEL MDE ===\n")
print(MDE_TypeI_Check)

if (any(
  abs(
    MDE_TypeI_Check$Difference) >
    0.03)) {

  warning(
    paste(
      "La tasa de rechazo bajo H0 se aleja",
      "más de 0.03 del 5 % nominal."))
}

write_csv(
  MDE_Results,
  file.path(
    power_dir,
    "Alpha_Shannon_PowerCurves.csv"))

write_csv(
  MDE_Summary,
  file.path(
    power_dir,
    "Alpha_Shannon_MDE80.csv"))

write_csv(
  MDE_TypeI_Check,
  file.path(
    power_dir,
    "Alpha_Shannon_TypeI_Check.csv"))

p_power <- ggplot(
  MDE_Results,
  aes(
    x = Beta,
    y = Power,
    colour = Analysis)) +
  geom_hline(
    yintercept =
      POWER_TARGET,
    linetype = 2) +
  geom_line() +
  geom_point() +
  scale_y_continuous(
    limits = c(
      0,1)) +
  labs(
    title =
      "Potencia para Shannon",
    subtitle =
      "MDE con potencia objetivo del 80 %",
    x =
      "Tamaño de efecto beta simulado",
    y = "Potencia") +
  theme_bw()

ggsave(
  file.path(
    power_images,
    "Alpha_Shannon_PowerCurves.pdf"),
  p_power,
  width = 7,
  height = 5)

###############################################################################
# 9. Sensibilidad 1: ajuste adicional por edad materna
###############################################################################

age_results <- list()
age_diagnostics <- list()

k <- 1
j <- 1

for (idx in indices) {

  for (a in analyses) {

    fit <- fit_lmm(
      alpha_window,
      idx,
      a,
      age_adjusted = TRUE)

    age_results[[k]] <-
      extract_main_effect(
        fit$model,
        idx,
        a,
        "Age_adjusted")

    age_diagnostics[[j]] <-
      save_diagnostics(
        fit$model,
        idx,
        a,
        "Age_adjusted")

    k <- k + 1
    j <- j + 1
  }
}

AgeAdjusted_Results <- bind_rows(
  age_results) %>%
  mutate(
    Significant =
      Pvalue < 0.05)

AgeAdjusted_Diagnostics <- bind_rows(
  age_diagnostics)

cat(
  "\n=== SENSIBILIDAD AJUSTADA POR EDAD ===\n")

print(
  AgeAdjusted_Results)

###############################################################################
# 10. Sensibilidad 2: promedio por participante
###############################################################################

participant_mean <- alpha_window %>%
  group_by(
    host_subject_id,
    cohort,
    Delivery) %>%
  summarise(
    Mean_Shannon =
      mean(
        Shannon,
        na.rm = TRUE),
    Mean_Simpson =
      mean(
        Simpson,
        na.rm = TRUE),
    Mean_GestWeek =
      mean(
        GestWeek,
        na.rm = TRUE),
    Mean_LogDepth =
      mean(
        Log_Sequencing_Depth,
        na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    Mean_Shannon_log1p =
      log1p(
        Mean_Shannon))

fit_subject_model <- function(
  dat,
  index,
  analysis) {

  response <- if (
    index == "Shannon") {
    "Mean_Shannon_log1p"
  } else {
    "Mean_Simpson"
  }

  d <- dat

  if (analysis == "Global") {

    rhs <- paste0(
      "Delivery * cohort + Mean_GestWeek + ",
      "Mean_LogDepth")

  } else {

    d <- d %>%
      filter(
        cohort == analysis) %>%
      droplevels()

    rhs <- paste0(
      "Delivery + Mean_GestWeek + ",
      "Mean_LogDepth")
  }

  lm(
    as.formula(
      paste(
        response,
        "~",
        rhs)),
    data = d)
}

extract_subject_effect <- function(
  model,
  index,
  analysis) {

  tt <- summary(
    model)$coefficients

  term <- if (
    analysis == "Global") {
    "DeliveryPreterm:cohortUAB"
  } else {
    "DeliveryPreterm"
  }

  if (!term %in%
    rownames(tt))
    stop(
      paste(
        "No se encontró",
        term,
        "en sensibilidad por participante."))

  ci <- confint(
    model,
    parm = term,
    level = .95)

  data.frame(
    Analysis = analysis,
    Index = index,
    Adjustment =
      "Participant_mean",
    Contrast = ifelse(
      analysis == "Global",
      "Delivery_x_Cohort",
      "Preterm_vs_Term"),
    Beta =
      tt[
        term,
        "Estimate"],
    SE =
      tt[
        term,
        "Std. Error"],
    DF =
      df.residual(model),
    T_value =
      tt[
        term,
        "t value"],
    Pvalue =
      tt[
        term,
        "Pr(>|t|)"],
    CI_low = ci[1],
    CI_high = ci[2],
    Direction = if (
      analysis == "Global") {
      NA_character_
    } else if (
      tt[
        term,
        "Estimate"] > 0) {
      "Higher_in_Preterm"
    } else {
      "Higher_in_Term"
    },
    N_samples =
      nobs(model),
    N_subjects =
      nobs(model),
    stringsAsFactors = FALSE)
}

participant_results <- list()

k <- 1

for (idx in indices) {

  for (a in analyses) {

    model <- fit_subject_model(
      participant_mean,
      idx,
      a)

    participant_results[[k]] <-
      extract_subject_effect(
        model,
        idx,
        a)

    k <- k + 1
  }
}

ParticipantMean_Results <- bind_rows(
  participant_results) %>%
  mutate(
    Significant =
      Pvalue < 0.05)

cat(
  "\n=== SENSIBILIDAD: PROMEDIO POR PARTICIPANTE ===\n")

print(
  ParticipantMean_Results)

###############################################################################
# 11. Comparación de sensibilidades
###############################################################################

Sensitivity_Comparison <- Primary_Results %>%
  select(
    Analysis,
    Index,
    Contrast,
    Beta_Primary = Beta,
    P_Primary = Pvalue) %>%
  left_join(
    AgeAdjusted_Results %>%
      select(
        Analysis,
        Index,
        Contrast,
        Beta_Age = Beta,
        P_Age = Pvalue),
    by = c(
      "Analysis",
      "Index",
      "Contrast")) %>%
  left_join(
    ParticipantMean_Results %>%
      select(
        Analysis,
        Index,
        Contrast,
        Beta_ParticipantMean =
          Beta,
        P_ParticipantMean =
          Pvalue),
    by = c(
      "Analysis",
      "Index",
      "Contrast"))

###############################################################################
# 12. Figuras descriptivas
###############################################################################

for (coh in c(
  "Stanford","UAB")) {

  for (idx in c(
    "Shannon","Gini-Simpson")) {

    dat_p <- alpha_long %>%
      filter(
        cohort == coh,
        Index == idx)

    p <- ggplot(
      dat_p,
      aes(
        x = Delivery,
        y = Value,
        fill = Delivery)) +
      geom_violin(
        trim = FALSE,
        alpha = .35) +
      geom_boxplot(
        width = .15,
        outlier.shape = NA) +
      labs(
        title =
          paste(
            coh,
            "-",
            idx),
        subtitle = paste0(
          "Ventana gestacional: ",
          GW_MIN,"-",GW_MAX,
          " semanas"),
        x =
          "Desenlace obstétrico",
        y = ifelse(
          idx ==
            "Gini-Simpson",
          "Gini-Simpson (1-D)",
          "Shannon")) +
      theme_bw() +
      theme(
        legend.position =
          "none")

    ggsave(
      file.path(
        images_dir,
        paste0(
          coh,
          "_",
          gsub(
            "-",
            "",
            idx),
          "_ViolinBoxplot.pdf")),
      p,
      width = 6,
      height = 5)
  }
}

###############################################################################
# 13. Guardar resultados
###############################################################################

write_csv(
  Primary_Results,
  file.path(
    results_dir,
    "AlphaDiversity_Final_PrimaryResults.csv"))

write_csv(
  AgeAdjusted_Results,
  file.path(
    results_dir,
    "AlphaDiversity_Final_AgeAdjustedResults.csv"))

write_csv(
  ParticipantMean_Results,
  file.path(
    results_dir,
    "AlphaDiversity_Final_ParticipantMeanSensitivity.csv"))

write_csv(
  Sensitivity_Comparison,
  file.path(
    results_dir,
    "AlphaDiversity_Final_SensitivityComparison.csv"))

write_csv(
  participant_mean,
  file.path(
    results_dir,
    "AlphaDiversity_Final_ParticipantMeans.csv"))

write_csv(
  Primary_Diagnostics,
  file.path(
    results_dir,
    "AlphaDiversity_Final_PrimaryDiagnostics.csv"))

write_csv(
  AgeAdjusted_Diagnostics,
  file.path(
    results_dir,
    "AlphaDiversity_Final_AgeAdjustedDiagnostics.csv"))

###############################################################################
# 14. Reproducibilidad
###############################################################################

software_versions <- data.frame(
  Package = c(
    "R",
    "phyloseq",
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
        "nlme"))),
  stringsAsFactors = FALSE)

write_csv(
  software_versions,
  file.path(
    results_dir,
    "AlphaDiversity_SoftwareVersions.csv"))

cat("\n=== VERSIONES ===\n")
print(software_versions)

cat("\n=== FIN ===\n")
print(sessionInfo())

###############################################################################
# FIN
###############################################################################
