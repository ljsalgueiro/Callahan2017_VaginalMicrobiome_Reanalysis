###############################################################################
#                                    TFM
###############################################################################
#              Script 10: Diversidad beta, PERMANOVA y potencia
###############################################################################

# Objetivo:
# Evaluar diferencias en la composición de la microbiota vaginal entre parto
# a término y pretérmino dentro de la ventana gestacional 15-33 semanas.
#
# Análisis principal:
# - unidad de análisis: participante;
# - perfil por participante: media de abundancias relativas de sus muestras;
# - distancia de Bray-Curtis;
# - global aditivo: cohort + Delivery, efectos marginales;
# - global con interacción: cohort * Delivery, reportando la interacción;
# - Stanford/UAB: Delivery por separado;
# - sensibilidad por edad materna;
# - BETADISPER y Tukey cuando procede.
#
# Sensibilidad longitudinal:
# - todas las muestras de 15-33 semanas;
# - Delivery se permuta a nivel de participante y se expande a sus muestras;
# - se evita pseudorreplicación por tamaños de cluster desiguales.
#
# Potencia/MDE:
# - simulación a nivel de participante;
# - Global, Stanford y UAB;
# - MDE expresado en R2;
# - comprobación obligatoria de tasa de rechazo bajo H0 (~5 %).

###############################################################################

rm(list = ls())

library(phyloseq)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)
library(ape)
library(tibble)

###############################################################################
# 1. Directorios y parámetros
###############################################################################

results_dir <- "VaginalMicrobiome/results/Beta_diversity/FINAL"
images_dir <- "VaginalMicrobiome/images/Beta_diversity/FINAL"

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
  power_dir,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  power_images,
  recursive = TRUE,
  showWarnings = FALSE)

SEED <- 20260911
NPERM <- 9999
GW_MIN <- 15
GW_MAX <- 33

N_SIM <- 300
NPERM_POWER <- 499
POWER_TARGET <- 0.80
MAX_K <- 128

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
  anyNA(meta$cohort) ||
    anyNA(meta$Delivery))
  stop(
    "Existen muestras sin cohorte o desenlace.")

###############################################################################
# 3. Ventana gestacional y abundancias relativas
###############################################################################

meta_win <- meta %>%
  filter(
    !Temporal_inconsistency,
    !is.na(GestWeek),
    GestWeek >= GW_MIN,
    GestWeek <= GW_MAX)

cat("\n=== VENTANA GESTACIONAL ===\n")
cat("Semanas:",GW_MIN,"-",GW_MAX,"\n")
cat("Muestras:",nrow(meta_win),"\n")
cat(
  "Mujeres:",
  n_distinct(
    meta_win$host_subject_id),
  "\n")

print(
  meta_win %>%
    distinct(
      host_subject_id,
      cohort,
      Delivery) %>%
    count(
      cohort,
      Delivery))

ps_win <- prune_samples(
  meta_win$SampleID,
  ps)

ps_win <- prune_taxa(
  taxa_sums(ps_win) > 0,
  ps_win)

ps_rel <- transform_sample_counts(
  ps_win,
  function(x)
    x / sum(x))

otu_rel <- as(
  otu_table(ps_rel),
  "matrix")

if (taxa_are_rows(ps_rel))
  otu_rel <- t(otu_rel)

otu_rel <- otu_rel[
  meta_win$SampleID,
  ,
  drop = FALSE]

###############################################################################
# 4. Perfil medio por participante
###############################################################################

subj <- as.character(
  meta_win$host_subject_id)

sums_subj <- rowsum(
  otu_rel,
  group = subj)

n_subj <- as.vector(
  table(subj)[
    rownames(sums_subj)])

prof_subj <- sums_subj /
  n_subj

meta_subj <- meta_win %>%
  group_by(
    host_subject_id) %>%
  summarise(
    cohort =
      first(cohort),
    Delivery =
      first(Delivery),
    AGE =
      first(AGE),
    N_samples = n(),
    GestWeek_mean =
      mean(GestWeek),
    .groups = "drop") %>%
  as.data.frame()

rownames(meta_subj) <-
  as.character(
    meta_subj$host_subject_id)

meta_subj <- meta_subj[
  rownames(prof_subj),
  ,
  drop = FALSE]

stopifnot(
  identical(
    rownames(prof_subj),
    rownames(meta_subj)))

cat("\n=== NIVEL PARTICIPANTE ===\n")
print(
  table(
    meta_subj$cohort,
    meta_subj$Delivery))

cat(
  "Muestras por mujer: mediana =",
  median(
    meta_subj$N_samples),
  "| rango =",
  min(
    meta_subj$N_samples),
  "-",
  max(
    meta_subj$N_samples),
  "\n")

D_subj <- vegdist(
  prof_subj,
  method = "bray")

###############################################################################
# 5. PERMANOVA principal
###############################################################################

set.seed(SEED)

pm_global <- adonis2(
  D_subj ~ cohort + Delivery,
  data = meta_subj,
  permutations = NPERM,
  by = "margin")

cat(
  "\n--- PERMANOVA global: efectos marginales ---\n")

print(pm_global)

set.seed(SEED)

pm_inter <- adonis2(
  D_subj ~ cohort * Delivery,
  data = meta_subj,
  permutations = NPERM,
  by = "terms")

cat(
  "\n--- PERMANOVA global: interacción ---\n")

print(pm_inter)

pm_strat <- list()

for (co in levels(
  meta_subj$cohort)) {

  idx <- which(
    meta_subj$cohort == co)

  Dc <- as.dist(
    as.matrix(D_subj)[
      idx,
      idx])

  mc <- droplevels(
    meta_subj[
      idx,
      ,
      drop = FALSE])

  set.seed(SEED)

  pm_strat[[co]] <- adonis2(
    Dc ~ Delivery,
    data = mc,
    permutations = NPERM,
    by = "margin")

  cat(
    "\n--- PERMANOVA",
    co,
    "---\n")

  print(
    pm_strat[[co]])
}

###############################################################################
# 6. Sensibilidad por edad materna
###############################################################################

ok_age <- !is.na(
  meta_subj$AGE)

cat(
  "\nMujeres con edad materna disponible:",
  sum(ok_age),
  "de",
  nrow(meta_subj),
  "\n")

D_age <- as.dist(
  as.matrix(D_subj)[
    ok_age,
    ok_age])

meta_age <- droplevels(
  meta_subj[
    ok_age,
    ,
    drop = FALSE])

set.seed(SEED)

pm_age <- adonis2(
  D_age ~ cohort + AGE + Delivery,
  data = meta_age,
  permutations = NPERM,
  by = "margin")

cat(
  "\n--- PERMANOVA global ajustada por edad ---\n")

print(pm_age)

pm_age_strat <- list()

for (co in levels(
  meta_subj$cohort)) {

  idx <- which(
    meta_subj$cohort == co &
      !is.na(
        meta_subj$AGE))

  Dc <- as.dist(
    as.matrix(D_subj)[
      idx,
      idx])

  mc <- droplevels(
    meta_subj[
      idx,
      ,
      drop = FALSE])

  set.seed(SEED)

  pm_age_strat[[co]] <- adonis2(
    Dc ~ AGE + Delivery,
    data = mc,
    permutations = NPERM,
    by = "margin")

  cat(
    "\n--- PERMANOVA",
    co,
    "ajustada por edad ---\n")

  print(
    pm_age_strat[[co]])
}

###############################################################################
# 7. BETADISPER y Tukey
###############################################################################

bd_res <- list()

for (fac in c(
  "Delivery","cohort")) {

  bd <- betadisper(
    D_subj,
    meta_subj[[fac]])

  set.seed(SEED)

  bd_res[[fac]] <- permutest(
    bd,
    permutations = NPERM)

  cat(
    "\n--- BETADISPER por",
    fac,
    "---\n")

  print(
    bd_res[[fac]])
}

grp4 <- droplevels(
  interaction(
    meta_subj$cohort,
    meta_subj$Delivery,
    sep = "-"))

bd4 <- betadisper(
  D_subj,
  grp4)

set.seed(SEED)

bd_res[["Grupos4"]] <- permutest(
  bd4,
  permutations = NPERM)

cat(
  "\n--- BETADISPER por los cuatro grupos ---\n")

print(
  bd_res[["Grupos4"]])

tukey4 <- NULL

if (
  bd_res[["Grupos4"]]$tab[
    1,
    "Pr(>F)"] <
    0.05) {

  tukey4 <- TukeyHSD(
    bd4)

  tukey_tab <- as.data.frame(
    tukey4$group) %>%
    rownames_to_column(
      "Comparacion")

  write_csv(
    tukey_tab,
    file.path(
      results_dir,
      "Beta_BETADISPER_Tukey_4grupos.csv"))
}

# BETADISPER estratificado por cohorte: Term vs Preterm.
bd_strat <- list()

for (co in levels(
  meta_subj$cohort)) {

  idx <- which(
    meta_subj$cohort == co)

  Dc <- as.dist(
    as.matrix(D_subj)[
      idx,
      idx])

  mc <- droplevels(
    meta_subj[
      idx,
      ,
      drop = FALSE])

  bd <- betadisper(
    Dc,
    mc$Delivery)

  set.seed(SEED)

  bd_strat[[co]] <- permutest(
    bd,
    permutations = NPERM)

  cat(
    "\n--- BETADISPER",
    co,
    "por Delivery ---\n")

  print(
    bd_strat[[co]])
}

###############################################################################
# 8. Sensibilidad longitudinal con permutación por participante
###############################################################################

gower_center <- function(D) {

  A <- -0.5 *
    as.matrix(D)^2

  n <- nrow(A)

  Cn <- diag(n) -
    1 / n

  Cn %*% A %*% Cn
}

tr_GH <- function(
  G,
  X) {

  q <- qr(X)

  Q <- qr.Q(q)[
    ,
    seq_len(
      q$rank),
    drop = FALSE]

  sum(
    Q *
      (G %*% Q))
}

pseudo_F <- function(
  G,
  X_full,
  X_red) {

  n <- nrow(G)

  s_f <- tr_GH(
    G,
    X_full)

  r_f <- qr(
    X_full)$rank

  s_r <- tr_GH(
    G,
    X_red)

  r_r <- qr(
    X_red)$rank

  ss_t <- s_f -
    s_r

  df_t <- r_f -
    r_r

  ss_e <- sum(
    diag(G)) -
    s_f

  df_e <- n -
    r_f

  c(
    F =
      (ss_t / df_t) /
      (ss_e / df_e),
    R2 =
      ss_t /
      sum(
        diag(G)))
}

permanova_cluster <- function(
  D,
  meta,
  cluster,
  label,
  covars = NULL,
  strata = NULL,
  nperm = 9999) {

  G <- gower_center(D)

  clus <- as.character(
    meta[[cluster]])

  rhs_r <- if (
    is.null(covars)) {
    "1"
  } else {
    paste(
      covars,
      collapse = " + ")
  }

  rhs_f <- paste(
    c(
      if (!is.null(
        covars)) {
        covars
      },
      label),
    collapse = " + ")

  X_red <- model.matrix(
    as.formula(
      paste(
        "~",
        rhs_r)),
    meta)

  X_full <- model.matrix(
    as.formula(
      paste(
        "~",
        rhs_f)),
    meta)

  obs <- pseudo_F(
    G,
    X_full,
    X_red)

  key_vars <- c(
    cluster,
    label,
    strata)

  key_vars <- key_vars[
    !is.na(key_vars)]

  key <- unique(
    meta[
      ,
      key_vars,
      drop = FALSE])

  key[[cluster]] <- as.character(
    key[[cluster]])

  stopifnot(
    !anyDuplicated(
      key[[cluster]]))

  Fp <- numeric(
    nperm)

  for (b in seq_len(
    nperm)) {

    kp <- key

    if (is.null(strata)) {

      kp[[label]] <- sample(
        kp[[label]])

    } else {

      kp[[label]] <- unsplit(
        lapply(
          split(
            kp[[label]],
            kp[[strata]]),
          sample),
        kp[[strata]])
    }

    mp <- meta

    mp[[label]] <- kp[[label]][
      match(
        clus,
        kp[[cluster]])]

    Xf <- model.matrix(
      as.formula(
        paste(
          "~",
          rhs_f)),
      mp)

    Fp[b] <- pseudo_F(
      G,
      Xf,
      X_red)["F"]
  }

  list(
    F =
      unname(
        obs["F"]),
    R2 =
      unname(
        obs["R2"]),
    p =
      (
        sum(
          Fp >=
            obs["F"]) +
          1
      ) /
      (
        nperm + 1
      ),
    nperm =
      nperm)
}

D_all <- vegdist(
  otu_rel,
  method = "bray")

cat(
  "\n=== SENSIBILIDAD LONGITUDINAL: ",
  "PERMUTACIÓN POR MUJER ===\n",
  sep = "")

set.seed(SEED)

sens_global <- permanova_cluster(
  D_all,
  meta_win,
  cluster =
    "host_subject_id",
  label = "Delivery",
  covars = "cohort",
  strata = "cohort",
  nperm = NPERM)

sens_strat <- list()

for (co in levels(
  meta_win$cohort)) {

  idx <- which(
    meta_win$cohort == co)

  Dc <- as.dist(
    as.matrix(D_all)[
      idx,
      idx])

  mc <- droplevels(
    meta_win[
      idx,
      ,
      drop = FALSE])

  set.seed(SEED)

  sens_strat[[co]] <- permanova_cluster(
    Dc,
    mc,
    cluster =
      "host_subject_id",
    label = "Delivery",
    nperm = NPERM)
}

cat(
  sprintf(
    "Global: F = %.3f | R2 = %.4f | p = %.4f\n",
    sens_global$F,
    sens_global$R2,
    sens_global$p))

for (co in names(
  sens_strat)) {

  cat(
    sprintf(
      "%-9s: F = %.3f | R2 = %.4f | p = %.4f\n",
      co,
      sens_strat[[co]]$F,
      sens_strat[[co]]$R2,
      sens_strat[[co]]$p))
}

###############################################################################
# 9. Potencia y MDE de PERMANOVA
###############################################################################

clr <- function(
  P,
  eps = NULL) {

  if (is.null(eps))
    eps <- min(
      P[
        P > 0]) /
      2

  L <- log(
    P + eps)

  L -
    rowMeans(L)
}

iclr <- function(L) {

  E <- exp(
    L -
      apply(
        L,
        1,
        max))

  E /
    rowSums(E)
}

adonis_row <- function(
  res,
  term) {

  i <- match(
    term,
    rownames(res))

  if (is.na(i))
    stop(
      paste(
        "No se encontró el término",
        term))

  c(
    p =
      res[["Pr(>F)"]][i],
    R2 =
      res[["R2"]][i])
}

one_power_rep <- function(
  k,
  P,
  md,
  use_cohort,
  direction,
  n_pre) {

  lab <- character(
    nrow(md))

  if (use_cohort) {

    for (co in levels(
      droplevels(
        md$cohort))) {

      idx <- which(
        md$cohort == co)

      lab[idx] <- sample(
        rep(
          c(
            "Preterm",
            "Term"),
          c(
            n_pre[[co]],
            length(idx) -
              n_pre[[co]])))
    }

  } else {

    lab <- sample(
      rep(
        c(
          "Preterm",
          "Term"),
        c(
          n_pre,
          nrow(md) -
            n_pre)))
  }

  lab <- factor(
    lab,
    levels = c(
      "Term","Preterm"))

  L <- clr(P)

  if (k > 0) {

    idx_pre <- lab ==
      "Preterm"

    L[
      idx_pre,
      ] <- L[
        idx_pre,
        ,
        drop = FALSE] +
      matrix(
        k * direction,
        nrow =
          sum(idx_pre),
        ncol =
          ncol(L),
        byrow = TRUE)
  }

  P_sim <- iclr(L)

  D <- vegdist(
    P_sim,
    method = "bray")

  d <- md
  d$lab <- lab


  if (use_cohort) {

    res <- adonis2(
      D ~ cohort + lab,
      data = d,
      permutations =
        NPERM_POWER,
      by = "terms")

  } else {

    res <- adonis2(
      D ~ lab,
      data = d,
      permutations =
        NPERM_POWER,
      by = "terms")
  }

  adonis_row(
    res,
    "lab")
}

evaluate_k <- function(
  k,
  P,
  md,
  use_cohort,
  direction,
  n_pre,
  analysis) {

  r <- t(
    replicate(
      N_SIM,
      one_power_rep(
        k,
        P,
        md,
        use_cohort,
        direction,
        n_pre)))

  out <- data.frame(
    Analisis = analysis,
    k = k,
    Potencia =
      mean(
        r[,"p"] <
          0.05),
    R2_medio =
      mean(
        r[,"R2"]),
    R2_median =
      median(
        r[,"R2"]),
    stringsAsFactors = FALSE)

  cat(
    sprintf(
      "%-9s | k = %6.2f | potencia = %.3f | R2 = %.4f\n",
      analysis,
      k,
      out$Potencia,
      out$R2_medio))

  out
}

power_curve <- function(
  P,
  md,
  analysis,
  use_cohort) {

  md <- droplevels(md)

  L <- clr(P)

  direction <- colMeans(
    L[
      md$Delivery ==
        "Preterm",
      ,
      drop = FALSE]) -
    colMeans(
      L[
        md$Delivery ==
          "Term",
        ,
        drop = FALSE])

  norm_dir <- sqrt(
    sum(
      direction^2))

  if (
    !is.finite(
      norm_dir) ||
      norm_dir == 0)
    stop(
      paste(
        "No se pudo definir la dirección para",
        analysis))

  direction <- direction /
    norm_dir

  if (use_cohort) {

    n_pre <- tapply(
      md$Delivery,
      md$cohort,
      function(x)
        sum(
          x == "Preterm"))

  } else {

    n_pre <- sum(
      md$Delivery ==
        "Preterm")
  }

  cat(
    "\n--- Curva de potencia:",
    analysis,
    "---\n")

  k_grid <- c(
    0,
    0.5,
    1,
    2,
    4)

  results <- lapply(
    k_grid,
    evaluate_k,
    P = P,
    md = md,
    use_cohort =
      use_cohort,
    direction =
      direction,
    n_pre = n_pre,
    analysis =
      analysis)

  out <- bind_rows(
    results)

  while (
    max(
      out$Potencia) <
      POWER_TARGET &&
      max(out$k) <
        MAX_K) {

    new_k <- max(
      out$k) *
      2

    new_result <- evaluate_k(
      new_k,
      P,
      md,
      use_cohort,
      direction,
      n_pre,
      analysis)

    out <- bind_rows(
      out,
      new_result)
  }

  out %>%
    arrange(k)
}

set.seed(SEED)

pow_global <- power_curve(
  prof_subj,
  meta_subj,
  "Global",
  TRUE)

idx_s <- which(
  meta_subj$cohort ==
    "Stanford")

set.seed(SEED)

pow_stanford <- power_curve(
  prof_subj[
    idx_s,
    ,
    drop = FALSE],
  meta_subj[
    idx_s,
    ,
    drop = FALSE],
  "Stanford",
  FALSE)

idx_u <- which(
  meta_subj$cohort ==
    "UAB")

set.seed(SEED)

pow_uab <- power_curve(
  prof_subj[
    idx_u,
    ,
    drop = FALSE],
  meta_subj[
    idx_u,
    ,
    drop = FALSE],
  "UAB",
  FALSE)

pow_all <- bind_rows(
  pow_global,
  pow_stanford,
  pow_uab)

type1_check <- pow_all %>%
  filter(
    k == 0) %>%
  transmute(
    Analisis,
    Rejection_rate =
      Potencia,
    Expected = 0.05,
    Difference =
      Rejection_rate -
      Expected)

cat("\n=== COMPROBACIÓN k = 0 ===\n")
print(type1_check)

if (any(
  abs(
    type1_check$Difference) >
    0.03))
  warning(
    paste(
      "La tasa de rechazo bajo H0 se aleja",
      "más de 0.03 del 5 % nominal."))

calculate_mde <- function(dat) {

  dat <- dat %>%
    arrange(k)

  i <- which(
    dat$Potencia >=
      POWER_TARGET)[1]

  if (is.na(i))
    return(
      data.frame(
        R2_MDE80 =
          NA_real_,
        k_MDE80 =
          NA_real_,
        Potencia_max =
          max(
            dat$Potencia)))

  if (i == 1)
    return(
      data.frame(
        R2_MDE80 =
          dat$R2_medio[i],
        k_MDE80 =
          dat$k[i],
        Potencia_max =
          max(
            dat$Potencia)))

  idx <- c(
    i - 1,
    i)

  data.frame(
    R2_MDE80 =
      approx(
        x =
          dat$Potencia[idx],
        y =
          dat$R2_medio[idx],
        xout =
          POWER_TARGET)$y,
    k_MDE80 =
      approx(
        x =
          dat$Potencia[idx],
        y =
          dat$k[idx],
        xout =
          POWER_TARGET)$y,
    Potencia_max =
      max(
        dat$Potencia))
}

mde_80 <- pow_all %>%
  group_split(
    Analisis) %>%
  lapply(
    function(x) {

      result <- calculate_mde(
        x)

      result$Analisis <-
        unique(
          x$Analisis)

      result
    }) %>%
  bind_rows() %>%
  select(
    Analisis,
    R2_MDE80,
    k_MDE80,
    Potencia_max)

r2_observed <- data.frame(
  Analisis = c(
    "Global",
    "Stanford",
    "UAB"),
  R2_observado = c(
    pm_global[
      "Delivery",
      "R2"],
    pm_strat[[
      "Stanford"]][
        "Delivery",
        "R2"],
    pm_strat[[
      "UAB"]][
        "Delivery",
        "R2"]))

MDE_Comparison <- left_join(
  r2_observed,
  mde_80,
  by = "Analisis") %>%
  mutate(
    Fraction_of_MDE =
      R2_observado /
      R2_MDE80)

cat(
  "\n=== MDE EN R2: POTENCIA 80 % ===\n")

print(mde_80)

cat(
  "\n=== R2 OBSERVADO VS MDE80 ===\n")

print(
  MDE_Comparison)

write_csv(
  pow_all,
  file.path(
    power_dir,
    "Power_Beta_Curves.csv"))

write_csv(
  type1_check,
  file.path(
    power_dir,
    "Power_Beta_TypeI_Check.csv"))

write_csv(
  mde_80,
  file.path(
    power_dir,
    "Power_Beta_MDE80.csv"))

write_csv(
  MDE_Comparison,
  file.path(
    power_dir,
    "Power_Beta_Observed_vs_MDE.csv"))

p_pow <- ggplot(
  pow_all,
  aes(
    x = R2_medio,
    y = Potencia,
    colour = Analisis)) +
  geom_hline(
    yintercept =
      POWER_TARGET,
    linetype = 2) +
  geom_line(
    linewidth = .7) +
  geom_point(
    size = 2) +
  scale_y_continuous(
    limits = c(
      0,1)) +
  labs(
    title =
      "Potencia de la PERMANOVA",
    subtitle = paste0(
      "Ventana gestacional ",
      GW_MIN,"-",GW_MAX,
      " semanas"),
    x =
      expression(
        R^2~"inducido"),
    y = "Potencia") +
  theme_bw()

ggsave(
  file.path(
    power_images,
    "Power_Beta_Curves.pdf"),
  p_pow,
  width = 7,
  height = 5)

###############################################################################
# 10. Ordenación a nivel de participante
###############################################################################

pcoa_res <- ape::pcoa(
  D_subj,
  correction = "cailliez")

ev <- if (!is.null(
  pcoa_res$values$Rel_corr_eig)) {
  pcoa_res$values$Rel_corr_eig
} else {
  pcoa_res$values$Relative_eig
}

pc_df <- data.frame(
  PCoA1 =
    pcoa_res$vectors[,1],
  PCoA2 =
    pcoa_res$vectors[,2],
  meta_subj)

cat("\n=== PCoA ===\n")
cat(
  sprintf(
    "Ejes 1-2: %.1f %% + %.1f %% = %.1f %%\n",
    100 * ev[1],
    100 * ev[2],
    100 * (
      ev[1] +
        ev[2])))

p_pcoa <- ggplot(
  pc_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = Delivery,
    shape = cohort)) +
  geom_point(
    size = 2.4,
    alpha = .85) +
  stat_ellipse(
    aes(
      group =
        interaction(
          cohort,
          Delivery)),
    type = "norm",
    linewidth = .4,
    show.legend = FALSE) +
  labs(
    title =
      "PCoA - Bray-Curtis",
    subtitle = paste0(
      "Un punto = una mujer; perfil medio, semanas ",
      GW_MIN,"-",GW_MAX),
    x = sprintf(
      "PCoA1 (%.1f %%)",
      100 * ev[1]),
    y = sprintf(
      "PCoA2 (%.1f %%)",
      100 * ev[2])) +
  theme_bw()

ggsave(
  file.path(
    images_dir,
    "PCoA_participant_level.pdf"),
  p_pcoa,
  width = 7,
  height = 5.5)

set.seed(SEED)

nmds <- metaMDS(
  D_subj,
  k = 2,
  trymax = 100,
  trace = 0)

cat(
  "NMDS stress:",
  round(
    nmds$stress,
    4),
  "\n")

###############################################################################
# 11. Tablas y objetos finales
###############################################################################

tidy_adonis <- function(
  x,
  label,
  terms = NULL) {

  out <- data.frame(
    Analisis = label,
    Termino = rownames(x),
    Df = x$Df,
    SumOfSqs =
      x$SumOfSqs,
    R2 = x$R2,
    F = x$F,
    Pvalue =
      x$`Pr(>F)`,
    stringsAsFactors = FALSE) %>%
    filter(
      !Termino %in%
        c(
          "Residual",
          "Total"))

  if (!is.null(terms))
    out <- out %>%
      filter(
        Termino %in%
          terms)

  out
}

permanova_tab <- bind_rows(
  tidy_adonis(
    pm_global,
    "Global (marginal)",
    c(
      "cohort",
      "Delivery")),
  tidy_adonis(
    pm_inter,
    "Global (interaccion)",
    "cohort:Delivery"),
  tidy_adonis(
    pm_strat[["Stanford"]],
    "Stanford",
    "Delivery"),
  tidy_adonis(
    pm_strat[["UAB"]],
    "UAB",
    "Delivery"),
  tidy_adonis(
    pm_age,
    "Global ajustado por edad",
    c(
      "cohort",
      "AGE",
      "Delivery")),
  tidy_adonis(
    pm_age_strat[["Stanford"]],
    "Stanford ajustado por edad",
    c(
      "AGE",
      "Delivery")),
  tidy_adonis(
    pm_age_strat[["UAB"]],
    "UAB ajustado por edad",
    c(
      "AGE",
      "Delivery")))

sens_tab <- bind_rows(
  data.frame(
    Analisis =
      "Global (Delivery | cohorte)",
    F = sens_global$F,
    R2 = sens_global$R2,
    Pvalue =
      sens_global$p),
  data.frame(
    Analisis = "Stanford",
    F =
      sens_strat[["Stanford"]]$F,
    R2 =
      sens_strat[["Stanford"]]$R2,
    Pvalue =
      sens_strat[["Stanford"]]$p),
  data.frame(
    Analisis = "UAB",
    F =
      sens_strat[["UAB"]]$F,
    R2 =
      sens_strat[["UAB"]]$R2,
    Pvalue =
      sens_strat[["UAB"]]$p))

tidy_betadisper <- function(
  x,
  label) {

  tab <- as.data.frame(
    x$tab)

  data.frame(
    Analisis = label,
    Df = tab[1,"Df"],
    SumSq = tab[1,"Sum Sq"],
    F = tab[1,"F"],
    Pvalue = tab[1,"Pr(>F)"],
    stringsAsFactors = FALSE)
}

betadisper_tab <- bind_rows(
  tidy_betadisper(
    bd_res[["Delivery"]],
    "Global: Delivery"),
  tidy_betadisper(
    bd_res[["cohort"]],
    "Global: cohort"),
  tidy_betadisper(
    bd_res[["Grupos4"]],
    "Global: 4 grupos"),
  tidy_betadisper(
    bd_strat[["Stanford"]],
    "Stanford: Delivery"),
  tidy_betadisper(
    bd_strat[["UAB"]],
    "UAB: Delivery"))

write_csv(
  permanova_tab,
  file.path(
    results_dir,
    "Beta_PERMANOVA_participant.csv"))

write_csv(
  sens_tab,
  file.path(
    results_dir,
    "Beta_PERMANOVA_clusterperm.csv"))

write_csv(
  betadisper_tab,
  file.path(
    results_dir,
    "Beta_BETADISPER.csv"))

write_csv(
  pc_df,
  file.path(
    results_dir,
    "Beta_PCoA_participant_scores.csv"))

saveRDS(
  list(
    D_subj = D_subj,
    meta_subj = meta_subj,
    prof_subj = prof_subj,
    pcoa = pcoa_res,
    nmds = nmds,
    power = pow_all,
    mde80 = mde_80),
  file.path(
    results_dir,
    "Beta_objects.rds"))

###############################################################################
# 12. Reproducibilidad
###############################################################################

versions <- data.frame(
  Package = c(
    "R",
    "phyloseq",
    "vegan",
    "ape"),
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
        "ape"))),
  stringsAsFactors = FALSE)

write_csv(
  versions,
  file.path(
    results_dir,
    "Beta_SoftwareVersions.csv"))

cat("\n=== VERSIONES ===\n")
print(versions)

cat("\n=== FIN ===\n")
print(sessionInfo())

###############################################################################
# FIN
###############################################################################
