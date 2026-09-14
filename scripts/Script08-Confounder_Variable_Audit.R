###############################################################################
#                                    TFM
###############################################################################
#         Script 08: Auditoría de potenciales variables de confusión
###############################################################################

# Objetivo:
#   Evaluar qué variables disponibles en sample_data (ps_prevotella_refined.rds)
#   podrían utilizarse como posibles confusores

# La auditoría evalúa:
#   1. Variables disponibles.
#   2. Tipo de variable.
#   3. Valores ausentes.
#   4. Número de valores/niveles únicos.
#   5. Constancia de cada variable dentro de cada mujer.
#   6. Distribución por cohorte.
#   7. Distribución por Term/Preterm.
#   8. Asociación preliminar con el desenlace.
#   9. Posible redundancia con cohorte y otras variables.

# IMPORTANTE:
#   Este script NO selecciona automáticamente los confusores definitivos.
#   Genera una auditoría para decidir posteriormente qué variables tienen
#   justificación clínica y estadística para entrar en el modelo.

# Entrada:
#   - ps_prevotella_refined.rds

# Salidas:
#   - Confounder_Variable_Audit.csv
#   - Confounder_SubjectConsistency.csv
#   - Confounder_Numeric_Associations.csv
#   - Confounder_Categorical_Associations.csv
#   - Confounder_Correlation_Matrix.csv
#   - Confounder_Audit.xlsx

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
# 1. Librerías
###############################################################################

library(phyloseq)
library(dplyr)
library(tidyr)
library(readr)

###############################################################################
# 2. Directorios
###############################################################################

results_dir <- file.path(project_dir, "results", "Confounder_audit")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE)

###############################################################################
# 3. Cargar objeto phyloseq definitivo
###############################################################################

ps <- readRDS(file.path(project_dir, "cleaned", "ps_prevotella_refined.rds"))

meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

n_metadata_original <- ncol(meta)
meta$SampleID <- rownames(meta)

# Corregir tipos de variables clínicas:
meta$AGE <- suppressWarnings(as.numeric(as.character(meta$AGE)))

# Conservar el valor original y extraer de forma robusta la semana gestacional.
# parse_number() recupera el componente numérico aunque el campo contenga texto.
if ("gest_week_collection" %in% colnames(meta)) {
  meta$gest_week_collection_raw <- as.character(meta$gest_week_collection)
  meta$gest_week_collection <- suppressWarnings(
    readr::parse_number(meta$gest_week_collection_raw))}

if ("gest_wk_delivery" %in% colnames(meta)) {
  meta$gest_wk_delivery_raw <- as.character(meta$gest_wk_delivery)
  meta$gest_wk_delivery_num <- suppressWarnings(
    readr::parse_number(meta$gest_wk_delivery_raw))}

cat("\nNúmero de muestras:", nrow(meta), "\n")
# Número de muestras: 2179 

cat("Número de variables:", ncol(meta), "\n")
# Número de variables: 58

###############################################################################
# 4. Conservar las muestras clínicas utilizadas en los análisis principales
###############################################################################

if (!"Groups" %in% colnames(meta)) {
  stop("ERROR: no se encontró la variable Groups.")}

meta <- meta %>% filter(Groups %in% c("ST", "SP", "UT", "UP"))

###############################################################################
# 5. Crear variables clínicas estandarizadas
###############################################################################

meta <- meta %>%
  mutate(
    Cohort_analysis = case_when(
      Groups %in% c("ST", "SP") ~ "Stanford",
      Groups %in% c("UT", "UP") ~ "UAB",
      TRUE ~ NA_character_),
    Delivery_analysis = case_when(
      Groups %in% c("ST", "UT") ~ "Term",
      Groups %in% c("SP", "UP") ~ "Preterm",
      TRUE ~ NA_character_))

###############################################################################
# 6. Identificar variable de sujeto
###############################################################################

if (!"host_subject_id" %in% colnames(meta)) {
  stop("ERROR: no se encontró host_subject_id.")}

###############################################################################
# 7. Variables que NO deben considerarse confusores
###############################################################################

# Se excluyen identificadores, variables técnicas obvias, outcome y variables
# derivadas directamente del outcome/cohorte.
# gest_wk_delivery y sus versiones derivadas se excluyen porque la semana de
# parto está directamente ligada al desenlace obstétrico y no es un confusor basal.
# gest_week_collection_raw se excluye porque es una copia auxiliar para el parseo;
# gest_week_collection se audita por separado como variable temporal.

exclude_variables <- c(
  "SampleID",
  "host_subject_id",
  "Groups",
  "cohort",
  "Cohort_analysis",
  "term_vs_preterm_delivery",
  "Delivery_analysis",
  "Experiment",
  "Library.Name",
  "Sample.Name",
  "BioSample",
  "BioProject",
  "SRA.Study",
  "Run",
  "IL_Run",
  "NumberInRun",
  "Platform",
  "Instrument",
  "Assay.Type",
  "LibraryLayout",
  "LibrarySelection",
  "LibrarySource",
  "Bases",
  "Bytes",
  "AvgSpotLen",
  "DATASTORE.filetype",
  "DATASTORE.provider",
  "DATASTORE.region",
  "create_date",
  "ReleaseDate",
  "version",
  "indication_for_PTB",
  "gest_day_delivery",
  "gest_wk_delivery",
  "gest_wk_delivery_raw",
  "gest_wk_delivery_num",
  "gest_week_collection_raw")

candidate_variables <- setdiff(
  colnames(meta),
  exclude_variables)

###############################################################################
# 8. Función para clasificar tipo de variable
###############################################################################

classify_variable <- function(x) {
  x_non_na <- x[!is.na(x)]
  if (length(x_non_na) == 0) {
    return("All_NA")}
  if (is.numeric(x) || is.integer(x)) {
    return("Numeric")}
  n_unique <- length(unique(x_non_na))
  if (n_unique == 2) {return("Binary")}
  if (n_unique <= 10) {return("Categorical")}
  return("Character_high_cardinality")}

###############################################################################
# 9. Auditoría general por variable
###############################################################################

variable_audit <- lapply(
  candidate_variables,
  function(v) {
    x <- meta[[v]]
    non_na <- !is.na(x)
    data.frame(
      Variable = v,
      Type = classify_variable(x),
      N_samples = length(x),
      N_non_missing = sum(non_na),
      N_missing = sum(!non_na),
      Missing_percent = 100 * mean(!non_na),
      N_unique = length(unique(x[non_na])), stringsAsFactors = FALSE)})

variable_audit <- bind_rows(variable_audit)

###############################################################################
# 10. Evaluar constancia dentro de cada mujer
###############################################################################

# Un confusor basal debería ser constante entre las muestras longitudinales
# de una misma mujer.

subject_consistency <- lapply(
  candidate_variables,
  function(v) {
    tmp <- meta %>%
      group_by(host_subject_id) %>%
      summarise(
        N_unique = n_distinct(.data[[v]], na.rm = TRUE),
        .groups = "drop")
    data.frame(
      Variable = v,
      N_subjects = nrow(tmp),
      Subjects_with_0_values = sum(tmp$N_unique == 0),
      Subjects_with_1_value = sum(tmp$N_unique == 1),
      Subjects_with_multiple_values = sum(tmp$N_unique > 1),
      Constant_within_subject = all(tmp$N_unique <= 1),
      stringsAsFactors = FALSE)})

subject_consistency <- bind_rows(subject_consistency)

###############################################################################
# 11. Unir auditoría general y consistencia
###############################################################################

variable_audit <- variable_audit %>%
  left_join(subject_consistency, by = "Variable")

###############################################################################
# 12. Crear metadata a nivel de mujer
###############################################################################

# Para evaluar confusores basales no conviene contar cada muestra longitudinal
# como una observación independiente.
# Por eso se crea una tabla con una fila por mujer.

subject_meta <- meta %>%
  group_by(host_subject_id) %>%
  summarise(
    Cohort = first(Cohort_analysis),
    Delivery = first(Delivery_analysis),
    across(all_of(candidate_variables), ~ {
        vals <- unique(na.omit(.x))
        if (length(vals) == 1) {vals[1]
        } else {NA}}),
    .groups = "drop")

###############################################################################
# 13. Variables candidatas constantes dentro de cada sujeto
###############################################################################

candidate_baseline <- variable_audit %>%
  filter(
    Constant_within_subject,
    Type %in% c("Numeric", "Binary", "Categorical"),
    Missing_percent < 50,
    N_unique > 1) %>%
  pull(Variable)

cat("\nVariables basales potencialmente analizables:\n")
print(candidate_baseline)
# AGE, ethnicity, race

###############################################################################
# 14. Asociación preliminar de variables numéricas con Delivery
###############################################################################

# Se utiliza Wilcoxon únicamente como auditoría descriptiva.
# NO se utiliza para decidir automáticamente qué variable debe entrar
# en el modelo.

numeric_candidates <- candidate_baseline[
  sapply(subject_meta[candidate_baseline], is.numeric)]

numeric_results <- list()

for (v in numeric_candidates) {
  dat <- subject_meta %>%
    select(Delivery, value = all_of(v)) %>%
    filter(!is.na(Delivery), !is.na(value))
  if (n_distinct(dat$Delivery) == 2 && nrow(dat) >= 3) {
    wt <- wilcox.test(value ~ Delivery, data = dat, exact = FALSE)
    descriptive <- dat %>%
      group_by(Delivery) %>%
      summarise(
        N = n(),
        Mean = mean(value),
        Median = median(value),
        SD = sd(value),
        IQR = IQR(value),
        .groups = "drop")
    numeric_results[[v]] <- data.frame(
      Variable = v,
      N = nrow(dat),
      Pvalue_Delivery = wt$p.value,
      Term_median = descriptive$Median[descriptive$Delivery == "Term"],
      Preterm_median = descriptive$Median[descriptive$Delivery == "Preterm"],
      stringsAsFactors = FALSE)}}

Numeric_Associations <- bind_rows(numeric_results)

###############################################################################
# 15. Asociación preliminar de variables categóricas con Delivery
###############################################################################

# Se utiliza Fisher cuando las tablas son pequeñas.
# Para tablas mayores se utiliza chi-cuadrado si es posible.

categorical_candidates <- setdiff(candidate_baseline, numeric_candidates)

categorical_results <- list()

for (v in categorical_candidates) {
  dat <- subject_meta %>%
    select(Delivery, value = all_of(v)) %>%
    filter(!is.na(Delivery), !is.na(value))
  tab <- table(dat$value, dat$Delivery)
  if (nrow(tab) >= 2 && ncol(tab) == 2) {
    expected <- suppressWarnings(chisq.test(tab)$expected)
    if (any(expected < 5)) {
      test <- fisher.test(
        tab,
        simulate.p.value = (nrow(tab) > 2),
        B = 10000)
      method <- "Fisher"
    } else {
      test <- chisq.test(tab, correct = FALSE)
      method <- "Chi-square"}
    categorical_results[[v]] <- data.frame(
      Variable = v,
      N = sum(tab),
      Levels = nrow(tab),
      Test = method,
      Pvalue_Delivery = test$p.value,
      stringsAsFactors = FALSE)}}

Categorical_Associations <- bind_rows(categorical_results)

###############################################################################
# 16. Evaluar asociación de cada candidato con cohorte
###############################################################################

# Una variable fuertemente asociada con Stanford/UAB puede reflejar
# diferencias estructurales entre las cohortes y debe revisarse antes
# de incorporarla junto con cohort en un modelo multivariable.

cohort_associations <- list()

for (v in candidate_baseline) {
  x <- subject_meta[[v]]
  if (is.numeric(x)) {
    dat <- subject_meta %>%
      dplyr::select(Cohort, value = all_of(v)) %>%
      dplyr::filter(!is.na(Cohort), !is.na(value))
    if (
      nrow(dat) > 0 &&
      dplyr::n_distinct(dat$Cohort) == 2 &&
      dplyr::n_distinct(dat$value) > 1) {
      wt <- wilcox.test(value ~ Cohort, data = dat, exact = FALSE)
      cohort_associations[[v]] <- data.frame(
        Variable = v,
        Test_Cohort = "Wilcoxon",
        Pvalue_Cohort = wt$p.value,
        stringsAsFactors = FALSE)}
  } else {
    dat <- subject_meta %>%
      dplyr::select(Cohort, value = all_of(v)) %>%
      dplyr::filter(!is.na(Cohort), !is.na(value))
    tab <- table(dat$value, dat$Cohort)
    if (nrow(tab) >= 2 && ncol(tab) == 2) {
      expected <- suppressWarnings(chisq.test(tab)$expected)
      if (any(expected < 5)) {
        test <- fisher.test(tab, simulate.p.value = nrow(tab) > 2, B = 10000)
        method <- "Fisher"
      } else {
        test <- chisq.test(tab, correct = FALSE)
        method <- "Chi-square"}
      cohort_associations[[v]] <- data.frame(
        Variable = v,
        Test_Cohort = method,
        Pvalue_Cohort = test$p.value,
        stringsAsFactors = FALSE)}}}

# Crear siempre una tabla con la estructura esperada
if (length(cohort_associations) > 0) {
  Cohort_Associations <- bind_rows(cohort_associations)
} else {
  Cohort_Associations <- data.frame(
    Variable = character(),
    Test_Cohort = character(),
    Pvalue_Cohort = numeric(),
    stringsAsFactors = FALSE)}

cat("\nAsociaciones evaluables con cohorte:", nrow(Cohort_Associations), "\n")
# Asociaciones evaluables con cohorte: 3

###############################################################################
# 17. Integrar asociaciones con Delivery y cohorte
###############################################################################

# Crear tablas vacías con las columnas necesarias si no hubo resultados
if (
  nrow(Numeric_Associations) == 0 ||
  !"Variable" %in% colnames(Numeric_Associations)) {
  Numeric_Associations <- data.frame(
    Variable = character(),
    Pvalue_Delivery = numeric(),
    stringsAsFactors = FALSE)}

if (
  nrow(Categorical_Associations) == 0 ||
  !"Variable" %in% colnames(Categorical_Associations)) {
  Categorical_Associations <- data.frame(
    Variable = character(),
    Pvalue_Delivery = numeric(),
    stringsAsFactors = FALSE)}

association_summary <- variable_audit %>%
  left_join(Numeric_Associations %>%
      dplyr::select(
        Variable,
        Pvalue_Delivery),
    by = "Variable") %>%
  rename(Pvalue_Delivery_numeric = Pvalue_Delivery) %>%
  left_join(Categorical_Associations %>%
      dplyr::select(Variable, Pvalue_Delivery), by = "Variable") %>%
  rename(Pvalue_Delivery_categorical = Pvalue_Delivery) %>%
  mutate(
    Pvalue_Delivery = coalesce(
      Pvalue_Delivery_numeric,
      Pvalue_Delivery_categorical)) %>%
  select(-Pvalue_Delivery_numeric, -Pvalue_Delivery_categorical) %>%
  left_join(Cohort_Associations, by = "Variable")

cat("\nNúmero de candidate_variables:", length(candidate_variables), "\n")
# Número de candidate_variables: 24

cat("Número de candidate_baseline:", length(candidate_baseline), "\n")
# Número de candidate_baseline: 3 

cat("\nVariables candidate_baseline:\n")
cat(candidate_baseline, sep = "\n")
# AGE, ethnicity, race

cat("\nNúmero de asociaciones con Delivery numéricas:", nrow(Numeric_Associations), "\n")
# Número de asociaciones con Delivery numéricas: 1

cat("Número de asociaciones con Delivery categóricas:", nrow(Categorical_Associations), "\n")
# Número de asociaciones con Delivery categóricas: 2

cat("Número de asociaciones con Cohort:", nrow(Cohort_Associations), "\n")
# Número de asociaciones con Cohort: 3

###############################################################################
# 18. Crear clasificación orientativa
###############################################################################

# Esta clasificación NO determina automáticamente qué variables deben usarse.
# Sirve para priorizar la revisión.

association_summary <- association_summary %>%
  mutate(
    Audit_status = case_when(
      !Constant_within_subject ~ "Varía dentro de la mujer",
      Missing_percent >= 50 ~ "Demasiados valores ausentes",
      N_unique <= 1 ~ "Sin variabilidad",
      Type == "Character_high_cardinality" ~ "Identificador o alta cardinalidad",
      TRUE ~ "Candidata para revisión"))

###############################################################################
# 19. Correlaciones entre variables numéricas
###############################################################################

numeric_baseline <- candidate_baseline[
  sapply(subject_meta[candidate_baseline], is.numeric)]

if (length(numeric_baseline) >= 2) {
  numeric_matrix <- subject_meta %>% select(all_of(numeric_baseline))
  Correlation_Matrix <- cor(
    numeric_matrix,
    use = "pairwise.complete.obs",
    method = "spearman")
  
  write_csv(
    data.frame(
      Variable = rownames(Correlation_Matrix),
      Correlation_Matrix,
      check.names = FALSE),
    file.path(
      results_dir,
      "Confounder_Correlation_Matrix.csv"))
} else {Correlation_Matrix <- NULL}

###############################################################################
# 20. Mostrar variables prioritarias para revisión manual
###############################################################################

review_table <- association_summary %>%
  filter(Audit_status == "Candidata para revisión") %>%
  arrange(Missing_percent, Pvalue_Delivery)

cat("\n========================================\n")
cat("VARIABLES CANDIDATAS PARA REVISIÓN\n")
cat("========================================\n")

print(as.data.frame(review_table))

###############################################################################
# 21. Guardar tablas
###############################################################################

write_csv(
  variable_audit,
  file.path(
    results_dir,
    "Confounder_Variable_Audit.csv"))

write_csv(
  subject_consistency,
  file.path(
    results_dir,
    "Confounder_SubjectConsistency.csv"))

write_csv(
  Numeric_Associations,
  file.path(
    results_dir,
    "Confounder_Numeric_Associations.csv"))

write_csv(
  Categorical_Associations,
  file.path(
    results_dir,
    "Confounder_Categorical_Associations.csv"))

write_csv(
  Cohort_Associations,
  file.path(
    results_dir,
    "Confounder_Cohort_Associations.csv"))

write_csv(
  association_summary,
  file.path(
    results_dir,
    "Confounder_Audit_Summary.csv"))

write_csv(
  review_table,
  file.path(
    results_dir,
    "Confounder_Candidates_For_Review.csv"))

###############################################################################
# 22. Guardar Excel
###############################################################################

excel_sheets <- list(
  Variable_Audit = variable_audit,
  Subject_Consistency = subject_consistency,
  Candidate_Review = review_table,
  Numeric_Delivery = Numeric_Associations,
  Categorical_Delivery = Categorical_Associations,
  Cohort_Associations = Cohort_Associations,
  Subject_Metadata = subject_meta)

if (!is.null(Correlation_Matrix)) {
  excel_sheets$Numeric_Correlations <- data.frame(
    Variable = rownames(Correlation_Matrix),
    Correlation_Matrix,
    check.names = FALSE)}

if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(
    excel_sheets,
    file = file.path(
      results_dir,
      "Confounder_Audit.xlsx"),
    rowNames = FALSE,
    overwrite = TRUE)
} else {
  warning(
    "openxlsx no está instalado; se omite la exportación XLSX.")
}

###############################################################################
# 23. Resumen final
###############################################################################

cat("\n========================================\n")
cat("AUDITORÍA DE CONFUSORES FINALIZADA\n")
cat("========================================\n")

cat("\nVariables originales en sample_data:", n_metadata_original, "\n")
# Variables originales en sample_data: 54

cat("Variables presentes en meta al final de la auditoría:", ncol(meta), "\n")
# Variables presentes en meta al final de la auditoría: 60 

cat("(La diferencia corresponde a variables auxiliares derivadas creadas por el script.)\n")

cat("Variables evaluadas:", length(candidate_variables), "\n")
# Variables evaluadas: 24

cat("Variables basales potencialmente analizables:", length(candidate_baseline), "\n")
# Variables basales potencialmente analizables: 3 

cat("Variables candidatas para revisión manual:", nrow(review_table), "\n")
# Variables candidatas para revisión manual: 3

cat("\nResultados guardados en:\n",
  results_dir,
  "\n")
# Resultados guardados en:
# VaginalMicrobiome/results/Confounder_audit

cat("\n========================================\n")

###############################################################################
# Distribución de potenciales confusores a nivel de mujer
###############################################################################

cat("\n========================================\n")
cat("EDAD POR COHORTE Y OUTCOME\n")
cat("========================================\n")

age_summary <- subject_meta %>%
  group_by(Cohort, Delivery) %>%
  summarise(
    N = sum(!is.na(AGE)),
    Mean = mean(AGE, na.rm = TRUE),
    SD = sd(AGE, na.rm = TRUE),
    Median = median(AGE, na.rm = TRUE),
    IQR = IQR(AGE, na.rm = TRUE),
    Min = min(AGE, na.rm = TRUE),
    Max = max(AGE, na.rm = TRUE),
    .groups = "drop")

print(age_summary)
#   Cohort   Delivery     N  Mean    SD Median   IQR   Min   Max
#    <chr>    <chr>    <int> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>
# 1 Stanford Preterm      9  34.9  3.55     34  5       31    42
# 2 Stanford Term        30  31.8  3.78     31  5       25    41
# 3 UAB      Preterm     40  26.2  4.35     25  7.25    20    34
# 4 UAB      Term        55  26.8  4.89     26  7       17    38

###############################################################################
# Race
###############################################################################

cat("\n========================================\n")
cat("RACE POR COHORTE Y OUTCOME\n")
cat("========================================\n")

race_table <- subject_meta %>%
  count(Cohort, Delivery, race, name = "N") %>%
  group_by(Cohort, Delivery) %>%
  mutate(Percent = 100 * N / sum(N)) %>%
  ungroup()

print(as.data.frame(race_table))
#      Cohort Delivery           race  N   Percent
# 1  Stanford  Preterm AmericanIndian  1 11.111111
# 2  Stanford  Preterm          Asian  1 11.111111
# 3  Stanford  Preterm          Other  4 44.444444
# 4  Stanford  Preterm          White  3 33.333333
# 5  Stanford     Term          Asian  3 10.000000
# 6  Stanford     Term          Black  1  3.333333
# 7  Stanford     Term        Decline  1  3.333333
# 8  Stanford     Term          Other  6 20.000000
# 9  Stanford     Term          White 19 63.333333
# 10      UAB  Preterm          Black 36 87.804878
# 11      UAB  Preterm          Other  2  4.878049
# 12      UAB  Preterm          White  3  7.317073
# 13      UAB     Term                 1  1.818182
# 14      UAB     Term          Asian  2  3.636364
# 15      UAB     Term          Black 43 78.181818
# 16      UAB     Term          Other  3  5.454545
# 17      UAB     Term          White  6 10.909091

###############################################################################
# Ethnicity
###############################################################################

cat("\n========================================\n")
cat("ETHNICITY POR COHORTE Y OUTCOME\n")
cat("========================================\n")

ethnicity_table <- subject_meta %>%
  count(Cohort, Delivery, ethnicity, name = "N") %>%
  group_by(Cohort, Delivery) %>%
  mutate(Percent = 100 * N / sum(N)) %>%
  ungroup()

print(as.data.frame(ethnicity_table))
#      Cohort Delivery      ethnicity  N   Percent
# 1  Stanford  Preterm       Hispanic  5 55.555556
# 2  Stanford  Preterm    NonHispanic  4 44.444444
# 3  Stanford     Term       Hispanic  5 16.666667
# 4  Stanford     Term    NonHispanic 25 83.333333
# 5       UAB  Preterm       Hispanic  1  2.439024
# 6       UAB  Preterm    NonHispanic 39 95.121951
# 7       UAB  Preterm not applicable  1  2.439024
# 8       UAB     Term        Decline  1  1.818182
# 9       UAB     Term       Hispanic  4  7.272727
# 10      UAB     Term    NonHispanic 50 90.909091

###############################################################################
# Race × Cohort
###############################################################################

cat("\n========================================\n")
cat("RACE × COHORT\n")
cat("========================================\n")

print(table(subject_meta$race, subject_meta$Cohort, useNA = "ifany"))
#                  Stanford UAB
#                       0   1
# AmericanIndian        1   0
# Asian                 4   2
# Black                 1  79
# Decline               1   0
# Other                10   5
# White                22   9

###############################################################################
# Ethnicity × Cohort
###############################################################################

cat("\n========================================\n")
cat("ETHNICITY × COHORT\n")
cat("========================================\n")

print(table(subject_meta$ethnicity, subject_meta$Cohort, useNA = "ifany"))
#                  Stanford UAB
# Decline               0   1
# Hispanic             10   5
# NonHispanic          29  89
# not applicable        0   1

###############################################################################
# Race × Ethnicity
###############################################################################

cat("\n========================================\n")
cat("RACE × ETHNICITY\n")
cat("========================================\n")

print(table(subject_meta$race, subject_meta$ethnicity, useNA = "ifany"))
#                  Decline Hispanic NonHispanic not applicable
#                      0        0           1              0
# AmericanIndian       0        1           0              0
# Asian                1        0           5              0
# Black                0        0          80              0
# Decline              0        1           0              0
# Other                0       12           2              1
# White                0        1          30              0

###############################################################################
# 24. Auditoría específica de la semana gestacional de toma de muestra
###############################################################################

# gest_week_collection es una variable dependiente del tiempo y NO un confusor
# basal. Se audita por separado para definir y comprobar la ventana gestacional
# de los análisis longitudinales de diversidad alfa y beta.

cat("\n========================================\n")
cat("AUDITORÍA DE gest_week_collection\n")
cat("========================================\n")

if (!"gest_week_collection" %in% colnames(meta)) {
  stop("ERROR: no se encontró gest_week_collection.")}

###############################################################################
# 24.1. Calidad del parseo
###############################################################################

raw_gw <- meta$gest_week_collection_raw
raw_nonmissing <- !is.na(raw_gw) & trimws(raw_gw) != ""
parsed_nonmissing <- !is.na(meta$gest_week_collection)

parse_summary <- data.frame(
  N_samples = nrow(meta),
  Raw_nonmissing = sum(raw_nonmissing),
  Parsed_nonmissing = sum(parsed_nonmissing),
  Parse_failed = sum(raw_nonmissing & !parsed_nonmissing),
  Missing_after_parse = sum(!parsed_nonmissing),
  Missing_percent = 100 * mean(!parsed_nonmissing),
  stringsAsFactors = FALSE)

cat("\nResumen de parseo:\n")
print(parse_summary)
#   N_samples Raw_nonmissing Parsed_nonmissing Parse_failed Missing_after_parse Missing_percent
# 1      2179           2179              2177            2                   2      0.09178522

cat("\nResumen global de semana gestacional:\n")
print(summary(meta$gest_week_collection))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.00   19.00   25.00   24.91   31.00   40.00       2

###############################################################################
# 24.2. Comprobar plausibilidad temporal
###############################################################################

meta <- meta %>%
  mutate(
    GestWeek_out_of_range = !is.na(gest_week_collection) &
      (gest_week_collection < 0 | gest_week_collection > 42),
    Collection_after_delivery = if ("gest_wk_delivery_num" %in% colnames(.)) {
      !is.na(gest_week_collection) & !is.na(gest_wk_delivery_num) &
        gest_week_collection > gest_wk_delivery_num
    } else FALSE)

timing_flags <- meta %>%
  filter(GestWeek_out_of_range | Collection_after_delivery) %>%
  dplyr::select(
    SampleID, host_subject_id, Cohort_analysis, Delivery_analysis,
    gest_week_collection_raw, gest_week_collection,
    any_of(c("gest_wk_delivery_raw", "gest_wk_delivery_num")),
    GestWeek_out_of_range, Collection_after_delivery)

cat("\nRegistros con posible inconsistencia temporal:", nrow(timing_flags), "\n")
# Registros con posible inconsistencia temporal: 1 

if (nrow(timing_flags) > 0) print(as.data.frame(timing_flags))

###############################################################################
# 24.3. Distribución por cohorte y outcome
###############################################################################

gestweek_sample_summary <- meta %>%
  group_by(Cohort_analysis, Delivery_analysis) %>%
  summarise(
    N_samples = sum(!is.na(gest_week_collection)),
    N_subjects = n_distinct(host_subject_id[!is.na(gest_week_collection)]),
    Mean = mean(gest_week_collection, na.rm = TRUE),
    SD = sd(gest_week_collection, na.rm = TRUE),
    Median = median(gest_week_collection, na.rm = TRUE),
    IQR = IQR(gest_week_collection, na.rm = TRUE),
    Min = min(gest_week_collection, na.rm = TRUE),
    Max = max(gest_week_collection, na.rm = TRUE),
    .groups = "drop")

cat("\nSemana gestacional por cohorte y outcome - nivel muestra:\n")
print(as.data.frame(gestweek_sample_summary))
#   Cohort_analysis Delivery_analysis N_samples N_subjects     Mean       SD Median   IQR Min Max
# 1        Stanford           Preterm       156          9 21.42949 7.348568     21 12.25   6  38
# 2        Stanford              Term       739         30 23.73748 8.689703     24 14.50   1  40
# 3             UAB           Preterm       392         41 23.83163 5.638488     23  9.00  10  36
# 4             UAB              Term       890         55 26.96067 6.707585     27 12.00  11  40

###############################################################################
# 24.4. Ventana gestacional de soporte común
###############################################################################

# La intersección de los rangos observados define una zona en la que los cuatro
# grupos tienen alguna representación. NO exige que todas las mujeres tengan
# muestras en cada semana; por ello se audita también el soporte semanal.
gw_group_ranges <- gestweek_sample_summary %>%
  dplyr::select(
    Cohort_analysis, Delivery_analysis, N_samples, N_subjects,
    Min_GW = Min, Max_GW = Max)

GW_common_min <- max(gw_group_ranges$Min_GW, na.rm = TRUE)
GW_common_max <- min(gw_group_ranges$Max_GW, na.rm = TRUE)

common_window <- data.frame(
  Common_GW_min = GW_common_min,
  Common_GW_max = GW_common_max,
  stringsAsFactors = FALSE)

cat("\nVentana gestacional de soporte común:",
    GW_common_min, "a", GW_common_max, "semanas\n")
# Ventana gestacional de soporte común: 11 a 36 semanas

###############################################################################
# 24.5. Cobertura semanal dentro de la ventana común
###############################################################################

gestweek_weekly_coverage <- meta %>%
  filter(!is.na(gest_week_collection)) %>%
  mutate(GestWeek_integer = floor(gest_week_collection)) %>%
  distinct(
    SampleID, host_subject_id, Cohort_analysis, Delivery_analysis,
    GestWeek_integer) %>%
  group_by(Cohort_analysis, Delivery_analysis, GestWeek_integer) %>%
  summarise(
    N_samples = n(),
    N_subjects = n_distinct(host_subject_id),
    .groups = "drop")

weekly_support <- tidyr::expand_grid(
  Cohort_analysis = c("Stanford", "UAB"),
  Delivery_analysis = c("Term", "Preterm"),
  GestWeek_integer = seq(floor(GW_common_min), floor(GW_common_max))) %>%
  left_join(
    gestweek_weekly_coverage,
    by = c("Cohort_analysis", "Delivery_analysis", "GestWeek_integer")) %>%
  mutate(
    N_samples = replace_na(N_samples, 0L),
    N_subjects = replace_na(N_subjects, 0L)) %>%
  group_by(GestWeek_integer) %>%
  summarise(
    Groups_with_samples = sum(N_samples > 0),
    Min_samples_per_group = min(N_samples),
    Min_subjects_per_group = min(N_subjects),
    .groups = "drop")

cat("\nCobertura semanal dentro de la ventana común:\n")
print(as.data.frame(weekly_support))

###############################################################################
# 24.6. Resumir calendario de muestreo por mujer
###############################################################################

gestweek_subject <- meta %>%
  group_by(host_subject_id, Cohort_analysis, Delivery_analysis) %>%
  summarise(
    N_samples = sum(!is.na(gest_week_collection)),
    First_GW = ifelse(all(is.na(gest_week_collection)), NA_real_,
                      min(gest_week_collection, na.rm = TRUE)),
    Last_GW = ifelse(all(is.na(gest_week_collection)), NA_real_,
                     max(gest_week_collection, na.rm = TRUE)),
    Mean_GW = ifelse(all(is.na(gest_week_collection)), NA_real_,
                     mean(gest_week_collection, na.rm = TRUE)),
    Median_GW = ifelse(all(is.na(gest_week_collection)), NA_real_,
                       median(gest_week_collection, na.rm = TRUE)),
    Followup_span = ifelse(all(is.na(gest_week_collection)), NA_real_,
                           max(gest_week_collection, na.rm = TRUE) -
                             min(gest_week_collection, na.rm = TRUE)),
    .groups = "drop")

gestweek_subject_summary <- gestweek_subject %>%
  group_by(Cohort_analysis, Delivery_analysis) %>%
  summarise(
    N_subjects = n(),
    Subjects_with_GW = sum(N_samples > 0),
    Median_N_samples = median(N_samples, na.rm = TRUE),
    Median_First_GW = median(First_GW, na.rm = TRUE),
    Median_Last_GW = median(Last_GW, na.rm = TRUE),
    Median_Mean_GW = median(Mean_GW, na.rm = TRUE),
    Median_Followup_span = median(Followup_span, na.rm = TRUE),
    .groups = "drop")

cat("\nCalendario de muestreo por mujer:\n")
print(as.data.frame(gestweek_subject_summary))

###############################################################################
# 24.7. Comparaciones descriptivas Term vs Preterm
###############################################################################

gestweek_tests <- list()
variables_gw <- c("First_GW", "Last_GW", "Mean_GW", "Followup_span", "N_samples")

for (coh in c("Stanford", "UAB")) {
  dat_cohort <- gestweek_subject %>% filter(Cohort_analysis == coh)
  for (v in variables_gw) {
    dat_test <- dat_cohort %>%
      dplyr::select(Delivery_analysis, value = all_of(v)) %>%
      filter(!is.na(Delivery_analysis), !is.na(value))
    if (n_distinct(dat_test$Delivery_analysis) == 2) {
      wt <- wilcox.test(value ~ Delivery_analysis, data = dat_test, exact = FALSE)
      gestweek_tests[[paste(coh, v, sep = "_")]] <- data.frame(
        Cohort = coh, Variable = v, N = nrow(dat_test),
        Pvalue = wt$p.value, stringsAsFactors = FALSE)}}}

GestWeek_Tests <- bind_rows(gestweek_tests)

cat("\nComparaciones Term vs Preterm a nivel de mujer:\n")
print(as.data.frame(GestWeek_Tests))

###############################################################################
# 24.8. Distribución por intervalos gestacionales
###############################################################################

meta <- meta %>%
  mutate(Gestational_interval = cut(
    gest_week_collection,
    breaks = c(-Inf, 12, 18, 24, 30, 36, Inf),
    labels = c("<=12", "12-18", "18-24", "24-30", "30-36", ">36"),
    right = TRUE))

gestweek_interval_table <- meta %>%
  filter(!is.na(Gestational_interval)) %>%
  count(
    Cohort_analysis, Delivery_analysis,
    Gestational_interval, name = "N_samples") %>%
  group_by(Cohort_analysis, Delivery_analysis) %>%
  mutate(Percent_samples = 100 * N_samples / sum(N_samples)) %>%
  ungroup()

cat("\nDistribución de muestras por intervalo gestacional:\n")
print(as.data.frame(gestweek_interval_table))

###############################################################################
# 24.9. Guardar resultados
###############################################################################

write_csv(parse_summary,
  file.path(results_dir, "GestWeek_Parse_Summary.csv"))
write_csv(timing_flags,
  file.path(results_dir, "GestWeek_Timing_Flags.csv"))
write_csv(gestweek_sample_summary,
  file.path(results_dir, "GestWeek_Sample_Summary.csv"))
write_csv(gw_group_ranges,
  file.path(results_dir, "GestWeek_Group_Ranges.csv"))
write_csv(common_window,
  file.path(results_dir, "GestWeek_Common_Window.csv"))
write_csv(gestweek_weekly_coverage,
  file.path(results_dir, "GestWeek_Weekly_Coverage.csv"))
write_csv(weekly_support,
  file.path(results_dir, "GestWeek_Weekly_CommonSupport.csv"))
write_csv(gestweek_subject,
  file.path(results_dir, "GestWeek_Subject_Data.csv"))
write_csv(gestweek_subject_summary,
  file.path(results_dir, "GestWeek_Subject_Summary.csv"))
write_csv(GestWeek_Tests,
  file.path(results_dir, "GestWeek_Term_Preterm_Tests.csv"))
write_csv(gestweek_interval_table,
  file.path(results_dir, "GestWeek_Intervals.csv"))

# Actualizar el Excel para incluir también la auditoría gestacional.
excel_sheets$GestWeek_Parse <- parse_summary
excel_sheets$GestWeek_Timing_Flags <- timing_flags
excel_sheets$GestWeek_Sample <- gestweek_sample_summary
excel_sheets$GestWeek_Group_Ranges <- gw_group_ranges
excel_sheets$GestWeek_Common_Window <- common_window
excel_sheets$GestWeek_Weekly <- gestweek_weekly_coverage
excel_sheets$GestWeek_WeeklySupport <- weekly_support
excel_sheets$GestWeek_Subject <- gestweek_subject
excel_sheets$GestWeek_Tests <- GestWeek_Tests
excel_sheets$GestWeek_Intervals <- gestweek_interval_table

if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(
    excel_sheets,
    file = file.path(
      results_dir,
      "Confounder_Audit.xlsx"),
    rowNames = FALSE,
    overwrite = TRUE)
}

cat("\nAuditoría de gest_week_collection finalizada.\n")


###############################################################################
# Reproducibilidad
###############################################################################

software_versions <- data.frame(
  Package = c("R", "phyloseq", "dplyr", "tidyr", "readr"),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("phyloseq")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("tidyr")),
    as.character(packageVersion("readr"))),
  stringsAsFactors = FALSE)

write_csv(
  software_versions,
  file.path(dirname(results_dir), "software_versions.csv"))

cat("\n=== VERSIONES ===\n")
print(software_versions)
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

###############################################################################
# FIN
###############################################################################
