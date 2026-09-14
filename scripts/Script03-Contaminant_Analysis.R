###############################################################################
#                                    TFM
###############################################################################
#                   Script 03: Evaluación de contaminación
###############################################################################

# Objetivo:
#   Identificar y eliminar posibles ASVs contaminantes mediante decontam,
#   utilizando los controles de extracción y considerando la corrida de
#   secuenciación. Posteriormente, eliminar los controles y conservar
#   exclusivamente las muestras vaginales.

# Método:
#   - Controles negativos: env_biome == "laboratory"
#   - decontam: método de prevalencia
#   - Lote: IL_Run
#   - Combinación entre corridas: Fisher
#   - Análisis principal: umbral 0.1 (valor por defecto de decontam)
#   - Sensibilidad: umbral 0.5
#   - Sensibilidad adicional con método combinado si existe una variable
#     válida de concentración de ADN en los metadatos

# Entrada:
#   VaginalMicrobiome/cleaned/phyloseqObject_sequences.rds

# Salidas:
#   VaginalMicrobiome/cleaned/ps_decontam.rds
#   VaginalMicrobiome/results/Contamination/Decontam_Results.csv
#   VaginalMicrobiome/results/Contamination/Decontam_Summary.csv
#   VaginalMicrobiome/results/Contamination/Decontam_Threshold_Sensitivity.csv
#   VaginalMicrobiome/results/Contamination/Decontam_Combined_Sensitivity.csv

###############################################################################

rm(list = ls())

###############################################################################
# Configuración de rutas
###############################################################################

# El script funciona tanto si R se inicia en la raíz del repositorio como si
# se inicia en la carpeta inmediatamente superior a "VaginalMicrobiome".
project_dir <- if (dir.exists("cleaned")) {
  "."
} else if (dir.exists(file.path("VaginalMicrobiome", "cleaned"))) {
  "VaginalMicrobiome"
} else {
  stop(
    "ERROR: no se encontró la carpeta 'cleaned'. ",
    "Ejecute el script desde la raíz del repositorio VaginalMicrobiome ",
    "o desde su carpeta padre."
  )
}

project_path <- function(...) file.path(project_dir, ...)

###############################################################################
# 1. Librerías
###############################################################################

library(phyloseq)
library(decontam)
library(dplyr)
library(readr)

###############################################################################
# 2. Carga del objeto phyloseq
###############################################################################

ps <- readRDS(project_path("cleaned", "phyloseqObject_sequences.rds"))

cat("\n=== OBJETO INICIAL ===\n")
cat("Registros:", nsamples(ps), "\n")
# Registros: 2367 

cat("ASVs:", ntaxa(ps), "\n")
# ASVs: 8298 

###############################################################################
# 3. Metadata e identificación de controles
###############################################################################

meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
meta$SampleID <- rownames(meta)

required_vars <- c("env_biome", "env_feature", "IL_Run")
missing_vars <- setdiff(required_vars, names(meta))

if (length(missing_vars) > 0) stop(paste("ERROR: faltan variables en metadata:", paste(missing_vars, collapse = ", ")))

meta$is.neg <- meta$env_biome == "laboratory"
sample_data(ps)$is.neg <- meta$is.neg

###############################################################################
# 4. Comprobación de muestras vaginales y controles
###############################################################################

cat("\n=== TIPO DE REGISTRO ===\n")
print(table(meta$env_biome, useNA = "ifany"))
# laboratory     vagina 
# 188       2179

cat("\n=== env_feature ===\n")
print(table(meta$env_feature, useNA = "ifany"))
#         mid-vaginal wall sham extraction controls 
#                    2179                      188 

cat("\n=== CONTROLES NEGATIVOS ===\n")
print(table(meta$is.neg, useNA = "ifany"))
# FALSE  TRUE 
# 2179   188 

n_controls <- sum(meta$is.neg, na.rm = TRUE)
n_vaginal <- sum(!meta$is.neg, na.rm = TRUE)

cat("\nControles de extracción:", n_controls, "\n")
# Controles de extracción: 188 

cat("Muestras vaginales:", n_vaginal, "\n")
# Muestras vaginales: 2179 

stopifnot(n_controls == 188, n_vaginal == 2179)

###############################################################################
# 5. Profundidad de secuenciación de los controles
###############################################################################

reads <- sample_sums(ps)

cat("\n=== LECTURAS EN CONTROLES ===\n")
print(summary(reads[meta$is.neg]))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 4.0    182.8    423.5   7360.5   2080.8 229960.0 

cat("Controles con >0 lecturas:", sum(reads[meta$is.neg] > 0), "\n")
# Controles con >0 lecturas: 188 

cat("Controles con 0 lecturas:", sum(reads[meta$is.neg] == 0), "\n")
# Controles con 0 lecturas: 0 

###############################################################################
# 6. Distribución por corrida
###############################################################################

if (anyNA(meta$IL_Run)) stop("ERROR: existen registros sin información de IL_Run.")

batch_table <- table(meta$is.neg, meta$IL_Run)

cat("\n=== TIPO DE REGISTRO POR CORRIDA ===\n")
print(batch_table)

if (any(batch_table == 0)) stop("ERROR: existe al menos una corrida sin controles o sin muestras vaginales.")

cat("\n=== CONTROLES POR CORRIDA ===\n")
print(table(meta$IL_Run[meta$is.neg]))
# IL02 IL03 IL05 IL06 IL07 IL08 IL09 
#   26   32   29   26   24   22   29 

cat("\n=== MUESTRAS VAGINALES POR CORRIDA ===\n")
print(table(meta$IL_Run[!meta$is.neg]))
# IL02 IL03 IL05 IL06 IL07 IL08 IL09 
# 156  219  387  682  114   21  600 

###############################################################################
# 7. Identificación de contaminantes mediante decontam
###############################################################################

# Análisis principal.
# threshold = 0.1 es el valor por defecto de isContaminant() y se mantiene
# como criterio principal para el filtrado.
contam <- isContaminant(
  ps,
  method = "prevalence",
  neg = "is.neg",
  batch = "IL_Run",
  batch.combine = "fisher",
  threshold = 0.1
)

cat("\n=== DECONTAM: ANÁLISIS PRINCIPAL (threshold = 0.1) ===\n")
print(table(contam$contaminant, useNA = "ifany"))

n_contaminants <- sum(contam$contaminant, na.rm = TRUE)
percent_contaminants <- round(100 * n_contaminants / ntaxa(ps), 2)

cat("\nASVs identificados como posibles contaminantes:", n_contaminants, "\n")
cat("Porcentaje de ASVs identificados:", percent_contaminants, "%\n")

###############################################################################
# 7A. Sensibilidad al umbral: threshold = 0.5
###############################################################################

contam_05 <- isContaminant(
  ps,
  method = "prevalence",
  neg = "is.neg",
  batch = "IL_Run",
  batch.combine = "fisher",
  threshold = 0.5
)

n_contaminants_05 <- sum(contam_05$contaminant, na.rm = TRUE)

threshold_sensitivity <- data.frame(
  Threshold = c(0.1, 0.5),
  Contaminant_ASVs = c(n_contaminants, n_contaminants_05),
  Percent_of_initial_ASVs = round(
    100 * c(n_contaminants, n_contaminants_05) / ntaxa(ps),
    2
  )
)

added_at_05 <- sum(
  contam_05$contaminant & !contam$contaminant,
  na.rm = TRUE
)

lost_at_05 <- sum(
  contam$contaminant & !contam_05$contaminant,
  na.rm = TRUE
)

cat("\n=== SENSIBILIDAD AL UMBRAL ===\n")
print(threshold_sensitivity)
cat(
  "ASVs adicionales identificados con threshold = 0.5:",
  added_at_05,
  "\n"
)
cat(
  "ASVs identificados con 0.1 pero no con 0.5:",
  lost_at_05,
  "\n"
)

###############################################################################
# 7B. Comprobación de concentración de ADN y método combinado
###############################################################################

# decontam requiere una medida cuantitativa de concentración de ADN para
# utilizar el componente de frecuencia. Se buscan nombres de variable
# compatibles de forma conservadora.
conc_candidates <- names(meta)[
  grepl(
    "dna.*conc|conc.*dna|concentration|ng[._ -]?(ul|µl)",
    names(meta),
    ignore.case = TRUE
  )
]

combined_status <- data.frame(
  Combined_method_run = FALSE,
  Concentration_variable = NA_character_,
  Reason = NA_character_,
  stringsAsFactors = FALSE
)

contam_combined <- NULL

if (length(conc_candidates) == 0) {

  combined_status$Reason <-
    "No se encontró una variable de concentración de ADN en los metadatos."

  cat(
    "\nMétodo combinado no ejecutado:",
    combined_status$Reason,
    "\n"
  )

} else if (length(conc_candidates) > 1) {

  combined_status$Reason <- paste0(
    "Se encontraron múltiples variables candidatas: ",
    paste(conc_candidates, collapse = ", "),
    ". Se requiere selección manual."
  )

  cat(
    "\nMétodo combinado no ejecutado:",
    combined_status$Reason,
    "\n"
  )

} else {

  conc_var <- conc_candidates[1]
  conc_values <- suppressWarnings(as.numeric(meta[[conc_var]]))

  if (anyNA(conc_values) || any(!is.finite(conc_values)) ||
      any(conc_values <= 0)) {

    combined_status$Concentration_variable <- conc_var
    combined_status$Reason <- paste0(
      "La variable ", conc_var,
      " contiene valores ausentes, no numéricos, no finitos o <= 0."
    )

    cat(
      "\nMétodo combinado no ejecutado:",
      combined_status$Reason,
      "\n"
    )

  } else {

    sample_data(ps)[[conc_var]] <- conc_values

    contam_combined <- isContaminant(
      ps,
      method = "combined",
      neg = "is.neg",
      conc = conc_var,
      batch = "IL_Run",
      batch.combine = "fisher",
      threshold = 0.1
    )

    combined_status$Combined_method_run <- TRUE
    combined_status$Concentration_variable <- conc_var
    combined_status$Reason <- "Método combinado ejecutado correctamente."

    cat("\n=== SENSIBILIDAD: MÉTODO COMBINADO ===\n")
    print(table(contam_combined$contaminant, useNA = "ifany"))
  }
}

###############################################################################
# 8. Tabla de resultados de decontam
###############################################################################

tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
tax$ASV <- rownames(tax)

otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)

stopifnot(identical(rownames(otu), taxa_names(ps)))
stopifnot(identical(colnames(otu), rownames(meta)))

decontam_results <- data.frame(
  ASV = taxa_names(ps),
  Contaminant = contam$contaminant,
  Contaminant_threshold_0.5 = contam_05$contaminant,
  Prev_control = rowMeans(otu[, meta$is.neg, drop = FALSE] > 0),
  Prev_vaginal = rowMeans(otu[, !meta$is.neg, drop = FALSE] > 0),
  Reads_control = rowSums(otu[, meta$is.neg, drop = FALSE]),
  Reads_vaginal = rowSums(otu[, !meta$is.neg, drop = FALSE]),
  stringsAsFactors = FALSE)

decontam_results <- left_join(decontam_results, tax, by = "ASV")

decontam_results$Contaminant_combined <- if (!is.null(contam_combined)) {
  contam_combined$contaminant
} else {
  NA
}

###############################################################################
# 9. Resumen taxonómico de los posibles contaminantes
###############################################################################

contaminant_taxa <- decontam_results %>% filter(Contaminant)

cat("\n=== POSIBLES CONTAMINANTES ===\n")
cat("Número de ASVs:", nrow(contaminant_taxa), "\n")

cat("\n=== GÉNEROS MÁS FRECUENTES ===\n")
genus_counts <- sort(table(contaminant_taxa$Genus), decreasing = TRUE)
print(head(genus_counts, 20))

cat("\n=== LACTOBACILLUS, GARDNERELLA Y PREVOTELLA ===\n")

interest_contaminants <- contaminant_taxa %>%
  filter(Genus %in% c("Lactobacillus", "Gardnerella", "Prevotella")) %>%
  arrange(Genus, desc(Reads_vaginal)) %>%
  select(ASV, Genus, Species, Prev_control, Prev_vaginal, Reads_control, Reads_vaginal)

print(interest_contaminants)
#           Genus        Species Prev_control Prev_vaginal Reads_control Reads_vaginal
#1 Lactobacillus           <NA>   0.02127660  0.001376778            72           838
#2    Prevotella melaninogenica   0.05851064  0.089490592          6381        385313
#3    Prevotella melaninogenica   0.04787234  0.004130335          2521           799
#4    Prevotella    nanceiensis   0.03191489  0.005048187           815           325

###############################################################################
# 10. Eliminación de ASVs contaminantes
###############################################################################

ps_decontam <- prune_taxa(!contam$contaminant, ps)

cat("\n=== DESPUÉS DE DECONTAM ===\n")
cat("ASVs iniciales:", ntaxa(ps), "\n")
# ASVs iniciales: 8298 
cat("ASVs contaminantes eliminados:", n_contaminants, "\n")
# ASVs contaminantes eliminados: 146 
cat("ASVs conservados:", ntaxa(ps_decontam), "\n")
# ASVs conservados: 8152 

###############################################################################
# 11. Eliminación de controles de extracción
###############################################################################

ps_clean <- subset_samples(ps_decontam, env_biome == "vagina")

cat("\n=== DESPUÉS DE RETIRAR CONTROLES ===\n")
cat("Muestras vaginales:", nsamples(ps_clean), "\n")
# Muestras vaginales: 2179 

###############################################################################
# 12. Eliminación de ASVs sin representación vaginal
###############################################################################

ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)

###############################################################################
# 13. Comprobaciones finales
###############################################################################

final_samples <- nsamples(ps_clean)
final_asvs <- ntaxa(ps_clean)
zero_samples <- sum(sample_sums(ps_clean) == 0)

cat("\n=== OBJETO FINAL ===\n")
cat("Muestras vaginales:", final_samples, "\n")
# Muestras vaginales: 2179 
cat("ASVs finales:", final_asvs, "\n")
# ASVs finales: 6539 
cat("Muestras con cero lecturas:", zero_samples, "\n")
# Muestras con cero lecturas: 0 
stopifnot(final_samples == 2179, zero_samples == 0)

###############################################################################
# 14. Tabla resumen
###############################################################################

decontam_summary <- data.frame(
  Initial_records = nsamples(ps),
  Vaginal_samples = final_samples,
  Extraction_controls = n_controls,
  Initial_ASVs = ntaxa(ps),
  Contaminant_ASVs = n_contaminants,
  ASVs_after_decontam = ntaxa(ps_decontam),
  Final_ASVs = final_asvs,
  Percent_contaminants = percent_contaminants,
  Contaminant_ASVs_threshold_0.5 = n_contaminants_05,
  Additional_ASVs_threshold_0.5 = added_at_05,
  Zero_read_vaginal_samples = zero_samples)

cat("\n=== RESUMEN FINAL ===\n")
print(decontam_summary)

###############################################################################
# 15. Guardado
###############################################################################

output_dir <- project_path("results", "Contamination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("cleaned"), recursive = TRUE, showWarnings = FALSE)

write_csv(
  decontam_results,
  file.path(output_dir, "Decontam_Results.csv")
)

write_csv(
  decontam_summary,
  file.path(output_dir, "Decontam_Summary.csv")
)

write_csv(
  threshold_sensitivity,
  file.path(output_dir, "Decontam_Threshold_Sensitivity.csv")
)

write_csv(
  combined_status,
  file.path(output_dir, "Decontam_Combined_Sensitivity.csv")
)

saveRDS(
  ps_clean,
  project_path("cleaned", "ps_decontam.rds")
)

cat("\n=== ARCHIVOS GENERADOS ===\n")
cat("- VaginalMicrobiome/cleaned/ps_decontam.rds\n")
cat("-", file.path(output_dir, "Decontam_Results.csv"), "\n")
cat("-", file.path(output_dir, "Decontam_Summary.csv"), "\n")

###############################################################################
# Fin del Script 03
###############################################################################

