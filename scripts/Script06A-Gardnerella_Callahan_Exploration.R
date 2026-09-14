###############################################################################
#                                    TFM
###############################################################################
#       Script 06A: Exploración de Gardnerella en Callahan et al. (2017)
###############################################################################

# Objetivo:
# Recuperar de processed.rda las variantes de Gardnerella utilizadas por
# Callahan et al. y exportar las tres primeras secuencias (G1, G2 y G3).
#
# Este script es exploratorio y no modifica ningún objeto phyloseq del TFM.
# El refinamiento definitivo de Gardnerella se realiza en Script 06B.

###############################################################################

rm(list = ls())

library(readr)

project_dir <- if (dir.exists("VaginalMicrobiome")) {
  "VaginalMicrobiome"
} else {
  "."
}

input_file <- file.path(
  project_dir,
  "Callahan2017",
  "RepRefine_Scripts",
  "input",
  "processed.rda")

results_dir <- file.path(
  project_dir,
  "results",
  "Gardnerella")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE)

if (!file.exists(input_file))
  stop(
    paste(
      "No se encontró:",
      input_file))

callahan_env <- new.env()

load(
  input_file,
  envir = callahan_env)

if (!exists(
  "tax",
  envir = callahan_env,
  inherits = FALSE))
  stop(
    "processed.rda no contiene el objeto 'tax'.")

tax_callahan <- get(
  "tax",
  envir = callahan_env)

if (!"Genus" %in%
  colnames(tax_callahan))
  stop(
    "La tabla tax de processed.rda no contiene la columna Genus.")

gard_idx <- which(
  tax_callahan[
    ,
    "Genus"] ==
    "Gardnerella")

if (length(gard_idx) < 3)
  stop(
    "processed.rda contiene menos de tres variantes Gardnerella.")

G_reference <- c(
  G1 =
    rownames(tax_callahan)[
      gard_idx[[1]]],
  G2 =
    rownames(tax_callahan)[
      gard_idx[[2]]],
  G3 =
    rownames(tax_callahan)[
      gard_idx[[3]]])

reference_table <- data.frame(
  GardnerellaVariant =
    names(G_reference),
  Callahan_index =
    gard_idx[1:3],
  Sequence =
    unname(G_reference),
  Length =
    nchar(
      unname(G_reference)),
  stringsAsFactors = FALSE)

cat(
  "\nVariantes Gardnerella en processed.rda:",
  length(gard_idx),
  "\n")

cat(
  "\nReferencias G1/G2/G3:\n")

print(
  reference_table)

write_csv(
  reference_table,
  file.path(
    results_dir,
    "Callahan_Gardnerella_G1_G2_G3_reference.csv"))

software_versions <- data.frame(
  Package = c(
    "R",
    "readr"),
  Version = c(
    paste(
      R.version$major,
      R.version$minor,
      sep = "."),
    as.character(
      packageVersion(
        "readr"))),
  stringsAsFactors = FALSE)

write_csv(
  software_versions,
  file.path(
    results_dir,
    "Script06A_SoftwareVersions.csv"))

cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

###############################################################################
# FIN
###############################################################################
