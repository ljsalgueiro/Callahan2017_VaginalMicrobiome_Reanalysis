###############################################################################
#                                    TFM
###############################################################################
# Script 07A: Prevotella - exportación, BLASTn y procesamiento
###############################################################################

# Objetivo:
#   1. Seleccionar dinámicamente los ASVs clasificados como Prevotella.
#   2. Exportar mapping y FASTA.
#   3. Ejecutar BLASTn remoto en NCBI.
#   4. Guardar RID, respuesta de envío y metadatos completos de la búsqueda.
#   5. Recuperar desde la salida remota el build de la base realmente usada.
#   6. Descargar XML2 y tabla tabular.
#   7. Calcular Query Cover y generar la tabla maestra de hits.
#
# BLASTn:
#   - program: blastn
#   - task: megablast
#   - database solicitada: nt
#   - filtro taxonómico: Prevotella[Organism]
#   - word_size: 28
#   - evalue: 10
#   - HITLIST_SIZE: 500
#
# Reproducibilidad de la base remota:
#   Además de la fecha de consulta, el script recupera de la propia salida
#   de NCBI BLAST la base reportada, su "Posted date" (build de la base),
#   el número de secuencias, el número de letras y la versión de BLAST.
#   Se guardan también el HTML de resultados y XML2 originales.
#
# IMPORTANTE:
#   HITLIST_SIZE limita el conjunto de hits devuelto por BLAST y no garantiza
#   que sean los 500 mejores alineamientos globales.
#
# Entrada:
#   - cleaned/ps_gard_refined.rds
#
# Salidas principales:
#   - Prevotella_ASV_mapping.csv
#   - Prevotella_ASVs.fasta
#   - Prevotella_BLAST_RID.txt
#   - Prevotella_BLAST_run_metadata.csv
#   - Prevotella_BLAST_database_metadata.csv
#   - Prevotella_BLAST_table_raw.csv
#   - Prevotella_BLAST_results_all.csv
#   - HTML y XML2 originales del RID

###############################################################################

rm(list = ls())

###############################################################################
# 1. Librerías, rutas y parámetros
###############################################################################

library(phyloseq)
library(Biostrings)
library(httr)
library(dplyr)
library(readr)
library(stringr)
library(xml2)
library(tibble)

project_dir <- if (
  dir.exists("VaginalMicrobiome")) {
  "VaginalMicrobiome"
} else {
  "."
}

input_file <- file.path(
  project_dir,
  "cleaned", "ps_gard_refined.rds")

base_dir <- file.path(
  project_dir,
  "results",
  "Prevotella")

blast_dir <- file.path(
  base_dir,
  "Blast")

dir.create(
  base_dir,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  blast_dir,
  recursive = TRUE,
  showWarnings = FALSE)

if (!file.exists(input_file))
  stop(
    paste(
      "No se encontró:",
      input_file))

blast_url <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

PROGRAM <- "blastn"
DATABASE <- "nt"
TASK <- "megablast"
WORD_SIZE <- 28
EVALUE <- 10
HITLIST_SIZE <- 500
ORGANISM_FILTER <- "Prevotella[Organism]"
TOOL <- "VaginalMicrobiome_TFM"
POLL_SECONDS <- 60

EMAIL <- Sys.getenv("NCBI_EMAIL")

if (EMAIL == "")
  stop(
    paste(
      "Defina NCBI_EMAIL antes de ejecutar BLAST.",
      "Ejemplo:",
      "Sys.setenv(NCBI_EMAIL = 'usuario@dominio.org')"))

###############################################################################
# 2. Funciones auxiliares
###############################################################################

clean_text <- function(x) {
  str_squish(
    as.character(x))
}

get_xml_text <- function(
  node,
  xpath) {

  x <- xml_find_first(
    node,
    xpath)

  if (inherits(
    x,
    "xml_missing"))
    return(NA_character_)

  xml_text(x)
}

collapse_unique <- function(x) {

  x <- unique(
    na.omit(
      clean_text(x)))

  x <- x[
    x != ""]

  if (length(x) == 0)
    return(NA_character_)

  paste(
    x,
    collapse = " | ")
}

get_taxon_sequences <- function(
  ps,
  taxa_ids) {

  rs <- tryCatch(
    refseq(ps),
    error = function(e)
      NULL)

  if (!is.null(rs)) {

    rs_chr <- as.character(rs)
    rs_names <- names(rs)

    if (
      !is.null(rs_names) &&
      all(
        taxa_ids %in%
          rs_names)) {

      seqs <- rs_chr[
        match(
          taxa_ids,
          rs_names)]

      names(seqs) <-
        taxa_ids

      return(seqs)
    }
  }

  seqs <- taxa_ids

  if (!all(
    grepl(
      "^[ACGTNacgtn]+$",
      seqs)))
    stop(
      paste(
        "No se pudieron recuperar secuencias desde refseq()",
        "y los taxa_names no parecen secuencias de ADN."))

  names(seqs) <-
    taxa_ids

  seqs
}

extract_html_field <- function(
  doc,
  label) {

  rows <- xml_find_all(
    doc,
    ".//*[local-name()='tr']")

  for (row in rows) {

    cells <- xml_find_all(
      row,
      "./*[local-name()='th' or local-name()='td']")

    values <- clean_text(
      xml_text(cells))

    if (length(values) < 2)
      next

    hit <- which(
      tolower(values) ==
        tolower(label))

    if (
      length(hit) > 0 &&
      hit[[1]] <
        length(values))
      return(
        values[
          hit[[1]] + 1])
  }

  page_lines <- str_split(
    xml_text(doc),
    "\n",
    simplify = FALSE)[[1]]

  page_lines <- clean_text(
    page_lines)

  page_lines <- page_lines[
    page_lines != ""]

  hit <- which(
    tolower(page_lines) ==
      tolower(label))

  if (
    length(hit) > 0 &&
    hit[[1]] <
      length(page_lines))
    return(
      page_lines[
        hit[[1]] + 1])

  NA_character_
}

extract_database_metadata <- function(
  xml_files,
  html_file) {

  xml_docs <- lapply(
    xml_files,
    read_xml)

  xml_values <- function(xpath) {
    collapse_unique(
      unlist(
        lapply(
          xml_docs,
          function(doc)
            xml_text(
              xml_find_all(
                doc,
                xpath)))))
  }

  database_reported <- xml_values(
    ".//*[local-name()='search-target']//*[local-name()='db']")

  if (is.na(database_reported))
    database_reported <- xml_values(
      ".//*[local-name()='Target']/*[local-name()='db']")

  blast_program_reported <- xml_values(
    ".//*[local-name()='Report']/*[local-name()='program']")

  blast_version_reported <- xml_values(
    ".//*[local-name()='Report']/*[local-name()='version']")

  db_num_xml <- xml_values(
    ".//*[local-name()='Statistics']/*[local-name()='db-num']")

  db_len_xml <- xml_values(
    ".//*[local-name()='Statistics']/*[local-name()='db-len']")

  html_doc <- read_html(
    html_file)

  db_posted_date <- extract_html_field(
    html_doc,
    "Posted date")

  db_sequences_html <- extract_html_field(
    html_doc,
    "Number of sequences")

  db_letters_html <- extract_html_field(
    html_doc,
    "Number of letters")

  if (is.na(db_posted_date))
    stop(
      paste(
        "No se pudo recuperar 'Posted date' de la salida BLAST.",
        "Revise el HTML guardado en:",
        html_file))

  data.frame(
    Database_requested =
      DATABASE,
    Database_reported =
      database_reported,
    Database_build_posted_date =
      db_posted_date,
    Database_sequences_reported =
      ifelse(
        is.na(db_sequences_html),
        db_num_xml,
        db_sequences_html),
    Database_letters_reported =
      ifelse(
        is.na(db_letters_html),
        db_len_xml,
        db_letters_html),
    BLAST_program_reported =
      blast_program_reported,
    BLAST_version_reported =
      blast_version_reported,
    stringsAsFactors = FALSE)
}

calculate_query_cover <- function(
  query_from,
  query_to,
  query_length) {

  valid <-
    !is.na(query_from) &
    !is.na(query_to)

  if (
    !any(valid) ||
    is.na(query_length) ||
    query_length <= 0)
    return(NA_real_)

  covered <- unique(
    unlist(
      Map(
        seq.int,
        query_from[valid],
        query_to[valid])))

  covered <- covered[
    covered >= 1 &
    covered <= query_length]

  length(covered) /
    query_length *
    100
}

###############################################################################
# 3. Selección y exportación de ASVs
###############################################################################

ps <- readRDS(
  input_file)

tax <- as.data.frame(
  tax_table(ps),
  stringsAsFactors = FALSE)

tax$Taxon_ID <- rownames(
  tax)

selected <- tax %>%
  filter(
    Genus == "Prevotella")

n_selected <- nrow(
  selected)

if (n_selected == 0)
  stop(
    "No se encontraron ASVs de Prevotella.")

taxa_ids <- selected$Taxon_ID

sequence_vector <- get_taxon_sequences(
  ps,
  taxa_ids)

sequence_vector <- sequence_vector[
  taxa_ids]

asv_ids <- paste0(
  "Prev_ASV_",
  seq_along(
    taxa_ids))


otu <- as(
  otu_table(ps),
  "matrix")

total_reads <- if (
  taxa_are_rows(ps)) {
  rowSums(otu)
} else {
  colSums(otu)
}

selected$Sequence <-
  unname(
    sequence_vector)

selected$Total_Reads <-
  total_reads[
    taxa_ids]

selected <- selected %>%
  arrange(
    desc(
      Total_Reads)) %>%
  mutate(
    Rank = row_number(),
    Relative_Prevotella =
      Total_Reads /
      sum(
        Total_Reads))

# El ranking por abundancia define los IDs Prev_ASV_1, Prev_ASV_2, etc.
taxa_ids <- selected$Taxon_ID
sequence_vector <- selected$Sequence
asv_ids <- paste0(
  "Prev_ASV_",
  seq_len(
    nrow(selected)))

mapping <- selected %>%
  transmute(
    ASV_ID = asv_ids,
    Taxon_ID = Taxon_ID,
    Sequence = Sequence,
    Species = Species,
    Total_Reads = Total_Reads,
    Relative_Prevotella =
      Relative_Prevotella)

write_csv(
  selected,
  file.path(
    base_dir,
    "Prevotella_ASV_audit.csv"))


write_csv(
  mapping,
  file.path(
    base_dir,
    "Prevotella_ASV_mapping.csv"))

seqs <- DNAStringSet(
  mapping$Sequence)

names(seqs) <-
  mapping$ASV_ID

writeXStringSet(
  seqs,
  filepath = file.path(
    base_dir,
    "Prevotella_ASVs.fasta"),
  format = "fasta")

cat(
  "\nASVs de Prevotella exportados:",
  length(seqs),
  "\n")

if ("Species" %in%
  names(selected)) {

  cat(
    "Clasificación SILVA a especie:\n")

  print(
    table(
      selected$Species,
      useNA = "always"))
}

###############################################################################
# 4. Enviar búsqueda BLASTn remota
###############################################################################

fasta_text <- paste(
  paste0(
    ">",
    names(seqs)),
  as.character(seqs),
  sep = "\n",
  collapse = "\n")

submission_time <- Sys.time()

response <- POST(
  blast_url,
  body = list(
    CMD = "Put",
    PROGRAM = PROGRAM,
    DATABASE = DATABASE,
    QUERY = fasta_text,
    ENTREZ_QUERY =
      ORGANISM_FILTER,
    MEGABLAST = "on",
    WORD_SIZE = WORD_SIZE,
    EXPECT = EVALUE,
    HITLIST_SIZE =
      HITLIST_SIZE,
    TOOL = TOOL,
    EMAIL = EMAIL),
  encode = "form")

stop_for_status(
  response)

response_text <- content(
  response,
  as = "text",
  encoding = "UTF-8")

writeLines(
  response_text,
  file.path(
    blast_dir,
    "Prevotella_BLAST_submission_response.txt"),
  useBytes = TRUE)

rid_line <- grep(
  "^\\s*RID\\s*=",
  strsplit(
    response_text,
    "\n")[[1]],
  value = TRUE)

if (length(rid_line) == 0)
  stop(
    "NCBI no devolvió un RID.")

RID <- trimws(
  sub(
    "^\\s*RID\\s*=\\s*",
    "",
    rid_line[[1]]))

if (
  RID == "" ||
  RID == "RTOE")
  stop(
    "RID no válido.")

writeLines(
  RID,
  file.path(
    blast_dir,
    "Prevotella_BLAST_RID.txt"))

submission_metadata <- data.frame(
  Submission_time =
    format(
      submission_time,
      "%Y-%m-%d %H:%M:%S %Z"),
  Program_requested =
    PROGRAM,
  Task_requested =
    TASK,
  Database_requested =
    DATABASE,
  Organism_filter =
    ORGANISM_FILTER,
  Word_size =
    WORD_SIZE,
  Evalue =
    EVALUE,
  Max_target_seqs =
    HITLIST_SIZE,
  RID = RID,
  N_ASVs =
    length(seqs),
  Min_sequence_length =
    min(
      width(seqs)),
  Max_sequence_length =
    max(
      width(seqs)),
  stringsAsFactors = FALSE)

write_csv(
  submission_metadata,
  file.path(
    blast_dir,
    "Prevotella_BLAST_submission_metadata.csv"))

cat(
  "RID:",
  RID,
  "\nEsperando finalización de BLAST...\n")

repeat {

  Sys.sleep(
    POLL_SECONDS)

  status_response <- GET(
    blast_url,
    query = list(
      CMD = "Get",
      RID = RID,
      FORMAT_OBJECT =
        "SearchInfo"))

  stop_for_status(
    status_response)

  status_text <- content(
    status_response,
    as = "text",
    encoding = "UTF-8")

  cat(
    format(
      Sys.time(),
      "%H:%M:%S"),
    "- comprobando estado...\n")

  if (grepl(
    "Status=READY",
    status_text,
    fixed = TRUE))
    break

  if (grepl(
    "Status=FAILED",
    status_text,
    fixed = TRUE))
    stop(
      "BLAST informó que la búsqueda falló.")

  if (grepl(
    "Status=UNKNOWN",
    status_text,
    fixed = TRUE))
    stop(
      "RID desconocido o expirado.")
}

retrieval_time <- Sys.time()

###############################################################################
# 5. Descargar XML2, HTML y tabla tabular
###############################################################################

zip_file <- file.path(
  blast_dir,
  paste0(
    "Prevotella_BLAST_XML2_RID_",
    RID,
    ".zip"))

xml_dir <- file.path(
  blast_dir,
  paste0(
    "XML2_RID_",
    RID))

dir.create(
  xml_dir,
  recursive = TRUE,
  showWarnings = FALSE)

xml_response <- GET(
  blast_url,
  query = list(
    CMD = "Get",
    RID = RID,
    FORMAT_TYPE = "XML2"),
  write_disk(
    zip_file,
    overwrite = TRUE))

stop_for_status(
  xml_response)

unzip(
  zip_file,
  exdir = xml_dir,
  overwrite = TRUE)

xml_files <- list.files(
  xml_dir,
  pattern = paste0(
    "^",
    RID,
    "_[0-9]+\\.xml$"),
  full.names = TRUE,
  recursive = TRUE)

if (length(xml_files) == 0)
  stop(
    "No se encontraron archivos XML2.")

html_file <- file.path(
  blast_dir,
  paste0(
    "Prevotella_BLAST_result_RID_",
    RID,
    ".html"))

html_response <- GET(
  blast_url,
  query = list(
    CMD = "Get",
    RID = RID,
    FORMAT_TYPE = "HTML"),
  write_disk(
    html_file,
    overwrite = TRUE))

stop_for_status(
  html_response)

csv_response <- GET(
  blast_url,
  query = list(
    CMD = "Get",
    RID = RID,
    FORMAT_TYPE = "CSV",
    FORMAT_OBJECT = "Alignment",
    ALIGNMENT_VIEW = "Tabular",
    DESCRIPTIONS =
      HITLIST_SIZE))

stop_for_status(
  csv_response)

raw_table_file <- file.path(
  blast_dir,
  "Prevotella_BLAST_table_raw.csv")

writeLines(
  content(
    csv_response,
    as = "text",
    encoding = "UTF-8"),
  raw_table_file,
  useBytes = TRUE)

###############################################################################
# 6. Recuperar build de la base y versión de BLAST desde la salida
###############################################################################

database_metadata <- extract_database_metadata(
  xml_files,
  html_file)

database_metadata <- database_metadata %>%
  mutate(
    Database_access_date =
      as.character(
        as.Date(
          retrieval_time)),
    Results_retrieval_time =
      format(
        retrieval_time,
        "%Y-%m-%d %H:%M:%S %Z"),
    RID = RID,
    .before = 1)

write_csv(
  database_metadata,
  file.path(
    blast_dir,
    "Prevotella_BLAST_database_metadata.csv"))

run_metadata <- bind_cols(
  submission_metadata,
  database_metadata %>%
    select(
      Database_access_date,
      Results_retrieval_time,
      Database_reported,
      Database_build_posted_date,
      Database_sequences_reported,
      Database_letters_reported,
      BLAST_program_reported,
      BLAST_version_reported))

write_csv(
  run_metadata,
  file.path(
    blast_dir,
    "Prevotella_BLAST_run_metadata.csv"))

cat(
  "\n=== METADATOS DE LA BASE BLAST ===\n")

print(
  database_metadata)

if (
  !is.na(
    database_metadata$Database_reported) &&
  database_metadata$Database_reported !=
    DATABASE) {

  warning(
    paste(
      "La base reportada por NCBI difiere de la solicitada:",
      DATABASE,
      "->",
      database_metadata$Database_reported))
}

###############################################################################
# 7. Extraer taxonomía desde XML2
###############################################################################

xml_results <- vector(
  "list",
  length(xml_files))

for (f in seq_along(
  xml_files)) {

  xml_doc <- read_xml(
    xml_files[[f]])

  searches <- xml_find_all(
    xml_doc,
    ".//*[local-name()='Search']")

  search_results <- vector(
    "list",
    length(searches))

  for (s in seq_along(
    searches)) {

    search <- searches[[s]]

    asv_id <- get_xml_text(
      search,
      "./*[local-name()='query-title']")

    query_length <- suppressWarnings(
      as.numeric(
        get_xml_text(
          search,
          "./*[local-name()='query-len']")))

    hits <- xml_find_all(
      search,
      "./*[local-name()='hits']/*[local-name()='Hit']")

    if (length(hits) == 0) {

      search_results[[s]] <- tibble(
        ASV_ID = asv_id,
        Query_Length =
          query_length,
        Hit_Rank =
          NA_integer_,
        Accession_XML =
          NA_character_,
        Description =
          NA_character_,
        Scientific_Name =
          NA_character_)

      next
    }

    hit_results <- vector(
      "list",
      length(hits))

    for (h in seq_along(
      hits)) {

      hit <- hits[[h]]

      hit_rank <- suppressWarnings(
        as.integer(
          get_xml_text(
            hit,
            "./*[local-name()='num']")))

      if (is.na(hit_rank))
        hit_rank <- h

      descr_nodes <- xml_find_all(
        hit,
        ".//*[local-name()='HitDescr']")

      if (length(descr_nodes) == 0) {

        hit_results[[h]] <- tibble(
          ASV_ID = asv_id,
          Query_Length =
            query_length,
          Hit_Rank =
            hit_rank,
          Accession_XML =
            NA_character_,
          Description =
            NA_character_,
          Scientific_Name =
            NA_character_)

      } else {

        hit_results[[h]] <- bind_rows(
          lapply(
            descr_nodes,
            function(descr)
              tibble(
                ASV_ID = asv_id,
                Query_Length =
                  query_length,
                Hit_Rank =
                  hit_rank,
                Accession_XML =
                  get_xml_text(
                    descr,
                    "./*[local-name()='accession']"),
                Description =
                  get_xml_text(
                    descr,
                    "./*[local-name()='title']"),
                Scientific_Name =
                  get_xml_text(
                    descr,
                    "./*[local-name()='sciname']"))))
      }
    }

    search_results[[s]] <- bind_rows(
      hit_results)
  }

  xml_results[[f]] <- bind_rows(
    search_results)
}

blast_taxonomy <- bind_rows(
  xml_results) %>%
  mutate(
    Accession_Key =
      sub(
        "\\.[0-9]+$",
        "",
        Accession_XML)) %>%
  filter(
    !is.na(ASV_ID),
    !is.na(
      Accession_Key)) %>%
  arrange(
    ASV_ID,
    Hit_Rank)

blast_taxonomy_unique <- blast_taxonomy %>%
  group_by(
    ASV_ID,
    Accession_Key) %>%
  slice_min(
    order_by =
      Hit_Rank,
    n = 1,
    with_ties = FALSE) %>%
  ungroup() %>%
  select(
    ASV_ID,
    Query_Length,
    Hit_Rank,
    Accession_Key,
    Description,
    Scientific_Name)

###############################################################################
# 8. Procesar HSPs y calcular Query Cover
###############################################################################

blast_hsps <- read_csv(
  raw_table_file,
  col_names = FALSE,
  show_col_types = FALSE)

if (ncol(blast_hsps) != 12)
  stop(
    "La tabla tabular BLAST no contiene 12 columnas.")

colnames(blast_hsps) <- c(
  "Query",
  "Accession",
  "PercentIdentity",
  "AlignmentLength",
  "Mismatches",
  "GapOpens",
  "QueryStart",
  "QueryEnd",
  "SubjectStart",
  "SubjectEnd",
  "Evalue",
  "BitScore")

query_lengths <- blast_taxonomy_unique %>%
  select(
    ASV_ID,
    Query_Length) %>%
  distinct()

blast_hsps <- blast_hsps %>%
  mutate(
    ASV_ID = Query,
    Accession_Key =
      sub(
        "\\.[0-9]+$",
        "",
        Accession),
    Query_From =
      pmin(
        QueryStart,
        QueryEnd),
    Query_To =
      pmax(
        QueryStart,
        QueryEnd)) %>%
  left_join(
    query_lengths,
    by = "ASV_ID")

query_coverage <- blast_hsps %>%
  group_by(
    ASV_ID,
    Accession) %>%
  summarise(
    Query_Cover =
      calculate_query_cover(
        Query_From,
        Query_To,
        first(
          Query_Length)),
    .groups = "drop")

best_hsp <- blast_hsps %>%
  arrange(
    ASV_ID,
    Accession,
    desc(
      BitScore),
    Evalue,
    desc(
      PercentIdentity)) %>%
  group_by(
    ASV_ID,
    Accession) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    ASV_ID,
    Accession,
    Accession_Key,
    Per_Ident =
      PercentIdentity,
    E_value =
      Evalue,
    Bit_Score =
      BitScore)

blast_results_all <- best_hsp %>%
  left_join(
    query_coverage,
    by = c(
      "ASV_ID",
      "Accession")) %>%
  left_join(
    blast_taxonomy_unique,
    by = c(
      "ASV_ID",
      "Accession_Key")) %>%
  transmute(
    ASV_ID,
    Hit_Rank,
    Accession,
    Description,
    Scientific_Name,
    Query_Cover,
    Per_Ident,
    E_value,
    Bit_Score) %>%
  filter(
    is.na(
      Hit_Rank) |
    Hit_Rank <=
      HITLIST_SIZE) %>%
  mutate(
    ASV_number =
      as.numeric(
        str_remove(
          ASV_ID,
          "^[^0-9]+"))) %>%
  arrange(
    ASV_number,
    Hit_Rank) %>%
  select(
    -ASV_number)

write_csv(
  blast_results_all,
  file.path(
    blast_dir,
    "Prevotella_BLAST_results_all.csv"))

cat(
  "\nBLAST Prevotella COMPLETADO\n")

cat(
  "ASVs con hits:",
  n_distinct(
    blast_results_all$ASV_ID),
  "\n")

cat(
  "Hits/accessions:",
  nrow(
    blast_results_all),
  "\n")

cat(
  "Máximo hits/ASV:",
  max(
    table(
      blast_results_all$ASV_ID)),
  "\n")

###############################################################################
# 9. Reproducibilidad
###############################################################################

software_versions <- data.frame(
  Package = c(
    "R",
    "phyloseq",
    "Biostrings",
    "httr",
    "dplyr",
    "readr",
    "stringr",
    "xml2",
    "tibble"),
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
        "Biostrings")),
    as.character(
      packageVersion(
        "httr")),
    as.character(
      packageVersion(
        "dplyr")),
    as.character(
      packageVersion(
        "readr")),
    as.character(
      packageVersion(
        "stringr")),
    as.character(
      packageVersion(
        "xml2")),
    as.character(
      packageVersion(
        "tibble"))),
  stringsAsFactors = FALSE)

write_csv(
  software_versions,
  file.path(
    blast_dir,
    "Prevotella_BLAST_software_versions.csv"))

session_info_file <- file.path(
  blast_dir,
  "Prevotella_BLAST_sessionInfo.txt")

writeLines(
  capture.output(
    sessionInfo()),
  session_info_file)

cat(
  "\n=== VERSIONES ===\n")

print(
  software_versions)

cat(
  "\n=== SESSION INFO ===\n")

print(
  sessionInfo())

###############################################################################
# FIN
###############################################################################
