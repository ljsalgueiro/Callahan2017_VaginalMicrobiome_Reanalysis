###############################################################################
#                                    TFM 
###############################################################################
#           Script 01 A: Carga e inspección de objeto phyloseq
###############################################################################

# Objetivo:
# Cargar el objeto phyloseq generado por DADA2 en Galaxy e inspeccionar
# su estructura antes de comenzar cualquier filtrado o análisis.
#
# No modifica el objeto original.

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
library(ggplot2)
library(dplyr)
library(tibble)

###############################################################################
# 2. Cargar objeto phyloseq
###############################################################################

ps_original <- readRDS(
  project_path("cleaned", "phyloseqObject.phyloseq"))

# Trabajaremos siempre sobre una copia
ps <- ps_original

##############################################################################
# 3. Inspección general de la estructura Phyloseq
##############################################################################

ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8298 taxa and 2367 samples ]
# sample_data() Sample Data:       [ 2367 samples by 53 sample variables ]
# tax_table()   Taxonomy Table:    [ 8298 taxa by 7 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 8298 reference sequences ]

nsamples(ps) # [1] 2367 muestras
ntaxa(ps) # [1] 8298 ASVs

sample_variables(ps) # Metadata
# [1] "AGE"                            "Assay.Type"                    
# [3] "AvgSpotLen"                     "BarcodeSequence"               
# [5] "Bases"                          "BioProject"                    
# [7] "BioSample"                      "BioSampleModel"                
# [9] "Bytes"                          "Center.Name"                   
# [11] "Collection_Date"                "Consent"                       
# [13] "DATASTORE.filetype"             "DATASTORE.provider"            
# [15] "DATASTORE.region"               "env_biome"                     
# [17] "env_feature"                    "env_material"                  
# [19] "ethnicity"                      "Experiment"                    
# [21] "geo_loc_name_country"           "geo_loc_name_country_continent"
# [23] "geo_loc_name"                   "gest_day_collection"           
# [25] "gest_day_delivery"              "gest_week_collection"          
# [27] "gest_wk_delivery"               "HOST"                          
# [29] "host_sex"                       "host_subject_id"               
# [31] "Groups"                         "cohort"                        
# [33] "host_tissue_sampled"            "IL_Run"                        
# [35] "indication_for_PTB"             "Instrument"                    
# [37] "lat_lon"                        "Library.Name"                  
# [39] "LibraryLayout"                  "LibrarySelection"              
# [41] "LibrarySource"                  "LinkerPrimerSequence"          
# [43] "NumberInRun"                    "Organism"                      
# [45] "Platform"                       "ReleaseDate"                   
# [47] "samp_collect_device"            "Sample.Name"                   
# [49] "SRA.Study"                      "term_vs_preterm_delivery"      
# [51] "create_date"                    "version"                       
# [53] "race"      

#############################################################################
# 4. Extraer tablas del objeto Phyloseq
#############################################################################

# Taxonomía (tax_table): Clasificación de los ASVs
#-------------------------------------
tax <- as.data.frame(tax_table(ps))

dim(tax) # [1] 8298 filas x 7 col
head(tax)
rank_names(ps) # columnas
# [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

# Abundancias (otu_table): Abundancia de taxones presentes en cada muestra
#---------------------------------------------------------------------------
# Es el conteo de leacturas ("reads") asignados a cada ASV por muestra
# Taxones = ASVs 
# Muestra = SRR

otu <- as(otu_table(ps), "matrix")

if(taxa_are_rows(ps)){otu <- t(otu)}
taxa_are_rows(ps)

dim(otu) # [1] 2367 filas (muestras = SRR) x 8298 col (ASVs)
head(otu)

# ASV (refseq):
#------
seqs <- as.character(refseq(ps))

names(seqs) <- taxa_names(ps)
names(seqs)
# [1] "ASV1"    "ASV2"    "ASV3"    "ASV4"    "ASV5"    "ASV6"    "ASV7"    "ASV8"  
#...
# Los ASVs no están nombrados con sus secuencias de nucleótidos.
# Estan únicamente numerados

# Metadata (sample_data): Datos clínicos
#--------------------------
meta <- as.data.frame(sample_data(ps))

dim(meta) # [1] 2367 muestras (SRRs) x 53 columnas (datos clínicos)
head(as.data.frame(meta))
names(meta)
# [1] "AGE"                            "Assay.Type"                    
# [3] "AvgSpotLen"                     "BarcodeSequence"               
# [5] "Bases"                          "BioProject"                    
# [7] "BioSample"                      "BioSampleModel"                
# [9] "Bytes"                          "Center.Name"                   
# [11] "Collection_Date"                "Consent"                       
# [13] "DATASTORE.filetype"             "DATASTORE.provider"            
# [15] "DATASTORE.region"               "env_biome"                     
# [17] "env_feature"                    "env_material"                  
# [19] "ethnicity"                      "Experiment"                    
# [21] "geo_loc_name_country"           "geo_loc_name_country_continent"
# [23] "geo_loc_name"                   "gest_day_collection"           
# [25] "gest_day_delivery"              "gest_week_collection"          
# [27] "gest_wk_delivery"               "HOST"                          
# [29] "host_sex"                       "host_subject_id"               
# [31] "Groups"                         "cohort"                        
# [33] "host_tissue_sampled"            "IL_Run"                        
# [35] "indication_for_PTB"             "Instrument"                    
# [37] "lat_lon"                        "Library.Name"                  
# [39] "LibraryLayout"                  "LibrarySelection"              
# [41] "LibrarySource"                  "LinkerPrimerSequence"          
# [43] "NumberInRun"                    "Organism"                      
# [45] "Platform"                       "ReleaseDate"                   
# [47] "samp_collect_device"            "Sample.Name"                   
# [49] "SRA.Study"                      "term_vs_preterm_delivery"      
# [51] "create_date"                    "version"                       
# [53] "race"      

###############################################################################
# 5. Profundidad de secuenciación (número total de lecturas por muestra)
###############################################################################

reads <- sample_sums(ps)

summary(reads)
#    Min. 1st Qu.  Median    Mean  3rd Qu.  Max. 
#     4  122273  153696    147014  180774  592866 

# La media y la mediana resumen la profundidad de secuenciación por muestra.
# La presencia de muestras con profundidades muy bajas o muy altas indica
# heterogeneidad en el número total de lecturas entre muestras.

###############################################################################
# 6. Distribución de los conteos totales por ASV
###############################################################################

# Número total de reads asignados a cada ASV sumando todas las muestras
asv_total_reads <- taxa_sums(ps)
# taxa_sums(ps) suma, para cada ASV, sus conteos en todas las muestras.

summary(asv_total_reads)
#      Min.   1st Qu.   Median   Mean   3rd Qu.      Max. 
#       1       8        25     41936     168    106945803 

################################################################################
# 7. Inspección taxonómica
################################################################################

# Distribución por género:
#--------------------------
table(tax$Genus, useNA = "ifany")
# Lactobacillus: 252 
# Gardnerella: 46
# Ureaplasma: 12
# NA: 1890
#...

# Distribución por especie:
#--------------------------
table(tax$Species, useNA="always")
# abscessus         acetatigenes          acidiphobus          acidophilus 
#         1                    1                    1                    2 
#aerilata          aerofaciens             aerolata            aerophila 
#       1                    2                    1                    1 
#         yunnanensis                 <NA> 
#                   1                 7513 
#...
# La presencia de asignaciones inesperadas a nivel de especie, junto con el
# elevado número de ASVs sin resolución a especie, justifica revisar de forma
# dirigida los principales taxones vaginales antes de los análisis finales.

# Taxonomía del género Lactobacillus:
#------------------------------------
lacto.tax <- tax %>% dplyr::filter(Genus == "Lactobacillus") 
nrow(lacto.tax) # [1] 252 
table(lacto.tax$Species, useNA="always") # Especies presentes
#  acidophilus       agilis coleohominis    crispatus       florum      gasseri        iners 
#           2            1            2            3            1            1            2 
#  kalixensis         oris      reuteri    rhamnosus      ruminis    vaginalis         <NA> 
#           1            2            1            1            2            1         232
sum(tax_table(ps)[,"Species"]=="crispatus", na.rm=TRUE)
# [1] 3
sum(tax_table(ps)[,"Species"]=="iners", na.rm=TRUE)
# [1] 2

# Taxonomía del género Gardnerella:
#----------------------------------
gard.tax <- tax %>% dplyr::filter(Genus == "Gardnerella")
nrow(gard.tax) # [1] 46
table(gard.tax$Species, useNA="always") # Especies presentes
#vaginalis      <NA> 
#       7        39 

# Taxonomía del género Ureaplasma:
#----------------------------------
urea.tax <- tax %>% dplyr::filter(Genus == "Ureaplasma")
nrow(urea.tax) # [1] 12
table(urea.tax$Species, useNA="always") # Especies presentes
# <NA> 
# 12 


# Lo que muestran los resultados:
# ASVs totales = 8298
# ASVs con especie asignada = 785
# ASVs sin especie = 7513

# y específicamente:
# Lactobacillus:
# acidophilus   2
# crispatus     3
# gasseri       1
# iners         2
# ...
# NA         232
# Gardnerella:
# vaginalis     7
# NA         39

# La asignación taxonómica de SILVA resulta conservadora a nivel de especie
# para varios taxones vaginales. Por ello, la comparación con Callahan y los
# refinamientos posteriores utilizan también información a nivel de secuencia
# para mejorar la resolución de Lactobacillus y Gardnerella.

###############################################################################
# 8. Guardar resumen
###############################################################################

dir.create(project_path("results"), recursive = TRUE, showWarnings = FALSE)

write.csv(
  data.frame(
    Sample=sample_names(ps),
    Reads=sample_sums(ps)),
  file=project_path("results", "Script01A_QC_Reads_per_sample.csv"),
  row.names=FALSE)

cat("\n=== FIN ===\n")
cat("sessionInfo() para la memoria:\n")
print(sessionInfo())
