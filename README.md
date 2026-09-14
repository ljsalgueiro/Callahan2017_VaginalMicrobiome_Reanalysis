# Vaginal Microbiome and Preterm Birth — Reanalysis of Callahan et al. (2017)

## Overview

This repository contains the bioinformatic and statistical workflow developed for a Master's Thesis investigating the vaginal microbiome during pregnancy and its association with preterm birth.

The project reanalyzes the publicly available 16S rRNA gene sequencing dataset described by Callahan et al. (2017):

> Callahan BJ, DiGiulio DB, Goltsman DSA, et al. Replication and refinement of a vaginal microbial signature of preterm birth in two racially distinct cohorts of US women. *Proceedings of the National Academy of Sciences*. 2017;114(37):9966–9971. doi:10.1073/pnas.1705899114.

The dataset contains longitudinal vaginal samples collected during pregnancy from two cohorts: Stanford and the University of Alabama at Birmingham (UAB), including women with term and preterm deliveries.

The workflow includes sequence processing and quality control, validation against the original study, contaminant removal, targeted taxonomic refinement, confounder assessment, alpha and beta diversity, longitudinal community stability, differential abundance, and Community State Type (CST) analyses using VALENCIA.

---

## Dataset

The study uses publicly available paired-end 16S rRNA gene sequencing data targeting the V4 region.

The final longitudinal dataset comprised:

| Cohort    |   Women |   Term | Preterm |  Samples |
| --------- | ------: | -----: | ------: | -------: |
| Stanford  |      39 |     30 |       9 |      897 |
| UAB       |      96 |     55 |      41 |     1282 |
| **Total** | **135** | **85** |  **50** | **2179** |

Multiple longitudinal samples were available for each participant.

Raw sequencing data are not redistributed through this repository and can be obtained from the public NCBI Sequence Read Archive (SRA) data associated with the original study.

---

## Bioinformatic processing

Because of the size of the sequencing dataset, raw-read processing was performed in **Galaxy** using a DADA2-based workflow.

The exported Galaxy workflow is available at:

`galaxy/Galaxy_DADA2_Workflow.ga`

It can be imported into a compatible Galaxy instance to inspect the sequence-processing steps, tools and parameters.

Downstream processing, quality control, taxonomic refinement and statistical analyses were performed in **R**.

---

## Repository structure

```text
.
├── README.md
├── LICENSE
├── .gitignore
│
├── galaxy/
│   └── Galaxy_DADA2_Workflow.ga
│
└── scripts/
    ├── Script01A-Load_and_Inspection.R
    ├── Script01B-Phyloseq_Reconstruction.R
    ├── Script02-Callahan_Validation.R
    ├── Script03-Contaminant_Analysis.R
    ├── Script04-Taxonomic_Filtering.R
    ├── Script05A-Lactobacillus_BLAST.R
    ├── Script05B-Lactobacillus_BLAST_Assignment_vSpeciate_Prepare.R
    ├── Script05C-Lactobacillus_Final_Refinement.R
    ├── Script06A-Gardnerella_Callahan_Exploration.R
    ├── Script06B-Gardnerella_Refined.R
    ├── Script07A-Prevotella_BLAST.R
    ├── Script07B-Prevotella_BLAST_Assignment_vSpeciate_Prepare.R
    ├── Script07C-Prevotella_Final_Refinement.R
    ├── Script08-Confounder_Variable_Audit.R
    ├── Script09-AlphaDiversity.R
    ├── Script10-BetaDiversity.R
    ├── Script11-TemporalStability.R
    ├── Script12-ANCOMBC2.R
    ├── Script13-VALENCIA_CST_Preparation.R
    └── Script14-VALENCIA_CST_Analysis.R
```

---

# Analysis workflow

The R scripts are numbered according to the analytical workflow.

### Script 01A — Load and initial inspection

`Script01A-Load_and_Inspection.R`

Loads the initial phyloseq object generated from the Galaxy/DADA2 output and performs quality-control checks of its structure, sample metadata, sequencing depth, ASV abundances, taxonomic assignments and reference sequences.

### Script 01B — Phyloseq reconstruction

`Script01B-Phyloseq_Reconstruction.R`

Reconstructs the phyloseq object using the ASV sequences and associated abundance and metadata information. The script verifies the correspondence between ASV sequences, abundance data, taxonomy and sample information to produce a consistent object for downstream analyses.

### Script 02 — Validation against Callahan et al. (2017)

`Script02-Callahan_Validation.R`

Compares the reprocessed dataset with the data and analysis objects associated with the original Callahan et al. study. The validation examines sample representation, ASV sequences, abundances and taxonomic information to assess the consistency of the independently reconstructed dataset with the original analysis.

### Script 03 — Contaminant analysis

`Script03-Contaminant_Analysis.R`

Identifies potential contaminants using the prevalence method implemented in `decontam`, defining laboratory extraction controls through `env_biome`, sequencing run (`IL_Run`) as the batch variable, and Fisher's method to combine evidence across batches. The primary analysis uses the conservative default threshold of 0.1, with a threshold of 0.5 evaluated as a sensitivity analysis. ASVs absent from all vaginal samples are subsequently removed together with the extraction controls.

### Script 04 — Taxonomic filtering

`Script04-Taxonomic_Filtering.R`

Performs post-denoising taxonomic quality control by removing ASVs explicitly assigned to non-target groups such as Eukaryota, Archaea, mitochondria and chloroplasts. Incompletely classified bacterial ASVs are retained unless there is explicit taxonomic evidence supporting their exclusion.

### Script 05A — *Lactobacillus* BLASTn analysis

`Script05A-Lactobacillus_BLAST.R`

Dynamically extracts *Lactobacillus* ASVs, exports their sequences and submits them to remote NCBI BLASTn searches. The script records the BLAST RID and search metadata, retrieves information on the database build and BLAST version, downloads the original XML2 and tabular results, calculates query coverage and generates a master table containing the sequence-level evidence used for subsequent species assignment.

### Script 05B — *Lactobacillus* BLAST assignment and vSpeciateDB preparation

`Script05B-Lactobacillus_BLAST_Assignment_vSpeciate_Prepare.R`

Filters *Lactobacillus* BLAST hits using minimum query coverage and sequence identity criteria and constructs a Best Hit Set using a hierarchical comparison of coverage, identity and E-value. Conservative rules are then applied to derive BLAST-based species assignments while preserving the original taxonomy. The script also evaluates sequence orientation and generates the FASTA input required for SpeciateIT/vSpeciateDB.

### Script 05C — Final *Lactobacillus* taxonomic refinement

`Script05C-Lactobacillus_Final_Refinement.R`

Integrates SpeciateIT/vSpeciateDB V4 classifications with the independent BLASTn results and applies conservative rules to obtain the final *Lactobacillus* assignments. Original SILVA taxonomy, BLAST-based classifications and vSpeciateDB classifications are preserved as separate fields, providing traceability for the final refined taxonomy.

### Script 06A — Exploration of *Gardnerella* variants from Callahan et al.

`Script06A-Gardnerella_Callahan_Exploration.R`

Recovers the *Gardnerella* variants used in the original Callahan et al. analysis from the published analysis object and extracts the reference sequences corresponding to G1, G2 and G3. This is an exploratory/reference-recovery step and does not modify the phyloseq object used in the TFM.

### Script 06B — *Gardnerella* taxonomic refinement

`Script06B-Gardnerella_Refined.R`

Compares *Gardnerella* ASVs from the reprocessed dataset with the G1, G2 and G3 reference sequences recovered from the original Callahan analysis. Exact sequence compatibility is used to classify ASVs as G1, G2, G3 or `G_other`, and the resulting classification is incorporated into a new phyloseq object without replacing the underlying taxonomic information.

### Script 07A — *Prevotella* BLASTn analysis

`Script07A-Prevotella_BLAST.R`

Dynamically extracts *Prevotella* ASVs and performs remote NCBI BLASTn searches. As in the *Lactobacillus* workflow, the script records search and database metadata, preserves the original BLAST output, calculates query coverage and produces a master table of sequence matches for reproducible downstream taxonomic assignment.

### Script 07B — *Prevotella* BLAST assignment and vSpeciateDB preparation

`Script07B-Prevotella_BLAST_Assignment_vSpeciate_Prepare.R`

Applies quality filters to the *Prevotella* BLAST results and constructs the Best Hit Set using query coverage, percent identity and E-value. It identifies true taxonomic ties, derives conservative BLAST-based assignments, incorporates them into the phyloseq object and prepares correctly oriented ASV sequences for independent classification with SpeciateIT/vSpeciateDB.

### Script 07C — Final *Prevotella* taxonomic refinement

`Script07C-Prevotella_Final_Refinement.R`

Integrates SpeciateIT/vSpeciateDB classifications with the BLASTn assignments and applies conservative decision rules to generate the final *Prevotella* taxonomy. BLAST-based, vSpeciateDB and final refined assignments are retained separately to preserve the evidence underlying each taxonomic decision.

### Script 08 — Audit of potential confounding variables

`Script08-Confounder_Variable_Audit.R`

Systematically audits the available clinical and study metadata to identify variables that could potentially act as confounders. The script evaluates variable type, missingness, number of levels, within-participant consistency, cohort and term/preterm distributions, preliminary associations with the outcome, correlations and potential redundancy. It provides an evidence base for covariate selection rather than automatically selecting the final adjustment variables.

### Script 09 — Alpha diversity

`Script09-AlphaDiversity.R`

Evaluates Shannon and Gini-Simpson diversity between term and preterm pregnancies within a common 15–33 gestational-week window while retaining repeated samples. Mixed-effects models include participant-specific random intercepts and adjustment for gestational week and sequencing depth, with global and cohort-specific analyses. Sensitivity analyses include maternal-age adjustment and participant-level averaging, and simulation-based minimum detectable effect analyses are performed for Shannon diversity.

### Script 10 — Beta diversity

`Script10-BetaDiversity.R`

Evaluates differences in overall vaginal microbial composition using Bray-Curtis dissimilarity and PERMANOVA. The primary analysis uses participant-level mean relative-abundance profiles to avoid pseudoreplication, with global additive and interaction models and separate Stanford and UAB analyses. The workflow also evaluates multivariate dispersion, maternal-age sensitivity, a longitudinal sample-level sensitivity analysis with participant-level permutations, and simulation-based minimum detectable effect analyses expressed as R².

### Script 11 — Temporal stability

`Script11-TemporalStability.R`

Evaluates whether vaginal microbial composition changes more between consecutive samples in women with preterm than term delivery. Bray-Curtis dissimilarity between consecutive within-participant samples is modeled using mixed-effects models with adjustment for the interval between samples and gestational timing, participant-specific random intercepts and a variance structure for residual heteroscedasticity. The primary 15–33-week analysis is complemented by full-pregnancy and minimum-pair sensitivity analyses.

### Script 12 — Differential abundance with ANCOM-BC2

`Script12-ANCOMBC2.R`

Tests for differentially abundant taxa between term and preterm pregnancies at genus and species levels using ANCOM-BC2. The primary analysis is restricted to 15–33 gestational weeks and incorporates gestational age, cohort where applicable, and participant-specific random effects to account for repeated measurements. Taxa are filtered using a 5% prevalence threshold, analyses are performed globally and separately for Stanford and UAB, and maternal-age sensitivity and effect-size precision assessments are included.

### Script 13 — VALENCIA CST preparation and integration

`Script13-VALENCIA_CST_Preparation.R`

Transforms the final refined phyloseq object into the taxonomic abundance format required by VALENCIA and audits its compatibility with the official VALENCIA centroid reference. Taxonomic labels are harmonized hierarchically while preserving the refined taxonomy, with specific handling of *Gardnerella* variants and other relevant vaginal genera. The script generates the input for external execution of `Valencia.py` and, once the output is available, integrates CST, sub-CST and similarity scores back into a phyloseq object.

### Script 14 — VALENCIA CST distribution and temporal dynamics

`Script14-VALENCIA_CST_Analysis.R`

Evaluates CST distribution and longitudinal CST transitions in relation to term/preterm delivery. The primary analysis uses the common 15–33-week window and mixed logistic one-vs-rest models for CSTs with sufficient prevalence, with global and Stanford/UAB-specific analyses and FDR correction. Temporal dynamics are assessed from consecutive within-participant sample pairs. Sensitivity analyses include the full pregnancy and restriction to VALENCIA assignments with similarity score ≥0.1, while unstable or non-convergent models are retained for audit but excluded from inference.

---

## Statistical analysis strategy

The main ecological and statistical analyses use a common gestational window of **15–33 weeks** to improve comparability between participants and cohorts.

Because the dataset contains repeated longitudinal observations from the same women, the analytical strategy explicitly accounts for within-participant dependence through mixed-effects models, participant-level aggregation or participant-level permutation schemes, depending on the analysis.

The principal analyses include:

* alpha diversity using Shannon and Gini-Simpson indices;
* beta diversity using Bray-Curtis dissimilarity and PERMANOVA;
* longitudinal within-participant microbial stability;
* differential abundance using ANCOM-BC2;
* VALENCIA Community State Type classification;
* CST distribution and transition analyses;
* sensitivity analyses for relevant covariates and analytical assumptions;
* power or minimum detectable effect analyses where applicable.

---

## Targeted taxonomic refinement

Standard reference-database classification provided limited species-level resolution for several biologically important vaginal taxa. Targeted sequence-based refinement was therefore performed for:

* *Lactobacillus*
* *Gardnerella*
* *Prevotella*

For *Lactobacillus* and *Prevotella*, independent BLASTn and SpeciateIT/vSpeciateDB evidence was integrated using conservative decision rules.

For *Gardnerella*, ASVs were compared with the G1, G2 and G3 sequence variants used in the original Callahan et al. analysis.

The original taxonomy and the results of the individual refinement methods are retained separately to preserve traceability.

---

## Reproducibility

The repository contains the code required to document and reproduce the analytical workflow, but does not redistribute large sequencing files, intermediate phyloseq objects, BLAST result archives or external reference databases.

The exported Galaxy workflow documents raw-read processing, while the numbered R scripts document the subsequent analytical workflow.

Remote NCBI BLAST scripts record information about the searches and database versions used. Individual scripts also contain their required inputs, generated outputs, parameters and software dependencies.

External resources such as the original sequencing dataset, vSpeciateDB and VALENCIA reference files must be obtained separately from their corresponding sources.

---

## Main software and packages

The workflow uses software and R packages including:

* Galaxy
* DADA2
* R
* phyloseq
* decontam
* vegan
* nlme
* lme4
* ANCOM-BC2
* NCBI BLASTn
* SpeciateIT / vSpeciateDB
* VALENCIA
* Biostrings
* tidyverse packages

Additional dependencies are documented within the individual scripts.

---

## Data availability

The sequencing data analyzed in this project are publicly available through the NCBI Sequence Read Archive and originate from the dataset analyzed by Callahan et al. (2017).

Raw sequencing data are not redistributed through this repository.

---

## References

**Callahan BJ, DiGiulio DB, Goltsman DSA, et al.** Replication and refinement of a vaginal microbial signature of preterm birth in two racially distinct cohorts of US women. *Proceedings of the National Academy of Sciences of the United States of America*. 2017;114(37):9966–9971. doi:10.1073/pnas.1705899114.

**Holm JB, Gajer P, Ravel J.** SpeciateIT and vSpeciateDB: novel, fast, and accurate per sequence 16S rRNA gene taxonomic classification of vaginal microbiota. *BMC Bioinformatics*. 2024;25:313. doi:10.1186/s12859-024-05930-3.

---

## Author

**Lucía Salgueiro**

Master's Thesis
Master's Degree in Bioinformatics for Health Sciences
Universidade da Coruña
