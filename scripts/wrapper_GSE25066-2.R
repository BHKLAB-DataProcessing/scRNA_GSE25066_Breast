library(MultiAssayExperiment)
library(SummarizedExperiment)
library(dplyr)
library(qs)
library(biomaRt)
library(purrr)
library(GEOquery)

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[[1]]
filename <- args[[2]]

data <- getGEO("GSE25066", GSEMatrix = TRUE)[[1]]
assay_data <- makeSummarizedExperimentFromExpressionSet(data)
coldata <- data@phenoData@data

sample_ids <- strsplit(coldata$characteristics_ch1, ": ") %>% map_chr(`[`, 2)
coldata$characteristics_ch1 <- sample_ids
drugs <- strsplit(coldata$characteristics_ch1.14, ": ") %>% map_chr(`[`, 2)
coldata$characteristics_ch1.14 <- drugs
coldata <- coldata[coldata$characteristics_ch1.14 == "Taxol", ]
# coldata <- coldata[c("characteristics_ch1","characteristics_ch1.14","characteristics_ch1.10")]
coldata$characteristics_ch1.14 <- "Paclitaxel"
coldata$characteristics_ch1.10 <- strsplit(coldata$characteristics_ch1.10, ": ") %>% map_chr(`[`, 2)
rownames(coldata) <- coldata$line
coldata <- coldata[!coldata$characteristics_ch1.10 == "NA", ]

new_coldata <- data.frame(patientid = coldata$characteristics_ch1)
new_coldata$treatmentid <- "Paclitaxel"

new_coldata$response <- coldata$characteristics_ch1.10
new_coldata[new_coldata$response == "RD", ]$response <- "NR"
new_coldata[new_coldata$response == "pCR", ]$response <- "R"
new_coldata$tissueid <- "Breast"
new_coldata$`survival_time_pfs/survival_time_os` <- strsplit(coldata$characteristics_ch1.13, ": ") %>% map_chr(`[`, 2)
new_coldata$survival_unit <- "years"
new_coldata$`event_occurred_pfs/event_occurred_os` <- strsplit(coldata$characteristics_ch1.12, ": ") %>% map_chr(`[`, 2)
new_coldata$age <- coldata$characteristics_ch1.2
new_coldata$stage <- coldata$characteristics_ch1.6

new_experiment <- MultiAssayExperiment(colData = new_coldata)

assays <- as.data.frame(assay(assay_data))
genes <- data@featureData@data
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
affyids <- genes$ID
id_version <- "affy_hg_u133a_2"

mapping <- getBM(
  attributes = c(
    id_version,
    "hgnc_symbol", "start_position", "end_position", "ensembl_gene_id_version"
  ), filters = id_version,
  values = affyids, mart = ensembl
)
mapping <- mapping[mapping$ensembl_gene_id_version != "", ]
assays <- assays[as.character(mapping[[id_version]]), ]
mapping <- mapping[!duplicated(mapping$ensembl_gene_id_version), ]
assays <- assays[as.character(mapping[[id_version]]), ]
rownames(assays) <- mapping$ensembl_gene_id_version
genes <- genes[!duplicated(genes$ID), ]
rownames(genes) <- genes$ID
genes <- genes[mapping$affy_hg_u133a_2, ]
rownames(genes) <- mapping$ensembl_gene_id_version
genes$ID <- mapping$ensembl_gene_id_version
colnames(assays) <- coldata$characteristics_ch1
assays <- assays[, new_coldata$patientid]
rownames(new_coldata) <- new_coldata$patientid
counts <- SummarizedExperiment(assays = list(assays), colData = new_coldata, rowData = genes)

assays <- assays[mapping$ensembl_gene_id_version, ]
rownames(assays) <- mapping$ensembl_gene_id_version
length <- mapping$end_position - mapping$start_position
x <- assays / length
tpm <- t(t(x) * 1e6 / colSums(x))
tpm <- as.data.frame(tpm)
tpm_counts <- SummarizedExperiment(assays = list(tpm), colData = new_coldata, rowData = genes)
experiment_list <- list(expr_gene_counts = counts, expr_gene_tpm = tpm)
experiment_list <- ExperimentList(experiment_list)
new_experiment@ExperimentList <- experiment_list

saveRDS(new_experiment, paste0(work_dir, filename))
