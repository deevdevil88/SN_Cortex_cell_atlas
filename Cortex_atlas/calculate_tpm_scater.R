library(scater)

library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)

cortex_raw <- read.table(file="Cortex_10x_raw_counts_gt_1_protein_codiing.txt", sep="\t", header=T)

cortex_metadata <- read.table(file="Cortex_10x_CCA_metadata.txt", sep="\t", header=T)

umi <- SingleCellExperiment(
  assays = list(counts = as.matrix(cortex_raw)), 
  colData = cortex_metadata
)

keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]

#tpm(umi) <- (calculateTPM(umi, effective_length = NULL))
tpm_cortex <- calculateTPM(umi, effective_length = NULL)
tpm_cortex <- as.data.frame(tpm_cortex)
assay(umi)
write.table(tpm_cortex, file="Cortex_10x_tpm_protein_coding.txt", sep="\t",quote=F)
