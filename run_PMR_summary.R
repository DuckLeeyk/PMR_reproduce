# library(devtools)
# install_github("yuanzhongshang/PMR")

library(PMR)

Gex_file <- "data/input/Gex_summary_ENSG00000186318.10.tsv"
GWAS_file <- "data/input/GWAS_summary_ENSG00000186318.10.tsv"
LD_file <- "data/input/LD_matrix_ENSG00000186318.10.txt"

# load the Zscore vector for the cis-SNP effect size from the eQTL data
Gex_zscore <- read.table(Gex_file, header = TRUE, stringsAsFactors = FALSE)
Zscore_1 <- as.vector(Gex_zscore$T)

# load the Zscore vector for the cis-SNP effect size from the GWAS data
GWAS_zscore <- read.table(GWAS_file, header = TRUE, stringsAsFactors = FALSE)
Zscore_2 <- as.vector(GWAS_zscore$tstat)

# load the LD matrix for the cis-SNPs in the eQTL data
LD_matrix <- read.table(LD_file)

# Please refer to https://github.com/yuanzhongshang/PMR/issues/6 about the problem of cholesky decomposition!
diag(LD_matrix) <- 1.5

# LDmatrix_1 and LDmatrix_2 are often from the same reference panel data
LDmatrix_1 <- as.matrix(LD_matrix)
LDmatrix_2 <- as.matrix(LD_matrix)

# load the sample size n1 from eQTL data and n2 from GWAS data
n1 <- mean(Gex_zscore$NMISS)
n2 <- mean(GWAS_zscore$n_complete_samples)

# run PMR with the summary data
result <- PMR_summary_Egger(Zscore_1, Zscore_2, LDmatrix_1, LDmatrix_2, samplen1 = n1, samplen2 = n2, lambda = 0, max_iterin = 1000, epsin = 1e-5, Heritability_geneexpression_threshold = 0)

print(result)
