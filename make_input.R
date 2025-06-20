# =============================================================================
# Description:
#   This script performs locus-specific data preparation by:
#     1. Loading cis-eQTL summary statistics for that gene and filtering SNPs
#        within ±100 kb of the gene interval on a specified chromosome.
#     2. Reading LD region definitions, selecting the best-matching LD block,
#        and loading both the variant list and full LD matrix for that block.
#     3. Reading GWAS variant information and summary statistics, filtering to
#        single-nucleotide variants, merging on the variant ID, and subsetting
#        to the same chromosome.
#     4. Identifying SNPs common to the eQTL data, LD variant list, and GWAS
#        summary statistics. Extracting the corresponding submatrix from the
#        LD matrix and subsetting all three datasets to these shared SNPs.
#     5. Writing out the filtered eQTL summary, LD matrix, and GWAS summary
#        statistics data for PMR analyses.
# =============================================================================

library(readr)
library(dplyr)

##############################################
# Section 1: Read and filter cis-SNP eQTL data
##############################################

# Define the target chromosome and path to the gene's eQTL summary file
Chr <- 11
Gex_file <- "data/Gex_summary/ENSG00000186318.10.tsv.gz"

# Read gene annotation file and extract start/end positions for the target gene
anno <- read_csv("data/Gex_summary/gencode_v12_gene_annotation.csv")

# Derive the ENSG gene ID from the eQTL filename using a regular expression
gene_ENSG_id <- sub(".*/(ENSG[^\\.]+\\.[^\\.]+)\\.tsv.gz$", "\\1", Gex_file)

# Select the annotation row for the target gene
target_row <- anno %>%
  filter(gene_id == gene_ENSG_id) %>%
  select(start, end)

# Extract numeric start and end positions
start_pos <- target_row$start
end_pos <- target_row$end

# Read the gene's eQTL data (TSV) and drop unneeded columns: BETA, SE, R2, P
Gex <- read_tsv(Gex_file) %>%
  select(-BETA, -SE, -R2, -P)

# Filter SNPs on the specified chromosome within ±100 kb of the gene interval
Gex <- Gex %>%
  filter(
    CHR == Chr,
    BP >= (start_pos - 100000),
    BP <= (end_pos + 100000)
  )

# Clean up intermediate objects and free memory
rm(anno, target_row)
gc()


########################################
# Section 2: Process LD region data
########################################

# Read LD region definitions
LD_region <- read.table(
  "data/LDetectregion1533.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Keep only the regions on the target chromosome
LD_region <- subset(LD_region, chr == paste0("chr", Chr))

# Find regions that fully contain the target ±100 kb interval
contain_idx <- which(
  LD_region$start <= (start_pos - 100000) &
    LD_region$stop >= (end_pos + 100000)
)

if (length(contain_idx) == 1) {
  # Exactly one region contains the interval
  selected <- LD_region[contain_idx, ]
} else if (length(contain_idx) > 1) {
  # Multiple regions found: warn and pick the first
  warning("Multiple LD regions fully contain the target interval; using the first match.")
  selected <- LD_region[contain_idx[1], ]
} else {
  # No region fully contains the interval: compute overlap lengths
  LD_region$overlap <- pmax(
    0,
    pmin(LD_region$stop, end_pos + 100000) -
      pmax(LD_region$start, start_pos - 100000)
  )
  # Select the region with the largest overlap (if ties, first one)
  best_idx <- which.max(LD_region$overlap)
  selected <- LD_region[best_idx, ]
  # (Optionally drop the overlap column if no longer needed)
  # LD_region$overlap <- NULL
}

# Construct filenames for the LD variant list and LD matrix based on the selected region
LD_Rvar_file <- paste0(
  "ukb_b37_0.1_chr", Chr,
  ".R_snp.", selected$start, "_", selected$stop, ".Rvar"
)
LD_RDS_file <- paste0(
  "ukb_b37_0.1_chr", Chr,
  ".R_snp.", selected$start, "_", selected$stop, ".RDS"
)

# Read the list of SNP IDs in the region and the full LD matrix
LD_Rvar <- read.table(
  file = file.path("data/LD", Chr, LD_Rvar_file),
  header = TRUE,
  stringsAsFactors = FALSE
)
LD_RDS <- readRDS(
  file = file.path("data/LD", Chr, LD_RDS_file)
)


###################################################
# Section 3: Process GWAS summary statistics data
###################################################

# Define file paths for GWAS variant info and summary statistics
GWAS_variants_file <- "data/GWAS_summary/variants.tsv.bgz"
GWAS_sumstats_file <- "data/GWAS_summary/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz"

# Read the GWAS variant reference file
variants_info <- read_tsv(
  file      = GWAS_variants_file,
  col_types = cols(), # let readr infer column types
  progress  = FALSE
)

# Read the GWAS summary statistics file
sumstats_data <- read_tsv(
  file      = GWAS_sumstats_file,
  col_types = cols(),
  progress  = FALSE
)

# Retain only necessary columns in each dataset
variants_info <- variants_info %>%
  select(variant, chr, pos, ref, alt, rsid)

sumstats_data <- sumstats_data %>%
  select(variant, tstat, n_complete_samples)

# Filter to SNPs with single-nucleotide alleles (exclude indels)
variants_info <- variants_info %>%
  filter(
    nchar(ref) == 1,
    nchar(alt) == 1
  )

# Merge summary stats with variant annotation by 'variant' identifier
merged_sumstats <- inner_join(sumstats_data, variants_info, by = "variant")

# Keep only SNPs on the target chromosome
merged_sumstats <- merged_sumstats %>%
  filter(chr == Chr)

# Remove duplicate rsids, retaining the first occurrence
merged_sumstats <- merged_sumstats %>%
  distinct(rsid, .keep_all = TRUE)

# Clean up intermediate objects and free memory
rm(variants_info, sumstats_data)
gc()


##########################################################
# Section 4: Identify common SNPs and assemble final data
##########################################################

# Determine SNPs common to eQTL data, LD variant list, and GWAS summary data
common_snps <- Reduce(
  intersect,
  list(Gex$SNP, LD_Rvar$id, merged_sumstats$rsid)
)

# Subset each dataset to the common SNPs
Gex_common <- Gex %>%
  filter(SNP %in% common_snps)

merged_sumstats_common <- merged_sumstats %>%
  filter(rsid %in% common_snps)

# Identify indices of common SNPs in the LD variant list
LD_common_indices <- which(LD_Rvar$id %in% common_snps)

# Subset the LD matrix to those common SNPs
LD_matrix <- LD_RDS[LD_common_indices, LD_common_indices]


###############################################
# Section 5: Write output for PMR analysis
###############################################

# Define output directory for input data
input_dir <- "data/input"

# Write the filtered eQTL summary to a TSV file
write_tsv(
  Gex_common,
  file.path(input_dir, paste0("Gex_summary_", gene_ENSG_id, ".tsv"))
)

# Write the LD matrix subset to a tab-delimited file
write.table(
  LD_matrix,
  file.path(input_dir, paste0("LD_matrix_", gene_ENSG_id, ".txt")),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# Write the filtered GWAS summary statistics to a TSV file
write_tsv(
  merged_sumstats_common,
  file.path(input_dir, paste0("GWAS_summary_", gene_ENSG_id, ".tsv"))
)
