# PMR Reproduce

This repository provides a tutorial on how to perform real data analysis using [PMR](https://github.com/yuanzhongshang/PMR) with summary-level gene expression data and GWAS data. 

PMR is a probabilistic Mendelian randomization (MR) method for TWAS applications. PMR relies on a MR likelihood framework that unifies many existing TWAS and MR methods, accommodates multiple correlated instruments, tests the causal effect of gene on trait in the presence of horizontal pleiotropy, and is scalable to hundreds of thousands of individuals. 

In this tutorial, we provide a real data example performing probabilistic Mendelian randomization analysis using **BMI** GWAS data from UK Biobank and gene expression data for **BACE1** from GEUVADIS. The BACE1 region is the same one used for simulations in the original paper. Users can extend the code provided in this repository to implement genome-wide PMR analysis.

## 1. Data Preparation

The analysis requires three main data components: gene expression summary statistics, GWAS summary statistics, and LD reference data. Below we describe how to obtain each of these datasets.

### Prepare Gene Annotation File (Optional)

> "In the expression data, we only focused on protein coding genes and lincRNAs that are annotated in GENCODE (release 12)."

Since gene expression data contains genome-wide SNPs, we need to extract cis-SNPs within and around the gene region for real data analysis. The original authors used gene position and functional annotations from GENCODE (release 12).

You can download the Comprehensive gene annotation from: [https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_12/gencode.v12.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_12/gencode.v12.annotation.gtf.gz)

After downloading, run the `make_gene_annotation.py` script in this repository to filter gene regions corresponding to protein coding genes and lincRNAs. 

(The output of this step is already provided in `data/gencode_v12_gene_annotation.csv`, so users can skip this step)

### Download Gene Expression Summary Statistics

For gene expression summary statistics, the original authors provided genome-wide summary-level data: [https://www.dropbox.com/scl/fo/4nqcmkblerspfmva5stwf/ANHZU_kX2AlveEEbx9DKbZU?rlkey=qjcxprlk83t7pw8ka2ne2v4w9&dl=0](https://www.dropbox.com/scl/fo/4nqcmkblerspfmva5stwf/ANHZU_kX2AlveEEbx9DKbZU?rlkey=qjcxprlk83t7pw8ka2ne2v4w9&dl=0)

Here we use gene BACE1 (ENSG00000186318.10) on chromosome 11 as an example. Users can download the `ENSG00000186318.10.tsv.gz` data from the chr11 folder in the Dropbox link.

### Download LD Reference Data

When performing PMR analysis with summary-level data, LD matrices from both gene expression data and GWAS data are required. As noted in the authors' example:

> "# LDmatrix_1 and LDmatrix_2 are often from the same reference panel data"

Since the authors recommend using LD reference from the same panel, and our GWAS data comes from UK Biobank, we use the UK Biobank LD reference for both gene expression and GWAS LD matrices in our analysis.

Here, we use the UK Biobank LD reference: [https://uchicago.app.box.com/s/jqocacd2fulskmhoqnasrknbt59x3xkn/folder/234629250877](https://uchicago.app.box.com/s/jqocacd2fulskmhoqnasrknbt59x3xkn/folder/234629250877)

Consistent with the summary statistics version, we select the **build 37** LD reference: UK Biobank LD reference/LD_matrices/b37/

For chromosome 11, users can download `11.tar.gz` from the chr11 folder in the Dropbox link and extract it.

This LD reference uses LDetect for block partitioning. The start position information for each LD block can be found in `data/LDetectregion1533.txt` in this repository (original source: [https://github.com/yuanzhongshang/GIFT/blob/main/reproduce/LDetectregion1533.txt](https://github.com/yuanzhongshang/GIFT/blob/main/reproduce/LDetectregion1533.txt)).

### Download GWAS Summary Statistics

UK Biobank GWAS summary data information is available in the online Google Docs: [https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=227859291#gid=227859291](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=227859291#gid=227859291)

First, download the variants reference file: [https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz](https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz)

Then, search for the desired traits in the "Manifest 201807" sheet and download the GWAS summary data.

For example, for the BMI trait example, download from: [https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz](https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz)

### Directory Structure Setup

After downloading the above data, please create the following directory structure and place the corresponding data files:

```
PMR_reproduce
├── data
│   ├── Gex_summary
│   │   └── ENSG00000186318.10.tsv.gz
│   ├── GWAS_summary
│   │   ├── 21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
│   │   └── variants.tsv.bgz
│   ├── LD
│   │   └── 11
│   │       ├── ukb_b37_0.1_chr11.R_snp.116383348_117747110.RDS
│   │       ├── ukb_b37_0.1_chr11.R_snp.116383348_117747110.Rvar
│   │       └── …
│   ├── input
│   └── ...
├── make_gene_annotation.py
├── make_input.R
└── run_PMR_summary.R
```

## 2. Data Filtering and Alignment

Users can refer to the `make_input.R` script provided in this repository to filter regions for PMR analysis and align gene expression summary statistics, GWAS summary statistics, and LD reference data to generate input data for PMR.

If the above file directory structure is ready, you can directly run the entire `make_input.R` script, and you will obtain 3 input data files in the `data/input` folder.

The `make_input.R` script includes the following processes:

1. Loading cis-eQTL summary statistics for the target gene and filtering SNPs within ±100 kb of the gene interval on the specified chromosome.
2. Reading LD region definitions, selecting the best-matching LD block, and loading both the variant list and full LD matrix for that block.
3. Reading GWAS variant information and summary statistics, filtering to single-nucleotide variants, merging on variant ID, and subsetting to the same chromosome.
4. Identifying common SNPs across eQTL data, LD variant list, and GWAS summary statistics. Extracting the corresponding submatrix from the LD matrix and subsetting all three datasets to these shared SNPs.
5. Writing out the filtered datasets: eQTL summary, LD matrix, and GWAS summary statistics data for PMR analyses.

## 3. Running PMR Analysis

It is easy to install the development version of PMR package using the 'devtools' package:

```r
# install.packages("devtools")
library(devtools)
install_github("yuanzhongshang/PMR")
```

Users can refer to the `run_PMR_summary.R` script provided in this repository to perform PMR analysis on the files in `data/input`.

If you have pre-downloaded the files in the `data/input/` directory from this repository, you can run the script directly.

**Important Note**: When running the PMR program, you may encounter issues with Cholesky decomposition. This problem has been reported in the original repository's issues ([https://github.com/yuanzhongshang/PMR/issues/6](https://github.com/yuanzhongshang/PMR/issues/6)).

In this example, we address this by setting the diagonal elements of the LD matrix to 1.5 to ensure PMR runs successfully. Users can refer to the author's suggestions in the issue for other adjustments.

## References

Yuan, Z., Zhu, H., Zeng, P. et al. Testing and controlling for horizontal pleiotropy with probabilistic Mendelian randomization in transcriptome-wide association studies. Nat Commun 11, 3861 (2020). [https://doi.org/10.1038/s41467-020-17668-6](https://doi.org/10.1038/s41467-020-17668-6)
