# Gut-Microbiota-in-correlation-with-the-risk-of-HNC
Gut Microbiota in correlation with the risk of HNC

# 1.Raw input data list
## The GWAS summary data of HNC patients and  Gut Microbiota data

| Variable | Consortium or study | Sample size | Year | Cohort | Nation of cohort | Number of samples | Data link |
|----------|---------------------|-------------|------|--------|------------------|-------------------|-----------|
| Gut microbiota | MiBioGen | 18,340 | 2021 | BSPSPC | Germany | 18,340 | https://mibiogen.gcc.rug.nl/menu/main/home |
| Head and neck cancer | UK Biobank | 373,122 | 2021 | ieu-b-4912 | European | 1106 cases / 372016 controls | https://opengwas.io/datasets/ieu-b-4912 |

# 2.Prerequisites

## R dependencies

### R version 4.4.3 

Users running Platform: x86_64-pc-linux-gnu and Running under: Ubuntu 20.04.6 LTS, to install the latest version of TwoSampleMR_0.5.8

```r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")

```
