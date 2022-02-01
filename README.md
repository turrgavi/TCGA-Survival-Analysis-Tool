# TCGA Survival Analysis Tool

This tool is designed to perform high throughput survival analyses on standard format TCGA mRNA expression data. 

## What goes in: 
- TCGA mRNA Expression data
- TCGA Phenotype Data
- List of genes you want to analyse

## What comes out:
- List of genes with significant survival analysis P-values (with their P-values and whether the upper of lower quartile of expression is associated with decreased survival)
- Survival plots for specified genes

## What you need:
- RStudio
- R
- *You can also use Microsoft Open R as an alternative to base R if you want to multithread. I recommend this if you are using Windows or Linux*

## Required Packages:
- stringR
- dplyr
- survival
- ggplot2
- survminer
- cowplot

# How to use:
1. Download required files from XenaBrowser.
2. Extract them and convert from TSV to CSV using Excel.
3. Convert ENSEMBL IDs to Gene symbols using online tool.
4. Compile a list of genes you are interested in checking and save as a CSV.
5. Set the working directory in the KMSA Script to the directory containing your CSV files.
6. Set the KMSA expression data name to the name of your mRNA expression data file (including the .csv).
7. Run the script!

# What does it do:
1. Loads the expression data, phenotype data and gene list.
2. Matches all genes on the gene list to those present on the TCGA expression data using Gene symbols.
3. Pairs matching genes with phenotype data by using the TCGA Barcodes.
4. Calculates and assigns upper and lower quartiles to gene expression.
5. Performs Kaplan-Meier Survival Analysis on all matched genes using the phenotype data to determine status and event variables. These are seperated into upper and lower quartile expression groups.
6. Writes all P-values (or only significant ones if you choose) to a CSV file alongside the Gene symbol and which quartile was associated with decreased survival to 5 years.

# What can you change to suit your needs:
- TCGA dataset
- Survival cutoff (default is 1825 days/5 years)
- Genes to analyse
- Write all or only singificant genes to output file
- Comparison groups (default is upper and lower expression quartiles)
- What sample type to analyse (default is Primary tumours only)

# Planned improvements:
1. Automatically convert ENSEMBL ID to Gene symbols
2. Automatically fetch TCGA data sets
3. Automatically convert TSV to CSV
4. Generate a report on completion of analysis
5. Create an option to generate and save survival plots

# About Me:
My name is Gavin Turrell and I am currently in my honours year. Currently working with liquid chromatography and mass spectroscopy. I have a really keen interest in cancer immunology and wrote this code myself as there was a need for high throughput KMSA that was more plug and play.
