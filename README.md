# Functional Profiling of Gut Microbiome of Children in Atlanta, USA and Maputo, Mozambique Using HUMAnN3

**Author:** Nina Sherman

**Institution:** Johns Hopkins University

**Course:** Practical Introduction to Metagenomics (Prof. Joshua Orvis)

**Date:** May 11, 2026

## Project Overview

This repository contains the complete bioinformatics pipeline and downstream analysis script used to compare the functional metabolic potential of gut microbiomes from young children (ages 2 and under) across distinct geographic and economic regions. Using publicly available shotgun metagenomic stool sequencing data, this project evaluates how differences in geography, diet, and environmental exposures influence the metabolic capabilities of the human gut.

### Key Biological Findings

#### Shared Core Functions: 
Essential pathways for cell survival—including L-rhamnose biosynthesis, L-valine biosynthesis, glycolysis IV, guanosine ribonucleotide biosynthesis, and L-isoleucine biosynthesis—were highly prevalent across all six samples.

#### Maputo Cohort Enrichment: 
Samples from Maputo, Mozambique showed a strong enrichment in biosynthetic potential. Key enriched pathways include queuosine biosynthesis (linked to bacterial stress response and translational accuracy), guanosine ribonucleotide de novo biosynthesis, and L-methionine biosynthesis. This aligns with an abundance of Prevotella species typically found in non-Western, high-fiber diets.

#### Atlanta Cohort Enrichment: 
Samples from Atlanta, USA demonstrated a relative enrichment in glycogen biosynthesis, suggesting variations in carbohydrate storage and energy metabolism linked to Western dietary patterns.

#### Statistical Power Note: 
Due to the small cohort size (\(n=3\) per group), no individual metabolic pathways achieved strict alpha-level statistical significance via the Mann–Whitney U test. The results instead serve as a clear indicator of functional trends.

## Dataset Metadata
All gut metagenomes were obtained from the NCBI Sequence Read Archive (SRA) under **BioProject: PRJNA747761**

Atlanta Samples: `SRR15202440`, `SRR15202441`, `SRR15202442`

Maputo Samples: `SRR15209203`, `SRR15209176`, `SRR15209171`

## Toolstack & Requirements

#### Bioinformatics Tools

**SRA Toolkit** — Raw sequence downloading (`prefetch`, `fasterq-dump`).

**fastp (v1.3.3)** — All-in-one FASTQ quality filtering and adapter trimming.

**MetaPhlAn 4** — Marker-gene-based taxonomic profiling using the `mpa_vJun23_CHOCOPhlAnSGB_202307` database.

**HUMAnN 3** — Metabolic gene family and pathway abundance reconstruction.

**Bowtie2 & DIAMOND** — Internal sequence alignment mapping utilities.
#### Downstream Python Packages

• `pandas` 

• `matplotlib` 

• `seaborn` 

• `scipy`

## Installation & Environment Setup

Replicate the exact runtime workspace by setting up this dedicated Conda environment:
```
# Create environment with Python 3.9
conda create -n humann3_env python=3.9 -y
conda activate humann3_env

# Install upstream command-line bioinformatic tools
conda install -c bioconda humann metaphlan fastp bowtie2 diamond -y

# Install downstream Python packages for statistics and data visualization
conda install -c conda-forge pandas matplotlib seaborn scipy -y
```
## Pipeline Execution Steps

#### 1. Data Retrieval

Download the raw paired-end sequence data from SRA. (Note: This step must be executed for all 6 target SRR accessions listed in the Metadata section).
```
# Fetch the raw SRA run file
prefetch SRR15202440

# Extract into split paired-end FASTQ data streams
fasterq-dump SRR15202440 --split-files
```

#### 2. Quality Control & Filtering
Filter out low-quality bases and trailing adapter sequences via `fastp`.
```
fastp \
  -i SRR15202440_1.fastq \
  -I SRR15202440_2.fastq \
  -o SRR15202440_1_clean.fastq \
  -O SRR15202440_2_clean.fastq
```

#### 3. Read Merging
Concatenate forward and reverse clean reads into a single master fastq file to maximize profiling alignment efficiency.

```
cat SRR15202440_1_clean.fastq SRR15202440_2_clean.fastq > SRR15202440_merged.fastq
```

#### 4. Taxonomic & Functional Profiling
Map community taxonomy first using `MetaPhlAn 4`. Use that output to guide `HUMAnN 3` through reconstructing localized functional pathway tables.
```
# Execute taxonomic classification
metaphlan SRR15202440_merged.fastq \
  --input_type fastq \
  --index mpa_vJun23_CHOCOPhlAnSGB_202307 \
  --bowtie2out SRR15202440.bowtie2out.txt \
  --nproc 8 \
  -o SRR15202440_profile.txt

# Reconstruct metabolic pathway abundances
humann \
  --input SRR15202440_merged.fastq \
  --output ~/humann_out/SRR15202440 \
  --threads 8 \
  --taxonomic-profile SRR15202440_profile.txt
```
#### 5. Table Aggregation & Normalization
Consolidate isolated sample pathway tables into a unified matrix, and normalize counts to Relative Abundance (`relab`) to account for varying sequencing depths.
```
# Collect individual pathabundance tables using symbolic links
mkdir -p ~/humann_pathabundance_files
find ~/humann_out -name "*pathabundance.tsv" -exec ln -sf {} ~/humann_pathabundance_files/ \;

# Merge matrices across all samples
humann_join_tables \
  --input ~/humann_pathabundance_files \
  --output ~/merged_pathabundance.tsv \
  --file_name pathabundance

# Normalize absolute values into decimal relative abundances (0.0 - 1.0)
humann_renorm_table \
  --input ~/merged_pathabundance.tsv \
  --output ~/merged_pathabundance_relab.tsv \
  --units relab
```

## Downstream Python Analysis

The `analyze_pathways.py` script automatically removes stratified taxonomic rows, cleans metadata tags, filters out unmapped entries, logs comparative metrics, checks for group divergence via the Mann-Whitney U test, and saves final plots.

```
# Setup output folder and stage the normalized data file
mkdir -p ~/metagenomics_analysis
cd ~/metagenomics_analysis
cp ~/merged_pathabundance_relab.tsv .

# Execute the script
python analyze_pathways.py
```

**Data Reproducibility Note**: The intermediate output spreadsheets (`.csv` and `.tsv` files) are not tracked directly in this repository's main branch due to server migration. However, the complete pipeline code to regenerate these files from scratch is provided in `analyze_pathways.py`, and the final visualized figures are preserved in the `figures/` directory.


#### Generated Analytical Outputs
| File Name | Output Classification | Scientific Purpose / Description |
| --------- | --------------------- | -------------------------------- |
| `top20_pathways.csv` | Processed Dataset |Top 20 overall highest abundance pathways across all cohorts |
| `top20_pathways_barplot.png` | Publication Figure | Horizontal bar plot displaying overall top pathway means |
| `top20_pathways_heatmap.png` | Publication Figure | Clustered heatmap visualizing functional profiles across individual samples |
| `atlanta_vs_maputo_pathway_means.csv` | Processed Dataset | Core comparative matrix showing fold changes and location differences |
| `top_differing_pathways.csv` | Processed Dataset | Top 20 pathways sorted by highest absolute divergence between locations |
| `top_differing_pathways_barplot.png` | Publication Figure | Horizontal bar plot explicitly tracking Maputo vs. Atlanta variance |
| `pathway_stats_mannwhitney.csv` | Statistical Report | Formal tracking document documenting output non-parametric p-values |

## Known Study Limitations

When interpreting these results or adapting this pipeline, please account for the following constraints noted in the paper:

• Cohort Size: The limited sample size (\(n=3\) per region) restricts overall statistical power to detect minor variations.

• DNA Potential vs. Activity: HUMAnN 3 maps metabolic pathway potential derived from metagenomic DNA content; it does not measure active gene transcription (metatranscriptomics) or actual running metabolic active processes (metabolomics).

• Sequencing Depth Variance: Despite scaling variables using `humann_renorm_table`, initial variations in sequencing depth between historical library prep runs can still skew fine-grain abundance estimates.
