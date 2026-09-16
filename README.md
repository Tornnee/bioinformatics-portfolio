# Bioinformatics Analysis Portfolio
**gh | Postgraduate Researcher — Bioinformatics & Molecular Biology | Nigeria**

This repository contains three documented analysis scripts covering the core
bioinformatics workflows I have performed during my postgraduate research.
Each script is fully commented, reproducible, and ready to run on your own data.

---

## Scripts

### 01 — Differential Gene Expression Analysis (`R`)
**File:** `01_differential_expression_analysis.R`
**Dataset:** GEO GSE120103 — Stage IV Ovarian Endometriosis vs Disease-Free Endometrium
**Platform:** Agilent Whole Human Genome Microarray 4x44K (GPL6480)

**What it does:**
- Downloads and loads GEO microarray data using `GEOquery`
- Performs quality control and quantile normalisation using `limma`
- Identifies differentially expressed genes (DEGs) using moderated t-tests
- Generates a volcano plot with labelled significant genes using `ggplot2`
- Performs GO Biological Process and KEGG pathway enrichment using `clusterProfiler`
- Exports results tables as CSV files

**Key findings from this analysis:**
- 137 DEGs identified (adj.P < 0.05): 130 upregulated, 7 downregulated
- Significant KEGG pathway: Neuroactive ligand-receptor interaction (FDR = 0.032)
- Hub genes identified: PRL, PRLHR, TSHB, IRS2, GHSR, APOA4, PITX1

**Tools used:** R, GEOquery, limma, clusterProfiler, ggplot2, ggrepel, org.Hs.eg.db

---

### 02 — ISSR Molecular Marker Analysis (`Python`)
**File:** `02_ISSR_marker_analysis.py`
**Data:** Binary gel band scoring from ISSR PCR — 9 samples, 8 loci

**What it does:**
- Combines replicate gel lanes per sample using logical OR
- Calculates Jaccard and Dice similarity coefficients for all sample pairs
- Calculates genetic distance matrix (1 − Jaccard)
- Performs UPGMA hierarchical clustering
- Generates a dendrogram figure using `scipy` and `matplotlib`
- Calculates band frequency and PIC (Polymorphism Information Content)
- Exports Newick tree file for MEGA 11 / FigTree visualisation
- Exports all matrices as CSV files

**Runs three analyses automatically:**
- Loci 1–4 only
- Loci 5–8 only
- All 8 loci combined

**Tools used:** Python 3, numpy, pandas, scipy, matplotlib

---

### 03 — ESBL Resistance Gene Pipeline (`Bash`)
**File:** `03_ESBL_pipeline.sh`
**Context:** In silico characterisation of ESBL resistance genes in
Enterobacteriaceae from borehole water sources — Lagos State, Nigeria

**What it does:**
- Step 1:  Raw read quality assessment (FastQC)
- Step 2:  Adapter trimming and quality filtering (Trimmomatic)
- Step 3:  De novo genome assembly (SPAdes)
- Step 4:  Assembly quality assessment with N50 threshold filter (QUAST)
- Step 5:  Resistance gene detection (CARD RGI)
- Step 6:  Supplementary HMM-based detection (AMRFinderPlus)
- Step 7:  BLAST confirmation of ESBL genes against NCBI nt database
- Step 8:  Genome annotation (Prokka)
- Step 9:  Multilocus sequence typing — ST131, ST307 detection (MLST)
- Step 10: Multiple sequence alignment with reference sequences (MAFFT)
- Step 11: Maximum likelihood phylogenetic tree construction (IQ-TREE 2)
- Step 12: Automated summary report generation

**Target genes:** blaTEM, blaSHV, blaCTX-M

**Usage:**
```bash
bash 03_ESBL_pipeline.sh ./raw_reads ./results 4
```

**Tools used:** FastQC, Trimmomatic, SPAdes, QUAST, CARD RGI,
AMRFinderPlus, BLAST+, Prokka, MLST 2.0, MAFFT, IQ-TREE 2

---

## Repository Structure

```
bioinformatics-portfolio/
│
├── 01_differential_expression_analysis.R   # DEG analysis (R)
├── 02_ISSR_marker_analysis.py              # Molecular marker analysis (Python)
├── 03_ESBL_pipeline.sh                     # ESBL resistance gene pipeline (Bash)
├── README.md                               # This file
│
├── data/                                   # Example input data (not uploaded)
│   ├── example_binary_scoring.csv          # Example ISSR scoring matrix
│   └── reference_seqs/                     # Reference FASTA sequences for alignment
│
└── results/                                # Example outputs (not uploaded — too large)
    ├── DEG_results_filtered.csv
    ├── ISSR_All_Loci_jaccard_similarity.csv
    ├── ISSR_All_Loci_dendrogram.png
    └── PIPELINE_SUMMARY.txt
```

---

## How to Run

### Script 01 — R
```r
# Install dependencies first (run once)
install.packages("BiocManager")
BiocManager::install(c("GEOquery","limma","clusterProfiler","org.Hs.eg.db","enrichplot"))
install.packages(c("ggplot2","ggrepel","dplyr","pheatmap"))

# Then run the script
source("01_differential_expression_analysis.R")
```

### Script 02 — Python
```bash
# Install dependencies
pip install numpy pandas scipy matplotlib

# Run the script
python3 02_ISSR_marker_analysis.py
```

### Script 03 — Bash (Linux / WSL2)
```bash
# Install dependencies via conda
conda install -c bioconda fastqc trimmomatic spades quast blast mafft prokka mlst iqtree
pip install rgi

# Make the script executable
chmod +x 03_ESBL_pipeline.sh

# Run with your data
bash 03_ESBL_pipeline.sh ./raw_reads ./results 4
```

---

## Research Context

These scripts were developed as part of my postgraduate research in
Bioinformatics and Molecular Biology. The three scripts together cover:

- **Transcriptomics** — microarray-based differential expression analysis
- **Population genetics** — ISSR marker scoring and genetic diversity analysis
- **Genomic epidemiology** — ESBL resistance gene detection and phylogenetics

Each script is designed to be reusable: replace the input data with your own
dataset and the pipeline runs end-to-end with no manual intervention.

---

## Contact

- **Institution:** Nigeria (Postgraduate, Bioinformatics & Molecular Biology)
- **Research areas:** Bioinformatics, Molecular Biology, Antimicrobial Resistance,
  Transcriptomics, Molecular Markers
- **Tools:** R, Python, Bash, MEGA 11, Cytoscape, STRING, DAVID, GEO2R,
  CARD, ResFinder, MAFFT, IQ-TREE 2

---

## Citation

If you use or adapt any of these scripts, please acknowledge this repository.
Reference datasets:
- GSE120103: NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120103)
- CARD: McArthur et al. (2023) Antimicrobial Agents and Chemotherapy
- MLST: Larsen et al. (2012) Journal of Clinical Microbiology
