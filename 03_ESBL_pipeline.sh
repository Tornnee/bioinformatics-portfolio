#!/bin/bash
# =============================================================================
# SCRIPT 03: Bioinformatics Pipeline Automation
# Author:    gh (Postgraduate Researcher, Bioinformatics & Molecular Biology)
# Purpose:   Automates sequence QC, assembly, resistance gene detection,
#            BLAST confirmation, alignment, and phylogenetic tree construction
#            for ESBL resistance gene characterisation in Enterobacteriaceae
# Usage:     bash 03_ESBL_pipeline.sh [INPUT_DIR] [OUTPUT_DIR] [THREADS]
# Example:   bash 03_ESBL_pipeline.sh ./raw_reads ./results 4
# =============================================================================
# DEPENDENCIES (must be installed and in PATH):
#   fastqc, trimmomatic, spades.py, quast.py
#   resfinder (via CGE API or local install)
#   rgi (CARD Resistance Gene Identifier)
#   amrfinder (NCBI AMRFinderPlus)
#   blastn (NCBI BLAST+)
#   mafft
#   iqtree2 or mega
#   mlst (from MLST 2.0 / CGE)
#   prokka
# Install via conda:
#   conda install -c bioconda fastqc trimmomatic spades quast blast mafft
#   conda install -c bioconda prokka mlst
#   pip install rgi
# =============================================================================

set -euo pipefail
# set -e : exit immediately if any command fails
# set -u : treat unset variables as errors
# set -o pipefail : catch errors in pipes

# ── COLOUR OUTPUT ─────────────────────────────────────────────────────────────
RED='\033[0;31m'
GRN='\033[0;32m'
YLW='\033[1;33m'
BLU='\033[0;34m'
NC='\033[0m' # No Colour

log_info()  { echo -e "${GRN}[INFO]${NC}  $(date '+%H:%M:%S')  $1"; }
log_warn()  { echo -e "${YLW}[WARN]${NC}  $(date '+%H:%M:%S')  $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%H:%M:%S')  $1"; }
log_step()  { echo -e "\n${BLU}══════ $1 ══════${NC}"; }

# ── ARGUMENTS ─────────────────────────────────────────────────────────────────
INPUT_DIR="${1:-./raw_reads}"   # Directory containing paired FASTQ files
OUTPUT_DIR="${2:-./results}"    # Where all output will be written
THREADS="${3:-4}"               # Number of CPU threads to use

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
ADAPTER_FILE="/opt/trimmomatic/adapters/TruSeq3-PE-2.fa"  # Adjust path if needed
MIN_READ_LEN=36          # Minimum read length after trimming
PHRED_CUTOFF=20          # Minimum Phred quality score (Q20)
MIN_N50=10000            # Minimum N50 for assembly quality (10 kb)
BLAST_DB="nt"            # NCBI BLAST database to search
BLAST_EVALUE="1e-5"      # E-value threshold for BLAST hits
BLAST_IDENTITY=90        # Minimum % identity to keep a BLAST hit
ORGANISM="Escherichia"   # Organism for AMRFinderPlus (Escherichia or Klebsiella)

# Target ESBL gene families to search
ESBL_GENES=("blaTEM" "blaSHV" "blaCTX-M")

# ── VALIDATE INPUTS ───────────────────────────────────────────────────────────
log_step "VALIDATING INPUTS"

if [[ ! -d "$INPUT_DIR" ]]; then
    log_error "Input directory not found: $INPUT_DIR"
    exit 1
fi

FASTQ_COUNT=$(ls "$INPUT_DIR"/*_R1*.fastq.gz 2>/dev/null | wc -l || echo 0)
if [[ "$FASTQ_COUNT" -eq 0 ]]; then
    log_error "No paired FASTQ files found in $INPUT_DIR (expected *_R1*.fastq.gz)"
    exit 1
fi
log_info "Found $FASTQ_COUNT paired FASTQ sample(s)"

# ── CREATE OUTPUT DIRECTORIES ─────────────────────────────────────────────────
mkdir -p "$OUTPUT_DIR"/{fastqc_raw,fastqc_trimmed,trimmed,assemblies,quast,\
resfinder,card_rgi,amrfinder,blast,alignment,phylogeny,mlst,prokka,logs}

LOG_FILE="$OUTPUT_DIR/logs/pipeline_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOG_FILE") 2>&1
log_info "All output will be saved to: $OUTPUT_DIR"
log_info "Log file: $LOG_FILE"

# ── TRACK QUALITY-FAILED SAMPLES ─────────────────────────────────────────────
FAILED_QC=()
PASSED_QC=()

# =============================================================================
# STEP 1: RAW READ QUALITY ASSESSMENT (FastQC)
# =============================================================================
log_step "STEP 1: RAW READ QUALITY ASSESSMENT"

for R1 in "$INPUT_DIR"/*_R1*.fastq.gz; do
    SAMPLE=$(basename "$R1" | sed 's/_R1.*//')
    R2="${R1/_R1/_R2}"

    if [[ ! -f "$R2" ]]; then
        log_warn "No R2 file found for $SAMPLE — skipping"
        continue
    fi

    log_info "FastQC: $SAMPLE"
    fastqc \
        "$R1" "$R2" \
        --outdir "$OUTPUT_DIR/fastqc_raw" \
        --threads "$THREADS" \
        --quiet \
        2>>"$OUTPUT_DIR/logs/fastqc_raw_${SAMPLE}.log"
done

log_info "Raw FastQC complete. Reports in: $OUTPUT_DIR/fastqc_raw/"

# =============================================================================
# STEP 2: ADAPTER TRIMMING & QUALITY FILTERING (Trimmomatic)
# =============================================================================
log_step "STEP 2: TRIMMING & QUALITY FILTERING"

for R1 in "$INPUT_DIR"/*_R1*.fastq.gz; do
    SAMPLE=$(basename "$R1" | sed 's/_R1.*//')
    R2="${R1/_R1/_R2}"
    [[ ! -f "$R2" ]] && continue

    log_info "Trimmomatic: $SAMPLE"

    OUT_R1P="$OUTPUT_DIR/trimmed/${SAMPLE}_R1_paired.fastq.gz"
    OUT_R1U="$OUTPUT_DIR/trimmed/${SAMPLE}_R1_unpaired.fastq.gz"
    OUT_R2P="$OUTPUT_DIR/trimmed/${SAMPLE}_R2_paired.fastq.gz"
    OUT_R2U="$OUTPUT_DIR/trimmed/${SAMPLE}_R2_unpaired.fastq.gz"

    trimmomatic PE \
        -threads "$THREADS" \
        -phred33 \
        "$R1" "$R2" \
        "$OUT_R1P" "$OUT_R1U" \
        "$OUT_R2P" "$OUT_R2U" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10:8:true \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:"$PHRED_CUTOFF" \
        MINLEN:"$MIN_READ_LEN" \
        2>"$OUTPUT_DIR/logs/trimmomatic_${SAMPLE}.log"

    log_info "Trimmomatic complete: $SAMPLE"
done

# Post-trim FastQC
log_info "Running FastQC on trimmed reads..."
fastqc \
    "$OUTPUT_DIR/trimmed"/*_paired.fastq.gz \
    --outdir "$OUTPUT_DIR/fastqc_trimmed" \
    --threads "$THREADS" \
    --quiet

log_info "Trimmed FastQC complete. Reports in: $OUTPUT_DIR/fastqc_trimmed/"

# =============================================================================
# STEP 3: DE NOVO GENOME ASSEMBLY (SPAdes)
# =============================================================================
log_step "STEP 3: DE NOVO ASSEMBLY"

for R1P in "$OUTPUT_DIR/trimmed"/*_R1_paired.fastq.gz; do
    SAMPLE=$(basename "$R1P" | sed 's/_R1_paired.*//')
    R2P="$OUTPUT_DIR/trimmed/${SAMPLE}_R2_paired.fastq.gz"
    [[ ! -f "$R2P" ]] && continue

    ASSEMBLY_DIR="$OUTPUT_DIR/assemblies/$SAMPLE"
    log_info "SPAdes assembly: $SAMPLE"

    spades.py \
        -1 "$R1P" \
        -2 "$R2P" \
        --careful \
        -k 21,33,55,77 \
        --threads "$THREADS" \
        --memory 16 \
        -o "$ASSEMBLY_DIR" \
        2>"$OUTPUT_DIR/logs/spades_${SAMPLE}.log"

    if [[ -f "$ASSEMBLY_DIR/contigs.fasta" ]]; then
        log_info "Assembly complete: $SAMPLE"
        # Copy contigs to a central location with sample name
        cp "$ASSEMBLY_DIR/contigs.fasta" \
           "$OUTPUT_DIR/assemblies/${SAMPLE}_contigs.fasta"
    else
        log_warn "Assembly failed for $SAMPLE — no contigs.fasta produced"
    fi
done

# =============================================================================
# STEP 4: ASSEMBLY QUALITY ASSESSMENT (QUAST)
# =============================================================================
log_step "STEP 4: ASSEMBLY QUALITY ASSESSMENT"

ALL_CONTIGS=("$OUTPUT_DIR/assemblies"/*_contigs.fasta)

if [[ ${#ALL_CONTIGS[@]} -gt 0 ]]; then
    quast.py \
        "${ALL_CONTIGS[@]}" \
        --output-dir "$OUTPUT_DIR/quast" \
        --threads "$THREADS" \
        --min-contig 500 \
        2>"$OUTPUT_DIR/logs/quast.log"
    log_info "QUAST report: $OUTPUT_DIR/quast/report.html"
fi

# Filter assemblies by N50 threshold
log_info "Checking assembly N50 (threshold: $MIN_N50 bp)..."
QUAST_REPORT="$OUTPUT_DIR/quast/transposed_report.tsv"

if [[ -f "$QUAST_REPORT" ]]; then
    while IFS=$'\t' read -r assembly n50 rest; do
        [[ "$assembly" == "Assembly" ]] && continue  # skip header
        SAMPLE=$(basename "$assembly" | sed 's/_contigs//')
        if (( $(echo "$n50 < $MIN_N50" | bc -l) )); then
            log_warn "FAILED QC (N50=$n50): $SAMPLE — excluded from downstream analysis"
            FAILED_QC+=("$SAMPLE")
        else
            log_info "PASSED QC (N50=$n50): $SAMPLE"
            PASSED_QC+=("$SAMPLE")
        fi
    done < <(awk 'NR>1 {print $1"\t"$18}' "$QUAST_REPORT" 2>/dev/null || \
             tail -n +2 "$QUAST_REPORT" | cut -f1,18)
fi

log_info "Samples passing QC: ${#PASSED_QC[@]}"
log_info "Samples failing QC: ${#FAILED_QC[@]}"
if [[ ${#FAILED_QC[@]} -gt 0 ]]; then
    log_warn "Excluded: ${FAILED_QC[*]}"
fi

# =============================================================================
# STEP 5: RESISTANCE GENE DETECTION (CARD RGI)
# =============================================================================
log_step "STEP 5: RESISTANCE GENE DETECTION — CARD RGI"

for SAMPLE in "${PASSED_QC[@]}"; do
    CONTIGS="$OUTPUT_DIR/assemblies/${SAMPLE}_contigs.fasta"
    [[ ! -f "$CONTIGS" ]] && continue

    log_info "CARD RGI: $SAMPLE"

    rgi main \
        --input_sequence "$CONTIGS" \
        --output_file "$OUTPUT_DIR/card_rgi/${SAMPLE}_rgi" \
        --input_type contig \
        --alignment_tool BLAST \
        --num_threads "$THREADS" \
        --clean \
        2>"$OUTPUT_DIR/logs/rgi_${SAMPLE}.log"

    log_info "RGI complete: $OUTPUT_DIR/card_rgi/${SAMPLE}_rgi.txt"
done

# =============================================================================
# STEP 6: SUPPLEMENTARY DETECTION — AMRFinderPlus
# =============================================================================
log_step "STEP 6: SUPPLEMENTARY DETECTION — AMRFinderPlus"

for SAMPLE in "${PASSED_QC[@]}"; do
    CONTIGS="$OUTPUT_DIR/assemblies/${SAMPLE}_contigs.fasta"
    [[ ! -f "$CONTIGS" ]] && continue

    log_info "AMRFinderPlus: $SAMPLE"

    amrfinder \
        --nucleotide "$CONTIGS" \
        --organism "$ORGANISM" \
        --output "$OUTPUT_DIR/amrfinder/${SAMPLE}_amrfinder.tsv" \
        --threads "$THREADS" \
        --plus \
        2>"$OUTPUT_DIR/logs/amrfinder_${SAMPLE}.log"

    log_info "AMRFinderPlus complete: $OUTPUT_DIR/amrfinder/${SAMPLE}_amrfinder.tsv"
done

# =============================================================================
# STEP 7: BLAST CONFIRMATION OF ESBL GENES
# =============================================================================
log_step "STEP 7: BLAST CONFIRMATION"

for SAMPLE in "${PASSED_QC[@]}"; do
    CONTIGS="$OUTPUT_DIR/assemblies/${SAMPLE}_contigs.fasta"
    [[ ! -f "$CONTIGS" ]] && continue

    log_info "BLASTn: $SAMPLE vs $BLAST_DB"

    blastn \
        -query "$CONTIGS" \
        -db "$BLAST_DB" \
        -out "$OUTPUT_DIR/blast/${SAMPLE}_blast_results.tsv" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        -evalue "$BLAST_EVALUE" \
        -perc_identity "$BLAST_IDENTITY" \
        -max_target_seqs 10 \
        -num_threads "$THREADS" \
        2>"$OUTPUT_DIR/logs/blast_${SAMPLE}.log"

    # Filter for ESBL gene hits only
    ESBL_PATTERN="blaTEM|blaSHV|blaCTX-M|ESBL|beta-lactamase"
    grep -E "$ESBL_PATTERN" "$OUTPUT_DIR/blast/${SAMPLE}_blast_results.tsv" \
        > "$OUTPUT_DIR/blast/${SAMPLE}_ESBL_hits.tsv" || true

    ESBL_HITS=$(wc -l < "$OUTPUT_DIR/blast/${SAMPLE}_ESBL_hits.tsv")
    log_info "ESBL BLAST hits for $SAMPLE: $ESBL_HITS"
done

# Add header to BLAST output files
BLAST_HEADER="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle"
for tsv in "$OUTPUT_DIR/blast"/*_ESBL_hits.tsv; do
    { echo -e "$BLAST_HEADER"; cat "$tsv"; } > "${tsv%.tsv}_with_header.tsv"
done

# =============================================================================
# STEP 8: GENOME ANNOTATION (Prokka)
# =============================================================================
log_step "STEP 8: GENOME ANNOTATION — Prokka"

for SAMPLE in "${PASSED_QC[@]}"; do
    CONTIGS="$OUTPUT_DIR/assemblies/${SAMPLE}_contigs.fasta"
    [[ ! -f "$CONTIGS" ]] && continue

    log_info "Prokka annotation: $SAMPLE"

    prokka \
        --outdir "$OUTPUT_DIR/prokka/$SAMPLE" \
        --prefix "$SAMPLE" \
        --genus Escherichia \
        --species coli \
        --cpus "$THREADS" \
        --force \
        "$CONTIGS" \
        2>"$OUTPUT_DIR/logs/prokka_${SAMPLE}.log"

    log_info "Annotation complete: $OUTPUT_DIR/prokka/$SAMPLE/$SAMPLE.gff"
done

# =============================================================================
# STEP 9: MULTILOCUS SEQUENCE TYPING (MLST)
# =============================================================================
log_step "STEP 9: SEQUENCE TYPING — MLST"

MLST_RESULTS="$OUTPUT_DIR/mlst/MLST_all_samples.tsv"
echo -e "Sample\tScheme\tST\tAlleles" > "$MLST_RESULTS"

for SAMPLE in "${PASSED_QC[@]}"; do
    CONTIGS="$OUTPUT_DIR/assemblies/${SAMPLE}_contigs.fasta"
    [[ ! -f "$CONTIGS" ]] && continue

    log_info "MLST: $SAMPLE"

    MLST_OUT=$(mlst \
        --scheme ecoli \
        --threads "$THREADS" \
        "$CONTIGS" 2>/dev/null || echo "MLST failed for $SAMPLE")

    echo -e "${SAMPLE}\t${MLST_OUT}" >> "$MLST_RESULTS"
done

log_info "MLST results: $MLST_RESULTS"

# =============================================================================
# STEP 10: MULTIPLE SEQUENCE ALIGNMENT (MAFFT)
# =============================================================================
log_step "STEP 10: MULTIPLE SEQUENCE ALIGNMENT — MAFFT"

# This step aligns extracted ESBL gene sequences with reference sequences
# Reference FASTA files should be placed in ./reference_seqs/
REF_DIR="./reference_seqs"

for GENE in "${ESBL_GENES[@]}"; do
    QUERY_FASTA="$OUTPUT_DIR/blast/${GENE}_query_sequences.fasta"
    REF_FASTA="$REF_DIR/${GENE}_references.fasta"

    if [[ ! -f "$QUERY_FASTA" ]]; then
        log_warn "No query sequences for $GENE — skipping alignment"
        continue
    fi

    if [[ ! -f "$REF_FASTA" ]]; then
        log_warn "No reference sequences for $GENE in $REF_DIR — skipping alignment"
        continue
    fi

    # Combine query + reference sequences
    COMBINED="$OUTPUT_DIR/alignment/${GENE}_combined.fasta"
    cat "$QUERY_FASTA" "$REF_FASTA" > "$COMBINED"

    log_info "MAFFT alignment: $GENE ($(grep -c '>' "$COMBINED") sequences)"

    mafft \
        --auto \
        --thread "$THREADS" \
        "$COMBINED" \
        > "$OUTPUT_DIR/alignment/${GENE}_aligned.fasta" \
        2>"$OUTPUT_DIR/logs/mafft_${GENE}.log"

    log_info "Alignment saved: $OUTPUT_DIR/alignment/${GENE}_aligned.fasta"
done

# =============================================================================
# STEP 11: PHYLOGENETIC TREE CONSTRUCTION (IQ-TREE 2)
# =============================================================================
log_step "STEP 11: PHYLOGENETIC TREE — IQ-TREE 2"

for GENE in "${ESBL_GENES[@]}"; do
    ALIGNED="$OUTPUT_DIR/alignment/${GENE}_aligned.fasta"
    [[ ! -f "$ALIGNED" ]] && continue

    log_info "IQ-TREE2 Maximum Likelihood tree: $GENE"

    iqtree2 \
        -s "$ALIGNED" \
        --prefix "$OUTPUT_DIR/phylogeny/${GENE}_ML_tree" \
        -m TEST \
        -B 1000 \
        --bnni \
        -T "$THREADS" \
        2>"$OUTPUT_DIR/logs/iqtree_${GENE}.log"

    log_info "ML tree saved: $OUTPUT_DIR/phylogeny/${GENE}_ML_tree.treefile"
done

# =============================================================================
# STEP 12: SUMMARY REPORT
# =============================================================================
log_step "STEP 12: GENERATING SUMMARY REPORT"

SUMMARY="$OUTPUT_DIR/PIPELINE_SUMMARY.txt"
{
    echo "============================================"
    echo "  ESBL Pipeline Summary Report"
    echo "  Generated: $(date)"
    echo "============================================"
    echo ""
    echo "Input directory:  $INPUT_DIR"
    echo "Output directory: $OUTPUT_DIR"
    echo "Threads used:     $THREADS"
    echo ""
    echo "SAMPLES PROCESSED"
    echo "-----------------"
    echo "  Total input samples: $FASTQ_COUNT"
    echo "  Passed QC (N50 > $MIN_N50 bp): ${#PASSED_QC[@]}"
    echo "  Failed QC: ${#FAILED_QC[@]}"
    if [[ ${#FAILED_QC[@]} -gt 0 ]]; then
        echo "  Failed samples: ${FAILED_QC[*]}"
    fi
    echo ""
    echo "OUTPUT FILES"
    echo "------------"
    echo "  FastQC reports:        $OUTPUT_DIR/fastqc_trimmed/"
    echo "  Trimmed reads:         $OUTPUT_DIR/trimmed/"
    echo "  Genome assemblies:     $OUTPUT_DIR/assemblies/"
    echo "  Assembly QC (QUAST):   $OUTPUT_DIR/quast/report.html"
    echo "  CARD RGI results:      $OUTPUT_DIR/card_rgi/"
    echo "  AMRFinderPlus results: $OUTPUT_DIR/amrfinder/"
    echo "  BLAST ESBL hits:       $OUTPUT_DIR/blast/"
    echo "  Prokka annotations:    $OUTPUT_DIR/prokka/"
    echo "  MLST sequence types:   $OUTPUT_DIR/mlst/MLST_all_samples.tsv"
    echo "  Sequence alignments:   $OUTPUT_DIR/alignment/"
    echo "  Phylogenetic trees:    $OUTPUT_DIR/phylogeny/"
    echo ""
    echo "ESBL GENE SUMMARY"
    echo "-----------------"
    for GENE in "${ESBL_GENES[@]}"; do
        COUNT=$(cat "$OUTPUT_DIR/blast"/*_ESBL_hits.tsv 2>/dev/null | \
                grep -c "$GENE" 2>/dev/null || echo 0)
        echo "  $GENE hits detected: $COUNT"
    done
    echo ""
    echo "Log file: $LOG_FILE"
    echo "============================================"
} > "$SUMMARY"

cat "$SUMMARY"
log_info "Pipeline complete! Summary: $SUMMARY"
