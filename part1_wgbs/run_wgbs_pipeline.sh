#!/usr/bin/env bash
# =============================================================================
# WGBS reproducibility script -- mirrors the Galaxy Training Network workflow
# used for Lin et al. (2015) breast-methylome reanalysis.
#
# Run on a Linux box with conda; tested on Ubuntu 22.04.  This script is
# documentation-quality reference: edit paths and resource flags before use.
#
# Author: Asim Ahmed (BSBI-2023, NUST SINES)
# =============================================================================

set -euo pipefail

# ----------------------- USER CONFIG -----------------------------------------
WORKDIR="${WORKDIR:-$PWD/wgbs_run}"
THREADS="${THREADS:-8}"
GENOME_FA="${GENOME_FA:-references/hg38.fa}"
CGI_BED="${CGI_BED:-references/CpGIslands.hg38.bed}"
CHROM_SIZES="${CHROM_SIZES:-references/hg38.chrom.sizes}"

# Two example samples; extend the list for the full Lin 2015 cohort.
declare -A R1=(
    [NB]="data/NB_1.fastq.gz"
    [MCF7]="data/MCF7_1.fastq.gz"
)
declare -A R2=(
    [NB]="data/NB_2.fastq.gz"
    [MCF7]="data/MCF7_2.fastq.gz"
)

mkdir -p "$WORKDIR"/{qc,trimmed,aligned,extracted,bigwig,profile,dmr}
cd "$WORKDIR"

# ----------------------- 0. ENVIRONMENT --------------------------------------
# Recommended conda env (run once, separately):
#   conda create -n wgbs -c bioconda -c conda-forge -y \
#       falco trim-galore bwameth samtools methyldackel \
#       deeptools ucsc-bedgraphtobigwig metilene
#   conda activate wgbs

# ----------------------- 1. QC (Falco) ---------------------------------------
for s in "${!R1[@]}"; do
    echo "[QC] $s"
    falco "${R1[$s]}" "${R2[$s]}" -o "qc/$s/"
done

# ----------------------- 2. ADAPTER TRIMMING ---------------------------------
for s in "${!R1[@]}"; do
    echo "[Trim] $s"
    trim_galore --paired --illumina --cores "$THREADS" \
                -o "trimmed/$s/" "${R1[$s]}" "${R2[$s]}"
done

# ----------------------- 3. INDEX REFERENCE (run once) -----------------------
if [ ! -f "${GENOME_FA}.bwameth.c2t" ]; then
    echo "[Index] bwameth on $GENOME_FA"
    bwameth.py index "$GENOME_FA"
fi

# ----------------------- 4. ALIGNMENT (bwameth) ------------------------------
for s in "${!R1[@]}"; do
    echo "[Align] $s"
    R1_TRIM="trimmed/$s/$(basename "${R1[$s]}" .fastq.gz)_val_1.fq.gz"
    R2_TRIM="trimmed/$s/$(basename "${R2[$s]}" .fastq.gz)_val_2.fq.gz"
    bwameth.py --reference "$GENOME_FA" -t "$THREADS" \
               "$R1_TRIM" "$R2_TRIM" \
        | samtools sort -@ "$THREADS" -o "aligned/${s}.bam" -
    samtools index "aligned/${s}.bam"
done

# ----------------------- 5. M-BIAS DIAGNOSTIC --------------------------------
for s in "${!R1[@]}"; do
    echo "[M-bias] $s"
    MethylDackel mbias -@ "$THREADS" "$GENOME_FA" \
                       "aligned/${s}.bam" "qc/${s}/${s}_mbias_"
done

# ----------------------- 6. METHYLATION EXTRACTION ---------------------------
for s in "${!R1[@]}"; do
    echo "[Extract] $s"
    MethylDackel extract --fraction --mergeContext -@ "$THREADS" \
                         "$GENOME_FA" "aligned/${s}.bam" -o "extracted/${s}"
done

# ----------------------- 7. BEDGRAPH -> BIGWIG --------------------------------
for s in "${!R1[@]}"; do
    BG="extracted/${s}_CpG.meth.bedGraph"
    BW="bigwig/${s}_CpG.meth.bw"
    sort -k1,1 -k2,2n "$BG" > "${BG}.sorted"
    bedGraphToBigWig "${BG}.sorted" "$CHROM_SIZES" "$BW"
done

# ----------------------- 8. METHYLATION PROFILE AROUND CGIs ------------------
echo "[Profile] computeMatrix + plotProfile"
BW_FILES=()
LABELS=()
for s in "${!R1[@]}"; do
    BW_FILES+=("bigwig/${s}_CpG.meth.bw")
    LABELS+=("$s")
done

computeMatrix reference-point \
    -S "${BW_FILES[@]}" \
    -R "$CGI_BED" \
    --referencePoint center -a 2000 -b 2000 \
    -p "$THREADS" \
    -o profile/cgi_matrix.gz

plotProfile -m profile/cgi_matrix.gz \
            --samplesLabel "${LABELS[@]}" \
            -o profile/cgi_profile.png \
            --plotTitle "Methylation around CGIs"

# ----------------------- 9. DMR DETECTION (Metilene) -------------------------
# Pairwise tumour-vs-normal example. For multi-sample studies, build the
# joint input with metilene_input.pl (shipped with metilene).
echo "[DMR] Metilene NB vs MCF7"
metilene -a "extracted/NB_CpG.meth.bedGraph" \
         -b "extracted/MCF7_CpG.meth.bedGraph" \
         -m 5 -d 0.1 -t "$THREADS" \
        > dmr/NB_vs_MCF7_dmrs.tsv

echo "Done. Top 5 DMRs:"
sort -k4,4n dmr/NB_vs_MCF7_dmrs.tsv | head -n 5

# ----------------------- 10. (Optional) HIERARCHICAL CLUSTERING --------------
# To replicate Lin 2015's HMR-based clustering, build an HMR x sample
# binary matrix using methylKit / DSS or the simple `methyldackel`
# below-threshold pipeline, then cluster with R's hclust(method="average")
# or scipy.cluster.hierarchy.linkage in Python. See
# part1_wgbs/notes_clustering.md for an outline.
