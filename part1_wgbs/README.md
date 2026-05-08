# Part 1 — WGBS: Hierarchical Clustering of Breast Cancer Methylomes

This part of the assignment summarises the methodology and findings of
Lin *et al.* (2015) — **"Hierarchical Clustering of Breast Cancer Methylomes
Revealed Differentially Methylated and Expressed Breast Cancer Genes"**
(*PLOS ONE* 10(2): e0118453) — and walks through how the same Whole-Genome
Bisulfite Sequencing (WGBS) workflow can be reproduced inside the Galaxy
training framework.

---

## 1. Background on WGBS

DNA methylation is the addition of a methyl group to the 5-position of
cytosine, almost always in a **CpG dinucleotide** context in mammalian
somatic tissue. Roughly 70-80 % of CpG sites are methylated in a typical
adult cell, with the dramatic exceptions being **CpG islands (CGIs)** —
clusters of CpGs commonly found at active promoters, where the cytosines
are kept unmethylated to permit transcription.

Aberrant DNA methylation is one of the hallmarks of cancer:

* **Promoter CGI hypermethylation** silences tumour-suppressor genes.
* **Genome-wide hypomethylation** destabilises chromosomal architecture,
  activates retrotransposons, and reawakens tissue-inappropriate
  enhancers.

WGBS measures methylation at **single-base resolution across the entire
genome**. The principle is:

1. Treat genomic DNA with **sodium bisulfite**, which deaminates
   *unmethylated* cytosine into uracil but leaves **5-methylcytosine
   (5mC) untouched**.
2. PCR-amplify and sequence — uracil now reads as thymine.
3. Map the converted reads to a reference genome with a **bisulfite-aware
   aligner** (which handles the C → T conversion ambiguity).
4. At each CpG, count the C reads (methylated) vs T reads (unmethylated)
   and compute a **β-value** = methylated / (methylated + unmethylated)
   in [0, 1].

Of the four major methylation profiling technologies (WGBS,
RRBS, MeDIP-seq, Infinium BeadChips) **WGBS is the gold standard** —
single-CpG resolution, no enrichment bias, full genome coverage. The
trade-off is cost: ~15-30× sequencing depth across ~3 Gb produces ~600 M
reads per sample.

---

## 2. The Lin *et al.* 2015 Study

### 2.1 Samples

The authors performed WGBS on **five new methylomes**:

| Sample  | Tissue / cell type               | Phenotype                                  |
|---------|----------------------------------|--------------------------------------------|
| NB      | Normal breast tissue             | Healthy reference                          |
| BT089   | Fibroadenoma                     | Benign tumour                              |
| BT126   | Invasive ductal carcinoma (IDC)  | Primary malignant tumour                   |
| BT198   | Invasive ductal carcinoma (IDC)  | Primary malignant tumour                   |
| MCF7    | Breast adenocarcinoma cell line  | Established malignant cell line            |

…and re-analysed two published methylomes for comparison:

| Sample  | Tissue / cell type                 | Source               |
|---------|------------------------------------|----------------------|
| HMEC    | Human mammary epithelial cell line | Hon *et al.* 2012    |
| HCC1954 | Ductal carcinoma cell line         | Hon *et al.* 2012    |

### 2.2 Sequencing & alignment statistics

* **Average reads per sample:** 405 M paired-end reads.
* **Aligned to hg18:** ~322 M pairs (79 %).
* **Mean coverage:** **18.8×**.
* **CpG sites covered:** ~26 M (91.3 % of the human genome's ~28 M CpGs).
* **Bisulfite conversion rate:** ≥ 99 % (validated against in-silico
  converted non-CpG cytosines).

### 2.3 Methodology pipeline (as described in the paper)

```
                FASTQ (paired-end, bisulfite-converted)
                                |
                                v
                   Quality control (FastQC / Falco)
                                |
                                v
                  Adapter & quality trimming
                                |
                                v
   Bisulfite-aware alignment to hg18 (BS Seeker / bwameth / Bismark)
                                |
                                v
              Methylation extraction at every CpG site
                                |
                                v
     β-value matrix (CpGs x samples)  +  per-CpG read depth
                                |
            +-------------------+-------------------+
            v                                       v
   HMR detection (per sample)                  PMD detection
   (e.g. methylKit, DSS)                  (Mb-scale low-meth blocks)
            |                                       |
            v                                       v
   Hierarchical clustering of                 Compare across samples
   HMRs across the 7 samples
            |
            v
   Joint analysis with RNA-seq (NB, MCF7, TCGA-BRCA)
            |
            v
   Identification of differentially methylated + expressed genes
```

### 2.4 Key biological findings

**(i) Genome-wide hypomethylation in cell lines.**
Mean methylation in HMEC, HCC1954, and MCF7 was significantly lower than
in primary tissue (t-test p = 0.0199). Cancer cell lines exhibited fewer
intermediately methylated CpGs and more strongly bimodal distributions.

**(ii) Inverted patterns at CGIs vs CpG-poor regions.**
In the malignant samples (BT126, BT198, HCC1954, MCF7):

* **CpG-rich regions** (CGIs and shores) became **hyper**-methylated.
* **CpG-poor regions** became **hypo**-methylated.

Promoter CGIs went from > 80 % unmethylated in NB / HMEC / BT089 down to
~70 % unmethylated in primary tumours and tumour cell lines.

**(iii) Hypomethylated regions (HMRs) reshape during tumourigenesis.**
HMR counts varied from 53 K (NB) up to ~116 K (HCC1954). Comparing each
tumour HMR with its matched NB HMR:

* Non-CGI HMRs **expanded** in tumours (often >8× wider in MCF7 / HCC1954).
* CGI-containing HMRs **contracted** as the underlying CGI was hypermethylated.

**(iv) Hierarchical clustering identified cancer-specific HMR clusters.**
Average-linkage clustering of HMRs across samples produced clusters
enriched for breast- and ovarian-cancer-associated genes (e.g. BRCA1,
ESR1, GATA3) and enhancer regions specific to either normal or
cancer-line cells.

**(v) Aberrant X-chromosome inactivation (XCI).**
In cancer samples and ~50 % of TCGA-BRCA tumours:

* The **XIST promoter was hypermethylated**, suppressing XIST RNA.
* X-linked genes that should be silenced were **hypomethylated and
  over-expressed**.
* Higher expression of the affected X-linked genes was associated with
  **poorer survival** in TCGA-BRCA.

**(vi) Methylation–expression coupling is weaker than the canonical model
implies.**
Many genes were differentially methylated *without* a matching change in
expression, suggesting that DNA methylation in human cancers acts in
concert with chromatin state, transcription-factor availability, and
replication timing rather than as a simple on/off switch.

---

## 3. Reproducing the analysis in Galaxy

The Galaxy Training Network tutorial *"DNA Methylation data analysis"*
(`https://training.galaxyproject.org/.../methylation-seq/tutorial.html`)
provides a hands-on subset of the Lin 2015 data hosted on Zenodo
(`https://zenodo.org/record/557099`). The tutorial uses small read
subsets and a precomputed BAM so the analysis runs in a teaching session,
but the workflow steps are identical to what one would run on the full
study data.

### 3.1 Workflow (Galaxy)

| # | Galaxy tool          | Purpose                                                            |
|---|----------------------|--------------------------------------------------------------------|
| 1 | **Falco** / FastQC   | Per-base QC. Bisulfite reads have a *characteristically* depleted C content and elevated T content because every unmethylated C reads as T after conversion — this is **expected**, not a quality failure. |
| 2 | **Trim Galore / Cutadapt** | Adapter / quality trimming. Optional: bias-aware trimming of the first ~5 bp at the 5' end where M-bias plots typically show non-uniform methylation. |
| 3 | **bwameth**          | Bisulfite-aware mapping to hg38. Uses a "three-letter alphabet" (C ≡ T on the forward strand, G ≡ A on the reverse) and resolves the ambiguity post-hoc. |
| 4 | **MethylDackel mbias** | Position-dependent methylation-bias diagnostic. The four panels (top/bottom × OT/OB) should be flat ± 5 %. Edges that drift get trimmed. |
| 5 | **MethylDackel extract** | Per-CpG methylation extraction (BedGraph of fractions). Use `--fraction` to write β-values directly. |
| 6 | **Wig/BedGraph-to-bigWig** | Convert BedGraph to bigWig for downstream visualisation. |
| 7 | **deepTools computeMatrix + plotProfile** | Aggregate β-values around CpG-island BED regions to generate the canonical "methylation around TSS / CGI" profile. |
| 8 | **Metilene**         | Detect Differentially Methylated Regions (DMRs) between conditions. Uses a circular binary segmentation + 2D-Kolmogorov–Smirnov test. |

### 3.2 Reference command-line equivalents

For users who prefer a command-line workflow, the equivalent steps using
standard tools (versions current as of 2024) are:

```bash
# 1. Quality control
falco subset_1.fastq.gz subset_2.fastq.gz -o qc/

# 2. Adapter / quality trim
trim_galore --paired --illumina subset_1.fastq.gz subset_2.fastq.gz -o trimmed/

# 3. Bisulfite-aware alignment (bwa-meth)
bwameth.py index hg38.fa
bwameth.py --reference hg38.fa \
           trimmed/subset_1_val_1.fq.gz trimmed/subset_2_val_2.fq.gz \
           -t 8 | samtools sort -@4 -o aligned.bam -
samtools index aligned.bam

# 4. M-bias diagnostic
MethylDackel mbias hg38.fa aligned.bam mbias_

# 5. Per-CpG extraction (β-values)
MethylDackel extract --fraction --mergeContext hg38.fa aligned.bam

# 6. Convert to bigWig
bedGraphToBigWig aligned_CpG.meth.bedGraph hg38.chrom.sizes aligned_meth.bw

# 7. Methylation profile around CGIs
computeMatrix reference-point -S aligned_meth.bw -R CpGIslands.bed \
              --referencePoint center -a 2000 -b 2000 -o matrix.gz
plotProfile -m matrix.gz -o cgi_profile.png

# 8. DMR detection between two conditions (example: tumour vs normal)
metilene -a normal_CpG.meth.bedGraph -b tumour_CpG.meth.bedGraph \
         -m 5 -d 0.1 -t 8 > dmrs.tsv
```

### 3.3 Interpreting the outputs

| Output                 | What you look for                                                     |
|------------------------|-----------------------------------------------------------------------|
| Falco "Per-base content" | Inverted C/T ratio (low C, high T) — *expected* in bisulfite reads. |
| MethylDackel M-bias plot | Flat ± 5 % across the read body. Drift at edges → trim those bases. |
| `*_CpG.meth.bedGraph` | One line per CpG: chrom, start, end, methylation %, methylated reads, unmethylated reads. |
| deepTools profile     | Sharp methylation **dip at TSS / CGI centre** in NB; **flatter / partially filled** dip in tumours (loss of CGI hypomethylation).            |
| Metilene DMR table    | CpG-rich DMRs with **q < 0.05** and **|Δβ| > 0.1** are candidates for downstream gene-enrichment or clustering analyses, mirroring the HMR-based clustering performed in Lin 2015. |

---

## 4. Hierarchical clustering of the methylomes (paper-specific result)

Lin *et al.* applied **average-linkage hierarchical clustering** on the
binarised HMR presence/absence matrix across all seven methylomes and
recovered three biologically meaningful super-clusters:

1. **Normal/benign cluster (NB, BT089, HMEC):** dominated by HMRs at
   active promoter CGIs and tissue-specific enhancers.
2. **Primary-tumour cluster (BT126, BT198):** retains many normal HMRs
   but acquires several thousand cancer-specific ones, especially at
   developmental TF binding sites.
3. **Cell-line cluster (MCF7, HCC1954):** exhibits the most extreme
   reorganisation — wholesale expansion of non-CGI HMRs and aggressive
   contraction of CGI HMRs. Many cell-line-specific HMRs overlap with
   X-inactivation escape genes.

The accompanying **joint analysis with RNA-seq (NB and MCF7) and the
TCGA-BRCA cohort** identified a curated set of differentially methylated
*and* expressed genes — including known breast-cancer drivers (e.g.,
*ESR1*, *FOXA1*, *GATA3*, *XIST*) — supporting the methylome as a
clinically informative axis of breast-cancer biology.

---

## 5. Take-aways

* WGBS is the gold-standard approach for genome-wide methylation
  profiling and was essential for resolving Lin 2015's findings on
  HMR dynamics and PMD hypomethylation, both of which require
  single-base resolution at low CpG density.
* Tumour methylomes diverge from normal in a **bidirectional** fashion:
  hyper-methylation at CpG-rich loci and hypo-methylation at CpG-poor
  loci.
* **Hierarchical clustering of HMRs** recovers biologically and
  clinically meaningful sample groupings — a useful template for
  unsupervised stratification of any new cohort of methylomes.
* The relationship between methylation and gene expression is **partial
  and context-dependent**, motivating multi-omic integration (e.g.,
  with RNA-seq, ATAC-seq and ChIP-seq) rather than reliance on
  methylation alone.

---

## References

* Lin, I-H. *et al.* (2015) *Hierarchical clustering of breast cancer
  methylomes revealed differentially methylated and expressed breast
  cancer genes.* PLOS ONE 10(2): e0118453.
  [DOI: 10.1371/journal.pone.0118453](https://doi.org/10.1371/journal.pone.0118453)
* Hon, G. C. *et al.* (2012) *Global DNA hypomethylation coupled to
  repressive chromatin domain formation and gene silencing in breast
  cancer.* Genome Res. 22:246-258.
* Galaxy Training Network — *DNA Methylation data analysis* tutorial:
  https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html
* Galaxy Training Network — *Introduction to DNA Methylation data
  analysis* slides:
  https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/introduction-dna-methylation/slides-plain.html
* Krueger F. & Andrews S. R. (2011) *Bismark: a flexible aligner and
  methylation caller for Bisulfite-Seq applications.* Bioinformatics
  27:1571-1572.
* Pedersen B. S. *et al.* (2014) *Fast and accurate alignment of long
  bisulfite-seq reads.* arXiv:1401.1129 (bwa-meth).
