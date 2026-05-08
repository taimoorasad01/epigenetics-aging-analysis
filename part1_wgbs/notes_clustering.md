# Notes — Replicating the Hierarchical Clustering of Methylomes

The flagship analysis in Lin *et al.* (2015) is a hierarchical clustering of
the seven WGBS methylomes by their **HMR (Hypomethylated Region) profile**.
This note distills the approach so it can be reproduced after running the
WGBS pipeline (`run_wgbs_pipeline.sh`) on a study cohort.

## Inputs

For each sample _i_, we need:

* **HMR call set** — a BED file of segments where the local mean β-value
  drops below a threshold (typically β < 0.1) over a minimum number of CpGs.
  Generate with `methylKit::methSeg`, `DSS`, or by writing a simple sliding-
  window caller on the `MethylDackel` BedGraph.

## Building the HMR x sample matrix

1. **Union of HMRs.** Merge all per-sample HMR BED files with
   `bedtools multiinter` to get a master list of HMR loci across the cohort.
2. **Presence/absence matrix.** For each master HMR _j_ and each sample _i_,
   set `M[j,i] = 1` if any HMR in sample _i_ overlaps locus _j_, else 0.
   This is the matrix used by Lin et al.

Optional: produce a continuous matrix of mean β-value per master HMR per
sample for quantitative clustering.

## Hierarchical clustering

```r
library(stats); library(pheatmap)

# Read binary HMR x sample matrix (rows = HMRs, cols = samples).
M <- read.table("HMR_matrix.tsv", header = TRUE, row.names = 1)

# Sample distance via Jaccard distance on binary HMR membership.
jaccard_dist <- function(x) {
    n <- ncol(x); D <- matrix(0, n, n)
    for (i in 1:n) for (j in 1:n) {
        a <- x[, i]; b <- x[, j]
        D[i, j] <- 1 - sum(a & b) / max(1, sum(a | b))
    }
    as.dist(D)
}

d  <- jaccard_dist(M)
hc <- hclust(d, method = "average")
plot(hc)

pheatmap(M,
         clustering_distance_cols = d,
         clustering_method        = "average",
         show_rownames            = FALSE)
```

Lin et al. recover three super-clusters:

* **Normal / benign:** NB, BT089, HMEC.
* **Primary tumours:** BT126, BT198.
* **Cancer cell lines:** MCF7, HCC1954.

## Joint methylation x expression analysis

For each *cluster-specific* HMR set (i.e., HMRs differential between two
super-clusters):

1. Annotate to nearest gene (e.g., `bedtools closest` against GENCODE
   TSS bed).
2. Cross-reference against an RNA-seq differential expression table
   (NB vs MCF7, or per-tumour TCGA-BRCA log-fold-change).
3. Retain genes with **both** |Δβ| > 0.2 across the cluster contrast and
   |log2FC| > 1 in expression.

The intersection in Lin 2015 highlighted *XIST*, *ESR1*, *FOXA1*, *GATA3*,
and several X-linked genes whose hypomethylation correlated with reduced
breast-cancer survival in TCGA-BRCA.

## Caveats

* The binary HMR-presence matrix is sensitive to HMR caller parameters
  (β threshold, minimum CpGs). Sensitivity-analyse by re-running at
  multiple thresholds.
* The seven-sample cohort is small; bootstrap stability of the
  dendrogram before drawing strong conclusions.
* Methylation–expression coupling is partial; mismatches between Δβ and
  ΔlogFC are biologically informative, not noise.
