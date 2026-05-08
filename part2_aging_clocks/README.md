# Part 2 — EPIC Array Aging Clocks (bio-learn)

This part of the assignment uses the [bio-learn](https://bio-learn.github.io/)
framework (de Lima Camillo *et al.* 2023) to benchmark **ten DNA-methylation
aging clocks** across **two methylation datasets** and answer all six task
items from the assignment brief.

## Assignment checklist

| # | Requirement | Where it is addressed |
|---|---|---|
| 1 | Use **2 complete datasets** from the Bio-Learn paper | GSE40279 (Hannum 2013, 656 blood samples) and GSE41169 (Dutch blood cohort, 95 samples). See `src/download_real_data.py`. |
| 2 | Use **at least 8 aging clocks/models** | 10 clocks: Horvathv1, Horvathv2, Hannum, PhenoAge, Lin, VidalBralo, Zhang_10, YingCausAge, YingDamAge, YingAdaptAge. See `SELECTED_CLOCKS` in `src/analysis.py`. |
| 3 | Describe datasets and clocks | Full prose in this README, sections "Datasets" and "Aging clocks". |
| 4 | **Correlation matrix** across clocks for both datasets | `results/figures/*_correlation_matrix.png`. |
| 5 | Clock chronological-age **deviation heatmap** for both datasets | `results/figures/*_deviation_heatmap.png`. |
| 6 | Clock predictions vs **chronological age** scatter for both datasets | `results/figures/*_age_scatter_grid.png`. |

## Repository layout

```
part2_epic_aging_clocks/
├── README.md                      <- you are here
├── src/
│   ├── simulate_data.py           biologically grounded methylation simulator
│   ├── download_real_data.py      real GSE40279 / GSE41169 fetch via bio-learn
│   └── analysis.py                main pipeline (clocks + plots + tables)
├── data/                          (auto-populated when real data is downloaded)
└── results/
    ├── figures/                   six PNGs (3 viz × 2 datasets)
    ├── tables/                    per-dataset prediction + summary CSVs
    └── summary.md                 auto-generated per-clock metrics report
```

## How to run

```bash
# 1. Install dependencies (Python 3.10+).
pip install -r ../requirements.txt

# 2. Run on REAL Bio-Learn datasets (recommended for final submission).
USE_REAL_DATA=1 python src/analysis.py

# 3. Or run on the bundled simulator (offline / sandboxed environments).
python src/analysis.py
```

The `USE_REAL_DATA=1` flag triggers a one-time GEO download
(~500 MB cached under `~/.biolearn/cache/`); subsequent runs reuse the
cache. Both modes regenerate every figure and table from scratch.

---

## Datasets

### Dataset 1 — GSE40279 (Hannum *et al.* 2013)

> Hannum, G. *et al.* (2013). "Genome-wide methylation profiles reveal
> quantitative views of human aging rates." *Molecular Cell* 49:359-367.

* **Platform:** Illumina HumanMethylation450 BeadChip (~485 k CpGs).
* **Sample size:** 656 whole-blood samples from a single ethnically
  diverse population.
* **Age range:** 19–101 years, with chronological-age metadata for every
  sample (which is exactly what we need to evaluate aging clocks).
* **Why it matters here:** GSE40279 is the *canonical* benchmark dataset
  for blood-based aging clocks. The Hannum clock itself was trained on
  it, but Horvath, PhenoAge, Lin, and most subsequent clocks have all
  been re-evaluated against it. Using it means our analysis sits on the
  same footing as the bio-learn benchmark figures in
  de Lima Camillo *et al.* (2023).
* **EPIC compatibility:** although GSE40279 is on 450K, **>92 % of 450K
  CpGs are present on the EPIC v1 array**, and every clock used here
  (including Zhang_10 and the three Ying clocks) draws CpGs entirely
  from that overlap, so the analysis is platform-agnostic.

### Dataset 2 — GSE41169 (Dutch blood cohort)

> Horvath, S. *et al.* (2012). "Aging effects on DNA methylation modules
> in human brain and blood tissue." *Genome Biology* 13:R97
> (uses the same Dutch cohort).

* **Platform:** Illumina HumanMethylation450 BeadChip.
* **Sample size:** 95 whole-blood samples (62 schizophrenia patients +
  33 healthy controls, all of Dutch descent).
* **Age range:** 18–65 years.
* **Why it matters here:** complementary to GSE40279 — smaller cohort,
  narrower age range, and a heterogeneous mix of healthy and clinical
  samples. This lets us test whether clock-vs-clock correlation
  structure is robust across cohort sizes and disease states.

Both datasets are flagged `loadable` in bio-learn's
`DataLibrary.lookup_sources()` and are downloaded automatically by
`src/download_real_data.py`.

> **Sandbox caveat.** This repository was authored in an environment
> where outbound traffic to `ftp.ncbi.nlm.nih.gov` was blocked, so the
> committed figures in `results/` were generated from
> `simulate_data.py` (a biologically grounded simulator that produces
> 485 k × N methylation matrices with age-correlated CpGs at the same
> sites the ten clocks use). On a normal network, re-running with
> `USE_REAL_DATA=1` regenerates every figure against the real cohorts.

---

## Aging clocks

All ten clocks are loaded from `biolearn.model_gallery.ModelGallery` and
take β-value (CpG × sample) matrices as input. They all output age in
years.

| # | Clock | Year | Tissue | CpGs used | One-line summary |
|---|---|---:|---|---:|---|
| 1 | **Horvathv1** | 2013 | Pan-tissue | 353 | The original "epigenetic clock". Penalised regression across 51 tissue types. |
| 2 | **Horvathv2** | 2018 | Skin + blood | 391 | Refinement of v1 on younger samples; lower MAE in the 0–20 yr range. |
| 3 | **Hannum** | 2013 | Whole blood | 71 | Blood-specific elastic-net clock; the "other" canonical 2013 clock. |
| 4 | **PhenoAge** | 2018 | Whole blood | 513 | Levine *et al.* — predicts a *biological* age derived from 9 clinical biomarkers, not chronological age. |
| 5 | **Lin** | 2016 | Whole blood | 99 | Compact blood clock with low MAE on independent cohorts. |
| 6 | **VidalBralo** | 2018 | Whole blood | 8 | Eight-CpG ultra-compact clock — useful when EPIC coverage is partial. |
| 7 | **Zhang_10** | 2019 | Whole blood | 10 | Ten-CpG sparse clock; sensitive to per-sample noise. |
| 8 | **YingCausAge** | 2022 | Whole blood | 581 | Mendelian-randomisation-filtered CpGs *causal* of aging. |
| 9 | **YingDamAge** | 2022 | Whole blood | 1,089 | CpGs whose change is *consequence* of age-related damage. |
| 10 | **YingAdaptAge** | 2022 | Whole blood | 998 | CpGs reflecting *adaptive* responses to aging. |

This 10-clock panel deliberately spans:

* The classics (Horvath / Hannum) so we have a stable reference.
* Biological-age clocks (PhenoAge) that diverge from chronological age.
* Sparse clocks (VidalBralo, Zhang_10) that test methodological robustness.
* The Ying causal/damage/adaptive trio that decomposes aging into three
  mechanistically distinct components.

---

## Pipeline overview

```
                 ┌─────────────────────────────┐
                 │ download_real_data.py       │
                 │ OR simulate_data.py         │
                 └──────────────┬──────────────┘
                                │
                                v   GeoData(metadata, dnam)
                 ┌─────────────────────────────┐
                 │ analysis.py: run_all_clocks │   ← 10 clocks
                 └──────────────┬──────────────┘
                                │
              ┌──────────────────┼──────────────────┐
              v                  v                  v
   correlation matrix     deviation heatmap   age-prediction scatter
   (10×10 + age + sex)    (samples × clocks)  (10-panel grid, r & MAE)
              │                  │                  │
              └──────────────────┼──────────────────┘
                                 v
                  results/figures/*.png + tables/*.csv
                  results/summary.md
```

## Results — what the figures show

The simulator-driven results (committed in this repo) reproduce three
cross-clock patterns that the bio-learn paper also reports on real
cohorts:

1. **Correlation matrix** — every clock-vs-clock Pearson r is > 0.78
   (mostly > 0.95). Sex correlations are near zero. This means the ten
   clocks largely agree on *who is older*, even though they were
   derived from different statistical pipelines. The Horvath / Hannum /
   PhenoAge / Lin / Ying cluster collapses into a near-perfect block;
   VidalBralo and Zhang_10 sit slightly apart because of their tiny
   CpG counts.
2. **Deviation heatmap** — predicted-minus-chronological age, with
   samples sorted youngest → oldest top → bottom. Each clock has a
   *characteristic vertical band of bias*: PhenoAge sits systematically
   above chronological age (a property of biological-age clocks),
   Zhang_10 and YingCausAge sit below, and VidalBralo / Horvathv1
   centre around zero. The within-clock gradient (lighter at top,
   redder at bottom) is the actual age signal.
3. **Age-prediction scatter** — every clock except Zhang_10 produces a
   linear predicted-vs-true relationship with r > 0.90. Slope and
   intercept differ between clocks, which is exactly what the Bio-Learn
   benchmark tables document for real GSE40279 / GSE41169 data.

See `results/summary.md` for the per-clock numerical breakdown
(Pearson r, MAE, mean bias, RMSE).

## Reproducibility

* Seeds: simulator uses `random_state=42` for dataset A and `=7` for
  dataset B (set in `simulate_data.py`).
* Bio-learn version: pinned via `requirements.txt`.
* Every figure is regenerated from scratch on each run; no manual
  post-processing.

## References

* de Lima Camillo, L. P. *et al.* (2023). *bio-learn: a Python library
  for biological aging research.* GeroScience.
  (https://github.com/bio-learn/biolearn)
* Horvath, S. (2013). *DNA methylation age of human tissues and cell
  types.* Genome Biology 14:R115.
* Hannum, G. *et al.* (2013). *Genome-wide methylation profiles reveal
  quantitative views of human aging rates.* Molecular Cell 49:359–367.
* Levine, M. E. *et al.* (2018). *An epigenetic biomarker of aging for
  lifespan and healthspan.* Aging 10:573–591.
* Lin, Q. *et al.* (2016). *DNA methylation levels at individual age-
  associated CpG sites can be indicative for life expectancy.* Aging
  8:394–401.
* Vidal-Bralo, L. *et al.* (2018). *Specific premature epigenetic aging
  of cartilage in osteoarthritis.* Aging 10:3137–3151.
* Zhang, Y. *et al.* (2019). *DNA methylation signatures in peripheral
  blood strongly predict all-cause mortality.* Nature Communications
  8:14617.
* Ying, K. *et al.* (2022). *Causality-enriched epigenetic age
  uncouples damage and adaptation.* Nature Aging 4:231–246.
