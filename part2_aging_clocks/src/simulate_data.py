"""
Realistic Methylation Data Simulator
=====================================
Generates EPIC/450K-style methylation matrices with biologically grounded
age-correlated CpG sites. Used as a fallback when GEO data cannot be
downloaded directly (e.g., sandboxed environments). On a normal machine,
prefer download_real_data.py to fetch true GSE40279 / GSE41169 datasets.

The simulator works by:
  1. Loading the bio-learn population-average beta values (~485k CpGs).
  2. Determining age direction (+/-) for each CpG from Horvath/Hannum
     coefficient signs.
  3. For each sample of true age A, perturbing the population mean by a
     small linear factor (A - 50) along the direction vector + Gaussian
     noise -- producing realistic beta values clipped to [0, 1].

Why this works: aging clocks are linear (or near-linear) over methylation,
so simulated samples produce predictions that strongly correlate with the
ground-truth age (Pearson r > 0.95 for the major clocks). Absolute
predictions may have systematic biases that real data would not, which we
flag explicitly in the report.

Author: Asim Ahmed (BSBI-2023, NUST SINES)
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


# Path to bio-learn's bundled population-average beta-value file.
_BIOLEARN_DATA_DIR = Path(
    os.path.dirname(os.path.dirname(__import__("biolearn").__file__))
) / "biolearn" / "data"


def _load_population_baseline() -> pd.Series:
    """Return CpG -> population-mean beta value (~485k CpGs)."""
    path = _BIOLEARN_DATA_DIR / "biolearn_averages_450k.csv"
    df = pd.read_csv(path).set_index("id")["average"]
    return df


_CLOCK_COEFFICIENT_FILES = [
    "Horvath1.csv",        # Horvathv1 (2013)
    "Horvath2.csv",        # Horvathv2 (2018)
    "Hannum.csv",          # Hannum (2013)
    "PhenoAge.csv",        # Levine PhenoAge (2018)
    "Lin.csv",             # Lin (2016)
    "VidalBralo.csv",      # Vidal-Bralo (2018)
    "Zhang_10.csv",        # Zhang sparse (2019)
    "YingCausAge.csv",     # Ying causal (2022)
    "YingDamAge.csv",      # Ying damage (2022)
    "YingAdaptAge.csv",    # Ying adaptive (2022)
]


def _build_cpg_age_directions() -> dict[str, float]:
    """
    Map each CpG to an aging-direction sign (+1 / -1) by majority vote across
    coefficient signs from ten major aging clocks. Voting (rather than just
    using one clock) means the simulated methylation carries enough signal
    for each downstream clock to recover age, including clocks like
    YingAdaptAge / Zhang_10 whose CpG sets are largely non-overlapping with
    Horvath/Hannum.
    """
    votes: dict[str, list[int]] = {}
    for fname in _CLOCK_COEFFICIENT_FILES:
        path = _BIOLEARN_DATA_DIR / fname
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if df.shape[1] < 2:
            continue
        cpg_col, coef_col = df.columns[0], df.columns[1]
        valid = df[df[cpg_col].astype(str).str.startswith("cg")]
        for cpg, coef in zip(valid[cpg_col], valid[coef_col]):
            cpg_s = str(cpg)
            votes.setdefault(cpg_s, []).append(int(np.sign(coef)))
    directions: dict[str, float] = {}
    for cpg, signs in votes.items():
        s = float(np.sign(np.sum(signs)))
        if s != 0.0:
            directions[cpg] = s
    return directions


def simulate_methylation_dataset(
    ages: np.ndarray,
    sex: Optional[np.ndarray] = None,
    sample_prefix: str = "S",
    slope_magnitude: float = 0.0015,
    noise_sd: float = 0.020,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Simulate a (CpG x sample) methylation beta-value matrix.

    Parameters
    ----------
    ages : np.ndarray
        Chronological age (years) for each sample.
    sex : np.ndarray, optional
        Sex codes (1=male, 2=female). Generated randomly if None.
    sample_prefix : str
        Prefix for synthetic GSM-like sample IDs.
    slope_magnitude : float
        Per-year change in beta value at age-correlated CpGs. Tuned so that
        major clocks reach r > 0.97 with simulated chronological age.
    noise_sd : float
        Standard deviation of Gaussian noise added to each (CpG, sample) cell.
    random_state : int
        RNG seed.

    Returns
    -------
    methylation_df : pd.DataFrame
        Index = CpG IDs (cgXXXXX), columns = sample IDs.
    metadata_df : pd.DataFrame
        Index = sample IDs, columns = ['age', 'sex'].
    """
    rng = np.random.default_rng(random_state)

    n_samples = len(ages)
    sample_ids = [f"{sample_prefix}_{i:04d}" for i in range(n_samples)]
    if sex is None:
        sex = rng.choice([1, 2], size=n_samples)

    baseline_series = _load_population_baseline()
    cpg_index = baseline_series.index.tolist()
    baseline = baseline_series.values
    cpg_idx = {c: i for i, c in enumerate(cpg_index)}

    directions = _build_cpg_age_directions()
    slopes = np.zeros(len(cpg_index), dtype=np.float32)
    for cpg, sign in directions.items():
        if cpg in cpg_idx:
            slopes[cpg_idx[cpg]] = slope_magnitude * sign

    # Build matrix column-by-column.
    M = np.empty((len(cpg_index), n_samples), dtype=np.float32)
    for j, age in enumerate(ages):
        age_centered = float(age) - 50.0
        noise = rng.normal(0.0, noise_sd, size=len(cpg_index)).astype(np.float32)
        M[:, j] = np.clip(baseline + slopes * age_centered + noise, 0.005, 0.995)

    methylation_df = pd.DataFrame(M, index=cpg_index, columns=sample_ids)
    metadata_df = pd.DataFrame(
        {"age": np.asarray(ages, dtype=float), "sex": np.asarray(sex, dtype=int)},
        index=sample_ids,
    )
    return methylation_df, metadata_df


def make_dataset_a(random_state: int = 42) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Synthetic stand-in for GSE40279 (Hannum et al. 2013, blood methylome,
    Illumina HumanMethylation450). 60 healthy adults, ages 25-75.
    """
    rng = np.random.default_rng(random_state)
    ages = np.sort(rng.uniform(25, 75, size=60))
    sex = rng.choice([1, 2], size=60)
    return simulate_methylation_dataset(
        ages=ages,
        sex=sex,
        sample_prefix="GSM_DA",
        random_state=random_state,
    )


def make_dataset_b(random_state: int = 7) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Synthetic stand-in for GSE41169 (Horvath/Dutch blood cohort, 450K).
    50 adults, age range 30-90 (older skew), mild noise increase to mimic
    technical heterogeneity in the original GEO cohort.
    """
    rng = np.random.default_rng(random_state)
    ages = np.sort(rng.uniform(30, 90, size=50))
    sex = rng.choice([1, 2], size=50)
    return simulate_methylation_dataset(
        ages=ages,
        sex=sex,
        sample_prefix="GSM_DB",
        noise_sd=0.025,
        random_state=random_state,
    )


if __name__ == "__main__":
    import time

    t0 = time.time()
    dnam_a, meta_a = make_dataset_a()
    print(f"Dataset A: {dnam_a.shape} CpGs x samples in {time.time()-t0:.1f}s")
    print(meta_a.describe())
    t0 = time.time()
    dnam_b, meta_b = make_dataset_b()
    print(f"Dataset B: {dnam_b.shape} CpGs x samples in {time.time()-t0:.1f}s")
    print(meta_b.describe())
