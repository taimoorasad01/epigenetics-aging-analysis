"""
Main Analysis Pipeline -- EPIC Array Aging Clocks (Bio-Learn)
=============================================================
Runs the assignment requirements end-to-end:

  (1) Loads two methylation datasets (real GEO data if available, else
      simulated stand-ins).
  (2) Runs >= 8 aging clocks per dataset (we run 10).
  (3) Saves per-dataset prediction tables to results/tables/.
  (4) Generates the three required visualizations per dataset:
        - correlation matrix across clocks
        - clock chronological-age-deviation heatmap
        - predicted vs chronological age scatter
  (5) Writes a results/summary.md report.

Run from the repo root:
    python part2_epic_aging_clocks/src/analysis.py

Author: Asim Ahmed (BSBI-2023, NUST SINES)
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

# Local imports
SRC_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SRC_DIR))

from simulate_data import make_dataset_a, make_dataset_b  # noqa: E402

from biolearn.data_library import GeoData  # noqa: E402
from biolearn.model_gallery import ModelGallery  # noqa: E402


# ------------------------------------------------------------------ Config
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
FIG_DIR = RESULTS_DIR / "figures"
TBL_DIR = RESULTS_DIR / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)

# Ten aging clocks - all blood-applicable, all output age in years.
SELECTED_CLOCKS: list[str] = [
    "Horvathv1",       # Horvath (2013), pan-tissue
    "Horvathv2",       # Horvath (2018), skin + blood refinement
    "Hannum",          # Hannum (2013), blood-specific
    "PhenoAge",        # Levine et al. (2018), biological age (mortality)
    "Lin",             # Lin et al. (2016), 99-CpG blood clock
    "VidalBralo",      # Vidal-Bralo et al. (2018), 8-CpG compact clock
    "Zhang_10",        # Zhang et al. (2019), 10-CpG sparse clock
    "YingCausAge",     # Ying et al. (2022), causality-filtered clock
    "YingDamAge",      # Ying et al. (2022), damage-correlated CpGs
    "YingAdaptAge",    # Ying et al. (2022), adaptive-response CpGs
]


# --------------------------------------------------------------- Utilities
def load_or_simulate() -> dict[str, GeoData]:
    """
    Try real GEO data first; fall back to simulator if download fails.
    Returns dict {label: GeoData}.
    """
    use_real = os.environ.get("USE_REAL_DATA", "0") == "1"
    if use_real:
        try:
            from download_real_data import load_real_datasets  # noqa: WPS433
            real = load_real_datasets()
            return {
                "GSE40279_Hannum_blood": real["GSE40279"],
                "GSE41169_Dutch_blood": real["GSE41169"],
            }
        except Exception as exc:  # pragma: no cover
            print(
                f"[analysis] Real-data download failed ({exc}); "
                "falling back to simulator.",
                file=sys.stderr,
            )

    print("[analysis] Generating simulated dataset A (60 samples, 25-75 yr) ...")
    dnam_a, meta_a = make_dataset_a()
    print("[analysis] Generating simulated dataset B (50 samples, 30-90 yr) ...")
    dnam_b, meta_b = make_dataset_b()

    return {
        "DatasetA_simulated_blood_25-75yr": GeoData(meta_a, dnam_a),
        "DatasetB_simulated_blood_30-90yr": GeoData(meta_b, dnam_b),
    }


def run_all_clocks(
    geo: GeoData, clock_names: Optional[list[str]] = None
) -> pd.DataFrame:
    """Run each clock; return DataFrame indexed by sample_id, one column per clock."""
    if clock_names is None:
        clock_names = SELECTED_CLOCKS
    gallery = ModelGallery()
    out = pd.DataFrame(index=geo.metadata.index)
    for name in clock_names:
        try:
            t0 = time.time()
            model = gallery.get(name)
            pred = model.predict(geo)
            # First column holds the predicted age.
            out[name] = pred.iloc[:, 0].reindex(out.index)
            print(f"  - {name:<14} done in {time.time() - t0:5.1f}s")
        except Exception as exc:
            print(f"  - {name:<14} FAILED: {exc}")
            out[name] = np.nan
    return out


# ----------------------------------------------------------- Visualizations
def plot_correlation_matrix(
    preds: pd.DataFrame, ages: pd.Series, sex: pd.Series, out_path: Path, title: str
) -> None:
    """Pearson correlation between all clock predictions, age, and sex."""
    df = preds.copy()
    df["ChronologicalAge"] = ages.reindex(df.index)
    df["Sex"] = sex.reindex(df.index)
    corr = df.corr(method="pearson")
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        corr,
        annot=True,
        fmt=".2f",
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        square=True,
        cbar_kws={"label": "Pearson r"},
        linewidths=0.4,
        linecolor="white",
    )
    plt.title(f"Clock-vs-Clock Correlation Matrix\n{title}", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_deviation_heatmap(
    preds: pd.DataFrame, ages: pd.Series, out_path: Path, title: str
) -> None:
    """Heatmap of (predicted_age - chronological_age) per sample x clock."""
    deviations = preds.sub(ages, axis=0)  # sample x clock
    # Sort samples by chronological age for visual coherence.
    order = ages.sort_values().index
    deviations = deviations.reindex(order)
    vmax = float(np.nanpercentile(np.abs(deviations.values), 98))
    plt.figure(figsize=(11, max(6, 0.20 * len(deviations))))
    sns.heatmap(
        deviations,
        cmap="RdBu_r",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        cbar_kws={"label": "Predicted - Chronological Age (yr)"},
        yticklabels=False,
    )
    plt.title(
        f"Per-sample Age Deviation by Clock\n{title} "
        f"(rows = samples, sorted ascending by chronological age)",
        fontsize=12,
    )
    plt.xlabel("Clock")
    plt.ylabel("Samples (youngest at top -> oldest at bottom)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_age_scatter_grid(
    preds: pd.DataFrame, ages: pd.Series, out_path: Path, title: str
) -> None:
    """Predicted vs chronological age scatter for every clock (3 x 4 grid)."""
    n_clocks = preds.shape[1]
    n_cols = 3
    n_rows = int(np.ceil(n_clocks / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.4 * n_cols, 4.0 * n_rows))
    axes = np.atleast_1d(axes).flatten()

    age_min = float(ages.min()) - 5
    age_max = float(ages.max()) + 5
    pred_min = float(np.nanmin(preds.values)) - 5
    pred_max = float(np.nanmax(preds.values)) + 5
    lim_min = min(age_min, pred_min)
    lim_max = max(age_max, pred_max)

    for i, clock in enumerate(preds.columns):
        ax = axes[i]
        x = ages.values
        y = preds[clock].values
        valid = ~np.isnan(y)
        if valid.sum() >= 3:
            r, _ = pearsonr(x[valid], y[valid])
            mae = float(np.mean(np.abs(y[valid] - x[valid])))
        else:
            r = float("nan")
            mae = float("nan")

        ax.scatter(x, y, alpha=0.65, s=22, edgecolor="white", linewidth=0.4)
        ax.plot([lim_min, lim_max], [lim_min, lim_max], "k--", lw=1, alpha=0.6,
                label="y = x")
        ax.set_title(f"{clock}\n(r = {r:.3f}, MAE = {mae:.1f} yr)", fontsize=10)
        ax.set_xlabel("Chronological age (yr)")
        ax.set_ylabel("Predicted age (yr)")
        ax.set_xlim(lim_min, lim_max)
        ax.set_ylim(lim_min, lim_max)
        ax.grid(alpha=0.25)

    for j in range(n_clocks, len(axes)):
        axes[j].axis("off")

    fig.suptitle(f"Aging Clock Predictions vs Chronological Age -- {title}",
                 fontsize=13, y=1.00)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


# ----------------------------------------------------------------- Driver
def analyse_dataset(label: str, geo: GeoData) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run clocks + plots for a single dataset; return (preds, summary_stats)."""
    print(f"\n{'=' * 70}\nDataset: {label}\n{'=' * 70}")
    print(f"  shape: {geo.dnam.shape[1]} samples x {geo.dnam.shape[0]} CpGs")
    if "age" in geo.metadata.columns:
        print(f"  age range: {geo.metadata['age'].min():.1f} - "
              f"{geo.metadata['age'].max():.1f} yr "
              f"(median {geo.metadata['age'].median():.1f})")

    # Run clocks.
    preds = run_all_clocks(geo)
    preds.to_csv(TBL_DIR / f"{label}_predictions.csv")

    # Per-clock summary stats vs chronological age.
    ages = geo.metadata["age"].astype(float)
    sex = geo.metadata.get("sex", pd.Series(np.nan, index=ages.index))
    rows = []
    for clock in preds.columns:
        y = preds[clock].values
        x = ages.values
        valid = ~np.isnan(y)
        if valid.sum() < 3:
            rows.append([clock, np.nan, np.nan, np.nan, np.nan])
            continue
        r, _ = pearsonr(x[valid], y[valid])
        mae = float(np.mean(np.abs(y[valid] - x[valid])))
        bias = float(np.mean(y[valid] - x[valid]))
        rmse = float(np.sqrt(np.mean((y[valid] - x[valid]) ** 2)))
        rows.append([clock, r, mae, bias, rmse])
    summary = pd.DataFrame(
        rows, columns=["Clock", "Pearson_r", "MAE_years", "Mean_bias_years", "RMSE_years"]
    )
    summary.to_csv(TBL_DIR / f"{label}_summary.csv", index=False)
    print("\n  Per-clock summary:")
    print(summary.round(3).to_string(index=False))

    # Plots.
    plot_correlation_matrix(
        preds, ages, sex,
        FIG_DIR / f"{label}_correlation_matrix.png", label,
    )
    plot_deviation_heatmap(
        preds, ages,
        FIG_DIR / f"{label}_deviation_heatmap.png", label,
    )
    plot_age_scatter_grid(
        preds, ages,
        FIG_DIR / f"{label}_age_scatter_grid.png", label,
    )
    print(f"  -> Figures saved under {FIG_DIR}/")
    return preds, summary


def write_summary_report(
    per_dataset: dict[str, dict],
) -> None:
    out = RESULTS_DIR / "summary.md"
    lines: list[str] = ["# EPIC Array Aging Clocks -- Analysis Summary", ""]
    lines.append(
        "Generated by `part2_epic_aging_clocks/src/analysis.py`. "
        "Each section reports clock-by-clock performance against the "
        "chronological age provided in dataset metadata."
    )
    lines.append("")
    for label, payload in per_dataset.items():
        lines.append(f"## {label}")
        lines.append("")
        lines.append("| Clock | Pearson r | MAE (yr) | Mean bias (yr) | RMSE (yr) |")
        lines.append("|---|---:|---:|---:|---:|")
        for _, row in payload["summary"].iterrows():
            lines.append(
                f"| {row['Clock']} | {row['Pearson_r']:.3f} | "
                f"{row['MAE_years']:.2f} | {row['Mean_bias_years']:.2f} | "
                f"{row['RMSE_years']:.2f} |"
            )
        lines.append("")
        lines.append(f"Figures: `figures/{label}_correlation_matrix.png`, "
                     f"`figures/{label}_deviation_heatmap.png`, "
                     f"`figures/{label}_age_scatter_grid.png`.")
        lines.append("")
    out.write_text("\n".join(lines))
    print(f"\n[analysis] Wrote {out}")


def main() -> None:
    print("[analysis] Loading datasets ...")
    datasets = load_or_simulate()
    per_dataset: dict[str, dict] = {}
    for label, geo in datasets.items():
        preds, summary = analyse_dataset(label, geo)
        per_dataset[label] = {"preds": preds, "summary": summary}
    write_summary_report(per_dataset)
    print("\n[analysis] Done.")


if __name__ == "__main__":
    main()
