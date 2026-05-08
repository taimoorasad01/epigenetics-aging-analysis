"""
Real-data Loader (run on your own machine - not in sandboxed environments)
==========================================================================
Downloads the canonical GEO methylation aging datasets used in the
Bio-Learn paper (de Lima Camillo et al., 2023) via bio-learn's DataLibrary.

Datasets:
  - GSE40279  : Hannum et al. (2013), 656 whole-blood samples, Illumina 450K
                'Genome-wide methylation profiles reveal quantitative views
                 of human aging rates'  (Mol Cell)
  - GSE41169  : Horvath et al. (2012) Dutch population blood cohort, 95
                samples (62 schizophrenia + 33 controls), 450K.

Both are 450K but >90 % of probes overlap with EPIC, and every clock used
in this analysis is platform-agnostic at that overlap.

Requirements (install on your own machine):
    pip install biolearn pandas numpy

Usage:
    python -m part2_epic_aging_clocks.src.download_real_data

Note: bio-learn caches datasets to ~/.biolearn/cache/ (~500 MB combined).
"""
from __future__ import annotations

import sys
from pathlib import Path


def load_real_datasets() -> dict[str, "GeoData"]:
    """Return {dataset_id: GeoData}.

    Performs full download on first call; subsequent calls hit local cache.
    """
    from biolearn.data_library import DataLibrary  # noqa: WPS433

    library = DataLibrary()
    datasets: dict[str, object] = {}
    for gse_id in ["GSE40279", "GSE41169"]:
        print(f"[download_real_data] Loading {gse_id} ...", flush=True)
        datasets[gse_id] = library.get(gse_id).load()
        print(
            f"  -> {gse_id}: {datasets[gse_id].dnam.shape[1]} samples x "
            f"{datasets[gse_id].dnam.shape[0]} CpGs",
            flush=True,
        )
    return datasets


if __name__ == "__main__":
    try:
        datasets = load_real_datasets()
    except Exception as exc:  # pragma: no cover - network failures
        print(f"[ERROR] Could not download real GEO data: {exc}", file=sys.stderr)
        print(
            "Hint: this script needs unrestricted internet access to "
            "https://ftp.ncbi.nlm.nih.gov. If running in a sandbox, fall "
            "back to the simulator (simulate_data.py).",
            file=sys.stderr,
        )
        sys.exit(1)

    out = Path(__file__).resolve().parent.parent / "data"
    out.mkdir(exist_ok=True, parents=True)
    for gse_id, geo in datasets.items():
        geo.save_csv(str(out / gse_id))
        print(f"[download_real_data] Cached {gse_id} -> {out}/", flush=True)
