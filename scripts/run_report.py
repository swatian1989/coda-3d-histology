#!/usr/bin/env python
"""Regenerate the full analysis report from cached artefacts.

    python scripts/run_report.py

Idempotent, and it recomputes nothing: every figure and table reads the CSVs
already written by the analysis scripts. Running it twice produces the same
report.

Outputs
    figures/F*.png, figures/F*.pdf     22 figures, 300 dpi raster and vector
    results/tables/T*.csv              14 tables
    reports/analysis_report.md         markdown, relative image links
    reports/analysis_report.html       self contained, images base64 embedded
    reports/analysis_report.docx       navy and steel blue, Calibri, justified

Every figure and table states REAL with its dataset and n, or SIMULATED, or
MISSING DATA naming the input required. The assembler refuses to emit a
report that silently drops a generated figure or table.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coda_my.reporting.report import build_report  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures-dir", default="figures")
    ap.add_argument("--tables-dir", default="results/tables")
    ap.add_argument("--reports-dir", default="reports")
    ap.add_argument("--config", default="config/coda_params.yaml")
    args = ap.parse_args()

    logging.basicConfig(level=logging.WARNING,
                        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger("report")

    out = build_report(figures_dir=args.figures_dir, tables_dir=args.tables_dir,
                       reports_dir=args.reports_dir, config_path=args.config)

    log.warning("figures: %d (%d MISSING-DATA placeholders)",
                len(out["figures"]), out["n_missing_figures"])
    log.warning("tables:  %d (%d MISSING-DATA placeholders)",
                len(out["tables"]), out["n_missing_tables"])
    for kind in ("md", "html", "docx"):
        p = Path(args.reports_dir) / f"analysis_report.{kind}"
        if p.exists():
            log.warning("wrote %s (%.1f KB)", p, p.stat().st_size / 1024)

    real = [f["id"] for f in out["figures"].values() if f["source"].startswith("REAL")]
    log.warning("figures from REAL data: %s", ", ".join(real) if real else "NONE")


if __name__ == "__main__":
    main()
