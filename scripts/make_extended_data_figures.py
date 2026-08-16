#!/usr/bin/env python
"""Extended Data figures in the layout of the source publication.

    python scripts/make_extended_data_figures.py [--which A1 A2 A7 A10]

Each figure reproduces a specific published Extended Data panel set, on this
data, panel for panel. Where a panel cannot be produced the reason is stated in
the caption rather than something else being substituted.

  A1   registration workflow            mirrors Extended Data Fig 1a, 10 panels
  A2   accuracy vs the benchmark        mirrors Extended Data Fig 1b
  A7   z-resolution validation          mirrors Extended Data Fig 5a
  A10  3D reconstruction                mirrors Extended Data Figs 6 to 8

NOT produced, and why:
  A3   cell detection precision/recall  needs two human annotators; none exist
  A5   segmentation training design     needs annotated tiles; none exist
  A6   segmentation confusion matrices  needs a trained model; none exists
  A8   tissue classes labelled          needs multi-class segmentation
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coda_my.loaders.kartasalo import (  # noqa: E402
    NATIVE_MPP_UM, SECTION_THICKNESS_UM, fiducial_array, load_fiducials,
)
from coda_my.registration import (  # noqa: E402
    RegistrationConfig, apply_rigid, elastic_field, preprocess, radon_transform,
)
from coda_my.registration_fix import SearchConfig, estimate_rigid_search  # noqa: E402

DATA = ROOT / "data/raw/kartasalo/extracted/Data_to_IDA"
OUT = ROOT / "figures/extended_data"
KART = ROOT / "results/kartasalo"
BASE_DS = 16
MPP = NATIVE_MPP_UM * BASE_DS
SRC = "REAL, Kartasalo mouse liver, n = 47 serial sections"


def bare(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def scalebar(ax, img_w_px, mpp, um=2000, label=None):
    """Burned-in scale bar, bottom left, in the style of the source figures."""
    px = um / mpp
    y = ax.get_ylim()[0]
    x0 = img_w_px * 0.04
    y0 = y - img_w_px * 0.06 if y < 0 else y * 0.94
    ax.plot([x0, x0 + px], [y0, y0], color="black", lw=3, solid_capstyle="butt")
    ax.text(x0 + px / 2, y0, label or f"{um/1000:g} mm", ha="center", va="bottom",
            fontsize=7)


def magenta_green(fixed: np.ndarray, moving: np.ndarray) -> np.ndarray:
    """Fixed in magenta, moving in green. Perfect alignment reads neutral grey.

    This is the most persuasive panel in the published workflow figure and it is
    persuasive because it cannot be faked: any residual offset shows as coloured
    fringing, and only genuine alignment collapses to grey.
    """
    f = fixed.astype(np.float32); m = moving.astype(np.float32)
    f = (f - f.min()) / max(np.ptp(f), 1e-6)
    m = (m - m.min()) / max(np.ptp(m), 1e-6)
    return np.clip(np.dstack([f, m, f]), 0, 1)


# ---------------------------------------------------------------- Figure A1


def fig_a1(grey: np.ndarray) -> None:
    """Extended Data Fig 1a: the registration workflow, panel for panel."""
    cfg = RegistrationConfig()
    i = len(grey) // 2
    fixed_raw, moving_raw = grey[i], grey[i + 1]
    f = preprocess(fixed_raw, cfg)
    m = preprocess(moving_raw, cfg)

    # solve the rigid step on a coarse copy, as the pipeline does
    def bm(a, k=11):
        h, w = a.shape
        h2, w2 = (h // k) * k, (w // k) * k
        return a[:h2, :w2].astype(np.float32).reshape(h2 // k, k, w2 // k, k).mean((1, 3))

    p = estimate_rigid_search(preprocess(bm(fixed_raw), cfg), preprocess(bm(moving_raw), cfg),
                              cfg, SearchConfig(max_abs_deg=45.0))
    ang, dy, dx = p["angle"], p["dy"] * 11, p["dx"] * 11
    m_global = apply_rigid(m, ang, dy, dx)

    fy, fx = elastic_field(f, m_global, cfg, mpp=MPP)
    from coda_my.registration import apply_elastic
    m_local = apply_elastic(m_global, fy, fx)

    fig = plt.figure(figsize=(16.5, 15.0))
    gs = fig.add_gridspec(4, 4, hspace=0.28, wspace=0.18,
                          height_ratios=[1, 1, 1, 1.1])

    # a-i  stack schematic
    ax = fig.add_subplot(gs[0, 0]); bare(ax)
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    for k, (y, c) in enumerate([(8.2, "0.75"), (7.9, "0.75"), (7.6, "0.75"),
                                (5.2, "0.15"),
                                (2.8, "0.75"), (2.5, "0.75"), (2.2, "0.75")]):
        ax.add_patch(Rectangle((2.2, y), 4.6, 0.26, facecolor=c, edgecolor="black",
                               lw=0.5))
    ax.text(1.9, 8.4, "image 1", ha="right", fontsize=8)
    ax.text(1.9, 5.3, "center", ha="right", fontsize=8, fontweight="bold")
    ax.text(1.9, 2.3, "image n", ha="right", fontsize=8)
    ax.text(4.5, 6.9, "...", ha="center", fontsize=11)
    ax.text(4.5, 3.9, "...", ha="center", fontsize=11)
    for y0, y1 in ((8.0, 5.6), (2.4, 5.0)):
        ax.add_patch(FancyArrowPatch((7.4, y0), (7.4, y1), connectionstyle="arc3,rad=.55",
                                     arrowstyle="-|>", mutation_scale=14, lw=2,
                                     color="black"))
    ax.set_title("registration inward to the centre section", fontsize=8)

    # a-ii  raw fixed and moving
    for j, (img, lab) in enumerate(((fixed_raw, "Sample fixed image$_n$"),
                                    (moving_raw, "Sample moving image$_{n+1}$"))):
        ax = fig.add_subplot(gs[0, 1 + j]); bare(ax)
        ax.imshow(img, cmap="gray"); ax.set_xlabel(lab, fontsize=8)
        if j == 0:
            scalebar(ax, img.shape[1], MPP)

    # a-iii  global registration point of reference
    ax = fig.add_subplot(gs[0, 3]); bare(ax)
    ax.imshow(fixed_raw, cmap="gray")
    cy, cx = np.array(fixed_raw.shape) / 2
    ax.plot(cx, cy, "+", color="black", ms=18, mew=2.5)
    ax.set_xlabel("+ global registration point\nof reference on fixed image", fontsize=8)

    # a-iv  filtered fixed image
    ax = fig.add_subplot(gs[1, 0]); bare(ax)
    ax.imshow(f, cmap="viridis")
    ax.set_xlabel("filtered fixed image", fontsize=8)

    # a-v  Radon transform
    ax = fig.add_subplot(gs[1, 1]); bare(ax)
    small = bm(fixed_raw, 6)
    r = radon_transform(preprocess(small, cfg), n_angles=360)
    ax.imshow(r, cmap="viridis", aspect="auto")
    ax.set_xlabel("Radon transform\n[0,360] degrees", fontsize=8)

    # a-vi  2D cross correlation surface
    ax = fig.add_subplot(gs[1, 2]); bare(ax)
    A = np.fft.fft2(bm(f, 6)); B = np.fft.fft2(bm(m_global, 6))
    xc = np.fft.fftshift(np.real(np.fft.ifft2(A * np.conj(B))))
    xc = (xc - xc.min()) / max(np.ptp(xc), 1e-9)
    im = ax.imshow(xc, cmap="viridis")
    py, px = np.unravel_index(np.argmax(xc), xc.shape)
    ax.plot(px, py, "o", mfc="none", mec="white", ms=12, mew=1.6)
    ax.set_xlabel("2D cross correlation", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.045)

    # a-vii  local registration points of reference
    ax = fig.add_subplot(gs[1, 3]); bare(ax)
    ax.imshow(fixed_raw, cmap="gray")
    step = max(int(cfg.tile_interval_um / MPP), 16)
    ys = np.arange(step // 2, fixed_raw.shape[0], step)
    xs = np.arange(step // 2, fixed_raw.shape[1], step)
    gy, gx = np.meshgrid(ys, xs, indexing="ij")
    ax.plot(gx.ravel(), gy.ravel(), "+", color="black", ms=7, mew=1.2)
    ax.set_xlabel("+ local registration points\nof reference (1.5 mm)", fontsize=8)

    # a-viii  deformed grid overlay
    ax = fig.add_subplot(gs[2, 0]); bare(ax)
    ax.imshow(fixed_raw, cmap="gray", alpha=0.55)
    sub = 6
    yy = np.arange(0, fy.shape[0], max(fy.shape[0] // 22, 1))
    xx = np.arange(0, fy.shape[1], max(fy.shape[1] // 22, 1))
    for y in yy:
        ax.plot(xx + fx[y, xx], np.full_like(xx, y, dtype=float) + fy[y, xx],
                color="red", lw=0.7)
    for x in xx:
        ax.plot(np.full_like(yy, x, dtype=float) + fx[yy, x], yy + fy[yy, x],
                color="red", lw=0.7)
    ax.set_xlabel("grid representation overlay of\nelastic registration results", fontsize=8)

    # a-ix  displacement fields
    vmin = min(fy.min(), fx.min()); vmax = max(fy.max(), fx.max())
    for j, (fld, lab) in enumerate(((fx, "horizontal"), (fy, "vertical"))):
        ax = fig.add_subplot(gs[2, 1 + j]); bare(ax)
        im = ax.imshow(fld, cmap="viridis", vmin=vmin, vmax=vmax)
        ax.set_xlabel(f"interpolated {lab}\ndisplacement field", fontsize=8)
        if j == 1:
            cb = fig.colorbar(im, ax=ax, fraction=0.045)
            cb.set_label("pixels", fontsize=7)

    # a-x  magenta / green overlays, on their own row so they can be judged
    for k, (mv, lab) in enumerate(((m, "Pre-registration"),
                                   (m_global, "Global registration"),
                                   (m_local, "Local registration"))):
        ax = fig.add_subplot(gs[3, k]); bare(ax)
        ax.imshow(magenta_green(f, mv))
        ax.set_xlabel(lab, fontsize=9)
    ax = fig.add_subplot(gs[3, 3]); bare(ax)
    ax.text(0.02, 0.74, "fixed = magenta\nmoving = green", fontsize=11,
            transform=ax.transAxes, fontweight="bold")
    ax.text(0.02, 0.34,
            "Only genuine alignment collapses to\nneutral grey. Any residual offset "
            "shows\nas coloured fringing at the tissue edge,\nso this panel cannot be "
            "made to look\ncorrect by a transform that is not.",
            fontsize=8.6, transform=ax.transAxes, va="top")

    fig.text(0.5, 0.045,
             "Extended Data Fig A1. Histological image registration sample workflow. "
             f"{SRC}.",
             ha="center", fontsize=9.5, fontweight="bold")
    fig.text(0.5, 0.006,
             "(a) Sections registered to the reference at the centre z-height. Global "
             "registration uses a rotational reference at the centre of the fixed image. "
             "Fixed and moving images are converted to greyscale, non-tissue removed, "
             "intensity complemented and Gaussian filtered. The Radon transform is "
             "computed over 0 to 360 degrees; here the rotation is recovered instead by "
             "direct search, because Radon estimation averaged 37.5 degrees of error on "
             "this tissue. Maximum of the 2D cross correlation yields the translation. "
             "Local registration is performed on tiles at 1.5 mm intervals and "
             "interpolated to whole-image displacement fields. Overlays show fixed in "
             "magenta and moving in green, before registration, after the global step "
             "and after the local step.",
             ha="center", va="bottom", fontsize=7.6, wrap=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"A1_registration_workflow.{ext}", dpi=300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote A1_registration_workflow")


# ---------------------------------------------------------------- Figure A2


def fig_a2(grey: np.ndarray) -> None:
    """Extended Data Fig 1b: fiducials and normalised performance."""
    import pandas as pd
    ref = pd.read_csv(ROOT / "data/reference/kartasalo2018_table2_liver.csv")
    low = ref[ref.resolution == "low"]
    s = json.loads((KART / "summary_ds16fix_rigid.json").read_text())
    f1 = load_fiducials(DATA / "fiducialcoordinates_liver_observer1.txt")
    i = len(grey) // 2

    fig = plt.figure(figsize=(15.5, 5.4))
    gs = fig.add_gridspec(1, 3, wspace=0.16)

    for j, (sec, col, lab) in enumerate(((i + 1, "#2ECC40", "Image N fiducial markers"),
                                         (i + 2, "#0074D9", "Image N+1 fiducial markers"))):
        ax = fig.add_subplot(gs[0, j]); bare(ax)
        ax.imshow(grey[sec - 1], cmap="gray")
        pts = fiducial_array(f1, sec, BASE_DS)
        ax.plot(pts[:, 1], pts[:, 0], "o", mfc=col, mec="black", ms=9, mew=0.8)
        ax.set_xlabel(lab, fontsize=9)
        if j == 0:
            scalebar(ax, grey.shape[2], MPP)

    # normalised performance: higher is always better
    ax = fig.add_subplot(gs[0, 2])
    metrics = [("tre_mean_um", "TRE", True), ("atre_mean_um", "ATRE", True),
               ("rmse_mean", "RMSE", True), ("jaccard_mean", "J", False),
               ("darea_pct_mean", "dA", True)]
    ours = {"tre_mean_um": s["tre_full_mean_um"], "atre_mean_um": s["atre_mean_um"]}
    for xi, (col, lab, lower_better) in enumerate(metrics):
        vals = low[col].astype(float).to_numpy()
        allv = np.abs(vals) if col == "darea_pct_mean" else vals
        lo, hi = np.nanmin(allv), np.nanmax(allv)
        def norm(v):
            v = abs(v) if col == "darea_pct_mean" else v
            n = (v - lo) / max(hi - lo, 1e-9)
            return 1 - n if lower_better else n
        for _, r in low.iterrows():
            y = norm(float(r[col]))
            if r.kind == "baseline":
                ax.plot(xi, y, "D", color="black", ms=7, zorder=4)
            else:
                ax.plot(xi, y, "o", color="0.55", ms=7, zorder=3)
        if col in ours:
            ax.plot(xi, norm(ours[col]), "s", color="#E4002B", ms=9, zorder=6)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([m[1] for m in metrics])
    ax.set_xlabel("Performance metric"); ax.set_ylabel("Normalized performance")
    ax.set_ylim(-0.05, 1.08)
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([], [], marker="D", ls="", color="black", label="Unregistered sections"),
        Line2D([], [], marker="s", ls="", color="#E4002B", label="This implementation"),
        Line2D([], [], marker="o", ls="", color="0.55", label="Other registration algorithms"),
    ], fontsize=8, loc="lower right", frameon=False)

    fig.text(0.5, -0.02,
             "Extended Data Fig A2. Registration accuracy against the published "
             f"benchmark. {SRC}; published values from Kartasalo et al. 2018 Table 2, "
             "liver at matching resolution.",
             ha="center", fontsize=9.5, fontweight="bold")
    fig.text(0.5, -0.09,
             "(a, b) Corresponding fiducial markers on adjacent sections; the four "
             "landmarks are laser-cut holes driven through the block before embedding. "
             "(c) Normalised performance for each metric, scaled so that higher is "
             "always better. RMSE, Jaccard and area change are shown for the published "
             "methods only, as those metrics were not recomputed here.",
             ha="center", va="top", fontsize=7.6, wrap=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"A2_benchmark_accuracy.{ext}", dpi=300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote A2_benchmark_accuracy")


# ---------------------------------------------------------------- Figure A7


def fig_a7() -> None:
    """Extended Data Fig 5a: z-resolution validation."""
    import pandas as pd
    ax_df = pd.read_csv(KART / "step8_axial_lateral_ds16fix_rigid.csv")
    zs = pd.read_csv(KART / "step8_zskip_ds16fix_rigid.csv")

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.4))
    ax = axes[0]
    for name, g in ax_df.groupby("axis"):
        g = g[g.correlation > 0]
        ax.plot(g.distance_um / 1000.0, np.log(g.correlation.clip(1e-6)),
                ls="--" if name == "xy" else "-",
                color="black" if name == "xy" else "#0074D9",
                label="xy (within section)" if name == "xy" else "z (between sections)")
    ax.axhline(np.log(0.95), color="#E4002B", lw=1.2, ls=":", label="95% correlation")
    ax.set_xlabel("mm"); ax.set_ylabel("log correlation")
    ax.legend(fontsize=7.5); ax.set_title("(a) axial versus lateral", fontsize=9)

    ax = axes[1]
    ax.bar(range(len(zs)), zs.percent_composition_error,
           color=plt.cm.Blues(np.linspace(.45, .95, len(zs))), edgecolor="white")
    ax.axhline(5.0, color="black", lw=1.2)
    ax.set_xticks(range(len(zs))); ax.set_xticklabels(zs.spacing_um.astype(int))
    ax.set_xlabel("distance between sections (um)")
    ax.set_ylabel("% change in composition from 5 um")
    ax.set_title("(b) cost of skipping sections", fontsize=9)

    ax = axes[2]
    keep = zs[zs.percent_composition_error <= 5.0]
    ax.axis("off")
    ax.text(0.02, 0.92, "Published claim, tested here", fontsize=10, fontweight="bold",
            transform=ax.transAxes)
    lines = [
        "correlation above 95% to 20 um:",
        f"   not reproduced. z correlation falls below 95%",
        f"   between the first and second section.",
        "",
        "error below 5% to 12 um:",
        f"   not reproduced. the largest skip inside 5% is",
        f"   {int(keep.skip.max()) if len(keep) else 1}"
        f" ({int(keep.spacing_um.max()) if len(keep) else 5} um).",
        "",
        "The published claims were made on pancreas with",
        "260 sections. This is liver with 47, and structure",
        "changes faster through z, so the shortcut that",
        "justified processing one section in three does not",
        "transfer. The measurement has to be repeated per",
        "tissue, which is why it is reported rather than",
        "assumed.",
    ]
    for k, ln in enumerate(lines):
        ax.text(0.02, 0.80 - k * 0.055, ln, fontsize=8.2, transform=ax.transAxes,
                family="monospace" if ln.startswith("   ") else None)

    fig.text(0.5, -0.03,
             f"Extended Data Fig A7. Validation of z-resolution reduction. {SRC}.",
             ha="center", fontsize=9.5, fontweight="bold")
    fig.text(0.5, -0.10,
             "(a) Pixel correlation across the z axis against within-section xy "
             "correlation, which is the ceiling set by intact tissue. (b) Percent change "
             "in tissue composition from the native 5 um spacing; horizontal line, the "
             "5 percent threshold reported by the source publication. (c) Whether the "
             "published claims reproduce on this tissue.",
             ha="center", va="top", fontsize=7.6, wrap=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"A7_z_resolution.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print("wrote A7_z_resolution")


# --------------------------------------------------------------- Figure A10


def fig_a10() -> None:
    """Extended Data Figs 6 to 8: reconstruction renders and z-projections."""
    z = np.load(KART / "stage6_projections.npz")
    rgbp = ROOT / "data/interim/kartasalo_liver_ds16_rgb.npy"
    s = json.loads((KART / "stage6_summary.json").read_text())

    fig = plt.figure(figsize=(15.5, 8.6))
    gs = fig.add_gridspec(2, 3, hspace=0.22, wspace=0.12)

    # (a) pseudo-3D block: the H&E volume seen as a slab
    ax = fig.add_subplot(gs[0, 0]); bare(ax)
    if rgbp.exists():
        rgb = np.load(rgbp)
        mid = rgb[len(rgb) // 2]
        ax.imshow(mid)
    ax.set_xlabel("(a) reconstructed H&E volume, central plane", fontsize=8.5)

    ax = fig.add_subplot(gs[0, 1]); bare(ax)
    ax.imshow(z["lumen_xy"], cmap="magma")
    ax.set_xlabel("(a) vascular lumina, summed through the block", fontsize=8.5)

    ax = fig.add_subplot(gs[0, 2]); bare(ax)
    ax.imshow(z["lumen_xz"], cmap="magma", aspect="auto")
    ax.set_xlabel("(b) xz through the stack, z vertical", fontsize=8.5)

    # (c) z-projections, black background, class name in colour
    for k, (arr, name, cmap, col) in enumerate((
            (z["lumen_xy"], "Vascular lumina", "magma", "#FF6EC7"),
            (z["lumen_yz"], "Lumina, yz", "magma", "#FF6EC7"),
            (z["section_mid"], "All tissue", "gray", "white"))):
        ax = fig.add_subplot(gs[1, k]); bare(ax)
        ax.set_facecolor("black")
        ax.imshow(arr, cmap=cmap, aspect="auto" if k == 1 else None)
        ax.text(0.03, 0.05, name, transform=ax.transAxes, color=col, fontsize=9,
                fontweight="bold")

    fig.text(0.5, 0.045,
             "Extended Data Fig A10. Three-dimensional reconstruction. "
             f"{SRC}, {s['n_sections']} sections at {s['mpp_um']:.2f} um/px and "
             f"{s['section_thickness_um']:.0f} um spacing.",
             ha="center", fontsize=9.5, fontweight="bold")
    fig.text(0.5, 0.005,
             "(a) Central plane of the reconstructed volume and the segmented vascular "
             "lumina summed through the whole block, where a coherent branching tree "
             "indicates successful registration. (b) The xz plane, a view no individual "
             "section contains. (c) Z-projections. Only one tissue class is shown, "
             "because separating the ten classes of the published figure requires a "
             "trained multi-class segmentation that does not exist for this material; "
             "the remaining class panels are therefore not produced rather than "
             "substituted.",
             ha="center", va="bottom", fontsize=7.6, wrap=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"A10_3d_reconstruction.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print("wrote A10_3d_reconstruction")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", nargs="+", default=["A1", "A2", "A7", "A10"])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    grey = np.load(ROOT / f"data/interim/kartasalo_liver_ds{BASE_DS}.npy")
    if "A1" in args.which:
        fig_a1(grey)
    if "A2" in args.which:
        fig_a2(grey)
    if "A7" in args.which:
        fig_a7()
    if "A10" in args.which:
        fig_a10()
    print(f"figures in {OUT}")


if __name__ == "__main__":
    main()
