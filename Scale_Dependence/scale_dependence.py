#!/usr/bin/env python3
"""
Scale-dependence plots from fnlo-tk-cppread outputs.

Per-bin 2D curves (kept), plus per observable & order:
  - 3D surface:  X = mu/mu_0, Y = bin center, Z = differential cross section
  - Heatmap

Usage:
    python scale_dependence.py /path/to/out/files [--out plots]
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict
import math

import numpy as np
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm


# ----------------------------- CLI -------------------------------------------
parser = argparse.ArgumentParser(
    description="Plot scale-factor dependence of LO/NLO/NNLO cross sections."
)
parser.add_argument("indir", type=Path, help="Directory containing *dscale1.txt.out files")
parser.add_argument("--out", type=Path, default=Path("plots"),
                    help="Output directory for the PNGs (default: ./plots)")
args = parser.parse_args()

# Use MathText (no full TeX install required)
mpl.rcParams["text.usetex"] = False

# ----------------------------- Config ----------------------------------------
# Scales to keep on X
S_MIN, S_MAX = 0.05, 20.0
# How many pT^W bins to drop at the low end for the ratio heatmap:
SKIP_PTW_BINS = 4

# Colorbar bottom annotations
CBAR_BOTTOM_ABS   = r"low"
CBAR_BOTTOM_RATIO = r"$<\!1$"

# --------------------------- Containers --------------------------------------
orders = ("LO", "NLO", "NNLO")
data  = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))  # data[obs][bin][order] -> [(scale, val)]
edges = defaultdict(dict)                                           # edges[obs][bin] -> (lo, hi)

OBS_LATEX = {
    "ptj1":    r"$p_T^{j1}$",
    "ptw":     r"$p_T^{W}$",
    "yj1": r"$y^{\,j1}$",
    "abs_yj1": r"$|y^{\,j1}|$",
    "ht_full": r"$H_T$",
    "ht_jets": r"$H_T^{jets}$",
}
OBS_UNIT = {"ptj1": "GeV", "ptw": "GeV", "yj1": "", "abs_yj1": "", "ht_full": "GeV","ht_jets": "GeV"}

# --------------------------- Parsing helpers ---------------------------------
fn_re    = re.compile(r"(?P<prefix>.+?)\.(?P<observable>ptj1|ht_full|ht_jets|yj1|abs_yj1|ptw)\.tab\.gz_[\d.]+dscale1\.txt\.out$")
scale_re = re.compile(r"scale factors xmur, xmuf chosen here are:\s*([0-9.]+)\s*,\s*[0-9.]+")
row_re   = re.compile(r"^\s*(\d+)\s")

def parse_file(fp: Path):
    m = fn_re.search(fp.name)
    if not m:
        raise ValueError(f"Cannot decode filename: {fp.name}")
    observable = m.group("observable")

    rows_current_block = []
    current_scale = None
    with fp.open() as fh:
        for line in fh:
            ms = scale_re.search(line)
            if ms:
                if current_scale is not None and rows_current_block:
                    yield observable, current_scale, rows_current_block
                    rows_current_block = []
                current_scale = float(ms.group(1))
                continue

            if row_re.match(line):
                cols = line.split(maxsplit=10)
                bin_idx   = int(cols[0])
                low_edge  = float(cols[3])
                high_edge = float(cols[4])
                lo, nlo, nnlo = map(float, cols[6:9])
                rows_current_block.append((bin_idx, lo, nlo, nnlo))
                if bin_idx not in edges[observable]:
                    edges[observable][bin_idx] = (low_edge, high_edge)

    if current_scale is not None and rows_current_block:
        yield observable, current_scale, rows_current_block

# --------------------------- Scan directory ----------------------------------
files = sorted(args.indir.glob("*dscale1.txt.out"))
if not files:
    raise SystemExit(f"No *dscale1.txt.out files found in {args.indir}")

for fp in files:
    for observable, scale, rows in parse_file(fp):
        for bin_idx, lo, nlo, nnlo in rows:
            for order, val in zip(orders, (lo, nlo, nnlo)):
                data[observable][bin_idx][order].append((scale, val))

# --------------------------- Utilities ---------------------------------------
def _strip_dollars(s: str) -> str:
    return s[1:-1] if len(s) >= 2 and s[0] == s[-1] == "$" else s

def _common_scales_for_order(bins_for_obs, order, smin=S_MIN, smax=S_MAX):
    """Scales present in ALL bins for this order (ensures rectangular grid)."""
    scale_sets = []
    for _, per_bin in bins_for_obs.items():
        if order not in per_bin:
            return []
        s_here = {s for (s, _) in per_bin[order] if smin <= s <= smax}
        if not s_here:
            return []
        scale_sets.append(s_here)
    return sorted(set.intersection(*scale_sets)) if scale_sets else []

def _surface_matrices(bins_for_obs, order, scales, obs):
    """Return (bin_centers, Z) for given common scales."""
    bin_indices = sorted(bins_for_obs.keys())
    bin_centers = np.array([0.5 * (edges[obs][i][0] + edges[obs][i][1]) for i in bin_indices])
    Z = np.zeros((len(bin_centers), len(scales)))
    for i, idx in enumerate(bin_indices):
        d = dict(bins_for_obs[idx][order])  # {scale: value}
        for j, s in enumerate(scales):
            Z[i, j] = d[s]
    return bin_centers, Z

def _nice_log_ticks(vmin, vmax):
    candidates = [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20]
    ticks = [v for v in candidates if vmin <= v <= vmax]
    return ticks if ticks else [vmin, vmax]

def _nice_log_ticks_z(vmin, vmax):
    """Decade ticks between vmin and vmax (both > 0)."""
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin <= 0 or vmax <= 0:
        return []
    e0 = int(math.floor(math.log10(vmin)))
    e1 = int(math.ceil (math.log10(vmax)))
    return [10.0**e for e in range(e0, e1 + 1)]


# --------------------------- Per-bin 2D plots --------------------------------
args.out.mkdir(parents=True, exist_ok=True)

for obs, bins in data.items():
    for bin_idx, odict in bins.items():
        plt.figure(figsize=(8, 6))
        for order in orders:
            if order not in odict:
                continue
            sx = sorted(odict[order])
            x, y = zip(*[(s, v) for s, v in sx if S_MIN <= s <= S_MAX])
            plt.plot(x, y, linewidth=2, label=order)

        unit = OBS_UNIT.get(obs, "")
        y_label = r"Differential cross section$\ [\mathrm{fb/GeV}]$" if unit == "GeV" \
                  else r"Differential cross section$\ [\mathrm{fb}]$"

        plt.xlabel(r"$\mu/\mu_{0}$", fontsize=20)
        plt.ylabel(y_label, fontsize=20)
        plt.xscale("log")

        ax = plt.gca()
        ax.yaxis.get_offset_text().set_fontsize(16)
        ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

        low, high = edges[obs][bin_idx]
        obs_core = _strip_dollars(OBS_LATEX.get(obs, obs))
        if unit:
            bin_text = fr"${low:g}\,\mathrm{{{unit}}} \leq {obs_core} \leq {high:g}\,\mathrm{{{unit}}}$"
        else:
            bin_text = fr"${low:g} \leq {obs_core} \leq {high:g}$"

        ax.text(0.06, 0.98, bin_text, transform=ax.transAxes,
                ha="left", va="top", fontsize=22, color="#444",
                bbox=dict(facecolor="white", alpha=0.0, edgecolor="none", boxstyle="round,pad=0.25"))

        plt.legend(fontsize=15, frameon=False)
        plt.tick_params(axis="both", which="major", labelsize=18)
        plt.xlim(min(x), max(x))
        plt.grid(True, ls="--", alpha=0.3)
        plt.tight_layout()

        outpng = args.out / f"{obs}_bin{bin_idx:02d}.png"
        plt.savefig(outpng, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"✓ {outpng}")

print(f"Done (per-bin). Plots in {args.out.resolve()}")

# --------------------------- Global color ranges ------------------------------
abs_vmins, abs_vmaxs = [], []
ratio_values = []

for obs, bins_for_obs in data.items():
    if obs == 'ht_jets':
        continue
    for order in orders:
        scales = _common_scales_for_order(bins_for_obs, order, S_MIN, S_MAX)
        if not scales:
            continue
        _, Z = _surface_matrices(bins_for_obs, order, scales, obs)
        Zpos = Z[np.isfinite(Z) & (Z > 0)]
        if Zpos.size:
            abs_vmins.append(Zpos.min())
            abs_vmaxs.append(Zpos.max())

        # Ratios: same selection we plot
        all_bin_indices = sorted(bins_for_obs.keys())
        sel_bin_indices = all_bin_indices[SKIP_PTW_BINS:] if obs == "ptw" else all_bin_indices
        if not sel_bin_indices:
            continue
        Zsel = np.zeros((len(sel_bin_indices), len(scales)))
        for i, idx in enumerate(sel_bin_indices):
            d = dict(bins_for_obs[idx][order])
            for j, s in enumerate(scales):
                Zsel[i, j] = d.get(s, np.nan)
        scales_arr = np.asarray(scales, float)
        cidx = int(np.argmin(np.abs(np.log10(scales_arr))))
        Z0 = Zsel[:, cidx]
        Znorm = np.divide(
            Zsel, Z0[:, None],
            out=np.full_like(Zsel, np.nan, dtype=float),
            where=np.isfinite(Z0[:, None]) & (Z0[:, None] != 0.0)
        )
        ratio_values.extend(Znorm[np.isfinite(Znorm)].ravel())

ABS_VMIN = float(np.min(abs_vmins)) if abs_vmins else 1e-12
ABS_VMAX = float(np.max(abs_vmaxs)) if abs_vmaxs else 1.0

if ratio_values:
    rmin, rmax = float(np.min(ratio_values)), float(np.max(ratio_values))
    delta = max(1 - rmin, rmax - 1)
    RATIO_VMIN = max(1 - delta, 1e-3)
    RATIO_VMAX = 1 + delta
else:
    RATIO_VMIN, RATIO_VMAX = 0.8, 1.2

print(f"Global absolute color range: [{ABS_VMIN:g}, {ABS_VMAX:g}] (log)")
print(f"Global ratio color range:    [{RATIO_VMIN:g}, {RATIO_VMAX:g}] (center=1)")


# --------------------------- 3D surfaces (log color) -------------------------
for obs, bins_for_obs in data.items():
    unit = OBS_UNIT.get(obs, "")
    obs_latex = OBS_LATEX.get(obs, obs)
    ylab_core = _strip_dollars(obs_latex)
    y_label = fr"${ylab_core}\ \mathrm{{bin\ center}}$" + (fr"$\ [\mathrm{{{unit}}}]$" if unit else "")

    zlabel_core = _strip_dollars(obs_latex)
    zlabel = (fr"$\mathrm{{d}}\sigma/\mathrm{{d}}{zlabel_core}\ [\mathrm{{fb/GeV}}]$"
              if unit == "GeV" else fr"$\mathrm{{d}}\sigma/\mathrm{{d}}{zlabel_core}\ [\mathrm{{fb}}]$")

    for order in orders:
        scales = _common_scales_for_order(bins_for_obs, order, S_MIN, S_MAX)
        if not scales:
            continue
        bin_centers, Z = _surface_matrices(bins_for_obs, order, scales, obs)

        Zpos = np.where(Z > 0, Z, np.nan)
        zmin = np.nanmin(Zpos)
        zmax = np.nanmax(Zpos)
        Z_log10 = np.log10(Zpos)

        fig = plt.figure(figsize=(12, 9), constrained_layout=False)
        gs  = fig.add_gridspec(ncols=1, nrows=1)
        fig.subplots_adjust(left=0, right=0.8, top=1, bottom=0.05)
        ax3 = fig.add_subplot(gs[0, 0], projection="3d")

        X_log = np.log10(np.array(scales))
        Xg, Yg = np.meshgrid(X_log, bin_centers)

        cmap_abs = plt.get_cmap("magma")
        norm_abs = LogNorm(vmin=ABS_VMIN, vmax=ABS_VMAX)
        facecolors = cmap_abs(norm_abs(Zpos))
        surf = ax3.plot_surface(
            Xg, Yg, Z_log10, rstride=1, cstride=1,
            facecolors=facecolors, linewidth=0, antialiased=True, shade=False
        )

        xtick_vals = _nice_log_ticks(min(scales), max(scales))
        xtick_pos  = [math.log10(v) for v in xtick_vals]
        ax3.set_xticks(xtick_pos)
        ax3.set_xticklabels([f"{v:g}" for v in xtick_vals], fontsize=14)

        # Show z ticks in original (linear) units
        zticks_vals = _nice_log_ticks_z(zmin, zmax)
        if zticks_vals:
            ax3.set_zticks([math.log10(v) for v in zticks_vals])
            ax3.set_zticklabels([f"{v:g}" for v in zticks_vals], fontsize=14)
            ax3.set_zlim(math.log10(zmin), math.log10(zmax))

        ax3.set_xlabel(r"$\mu/\mu_{0}$", fontsize=22, labelpad=16)
        ax3.set_ylabel(y_label,         fontsize=22, labelpad=16)
        ax3.set_zlabel(zlabel,          fontsize=22, labelpad=26)

        ax3.tick_params(axis='x', which='major', labelsize=14, pad=6)
        ax3.tick_params(axis='y', which='major', labelsize=14, pad=6)
        ax3.tick_params(axis='z', which='major', labelsize=14, pad=12)
        for pane in (ax3.xaxis.pane, ax3.yaxis.pane, ax3.zaxis.pane):
            pane.set_facecolor((1, 1, 1, 0.92))
        for axis in (ax3.xaxis, ax3.yaxis, ax3.zaxis):
            axis._axinfo["grid"]["linestyle"] = "--"
            axis._axinfo["grid"]["color"] = (0, 0, 0, 0.15)
        ax3.view_init(elev=28, azim=-60)
        try:
            ax3.set_box_aspect((1.2, 1.2, 1.2))
        except Exception:
            pass

        # Colorbar from a ScalarMappable using the same norm/cmap
        cax = fig.add_axes([0.85, 0.25, 0.02, 0.5])
        sm  = mpl.cm.ScalarMappable(norm=norm_abs, cmap=cmap_abs)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.locator   = mticker.LogLocator(base=10)
        cbar.formatter = mticker.LogFormatter(base=10)
        cbar.update_ticks()
        cbar.ax.set_ylabel(zlabel, rotation=90, va="center", fontsize=22, labelpad=20)
        cbar.ax.tick_params(labelsize=14)
        #_add_cbar_bottom_label(cbar, CBAR_BOTTOM_ABS)

        fig.suptitle(fr"{obs_latex} – {order}", fontsize=30, y=0.94)
        outpng = args.out / f"{obs}_surface_{order}_3d.png"
        fig.savefig(outpng, dpi=400)
        plt.close(fig)
        print(f"✓ {outpng} (3D surface, log color)")

# --------------------------- Heatmap (ratio to central, shared range) --------
for obs, bins_for_obs in data.items():
    unit = OBS_UNIT.get(obs, "")
    obs_latex = OBS_LATEX.get(obs, obs)
    ylab_core = _strip_dollars(obs_latex)
    y_label = fr"${ylab_core}\ \mathrm{{bin\ center}}$" + (fr"$\ [\mathrm{{{unit}}}]$" if unit else "")

    zlabel_core = _strip_dollars(obs_latex)
    zlabel_ratio = fr"$\left(\frac{{\mathrm{{d}}\sigma}}{{\mathrm{{d}}{zlabel_core}}}\right)\,/\,\left(\frac{{\mathrm{{d}}\sigma}}{{\mathrm{{d}}{zlabel_core}}}\right)_{{\mu=\mu_0}}$"

    for order in orders:
        all_bin_indices = sorted(bins_for_obs.keys())
        sel_bin_indices = all_bin_indices[SKIP_PTW_BINS:] if obs == "ptw" else all_bin_indices
        if not sel_bin_indices:
            continue

        scales = _common_scales_for_order(bins_for_obs, order, S_MIN, S_MAX)
        if not scales:
            continue

        bin_centers = np.array([0.5 * (edges[obs][i][0] + edges[obs][i][1]) for i in sel_bin_indices])
        Z = np.zeros((len(sel_bin_indices), len(scales)))
        for i, idx in enumerate(sel_bin_indices):
            d = dict(bins_for_obs[idx][order])
            for j, s in enumerate(scales):
                Z[i, j] = d.get(s, np.nan)

        scales_arr  = np.asarray(scales, dtype=float)
        central_idx = int(np.argmin(np.abs(np.log10(scales_arr))))
        Z0          = Z[:, central_idx]
        Z_norm = np.divide(
            Z, Z0[:, None],
            out=np.full_like(Z, np.nan, dtype=float),
            where=np.isfinite(Z0[:, None]) & (Z0[:, None] != 0.0)
        )

        fig = plt.figure(figsize=(6, 6), constrained_layout=True)
        gs  = fig.add_gridspec(ncols=2, nrows=1, width_ratios=[1.0, 0.05])

        ax  = fig.add_subplot(gs[0, 0])
        cax = fig.add_subplot(gs[0, 1])

        vmax_for_norm = max(RATIO_VMAX, 1.0001)      # ensure vcenter < vmax
        norm_ratio = TwoSlopeNorm(vmin=0.5, vcenter=1.0, vmax=vmax_for_norm)  # no 'clip' kwarg
        cmap_ratio = mpl.cm.get_cmap("RdBu_r")
        cmap_ratio.set_under(cmap_ratio(0.0))        # values <0.5 use the lowest displayed color

        # (swap axes):
        X, Y = np.meshgrid(bin_centers, scales_arr)      # x = bin centers, y = scales
        im = ax.pcolormesh(X, Y, Z_norm.T, shading="auto",
                        cmap=cmap_ratio, norm=norm_ratio)

        ax.set_yscale("log")
        ax.set_xlabel(y_label,              fontsize=22, labelpad=6)
        ax.set_ylabel(r"$\mu/\mu_{0}$",   fontsize=22, labelpad=6)
        ax.tick_params(axis="both", which="major", labelsize=17)  # adjust size 

        if obs in ["ht_full","ht_jets"]:
            # Limit how many x ticks are shown and hide the last label
            ax.set_xlim(bin_centers.min(), bin_centers.max())
            ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=3, prune="upper"))  # ~6 ticks total
            ax.xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:.0f}"))      # integers
            ax.tick_params(axis="x", which="major", labelsize=22, pad=2)


        # Colorbar: arrow (triangle) at the bottom and include a real 0.5 tick
        cbar = fig.colorbar(im, cax=cax, extend="min")
        cbar.ax.set_ylabel(zlabel_ratio, rotation=90, va="center", fontsize=22, labelpad=16)
        cbar.ax.tick_params(labelsize=17)

        # ---- ticks and labels: 0.5, 1, 2, 3, 4, 5 (show decimals only for 0.5) ----
        eps = 1e-9
        ticks = [0.5] + [i for i in range(1, 6) if i <= vmax_for_norm + eps]
        cbar.set_ticks(ticks)

        cbar.formatter = mticker.FuncFormatter(
            lambda x, pos: ("0.5" if abs(x - 0.5) < 1e-9
                            else str(int(round(x))) if abs(x - round(x)) < 1e-9
                            else f"{x:g}"))
        cbar.update_ticks()
        # --------------------------------------------------------------------------

        outpng = args.out / f"{obs}_surface_{order}_heatmap_ratio_to_central.png"
        fig.savefig(outpng, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"✓ {outpng} (heatmap ratio, shared global range; "
              f"{'skipped first '+str(SKIP_PTW_BINS)+' ptw bins' if obs=='ptw' else 'all bins'})")
