#!/usr/bin/env python3
# Visualizations for MAT2007 Project - Statistical Analysis of Root-Hormone: Effect of Phytohormones on Plant 
#                                      Root Lengths and the Correlations of Change


import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONDITIONS = ["Control", "IAA", "ABA", "CK"]
HORMONES = CONDITIONS[1:]
TRL_COL = {"Control": "TRL_control", "IAA": "TRL_IAA", "ABA": "TRL_ABA", "CK": "TRL_CK"}
RATIO_COL = {"IAA": "ratio_IAA", "ABA": "ratio_ABA", "CK": "ratio_CK"}

# Use the same colour for each condition throughout the figures
COLORS = {"Control": "#8C8C8C", "IAA": "#E69F00", "ABA": "#0072B2", "CK": "#009E73"}

SE_COLOR = "black"   # Standard error bars
SUB_COLOR = "red"    # Subsampling uncertainty bars
SIGMA_THRESHOLD = 3.0



# LOADING THE DATA
def load_data(folder):
    names = {"summary": "results_summary.csv", "acc": "ratios_per_accession.csv",
             "pair": "pairwise_comparisons.csv", "corr": "correlations.csv"}
    data = {}
    for key, name in names.items():
        path = os.path.join(folder, name)
        if not os.path.exists(path):
            sys.exit(f"Missing {path}. Run the C++ programme first (in the same folder).")
        data[key] = pd.read_csv(path)

    needed = {"summary": ["mean_se", "sigma_from_1", "sigma_from_1_se"],
              "pair": ["combined_se", "sigma", "sigma_se"],
              "corr": ["pearson_r_se", "pearson_r_full"]}
    for key, cols in needed.items():
        for col in cols:
            if col not in data[key].columns:
                sys.exit(f"Column '{col}' not found in {names[key]}. "
                         "Use hormone_analysis_se.cpp to create the CSV files.")
    return data


def get_row(summary, quantity, condition):
    sel = summary[(summary["quantity"] == quantity) & (summary["condition"] == condition)]
    if sel.empty:
        sys.exit(f"No row for {quantity} / {condition} in results_summary.csv")
    return sel.iloc[0]


def save(fig, outdir, name):
    """Save a figure and retry if Windows temporarily blocks the file"""
    fig.tight_layout()
    path = os.path.join(os.path.abspath(outdir), name + ".png")
    for attempt in range(1, 6):
        try:
            fig.savefig(path, dpi=150, bbox_inches="tight")
            break
        except OSError as err:
            if attempt == 5:
                sys.exit(f"Could not write {path}: {err}\n"
                         "Try running the script from a folder that is not synced by OneDrive (for example C:\\MAT2007).")
            print(f"  could not write {name}.png ({err}), retrying ({attempt}/5)...")
            time.sleep(1)
    plt.close(fig)
    print("  saved", name + ".png")


def sigma_text(sig, signed=True):
    """Format a sigma value and add a * when it is greater than the chosen threshold"""
    text = f"{sig:+.1f} σ" if signed else f"{abs(sig):.1f} σ"
    return text + (" *" if abs(sig) > SIGMA_THRESHOLD else "")


# PART 1 - make a histogram for each condition and show its mean with both error bars
def fig_condition(cond, data, outdir):
    trl = data["acc"][TRL_COL[cond]].dropna()
    row = get_row(data["summary"], "TRL_mm", cond)

    fig, ax = plt.subplots(figsize=(6, 4))
    counts, _, _ = ax.hist(trl, bins=15, color=COLORS[cond], edgecolor="black")
    top = counts.max()

    # Show the mean as a dot and place the two error bars at different heights to keep them separate
    ax.errorbar(row["mean_full"], top * 1.10, xerr=row["mean_se"], fmt="o", color=SE_COLOR,
                capsize=5, label=f"Standard error: {row['mean_full']:.2f} ± {row['mean_se']:.2f}")
    ax.errorbar(row["mean"], top * 1.25, xerr=row["mean_unc"], fmt="o", color=SUB_COLOR,
                capsize=5, label=f"Subsample uncertainty: {row['mean']:.2f} ± {row['mean_unc']:.2f}")
    ax.set_ylim(0, top * 1.40)

    ax.set_xlabel("Total root length (mm)")
    ax.set_ylabel("Number of accessions")
    ax.set_title(f"{cond} (N = {int(row['N'])})")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), fontsize=8, frameon=False)
    save(fig, outdir, f"fig1_{cond}")


# PART 2 - figures comparing the different conditions and hormone responses
def two_panel_bars(labels, se_vals, se_errs, sub_vals, sub_errs, ylabel, title, name, outdir, digits,
                   ref=None, se_sigmas=None, sub_sigmas=None):
    """Compare the means in two panels, using standard errors on the left and subsampling uncertainty on the right
    Each condition keeps its own colour, and optional sigma values are shown above the bars"""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.3), sharey=True)
    panels = [(axes[0], se_vals, se_errs, se_sigmas, SE_COLOR, "Mean ± standard error"),
              (axes[1], sub_vals, sub_errs, sub_sigmas, SUB_COLOR, "Mean ± subsample uncertainty")]
    top = max(max(v + e for v, e in zip(se_vals, se_errs)), max(v + e for v, e in zip(sub_vals, sub_errs)))
    colors = [COLORS[label] for label in labels]

    for ax, vals, errs, sigmas, color, ptitle in panels:
        ax.bar(labels, vals, color=colors, edgecolor="black")
        ax.errorbar(labels, vals, yerr=errs, fmt="none", ecolor=color, elinewidth=2, capsize=6)
        for i, (v, e) in enumerate(zip(vals, errs)):
            text = f"{v:.{digits}f} ± {e:.{digits}f}"
            if sigmas is not None:
                text += "\n" + sigma_text(sigmas[i])
            ax.annotate(text, (i, v + e), xytext=(0, 4), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))
        if ref is not None:
            ax.axhline(ref, color="grey", linestyle="--")
        ax.set_ylim(0, top * 1.25)
        ax.set_title(ptitle)
    axes[0].set_ylabel(ylabel)
    fig.suptitle(title)
    if se_sigmas is not None:
        fig.text(0.5, -0.02, "σ = distance from the dashed line in units of the uncertainty shown in that panel.   "
                 f"* = |σ| > {SIGMA_THRESHOLD:.0f}", ha="center", fontsize=8, color="#444444")
    save(fig, outdir, name)


def two_panel_dots(labels, face_colors, edge_colors, se_vals, se_errs, sub_vals, sub_errs, xlabel, title,
                   name, outdir, ref, se_sigmas, sub_sigmas, signed, digits=3):
    """Compare values and their horizontal error bars in two panels
    The left panel uses standard errors and the right panel uses subsampling uncertainty
    Each dot uses the first hormone's colour inside and the second hormone's colour around the edge"""
    n = len(labels)
    fig, axes = plt.subplots(1, 2, figsize=(10, 2.4 + 0.9 * n), sharex=True, sharey=True)
    panels = [(axes[0], se_vals, se_errs, se_sigmas, SE_COLOR, "Value ± standard error"),
              (axes[1], sub_vals, sub_errs, sub_sigmas, SUB_COLOR, "Value ± subsample uncertainty")]

    # Set shared x-axis limits with enough space for the labels
    lo = min(min(np.array(v) - np.array(e)) for v, e in [(se_vals, se_errs), (sub_vals, sub_errs)] + [([ref], [0])])
    hi = max(max(np.array(v) + np.array(e)) for v, e in [(se_vals, se_errs), (sub_vals, sub_errs)] + [([ref], [0])])
    pad = 0.25 * (hi - lo)

    for ax, vals, errs, sigmas, color, ptitle in panels:
        for i in range(n):
            ax.errorbar(vals[i], i, xerr=errs[i], fmt="o", ms=10, mfc=face_colors[i], mec=edge_colors[i], mew=2.5,
                        ecolor=color, elinewidth=2, capsize=6, zorder=3)
            ax.annotate(f"{vals[i]:.{digits}f} ± {errs[i]:.{digits}f}\n" + sigma_text(sigmas[i], signed),
                        (vals[i], i), xytext=(0, 10), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))
        ax.axvline(ref, color="grey", linestyle="--", zorder=1)
        ax.set_yticks(range(n))
        ax.set_yticklabels(labels)
        ax.set_ylim(n - 0.5, -0.9)  # Keep the first label at the top and leave room for its annotation
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_xlabel(xlabel)
        ax.set_title(ptitle)
    fig.suptitle(title)
    fig.text(0.5, -0.02, "Dot: first hormone = fill colour, second hormone = edge colour.   σ = distance from the dashed line "
             f"in units of the uncertainty shown in that panel.   * = |σ| > {SIGMA_THRESHOLD:.0f}",
             ha="center", fontsize=8, color="#444444")
    save(fig, outdir, name)


def fig_correlation_scatter(data, outdir):
    """Create one scatter plot for each hormone pair to compare their responses across accessions"""
    acc, corr = data["acc"], data["corr"]
    n = len(corr)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5))
    axes = np.atleast_1d(axes)

    for ax, (_, r) in zip(axes, corr.iterrows()):
        a, b = r["hormone_1"], r["hormone_2"]
        d = acc[[RATIO_COL[a], RATIO_COL[b]]].dropna()
        x, y = d.iloc[:, 0].to_numpy(), d.iloc[:, 1].to_numpy()

        # Use the x-axis hormone for the dot fill and the y-axis hormone for its outline
        ax.scatter(x, y, s=75, facecolor=COLORS[a], edgecolor=COLORS[b], linewidth=1.8, alpha=0.9, zorder=3)
        slope, intercept = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, slope * xs + intercept, color="black", linewidth=1.2, zorder=4)
        ax.axvline(1, color="grey", linestyle=":", linewidth=1.2)  # A ratio of 1 means no effect on the x-axis
        ax.axhline(1, color="grey", linestyle=":", linewidth=1.2)  # A ratio of 1 means no effect on the y-axis

        sig_se = r["pearson_r_full"] / r["pearson_r_se"]
        sig_sub = r["pearson_r"] / r["pearson_r_unc"]
        ax.text(0.03, 0.97,
                f"r = {r['pearson_r_full']:.2f} ± {r['pearson_r_se']:.2f} (SE)  {sig_se:.1f} σ\n"
                f"r = {r['pearson_r']:.2f} ± {r['pearson_r_unc']:.2f} (subsample)  {sig_sub:.1f} σ",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox=dict(boxstyle="round", fc="white", ec="#CCCCCC"))
        ax.set_xlabel(f"{a} / control")
        ax.set_ylabel(f"{b} / control")
        ax.set_title(f"{a} vs {b} (N = {int(r['N'])})")

    fig.suptitle("Do accessions respond similarly to different hormones?")
    fig.text(0.5, -0.02, "Dots: fill = x-axis hormone, edge = y-axis hormone.   Black line = linear fit.   "
             "Dotted lines = ratio of 1 (no effect).   σ = r / its uncertainty (distance from r = 0).   "
             "All ratios share the same control value, which adds some positive correlation.",
             ha="center", fontsize=8, color="#444444", wrap=True)
    save(fig, outdir, "fig6_correlation_scatter")


def main():
    parser = argparse.ArgumentParser(description="Plots for the plant hormone root-length analysis")
    parser.add_argument("--dir", default=".", help="folder with the CSV files written by the C++ programme")
    parser.add_argument("--out", default="figures", help="output folder for the figures")
    args = parser.parse_args()

    data = load_data(args.dir)
    os.makedirs(args.out, exist_ok=True)
    summary, pair, corr = data["summary"], data["pair"], data["corr"]

    print("PART 1 - figures per condition")
    for cond in CONDITIONS:
        fig_condition(cond, data, args.out)

    print("PART 2 - overall figures")

    # Compare the mean total root length for the control and all three hormones
    rows = [get_row(summary, "TRL_mm", c) for c in CONDITIONS]
    two_panel_bars(CONDITIONS,
                   [r["mean_full"] for r in rows], [r["mean_se"] for r in rows],
                   [r["mean"] for r in rows], [r["mean_unc"] for r in rows],
                   "Total root length (mm)", "Mean total root length", "fig2_mean_root_length", args.out, digits=2)

    # Compare relative root length under each hormone with the control, where a ratio of 1 means no effect
    rows = [get_row(summary, "ratio_to_control", h) for h in HORMONES]
    two_panel_bars(HORMONES,
                   [r["mean_full"] for r in rows], [r["mean_se"] for r in rows],
                   [r["mean"] for r in rows], [r["mean_unc"] for r in rows],
                   "Relative root length (hormone / control)",
                   "Relative root length (dashed line = no effect)", "fig3_relative_root_length", args.out,
                   digits=3, ref=1.0,
                   se_sigmas=[r["sigma_from_1_se"] for r in rows], sub_sigmas=[r["sigma_from_1"] for r in rows])

    # Make scatter plots to compare the responses to pairs of hormones
    fig_correlation_scatter(data, args.out)

    print(f"\nDone. Figures are in '{args.out}/'")


if __name__ == "__main__":
    main()
