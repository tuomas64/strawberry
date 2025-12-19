#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Direct percentile confidence intervals for MSMC-IM median curves
#
# 4×2 figure:
#   Rows 1–3: Ne (left) and m (right) with direct percentile CIs + median
#   Row 4:    Histograms (G,H) for isolation midpoint CSVs
#
# Bottom-row CSV globs via:
#   HIST_G_GLOB=".../*.csv"
#   HIST_H_GLOB=".../*.csv"
#
# Bottom-row axes & bins:
#   HIST_XMIN=10000 HIST_XMAX=420000 HIST_BINS=30 HIST_UNITS=auto
#   HIST_LOGX=0  # 0=linear, 1=log x-axis
#
# Save histogram bin counts:
#   HIST_SAVE_TABLE=1 HIST_TABLE_OUT="histogram_tables_row4.csv"
#
# Styling / content rules (per your latest requests):
#   - Smaller axis/tick fonts
#   - Keep only: panel letters + sample sizes + legend
#   - Legend labels simplified: "95% CI", "75% CI", "Median"
#   - Panel letters:
#       left column  A, C, E, G
#       right column B, D, F, H
#     placed outside top-left of each axes
#   - Bottom row: show MIS2, PGP, MIS8, MIS10 inside histogram panels (like earlier)
# -----------------------------------------------------------------------------

import os, hashlib
from glob import glob
from collections import Counter
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, MultipleLocator, FuncFormatter
import matplotlib.font_manager as fm

# ============================
# User-tunable global options
# ============================
VERBOSE = True
DISABLE_QC = (os.getenv("DISABLE_QC", "0") == "0")  # keep existing behavior

# Smoothing (publish-friendly visualization)
SMOOTH_CURVES  = True
SMOOTH_METHOD  = "pchip"        # "pchip" (preferred) or "moving_avg"
SMOOTH_POINTS  = 80
MOVAVG_WINDOW  = 7

# Optional: faint individual replicate lines (“spaghetti”)
SHOW_SPAGHETTI = False
SPAGHETTI_ALPHA = 0.12
SPAGHETTI_LW = 0.7
HSPACE = 0.12

# Direct percentile CI levels across replicates
PERCENTILE_CI_95 = (2.5, 97.5)
PERCENTILE_CI_75 = (12.5, 87.5)

# Separate CI spaces: Ne=linear (default), m=log (default)
NE_CI_SPACE = os.environ.get("NE_CI_SPACE", "linear").strip().lower()  # "linear" or "log"
M_CI_SPACE  = os.environ.get("M_CI_SPACE",  "log").strip().lower()     # "linear" or "log"

# Visual style
PERC_CI_COLOR  = "#059669"
CI95_STYLE, CI95_LW = "-", 1.2
CI75_STYLE, CI75_LW = ":", 1.0
MEDIAN_COLOR = "black"

# ── LR04 overlay config ───────────────────────────────────────────────────────
LR04_PATH   = os.environ.get(
    "CLIMATE_FILE",
    "/scratch/project_2001113/MSMC-IM/MSMC-IM-output/spain_italy/prop/Lisiecki_and_Raymo_2005-Data.txt"
)
LR04_LINE_KW = dict(lw=1.0, alpha=0.7, color="0.25")
LR04_LABEL   = r"$-\delta^{18}\\mathrm{O}$"
LR04_TICKS   = [-5.0, -4.5, -4.0, -3.5, -3.0]
LR04_YLIMS   = (-5.0, -3.0)

# ── paths & config ────────────────────────────────────────────────────────────
LOCAL_ROOT = os.environ.get("MSMC_LOCAL_ROOT", "").rstrip("/")
def remap(p):
    if LOCAL_ROOT and p and p.startswith("/scratch/"):
        return os.path.join(LOCAL_ROOT, p.split("/scratch/", 1)[1])
    return p

dirs_top  = [remap(p) for p in ["/scratch/project_2001113/MSMC-IM/CORE_collected"]]
dirs_mis6 = [remap(p) for p in ["/scratch/project_2001113/MSMC-IM/PGP_collected"]]
dirs_mis8 = [remap(p) for p in ["/scratch/project_2001113/MSMC-IM/MIS8_collected"]]

# ── global plotting ranges & constants ────────────────────────────────────────
GEN=2
X_MIN,X_MAX=1_000,375_000
M_Y_MIN,M_Y_MAX=1e-7,5e-3
N1_Y_MIN,N1_Y_MAX=7e2,1e5
major_x  = [k*1_000 for k in (1,4,10,20,50,100,150,250,350)]
major_xt = [str(k) for k in (1,4,10,20,50,100,150,250,350)]

# Custom N1 y-axis ticks shown as ×1000-style labels (1, 5, 10, 20, 50)
N1_YTICKS = [1_000, 5_000, 10_000, 20_000, 50_000, 80_000]
N1_YTICKLABELS = ["1", "5", "10", "20", "50", "80"]

# NE QC thresholds (applied only for times >= CUTOFF_YEARS)
LO_NE, HI_NE = 10, 200_000
CUTOFF_YEARS = 14_000

# ── climate background bands (all panels) ─────────────────────────────────────
rectangles=[(s*1_000,e*1_000,c,a) for s,e,c,a in [
 (1,4,"red",.05),(4,9,"red",.20),(9,14,"red",.05),
 (14,17,"blue",.05),(17,22,"blue",.20),(22,29,"blue",.05),
 (29,59,"red",.05),(59,71,"blue",.05),(71,84,"red",.05),(84,93,"blue",.05),
 (93,108,"red",.05),(108,115,"blue",.20),(115,130,"red",.05),
 (130,190,"blue",.20),(190,247,"red",.05),(247,300,"blue",.20),
 (300,337,"red",.05),(337,374,"blue",.20)]]

def add_climate_labels(ax):
    # Used for A–F (top 3 rows). Places labels above the axes (outside).
    for s,e,label,color in [
        (4_000, 9_000,"HTM","red"), (17_000,22_000,"LGM","blue"),
        (29_000,59_000,"3","red"), (59_000,71_000,"4","blue"),
        (108_000,115_000,"5d","blue"), (115_000,130_000,"5e","red"),
        (130_000,190_000,"PGP","blue"), (190_000,247_000,"7","red"),
        (247_000,300_000,"8","blue"), (300_000,337_000,"9","red"),
        (337_000,374_000,"10","blue"),
    ]:
        if e < X_MIN or s > X_MAX: continue
        xmid = max(X_MIN, min(X_MAX, 0.5*(s+e)))
        ax.text(xmid, 1.06, label, color=color, transform=ax.get_xaxis_transform(),
                fontsize=7, fontweight='bold', ha='center', va='bottom',
                bbox=dict(fc='white', ec='none', alpha=1.0, pad=0.8))

# ── style defaults (smaller fonts) ────────────────────────────────────────────
plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "Arial" if any("Arial" in f.name for f in fm.fontManager.ttflist) else "DejaVu Sans",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
})
plt.rc("text", usetex=False)

M_COLS  = ("m","M","mig","migration","migration_rate","m_rate")
N1_COLS = ("im_N1","N1","ne1","Ne1","N1_eff","N1_effective")
N2_COLS = ("im_N2","N2","ne2","Ne2","N2_eff","N2_effective")

# ── helpers ──────────────────────────────────────────────────────────────────
def pick_col(df, cands):
    for c in cands:
        if c in df.columns: return c
    return None

def qc_ne_reasons(df, lo=LO_NE, hi=HI_NE, cutoff=CUTOFF_YEARS):
    if DISABLE_QC: return []
    reasons=[]
    c1 = pick_col(df, N1_COLS); c2 = pick_col(df, N2_COLS)
    if c1 is None: reasons.append("missing_im_N1")
    if c2 is None: reasons.append("missing_im_N2")
    if reasons: return reasons

    t_years = df["left_time_boundary"].to_numpy(dtype=float) * GEN
    mask_old = np.isfinite(t_years) & (t_years >= cutoff)

    v1 = df[c1].to_numpy(dtype=float); v2 = df[c2].to_numpy(dtype=float)
    v1_old = v1[mask_old]; v2_old = v2[mask_old]

    if v1_old.size and not np.all(np.isfinite(v1_old)): reasons.append("nonfinite_im_N1_ge14ka")
    if v2_old.size and not np.all(np.isfinite(v2_old)): reasons.append("nonfinite_im_N2_ge14ka")

    if v1_old.size and np.nanmin(v1_old) < lo: reasons.append("im_N1_lt_10_ge14ka")
    if v1_old.size and np.nanmax(v1_old) > hi: reasons.append("im_N1_gt_200000_ge14ka")
    if v2_old.size and np.nanmin(v2_old) < lo: reasons.append("im_N2_lt_10_ge14ka")
    if v2_old.size and np.nanmax(v2_old) > hi: reasons.append("im_N2_gt_200000_ge14ka")
    return reasons

EXCLUDED_RUNS=[]
def record_exclusion(path, reasons, df):
    EXCLUDED_RUNS.append({"path": path, "reasons": reasons})
    if VERBOSE:
        print(f"[QC-EXCLUDE ≥14ka] {os.path.basename(path)} :: {', '.join(reasons)}")

def collect_series(roots, col_cands):
    out, all_t, seen = [], [], set()
    roots = roots if isinstance(roots,(list,tuple)) else [roots]
    total=used=qc_fail=miss_col=empty_mask=0
    for r in roots:
        if not os.path.exists(r):
            if VERBOSE: print(f"[collect] root not found: {r}")
            continue
        for f in glob(os.path.join(r,"**","*.estimates.txt"), recursive=True):
            total += 1
            try:
                df = pd.read_table(f, comment="#")
            except Exception as e:
                if VERBOSE: print(f"[collect] read fail: {f} :: {e}")
                continue
            if "left_time_boundary" not in df.columns:
                if VERBOSE: print(f"[collect] no left_time_boundary: {f}")
                continue

            reasons = qc_ne_reasons(df)
            if reasons:
                qc_fail += 1
                record_exclusion(f, reasons, df)
                continue

            col = pick_col(df, col_cands)
            if not col:
                miss_col += 1
                continue

            t = df["left_time_boundary"].to_numpy()*GEN
            v = df[col].to_numpy()
            msk = (t>=X_MIN) & (t<=X_MAX) & np.isfinite(v) & (v>0)
            if not msk.any():
                empty_mask += 1
                continue

            buf = np.column_stack((t[msk],v[msk])).tobytes()
            md5 = hashlib.md5(buf).hexdigest()
            if md5 in seen:
                continue
            seen.add(md5)
            s = pd.Series(v[msk], index=t[msk]).sort_index()
            out.append(s); all_t.extend(s.index.tolist()); used += 1

    if VERBOSE:
        print(f"[collect] roots={len(roots)} files={total} used={used} "
              f"qc_fail={qc_fail} miss_col={miss_col} empty_mask={empty_mask}")
    return out, all_t

# ── percentile summary ────────────────────────────────────────────────────────
def _union_grid(all_times):
    ut = np.array(sorted({float(t) for t in all_times
                          if np.isfinite(t) and X_MIN <= t <= X_MAX}))
    return ut

def _matrix_for_space(series_list, ut, space):
    mats = []
    for s in series_list:
        y = pd.Series(s).reindex(ut).ffill().bfill().to_numpy()
        y = np.where((~np.isfinite(y)) | (y <= 0), np.nan, y)
        y = pd.Series(y, index=ut).ffill().bfill().to_numpy()
        if space == "log":
            y = np.log10(y)
        mats.append(y)
    return np.vstack(mats) if mats else np.empty((0, 0))

def _summary(series_list, all_times, space):
    if not series_list or not all_times:
        if VERBOSE: print("[summary] no series or times")
        return None

    ut = _union_grid(all_times)
    if ut.size < 2:
        if VERBOSE: print(f"[summary] union grid too small: ut.size={ut.size}")
        return None

    mat = _matrix_for_space(series_list, ut, space=space)
    if mat.size == 0:
        return None

    def p(pct):
        arr = np.percentile(mat, pct, axis=0, method="linear")
        return (10 ** arr) if space == "log" else arr

    res = {
        "x": ut.astype(float),
        "med": p(50.0),
        "B": mat.shape[0],
        "_space": space,
        "perc_lo95": p(PERCENTILE_CI_95[0]),
        "perc_hi95": p(PERCENTILE_CI_95[1]),
        "perc_lo75": p(PERCENTILE_CI_75[0]),
        "perc_hi75": p(PERCENTILE_CI_75[1]),
    }
    return res

# ── smoothing helpers ────────────────────────────────────────────────────────
try:
    from scipy.interpolate import PchipInterpolator as _PCHIP
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False
    _PCHIP = None

def _logspace_dense(xmin, xmax, n=500):
    xmin=float(xmin); xmax=float(xmax)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmin<=0 or xmax<=xmin: return None
    return np.logspace(np.log10(xmin), np.log10(xmax), int(n))

def _smooth_y_on_xdense_loglog(x, y, x_dense, method="pchip"):
    x=np.asarray(x,float); y=np.asarray(y,float)
    m = np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0); x=x[m]; y=y[m]
    if x.size<2: return None
    if method=="pchip" and _HAVE_SCIPY:
        f=_PCHIP(np.log10(x), np.log10(y))
        ylog=f(np.log10(x_dense))
        return 10**ylog
    ylog=np.interp(np.log10(x_dense), np.log10(x), np.log10(y))
    k=int(MOVAVG_WINDOW) if int(MOVAVG_WINDOW)>0 else 5
    if k%2==0: k+=1
    if k>1:
        pad=k//2; ypad=np.pad(ylog, pad_width=pad, mode="edge"); ker=np.ones(k)/k
        ylog=np.convolve(ypad, ker, mode="valid")
    return 10**ylog

def smooth_summary_loglog(S, n_points=SMOOTH_POINTS, method=SMOOTH_METHOD):
    if S is None or "x" not in S or S["x"].size<2: return S
    x=np.asarray(S["x"],float); x=x[np.isfinite(x)&(x>0)]
    if x.size<2: return S
    xd=_logspace_dense(np.nanmin(x), np.nanmax(x), n=n_points)
    if xd is None: return S
    out={"x":xd, "B":S.get("B",0), "_space":S.get("_space","")}
    for key,val in S.items():
        if key in ("x","B","_space"): continue
        ys=_smooth_y_on_xdense_loglog(S["x"], val, xd, method=method)
        if ys is None: return S
        out[key]=ys
    return out

# ── LR04 loader & overlay ────────────────────────────────────────────────────
def load_lr2005(path):
    try:
        df = pd.read_table(path, comment="#")
    except Exception as e:
        print(f"[LR04] read fail: {e} :: {path}")
        return None
    age_col = next((c for c in df.columns if c.lower().startswith("age") and "ka" in c.lower()), None)
    d18_col = next((c for c in df.columns if "d18" in c.lower()), None)
    if age_col is None: age_col = "age_ka" if "age_ka" in df.columns else None
    if d18_col is None:
        d18_col = "bent_d18O" if "bent_d18O" in df.columns else ("d18O" if "d18O" in df.columns else None)
    if age_col is None or d18_col is None:
        print(f"[LR04] Could not detect columns. Headers: {list(df.columns)}")
        return None
    years = pd.to_numeric(df[age_col], errors="coerce").to_numpy(float) * 1_000.0
    d18   = pd.to_numeric(df[d18_col], errors="coerce").to_numpy(float)
    m = np.isfinite(years) & np.isfinite(d18)
    if not np.any(m): return None
    years, d18 = years[m], d18[m]
    o = np.argsort(years); years, d18 = years[o], d18[o]
    keep = (years >= X_MIN) & (years <= X_MAX)
    years, d18 = years[keep], d18[keep]
    if years.size == 0: return None
    return {"years": years, "bent_d18O": d18}

def overlay_lr04_on_ax_like_old(ax_left, iso_label=False):
    iso = load_lr2005(LR04_PATH)
    if iso is None:
        print("[LR04] no data; skip overlay")
        return None
    years = np.asarray(iso["years"], float)
    yvals = -np.asarray(iso["bent_d18O"], float)
    ax2 = ax_left.twinx()
    ax2.set_xscale('log'); ax2.set_yscale('linear'); ax2.set_xlim(ax_left.get_xlim())
    ax2.xaxis.set_major_locator(FixedLocator(major_x))
    ax2.xaxis.set_major_formatter(FixedFormatter(major_xt))
    ax2.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax2.patch.set_alpha(0.0); ax2.set_zorder(ax_left.get_zorder()+1)
    ax2.plot(years, yvals, **LR04_LINE_KW, zorder=10)
    ax2.set_autoscale_on(False)
    ax2.set_ylim(*LR04_YLIMS)
    ax2.yaxis.set_major_locator(FixedLocator(LR04_TICKS))
    ax2.yaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in LR04_TICKS]))
    ax2.tick_params(axis='y', labelsize=8, colors="0.25")
    ax2.set_ylabel(LR04_LABEL if iso_label else "", fontweight='bold', color="0.25", labelpad=10)
    for sp in ax2.spines.values(): sp.set_linewidth(1.2)
    return ax2

# ── gather data (top 3 rows) ─────────────────────────────────────────────────
def load_stage(name, dirs):
    m,tm=collect_series(dirs,M_COLS); n,tn=collect_series(dirs,N1_COLS)
    m_sum=_summary(m,tm, space=M_CI_SPACE); n_sum=_summary(n,tn, space=NE_CI_SPACE)
    return {"name":name, "m":m, "n":n, "m_sum":m_sum, "n_sum":n_sum}

core={}
core["m"], tm = collect_series(dirs_top, M_COLS)
core["n"], tn = collect_series(dirs_top, N1_COLS)
core["m_sum"] = _summary(core["m"], tm, space=M_CI_SPACE)
core["n_sum"] = _summary(core["n"], tn, space=NE_CI_SPACE)

stage_mis6 = load_stage("MIS6", dirs_mis6)
stage_mis8 = load_stage("MIS8", dirs_mis8)

# ── plotting helpers for top 3 rows ──────────────────────────────────────────
def prep_ax(ax, xlim, ylim, xlabel=False, ylabel=None, ytick_locs=None, ytick_labels=None):
    for s,e,c,a in rectangles:
        if e>=xlim[0] and s<=xlim[1]:
            ax.add_patch(Rectangle((s, ylim[0]), e - s, ylim[1] - ylim[0],
                       facecolor=c, alpha=a,
                       edgecolor="none", linewidth=0,
                       zorder=0))

    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(FixedLocator(major_x)); ax.xaxis.set_major_formatter(FixedFormatter(major_xt))
    if ytick_locs is not None:
        ax.yaxis.set_major_locator(FixedLocator(ytick_locs))
        if ytick_labels is not None:
            ax.yaxis.set_major_formatter(FixedFormatter(ytick_labels))
    if ylabel: ax.set_ylabel(ylabel, fontweight='bold')
    if xlabel: ax.set_xlabel('Years ago (×1000)', fontweight='bold')
    ax.grid(True, which='major', ls='--', alpha=.3)
    for sp in ax.spines.values(): sp.set_linewidth(1.4); sp.set_zorder(15)

def plot_stage_row_ci(ax_left_ne, ax_right_m, stage, add_legend=False, show_n_on_left=True):
    prep_ax(ax_left_ne,(X_MIN,X_MAX),(N1_Y_MIN,N1_Y_MAX),
            ylabel="Effective population size x 1000",
            ytick_locs=N1_YTICKS, ytick_labels=N1_YTICKLABELS)
    prep_ax(ax_right_m,(X_MIN,X_MAX),(M_Y_MIN,M_Y_MAX),
            ylabel="Migration rate (m)")

    n_files = max(len(stage.get("n", [])), len(stage.get("m", [])))

    def draw(ax, S, which):
        if S is None or S["x"].size < 2:
            if VERBOSE:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=9, color="0.4")
            return

        Sd = smooth_summary_loglog(S) if SMOOTH_CURVES else S

        if SHOW_SPAGHETTI:
            series_list = stage["n"] if which == "n" else stage["m"]
            for s in series_list:
                xk = np.asarray(s.index, float); yk = np.asarray(s.values, float)
                m = np.isfinite(xk) & np.isfinite(yk) & (xk>0) & (yk>0)
                xk, yk = xk[m], yk[m]
                if xk.size < 2: continue
                if SMOOTH_CURVES:
                    xd = np.logspace(np.log10(xk.min()), np.log10(xk.max()), 300)
                    ys = _smooth_y_on_xdense_loglog(xk, yk, xd, method=SMOOTH_METHOD)
                    if ys is None: continue
                    ax.plot(xd, ys, lw=SPAGHETTI_LW, alpha=SPAGHETTI_ALPHA, color="0.2", zorder=1)
                else:
                    ax.plot(xk, yk, lw=SPAGHETTI_LW, alpha=SPAGHETTI_ALPHA, color="0.2",
                            drawstyle="steps-post", zorder=1)

        x = Sd["x"]

        # 95% CI
        ax.plot(x, Sd["perc_lo95"], linestyle=CI95_STYLE, color=PERC_CI_COLOR, lw=CI95_LW, zorder=2.2)
        ax.plot(x, Sd["perc_hi95"], linestyle=CI95_STYLE, color=PERC_CI_COLOR, lw=CI95_LW, zorder=2.2)
        # 75% CI
        ax.plot(x, Sd["perc_lo75"], linestyle=CI75_STYLE, color=PERC_CI_COLOR, lw=CI75_LW, zorder=2.3)
        ax.plot(x, Sd["perc_hi75"], linestyle=CI75_STYLE, color=PERC_CI_COLOR, lw=CI75_LW, zorder=2.3)
        # Median
        ax.plot(x, Sd["med"], color=MEDIAN_COLOR, lw=1.8, zorder=4)

    draw(ax_left_ne, stage["n_sum"], which="n")
    draw(ax_right_m, stage["m_sum"], which="m")

    # Sample size ONLY (kept)
    if show_n_on_left:
        ax_left_ne.text(0.98, 0.98, f"n={n_files}",
                        transform=ax_left_ne.transAxes,
                        ha="right", va="top",
                        fontsize=8, color="0.25")

    # Simplified legend ONLY (kept)
    if add_legend:
        handles = [
            Line2D([0],[0], color=PERC_CI_COLOR, lw=CI95_LW, linestyle=CI95_STYLE, label="95% CI"),
            Line2D([0],[0], color=PERC_CI_COLOR, lw=CI75_LW, linestyle=CI75_STYLE, label="75% CI"),
            Line2D([0],[0], color=MEDIAN_COLOR, lw=1.8, label="Median"),
        ]
        ax_right_m.legend(handles=handles, loc="upper right",
                          frameon=False, fontsize=8)

# ── bottom-row histogram inputs (row 4: panels G, H) ──────────────────────────
HIST_G_GLOB = os.environ.get("HIST_G_GLOB", "").strip()
HIST_H_GLOB = os.environ.get("HIST_H_GLOB", "").strip()
HIST_XMIN = int(os.environ.get("HIST_XMIN", "10000"))
HIST_XMAX = int(os.environ.get("HIST_XMAX", "420000"))
HIST_BINS = int(os.environ.get("HIST_BINS", "30"))
HIST_UNITS = os.environ.get("HIST_UNITS", "auto").strip().lower()  # 'auto'|'kya'|'years'
HIST_G_COMBINE = bool(int(os.environ.get("HIST_G_COMBINE", "0")))
HIST_H_COMBINE = bool(int(os.environ.get("HIST_H_COMBINE", "1")))
HIST_LOGX = bool(int(os.environ.get("HIST_LOGX", "0")))

# Save histogram tables
HIST_SAVE_TABLE = bool(int(os.environ.get("HIST_SAVE_TABLE", "1")))
HIST_TABLE_OUT  = os.environ.get("HIST_TABLE_OUT", "histogram_tables_row4.csv")

# Bottom-row in-panel MIS labels (requested)
SELECTED_BOTTOM_LABELS = [
    ("MIS2",  14_000, 29_000, "blue"),
    ("PGP",  130_000,190_000, "blue"),
    ("MIS8", 247_000,300_000, "blue"),
    ("MIS10",337_000,374_000, "blue"),
]

def _place_hist_labels_inside(ax, xmin: int, xmax: int):
    for name, s, e, color in SELECTED_BOTTOM_LABELS:
        if e < xmin or s > xmax:
            continue
        xmid = max(xmin, min(xmax, 0.5 * (s + e)))
        ax.text(
            xmid, 0.985, name,
            transform=ax.get_xaxis_transform(),  # x=data, y=axes fraction
            ha="center", va="top",
            fontsize=9, fontweight="bold", color=color,
            bbox=dict(fc="white", ec="none", alpha=0.9, pad=0.6),
            zorder=10,
        )

def _detect_units_from_kya(series_kya: pd.Series) -> str:
    s = pd.to_numeric(series_kya, errors='coerce').dropna()
    if s.empty: return 'kya'
    return 'years' if float(s.quantile(0.95)) > 10000 else 'kya'

def _load_midpoints_years(csv_path: str, units: str = 'auto', gen_time: float | None = None, verbose=VERBOSE) -> pd.Series:
    df = pd.read_csv(csv_path)
    if 'midpoint_kya' in df.columns:
        col_units = units
        if units == 'auto':
            col_units = _detect_units_from_kya(df['midpoint_kya'])
            if verbose: print(f"[midpoints] {os.path.basename(csv_path)} midpoint_kya detected as {col_units}")
        vals = pd.to_numeric(df['midpoint_kya'], errors='coerce')
        return (vals * 1000.0) if col_units == 'kya' else vals
    if 'midpoint_gen' in df.columns:
        if not gen_time:
            raise RuntimeError("[midpoints] midpoint_gen present but no generation time provided")
        vals = pd.to_numeric(df['midpoint_gen'], errors='coerce')
        return vals * float(gen_time)
    raise RuntimeError(f"[midpoints] No midpoint_kya or midpoint_gen in {csv_path}")

def _shade_mis_like_main(ax, x_min: int, x_max: int):
    for s,e,c,a in rectangles:
        if e >= x_min and s <= x_max:
            ax.add_patch(Rectangle((s, 0), e - s, 1,
                       transform=ax.get_xaxis_transform(),
                       facecolor=c, alpha=a,
                       edgecolor="none", linewidth=0,
                       zorder=0))


def _hist_table(series_list, files, xmin, xmax, bins, logx, panel_label):
    if logx:
        edges = np.logspace(np.log10(max(xmin, 1)), np.log10(xmax), bins + 1)
    else:
        edges = np.linspace(xmin, xmax, bins + 1)

    out = pd.DataFrame({
        "panel": panel_label,
        "bin_left_years": edges[:-1],
        "bin_right_years": edges[1:],
        "bin_mid_years": 0.5 * (edges[:-1] + edges[1:]),
        "bin_left_kya": edges[:-1] / 1000.0,
        "bin_right_kya": edges[1:] / 1000.0,
        "bin_mid_kya": 0.5 * (edges[:-1] + edges[1:]) / 1000.0,
    })

    all_combined = (pd.concat(series_list, ignore_index=True)
                    if series_list else pd.Series(dtype=float))

    counts_total, _ = np.histogram(all_combined.to_numpy(dtype=float), bins=edges)
    out["count_total"] = counts_total

    for s_in, f in zip(series_list, files):
        name = os.path.splitext(os.path.basename(f))[0]
        c, _ = np.histogram(s_in.to_numpy(dtype=float), bins=edges)
        out[f"count__{name}"] = c

    tot = int(out["count_total"].sum())
    out["prop_total"] = out["count_total"] / tot if tot > 0 else 0.0
    return out

def draw_midpoint_histogram(ax, csv_glob: str, xmin: int, xmax: int, bins: int,
                           units: str = 'auto', panel_letter: str = "",
                           color_cycle=None, combine: bool=False, logx: bool = HIST_LOGX):
    if not csv_glob:
        ax.text(0.5, 0.5, "No CSV glob", transform=ax.transAxes,
                ha="center", va="center", color="0.4")
        return None
    files = sorted(glob(csv_glob))
    if not files:
        ax.text(0.5, 0.5, "No files matched", transform=ax.transAxes,
                ha="center", va="center", color="0.4")
        return None

    series_list, totals = [], []
    for p in files:
        s = _load_midpoints_years(p, units=units, gen_time=None, verbose=VERBOSE).dropna()
        totals.append(int(s.shape[0]))
        s_in = s[(s >= xmin) & (s <= xmax)]
        series_list.append(s_in)

    # bins for plotting
    if logx:
        bin_edges = np.logspace(np.log10(max(xmin, 1)), np.log10(xmax), bins + 1)
        use_range = False
    else:
        bin_edges = None
        use_range = True

    if combine:
        all_combined = (pd.concat(series_list, ignore_index=True)
                        if series_list else pd.Series(dtype=float))
        if use_range:
            counts, _, _ = ax.hist(
                all_combined, bins=bins, range=(xmin, xmax),
                alpha=0.55, color="#2E86AB", edgecolor="black"
            )
        else:
            counts, _, _ = ax.hist(
                all_combined, bins=bin_edges,
                alpha=0.55, color="#2E86AB", edgecolor="black"
            )
        max_count = counts.max() if counts.size else 1
        drops_in = int(all_combined.shape[0])
    else:
        if not color_cycle:
            color_cycle = ["#2E86AB", "#D1495B", "#7A5195", "#FFA600", "#003F5C"]
        max_count = 1
        drops_in = 0
        for i, s_in in enumerate(series_list):
            if s_in.empty: continue
            if use_range:
                counts, _, _ = ax.hist(
                    s_in, bins=bins, range=(xmin, xmax),
                    alpha=0.55, color=color_cycle[i % len(color_cycle)],
                    edgecolor="black"
                )
            else:
                counts, _, _ = ax.hist(
                    s_in, bins=bin_edges,
                    alpha=0.55, color=color_cycle[i % len(color_cycle)],
                    edgecolor="black"
                )
            drops_in += len(s_in)
            if counts.size:
                max_count = max(max_count, int(np.max(counts)))

    ax.set_ylim(0, max_count * 1.12)
    ax.set_xlim(xmin, xmax)

    if logx:
        ax.set_xscale('log')
        ax.xaxis.set_major_locator(FixedLocator(major_x))
        ax.xaxis.set_major_formatter(FixedFormatter(major_xt))
        ax.grid(axis="x", which="both", alpha=0.18)
    else:
        step = 50_000 if (xmax - xmin) <= 300_000 else 100_000
        ax.xaxis.set_major_locator(MultipleLocator(step))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda t, pos: f"{int(t/1000)}"))

    ax.set_xlabel("Time (ka)", fontweight='bold')
    ax.set_ylabel("Count", fontweight='bold')
    ax.grid(axis="y", alpha=0.28, linestyle="--")

    _shade_mis_like_main(ax, xmin, xmax)
    _place_hist_labels_inside(ax, xmin, xmax)  # requested (MIS2/PGP/MIS8/MIS10)

    # Sample size (drops in the plotted window) ONLY (kept)
    ax.text(0.98, 0.98, f"n={drops_in}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="0.25")

    # return table
    try:
        return _hist_table(series_list, files, xmin, xmax, bins, logx, panel_letter)
    except Exception as e:
        if VERBOSE:
            print(f"[hist-table] Could not build table for {panel_letter}: {e}")
        return None

# ── main 4×2 figure ─────────────────────────────────────────────────────────
from matplotlib.gridspec import GridSpec

# --- Figure size: width = 180 mm ---
W_MM = 180
W_IN = W_MM / 25.4

old_w, old_h = 14, 16.6      # your previous figsize
H_IN = W_IN * (old_h / old_w)

fig = plt.figure(figsize=(W_IN, H_IN))

gs = fig.add_gridspec(nrows=4, ncols=2, height_ratios=[1.0, 1.0, 1.0, 1.0], hspace=HSPACE)

# Axes
ax00 = fig.add_subplot(gs[0,0]); ax01 = fig.add_subplot(gs[0,1])
ax10 = fig.add_subplot(gs[1,0]); ax11 = fig.add_subplot(gs[1,1])
ax20 = fig.add_subplot(gs[2,0]); ax21 = fig.add_subplot(gs[2,1])
ax30 = fig.add_subplot(gs[3,0]); ax31 = fig.add_subplot(gs[3,1])

# Row 1: core (legend here only)
plot_stage_row_ci(
    ax00, ax01,
    {"n":core["n"],"m":core["m"],"n_sum":core["n_sum"],"m_sum":core["m_sum"]},
    add_legend=True,
    show_n_on_left=True
)
add_climate_labels(ax00); add_climate_labels(ax01)
overlay_lr04_on_ax_like_old(ax00, iso_label=False)

# Row 2: MIS6
plot_stage_row_ci(ax10, ax11, stage_mis6, add_legend=False, show_n_on_left=True)

# Row 3: MIS8
plot_stage_row_ci(ax20, ax21, stage_mis8, add_legend=False, show_n_on_left=True)

# Row 4: histograms (G,H)
tabG = draw_midpoint_histogram(ax30, HIST_G_GLOB, xmin=HIST_XMIN, xmax=HIST_XMAX,
                               bins=HIST_BINS, units=HIST_UNITS,
                               panel_letter="G", combine=HIST_G_COMBINE, logx=HIST_LOGX)

tabH = draw_midpoint_histogram(ax31, HIST_H_GLOB, xmin=HIST_XMIN, xmax=HIST_XMAX,
                               bins=HIST_BINS, units=HIST_UNITS,
                               panel_letter="H", combine=HIST_H_COMBINE, logx=HIST_LOGX)

# Save histogram tables (optional)
if HIST_SAVE_TABLE:
    tabs = [t for t in (tabG, tabH) if isinstance(t, pd.DataFrame) and not t.empty]
    if tabs:
        pd.concat(tabs, ignore_index=True).to_csv(HIST_TABLE_OUT, index=False)
        print("[saved]", HIST_TABLE_OUT)
    else:
        print("[hist-table] No tables to save (missing globs or no matched files).")

# Panel letters:
#   left  column:  A, C, E, G
#   right column:  B, D, F, H
# Placed outside top-left corners
left_axes  = [ax00, ax10, ax20, ax30]
right_axes = [ax01, ax11, ax21, ax31]
left_letters  = ["A", "C", "E", "G"]
right_letters = ["B", "D", "F", "H"]

for ax, lab in zip(left_axes, left_letters):
    ax.text(-0.10, 1.06, lab, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=12, fontweight="bold",
            clip_on=False)

for ax, lab in zip(right_axes, right_letters):
    ax.text(-0.10, 1.06, lab, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=12, fontweight="bold",
            clip_on=False)

# Layout + save
fig.subplots_adjust(hspace=HSPACE, wspace=0.20, top=0.95, bottom=0.06, left=0.08, right=0.99)

out = os.getenv("OUT_PREFIX", "Fig_percentile_CIs_4rows_clean")
plt.savefig(out + ".pdf", bbox_inches="tight", dpi=600)
plt.savefig(out + ".png", dpi=300, bbox_inches="tight")
print("[saved]", out + ".pdf", "/", out + ".png")
print(f"[info] Ne CI space: {NE_CI_SPACE} | m CI space: {M_CI_SPACE} | hist_logx={int(HIST_LOGX)}")

# ── excluded runs report ─────────────────────────────────────────────────────
def print_excluded_runs():
    if not EXCLUDED_RUNS:
        print("[QC ≥14ka] No runs excluded by im_N1/im_N2 thresholds.")
        return
    print("\n[QC ≥14ka] Excluded runs due to im_N thresholds (times ≥ 14,000 years):")
    rows=[]
    for rec in EXCLUDED_RUNS:
        rows.append({
            "file": os.path.basename(rec["path"]),
            "reasons": ";".join(rec["reasons"]),
        })
    df = pd.DataFrame(rows).sort_values(["reasons","file"])
    csv_path = "excluded_runs_qc_imN1_imN2_ge14ka.csv"
    df.to_csv(csv_path, index=False)
    flat = [r for rec in EXCLUDED_RUNS for r in rec["reasons"]]
    cnt = Counter(flat)
    print("Reason counts:", dict(cnt))
    with pd.option_context('display.max_rows', None, 'display.max_colwidth', 200):
        print(df.to_string(index=False))
    print(f"[QC] Saved CSV: {csv_path}")

print_excluded_runs()
