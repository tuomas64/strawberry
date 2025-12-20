#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MSMC-IM timeseries: MIS2, MIS6, MIS8, MIS10. Ne (left), migration rate m (right).

This variant enforces **ROOTS-ONLY** indexing and lets you choose how to
compute uncertainty bands for the median curve:
- "bootstrap" (default): resample runs with replacement
- "percentile": direct percentiles across runs (no resampling)

Select with env var: CI_METHOD=bootstrap|percentile

Usage (example):
  export ROOT_ARAB=/scratch/project_2001113/MSMC-IM/MSMC-IM-output/arabidopsis_mutation_rate
  export ROOT_FRAG=/scratch/project_2001113/MSMC-IM/MSMC-IM-output/fragaria_mutation_rate
  export CI_METHOD=percentile   # ← direct percentile CIs
  export N_POINTS=80           # ← interpolation grid size (default 100)
  python fig_S11-14_CI_percentile_or_boots.py
"""

import os, re, hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.lines import Line2D
import matplotlib.font_manager as fm

# ── Inputs (TSVs for the MIS windows) ──────────────────────────────────────────
MIS2_TSV  = os.environ.get("MIS2_TSV",  "DropPairs_MUTRATIO_range_0.20_0.30_1e-7.tsv.ARABwithin.MIS2.tsv")
MIS6_TSV  = os.environ.get("MIS6_TSV",  "DropPairs_MUTRATIO_range_0.20_0.30_1e-7.tsv.ARABwithin.MIS6.tsv")
MIS8_TSV  = os.environ.get("MIS8_TSV",  "DropPairs_MUTRATIO_range_0.20_0.30_1e-7.tsv.ARABwithin.MIS8.tsv")
MIS10_TSV = os.environ.get("MIS10_TSV", "mis10.tsv")

ROOT_FRAG = os.environ.get("ROOT_FRAG", "/scratch/project_2001113/MSMC-IM/MSMC-IM-output/fragaria_mutation_rate")
ROOT_ARAB = os.environ.get("ROOT_ARAB", "/scratch/project_2001113/MSMC-IM/MSMC-IM-output/arabidopsis_mutation_rate")
CLIMATE_FILE = os.environ.get("CLIMATE_FILE", "/scratch/project_2001113/MSMC-IM/MSMC-IM-output/spain_italy/prop/Lisiecki_and_Raymo_2005-Data.txt")

USE_PASSED_ONLY  = os.environ.get("USE_PASSED_ONLY", "1") not in ("0","false","False","no","NO")
STRICT_IMN_GUARD = os.environ.get("STRICT_IMN_GUARD", "1") not in ("0","false","False","no","NO")

# Selection controls (optional)
FILTER_KEYS = [k.strip() for k in os.environ.get("FILTER_KEYS","").split(",") if k.strip()]
FILTER_REGEX = os.environ.get("FILTER_REGEX","").strip() or None
FILTER_LIST_FILE = os.environ.get("FILTER_LIST_FILE","").strip() or None

# ===== CI controls ============================================================
# Method selector: "bootstrap" (default) or "percentile"
CI_METHOD = os.environ.get("CI_METHOD", "bootstrap").strip().lower()

# Bootstrap knobs (kept for backward compatibility)
BOOTSTRAP_ENABLED = os.environ.get("BOOTSTRAP_ENABLED", "1") not in ("0","false","False","no","NO")
BOOTSTRAP_B       = int(os.environ.get("BOOTSTRAP_B", "1000"))
BOOTSTRAP_SEED    = int(os.environ.get("BOOTSTRAP_SEED", "42"))
BOOT_CI_95 = (2.5, 97.5)
BOOT_CI_75 = (12.5, 87.5)
NE_CI_SPACE = os.environ.get("NE_CI_SPACE", "linear").strip().lower()  # "linear" or "log"
M_CI_SPACE  = os.environ.get("M_CI_SPACE",  "log").strip().lower()     # "linear" or "log"

# Hard diagnostics
HARD_ERROR_ON_IDENTICAL = os.environ.get("HARD_ERROR_ON_IDENTICAL", "1") not in ("0","false","False","no","NO")
DIAG_ROWS = int(os.environ.get("DIAG_ROWS", "8"))

# Scales
GEN   = 2.0
X_MIN = float(os.environ.get("X_MIN", 1_000))
X_MAX = float(os.environ.get("X_MAX", 400_000))
_m = os.environ.get("M_YLIMS", "").replace(","," ").split()
_n = os.environ.get("N1_YLIMS","").replace(","," ").split()
M_Y_MIN,  M_Y_MAX  = (float(_m[0]),  float(_m[1]))  if len(_m)==2 else (1e-7, 5e-3)
N1_Y_MIN, N1_Y_MAX = (float(_n[0]), float(_n[1]))  if len(_n)==2 else (1.9e2, 1.2e5)

MIS_LIST = ["MIS2", "MIS6", "MIS8", "MIS10"]

# Interpolation grid size (log-spaced)
N_POINTS = int(os.environ.get("N_POINTS", "100"))  # ← NEW: controls smooth_interpolate density

# x ticks (log years → ka labels)
major_x  = [1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000, 400_000]
major_xt = ["1","2","5","10","20","50","100","200","400"]

# ── Style ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({"pdf.fonttype":42, "ps.fonttype":42, "font.size":12})
plt.rc("text", usetex=False)
FONT_FAMILY = "Arial" if any("Arial" in f.name for f in fm.fontManager.ttflist) else "DejaVu Sans"
plt.rcParams["font.family"] = FONT_FAMILY

MEDIAN_COLOR = os.environ.get("MEDIAN_COLOR", "#111111")
CI95_COLOR   = os.environ.get("CI95_COLOR", "#fc0303")
CI75_COLOR   = os.environ.get("CI75_COLOR", "#fc0303")
CI95_LS = os.environ.get("CI95_LS", "-")
CI75_LS = os.environ.get("CI75_LS", ":")
CI95_LW = float(os.environ.get("CI95_LW", "1.6"))
CI75_LW = float(os.environ.get("CI75_LW", "1.4"))

LEGEND_X = float(os.environ.get("LEGEND_X", "0.01"))
LEGEND_Y = float(os.environ.get("LEGEND_Y", "0.82"))

# ── Climate rectangles & labels ────────────────────────────────────────────────
ka = lambda k: k*1_000
rectangles = [(ka(s),ka(e),c,a) for s,e,c,a in [
    (1, 4, "red", .05), (4, 9, "red", .20), (9, 14, "red", .05),
    (14, 17, "blue", .05), (17, 22, "blue", .20), (22, 29, "blue", .05),
    (29, 59, "red", .05),
    (59, 71, "blue", .05), (71, 84, "red", .05), (84, 93, "blue", .05),
    (93,108, "red", .05), (108,115,"blue", .20), (115,130,"red", .05),
    (130,190,"blue", .20), (190,247,"red", .05), (247,300,"blue", .20),
    (300,337,"red", .05), (337,374,"blue", .20)
]]
climate_labels = {
    "HTM":  (4_000, 9_000,  "HTM", "red"),
    "LGM":  (17_000, 22_000, "LGM", "blue"),
    "MIS 3":(29_000, 59_000, "3",   "red"),
    "MIS 4":(59_000, 71_000, "4",   "blue"),
    "MIS 5d":(108_000,115_000,"5d", "blue"),
    "MIS 5e":(115_000,130_000,"5e", "red"),
    "PGP":  (130_000,190_000,"PGP", "blue"),
    "MIS 7":(190_000,247_000,"7",   "red"),
    "MIS 8":(247_000,300_000,"8",   "blue"),
    "MIS 9":(300_000,337_000,"9",   "red"),
    "MIS 10":(337_000,374_000,"10", "blue"),
}

def add_climate_background(ax, ymin, ymax):
    for s,e,col,a in rectangles:
        if e>=X_MIN and s<=X_MAX:
            ax.add_patch(Rectangle((s, ymin), e-s, ymax-ymin, color=col, alpha=a, zorder=0))

def add_climate_labels_above(ax):
    for s,e,label,color in climate_labels.values():
        xm = max(X_MIN, min(X_MAX, (s+e)/2))
        ax.text(xm, 1.04, label, color=color,
                transform=ax.get_xaxis_transform(),
                fontsize=8, fontweight='bold',
                ha='center', va='bottom',
                bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.8))

# ── ROOTS-ONLY indexing ───────────────────────────────────────────────────────

def build_path_index_roots_only(root: str) -> dict:
    """
    Index only files that are DIRECT children of `root` (no subdirectories).
    Keys are basenames.
    """
    root_abs = os.path.abspath(root)
    try:
        entries = os.listdir(root_abs)
    except FileNotFoundError:
        return {}
    idx = {}
    for name in entries:
        ap = os.path.join(root_abs, name)
        if os.path.isfile(ap) and name.endswith("MSMC_IM.estimates.txt"):
            idx[name] = ap
    return idx

def build_path_indexes():
    path_A = build_path_index_roots_only(ROOT_ARAB)
    path_F = build_path_index_roots_only(ROOT_FRAG)
    return path_A, path_F

def _assert_roots_only(idx: dict, root: str, label: str):
    root_abs = os.path.abspath(root)
    bad = [p for p in idx.values() if os.path.dirname(os.path.abspath(p)) != root_abs]
    if bad:
        print(f"[fatal] {label}: found {len(bad)} file(s) not directly in root: {root_abs}")
        for p in bad[:5]:
            print("       ", p)
        raise SystemExit(f"{label} indexing violated 'roots-only' constraint.")

# ── Helpers ───────────────────────────────────────────────────────────────────

def _file_md5(path, chunk=1<<20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b: break
            h.update(b)
    return h.hexdigest()

def strict_imN_guard(path, nmin=10, nmax=200000):
    try:
        df = pd.read_csv(path, sep=r"\s+", engine="python", comment="#")
    except Exception:
        return False
    if not {"im_N1","im_N2"}.issubset(df.columns): return False
    n1 = pd.to_numeric(df["im_N1"], errors="coerce").to_numpy(float)
    n2 = pd.to_numeric(df["im_N2"], errors="coerce").to_numpy(float)
    ok = np.isfinite(n1) & np.isfinite(n2)
    if not np.any(ok): return False
    bad = (n1[ok] < nmin) | (n1[ok] > nmax) | (n2[ok] < nmin) | (n2[ok] > nmax)
    return not np.any(bad)

def load_series(path, col):
    try:
        df = pd.read_table(path, comment="#")
    except Exception:
        return None
    if "left_time_boundary" not in df.columns or col not in df.columns:
        return None
    t = pd.to_numeric(df["left_time_boundary"], errors="coerce").to_numpy(float) * GEN
    v = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
    m = (t >= X_MIN) & (t <= X_MAX) & np.isfinite(v) & (v > 0)
    if not np.any(m): return None
    t, v = t[m], v[m]
    o = np.argsort(t); t, v = t[o], v[o]
    if t[0] > X_MIN: t = np.insert(t, 0, X_MIN); v = np.insert(v, 0, v[0])
    if t[-1] < X_MAX: t = np.append(t, X_MAX); v = np.append(v, v[-1])
    return pd.Series(v, index=t)

def series_md5(s: pd.Series):
    arr = np.column_stack((s.index.values, s.values)).astype(np.float64)
    return hashlib.md5(arr.tobytes()).hexdigest()

def _union_times(series_list):
    times = np.array(sorted({float(t)
                             for s in series_list for t in s.index
                             if np.isfinite(t) and X_MIN <= t <= X_MAX}))
    return times if times.size else None

def _matrix_for_space(series_list, times, space):
    mats = []
    for s in series_list:
        y = s.reindex(times).ffill().bfill().to_numpy()
        y = np.where((~np.isfinite(y)) | (y <= 0), np.nan, y)
        y = pd.Series(y, index=times).ffill().bfill().to_numpy()
        if space == "log":
            y = np.log10(y)
        mats.append(y)
    return np.vstack(mats) if mats else np.empty((0, 0))

def smooth_interpolate(times, *ys, n_points=100):
    if times is None:
        return (None,) * (1 + len(ys))
    grid = np.logspace(np.log10(X_MIN), np.log10(X_MAX), n_points)
    lx  = np.log10(times); lg = np.log10(grid)
    outs = [grid]
    for y in ys:
        if y is None:
            outs.append(None); continue
        outs.append(np.interp(lg, lx, y))
    return tuple(outs)

# ===== Bootstrap helpers ======================================================

def _bootstrap_median_ci(series_list, space="linear", B=1000, ci_lo_hi=(2.5,97.5), seed=42):
    if not series_list:
        return None, None, None, None
    times = _union_times(series_list)
    if times is None or times.size < 2:
        return None, None, None, None
    mat = _matrix_for_space(series_list, times, space)
    if mat.size == 0:
        return None, None, None, None

    med_space = np.nanmedian(mat, axis=0)

    n, T = mat.shape
    rng = np.random.default_rng(seed)
    meds = np.empty((B, T), dtype=float)
    for b in range(B):
        idx = rng.integers(0, n, size=n)  # resample runs
        meds[b, :] = np.nanmedian(mat[idx, :], axis=0)
    lo_space = np.percentile(meds, ci_lo_hi[0], axis=0, method="linear")
    hi_space = np.percentile(meds, ci_lo_hi[1], axis=0, method="linear")

    if space == "log":
        med = 10.0 ** med_space
        lo = 10.0 ** lo_space
        hi = 10.0 ** hi_space
    else:
        med, lo, hi = med_space, lo_space, hi_space
    return times, med, lo, hi


def _bootstrap_both_levels(series_list, space, B, seed):
    t, med, lo95, hi95 = _bootstrap_median_ci(series_list, space, B, BOOT_CI_95, seed)
    if t is None:
        return dict(times=None, med=None, lo95=None, hi95=None, lo75=None, hi75=None)
    _, _, lo75, hi75 = _bootstrap_median_ci(series_list, space, B, BOOT_CI_75, seed+1)
    return dict(times=t, med=med, lo95=lo95, hi95=hi95, lo75=lo75, hi75=hi75)

# ===== Percentile helpers (no resampling) =====================================

def _percentile_median_ci(series_list, space="linear", ci_lo_hi=(2.5,97.5)):
    if not series_list:
        return None, None, None, None
    times = _union_times(series_list)
    if times is None or times.size < 2:
        return None, None, None, None
    mat = _matrix_for_space(series_list, times, space)
    if mat.size == 0:
        return None, None, None, None

    # Compute directly across runs at each time
    med_space = np.nanmedian(mat, axis=0)
    lo_space  = np.nanpercentile(mat, ci_lo_hi[0], axis=0, method="linear")
    hi_space  = np.nanpercentile(mat, ci_lo_hi[1], axis=0, method="linear")

    if space == "log":
        med = 10.0 ** med_space
        lo  = 10.0 ** lo_space
        hi  = 10.0 ** hi_space
    else:
        med, lo, hi = med_space, lo_space, hi_space
    return times, med, lo, hi


def _percentile_both_levels(series_list, space):
    t, med, lo95, hi95 = _percentile_median_ci(series_list, space, BOOT_CI_95)
    if t is None:
        return dict(times=None, med=None, lo95=None, hi95=None, lo75=None, hi75=None)
    _, _, lo75, hi75 = _percentile_median_ci(series_list, space, BOOT_CI_75)
    return dict(times=t, med=med, lo95=lo95, hi95=hi95, lo75=lo75, hi75=hi75)

# ── Selection helpers ─────────────────────────────────────────────────────────

def read_filter_list(path):
    if not path or not os.path.exists(path): return []
    vals = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s: vals.append(s)
    return vals

def filter_rows(df):
    if df is None or df.empty: return df
    want = set(FILTER_KEYS)
    want.update(read_filter_list(FILTER_LIST_FILE))
    regex = re.compile(FILTER_REGEX) if FILTER_REGEX else None
    if not want and not regex:
        return df
    mask = pd.Series(False, index=df.index)
    if want:
        mask = mask | df["pair_key"].astype(str).isin(want)
    if regex:
        mask = mask | df["pair_key"].astype(str).str.contains(regex)
    return df.loc[mask].copy()

def load_rows_from_tsv(tsv_path, use_passed_only=True):
    if not tsv_path or not os.path.exists(tsv_path):
        raise SystemExit(f"Missing TSV: {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t", usecols=lambda c: True)
    need = ["pair_key","arab_file","frag_file"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit(f"TSV {tsv_path} is missing required columns: {missing}")
    if use_passed_only and "passes_mutratio" in df.columns:
        df = df.loc[df["passes_mutratio"]==True]
    df = df.dropna(subset=["pair_key","arab_file","frag_file"]).copy()
    df["arab_file"] = df["arab_file"].astype(str).str.strip()
    df["frag_file"] = df["frag_file"].astype(str).str.strip()
    df = filter_rows(df)
    df = df.drop_duplicates(subset=["arab_file","frag_file"]).reset_index(drop=True)
    return df

# ── Collect curves (roots-only maps) ──────────────────────────────────────────

def collect_curves_for_subset(rows, path_A, path_F):
    A_m, A_n1, F_m, F_n1 = [], [], [], []
    seenA_m, seenA_n, seenF_m, seenF_n = set(), set(), set(), set()
    used_row_indices = []
    diag = []

    for idx, r in rows.iterrows():
        ba = r["arab_file"]; bf = r["frag_file"]
        pa = path_A.get(ba); pf = path_F.get(bf)
        if not pa or not os.path.exists(pa) or not pf or not os.path.exists(pf):
            continue

        # Enforce ROOTS-ONLY at use time as well
        if os.path.dirname(os.path.abspath(pa)) != os.path.abspath(ROOT_ARAB):
            continue
        if os.path.dirname(os.path.abspath(pf)) != os.path.abspath(ROOT_FRAG):
            continue

        # Guards
        md5a = _file_md5(pa)
        md5f = _file_md5(pf)
        same_bytes = (md5a == md5f)
        if len(diag) < DIAG_ROWS:
            diag.append((idx, ba, pa, md5a, bf, pf, md5f, same_bytes))
        if same_bytes and HARD_ERROR_ON_IDENTICAL:
            raise SystemExit(
                f"[IDENTICAL INPUTS] row={idx}\n"
                f"  ARAB: {pa}\n  FRAG: {pf}\n"
                f"  md5(arab)={md5a}  md5(frag)={md5f}\n"
                f"→ Inputs are byte-identical. Fix Fragaria folder contents."
            )

        if STRICT_IMN_GUARD and not (strict_imN_guard(pa) and strict_imN_guard(pf)):
            continue

        sAm = load_series(pa, "m");    sAn = load_series(pa, "im_N1")
        sFm = load_series(pf, "m");    sFn = load_series(pf, "im_N1")
        if sAm is None or sAn is None or sFm is None or sFn is None:
            continue

        hAm = series_md5(sAm); hAn = series_md5(sAn)
        hFm = series_md5(sFm); hFn = series_md5(sFn)
        added = False
        if hAm not in seenA_m:
            A_m.append(sAm); seenA_m.add(hAm); added = True
        if hAn not in seenA_n:
            A_n1.append(sAn); seenA_n.add(hAn); added = True
        if hFm not in seenF_m:
            F_m.append(sFm); seenF_m.add(hFm); added = True
        if hFn not in seenF_n:
            F_n1.append(sFn); seenF_n.add(hFn); added = True

        if added:
            used_row_indices.append(idx)

    if diag:
        df_diag = pd.DataFrame(diag, columns=[
            "row","arab_file","arab_path","arab_md5","frag_file","frag_path","frag_md5","same_bytes"
        ])
        print("[diag] First rows (Arab vs Frag path + MD5):")
        with pd.option_context('display.max_colwidth', None):
            print(df_diag.to_string(index=False))

    return A_m, A_n1, F_m, F_n1, used_row_indices

# ── δ18O loader ───────────────────────────────────────────────────────────────

def load_lr2005(path):
    try:
        df = pd.read_table(path, comment="#")
    except Exception:
        return None
    age_col = next((c for c in df.columns if c.lower().startswith("age") and "ka" in c.lower()), None)
    d18_col = next((c for c in df.columns if "d18" in c.lower()), None)
    if age_col is None or d18_col is None:
        age_col = "age_ka" if "age_ka" in df.columns else age_col
        d18_col = "bent_d18O" if "bent_d18O" in df.columns else d18_col
    if age_col is None or d18_col is None:
        return None
    years = pd.to_numeric(df[age_col], errors="coerce").to_numpy(float)*1_000.0
    d18   = pd.to_numeric(df[d18_col], errors="coerce").to_numpy(float)
    m = np.isfinite(years) & np.isfinite(d18)
    if not np.any(m): return None
    o = np.argsort(years); years, d18 = years[o], d18[o]
    return {"years": years, "bent_d18O": d18}

# ── Plotting (LEFT: Ne, RIGHT: m) with letters A,B,C,D ────────────────────────

def draw_2x2_for_mis(mis, rows, iso, path_A, path_F):
    if rows.empty:
        print(f"[{mis}] 0 rows after filters — skipping.")
        return None

    same_col = (rows["arab_file"] == rows["frag_file"]).mean()
    if same_col >= 0.8:
        print(f"[warn] TSV {mis}: {same_col:.0%} of rows share identical basenames across arab/frag — OK if roots differ.")

    A_m, A_n1, F_m, F_n1, used_row_indices = collect_curves_for_subset(rows, path_A, path_F)
    print(f"[{mis}] curves — Arab: m={len(A_m)}, N1={len(A_n1)}; Frag: m={len(F_m)}, N1={len(F_n1)}")

    # ===== CI computation (median curve) =====
    use_bootstrap = (CI_METHOD == "bootstrap") and BOOTSTRAP_ENABLED
    if use_bootstrap:
        A_N_ci = _bootstrap_both_levels(A_n1, space=NE_CI_SPACE, B=BOOTSTRAP_B, seed=BOOTSTRAP_SEED)
        F_N_ci = _bootstrap_both_levels(F_n1, space=NE_CI_SPACE, B=BOOTSTRAP_B, seed=BOOTSTRAP_SEED+10)
        A_m_ci = _bootstrap_both_levels(A_m,  space=M_CI_SPACE,  B=BOOTSTRAP_B, seed=BOOTSTRAP_SEED+1)
        F_m_ci = _bootstrap_both_levels(F_m,  space=M_CI_SPACE,  B=BOOTSTRAP_B, seed=BOOTSTRAP_SEED+11)
        CI_LABEL = "Bootstrap"
    else:
        A_N_ci = _percentile_both_levels(A_n1, space=NE_CI_SPACE)
        F_N_ci = _percentile_both_levels(F_n1, space=NE_CI_SPACE)
        A_m_ci = _percentile_both_levels(A_m,  space=M_CI_SPACE)
        F_m_ci = _percentile_both_levels(F_m,  space=M_CI_SPACE)
        CI_LABEL = "Percentile"

    # Smooth onto dense grid (uses global N_POINTS)
    A_gx_n, A_gmed_n, A_glo95_n, A_ghi95_n, A_glo75_n, A_ghi75_n = smooth_interpolate(
        A_N_ci["times"], A_N_ci["med"], A_N_ci["lo95"], A_N_ci["hi95"], A_N_ci["lo75"], A_N_ci["hi75"],
        n_points=N_POINTS)
    F_gx_n, F_gmed_n, F_glo95_n, F_ghi95_n, F_glo75_n, F_ghi75_n = smooth_interpolate(
        F_N_ci["times"], F_N_ci["med"], F_N_ci["lo95"], F_N_ci["hi95"], F_N_ci["lo75"], F_N_ci["hi75"],
        n_points=N_POINTS)
    A_gx_m, A_gmed_m, A_glo95_m, A_ghi95_m, A_glo75_m, A_ghi75_m = smooth_interpolate(
        A_m_ci["times"], A_m_ci["med"], A_m_ci["lo95"], A_m_ci["hi95"], A_m_ci["lo75"], A_m_ci["hi75"],
        n_points=N_POINTS)
    F_gx_m, F_gmed_m, F_glo95_m, F_ghi95_m, F_glo75_m, F_ghi75_m = smooth_interpolate(
        F_m_ci["times"], F_m_ci["med"], F_m_ci["lo95"], F_m_ci["hi95"], F_m_ci["lo75"], F_m_ci["hi75"],
        n_points=N_POINTS)

    # figure (Ne LEFT; m RIGHT)
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.4), sharex='all')
    axNeTop, axMTop = axes[0]
    axNeBot, axMBot = axes[1]

    # NEW: show x-axis tick labels under the first-row plots as well
    for ax in (axNeTop, axMTop):
        ax.tick_params(axis='x', which='both', labelbottom=True)

    def setup_axes(ax, ylabel, ylims, add_labels_above=False, add_iso=False, iso_label=False):
        ylo,yhi = ylims
        add_climate_background(ax, ylo, yhi)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(X_MIN, X_MAX); ax.set_ylim(ylo, yhi)
        ax.xaxis.set_major_locator(FixedLocator(major_x))
        ax.xaxis.set_major_formatter(FixedFormatter(major_xt))
        ax.set_ylabel(ylabel, fontweight='bold')
        ax.grid(axis='y', which='major', ls='--', alpha=.3)
        if add_labels_above: add_climate_labels_above(ax)
        if add_iso and iso is not None:
            iso_ax = ax.twinx()
            iso_ax.set_xscale('log')
            iso_ax.plot(iso["years"], -iso["bent_d18O"], lw=1.0, alpha=0.7, color="0.25")
            iso_ax.set_ylabel(r'$\delta^{18}\mathrm{O}$' if iso_label else "")
            iso_ax.tick_params(axis='y', labelsize=10)

    # Ne panels (A, B)
    setup_axes(axNeTop, "Effective population size", (N1_Y_MIN, N1_Y_MAX), add_labels_above=True, add_iso=True, iso_label=False)
    setup_axes(axNeBot, "Effective population size", (N1_Y_MIN, N1_Y_MAX), add_labels_above=False, add_iso=True, iso_label=True)

    # m panels (C, D)
    setup_axes(axMTop, "Migration rate (m)", (M_Y_MIN, M_Y_MAX), add_labels_above=True, add_iso=False, iso_label=False)
    setup_axes(axMBot, "Migration rate (m)", (M_Y_MIN, M_Y_MAX), add_labels_above=False, add_iso=False, iso_label=False)

    def draw_ci_lines(ax, gx, lo75, hi75, lo95, hi95):
        if gx is not None and lo75 is not None and hi75 is not None:
            ax.plot(gx, lo75, linestyle=CI75_LS, linewidth=CI75_LW, color=CI75_COLOR, zorder=2.1)
            ax.plot(gx, hi75, linestyle=CI75_LS, linewidth=CI75_LW, color=CI75_COLOR, zorder=2.1)
        if gx is not None and lo95 is not None and hi95 is not None:
            ax.plot(gx, lo95, linestyle=CI95_LS, linewidth=CI95_LW, color=CI95_COLOR, zorder=2.2)
            ax.plot(gx, hi95, linestyle=CI95_LS, linewidth=CI95_LW, color=CI95_COLOR, zorder=2.2)

    def draw_individuals(ax, series_list, lw=1.0, alpha=0.6, color="0.6", where="post"):
        for s in series_list:
            ax.step(s.index.values, s.values, where=where, lw=lw, alpha=alpha, color=color, zorder=1.4)

    simple_mode = False

    # ---- Plot Ne: Arab (A)
    if not simple_mode:
        draw_ci_lines(axNeTop, A_gx_n, A_glo75_n, A_ghi75_n, A_glo95_n, A_ghi95_n)
    else:
        draw_individuals(axNeTop, A_n1)
    if A_gx_n is not None and A_gmed_n is not None:
        axNeTop.plot(A_gx_n, A_gmed_n, lw=2.2, zorder=3, color=MEDIAN_COLOR)
    axNeTop.text(0.01, 0.98, "A", transform=axNeTop.transAxes, ha='left', va='top',
                 fontsize=12, fontweight='bold', bbox=dict(fc='white', ec='none', alpha=0.9, pad=0.4))
    axNeTop.text(0.02, 0.04, f"n={len(A_n1)}", transform=axNeTop.transAxes,
                 ha='left', va='bottom', fontsize=10, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.4))

    # ---- Plot Ne: Frag (B)
    if not simple_mode:
        draw_ci_lines(axNeBot, F_gx_n, F_glo75_n, F_ghi75_n, F_glo95_n, F_ghi95_n)
    else:
        draw_individuals(axNeBot, F_n1)
    if F_gx_n is not None and F_gmed_n is not None:
        axNeBot.plot(F_gx_n, F_gmed_n, lw=2.2, zorder=3, color=MEDIAN_COLOR)
    axNeBot.text(0.01, 0.98, "B", transform=axNeBot.transAxes, ha='left', va='top',
                 fontsize=12, fontweight='bold', bbox=dict(fc='white', ec='none', alpha=0.9, pad=0.4))

    # ---- Plot m: Arab (C)
    if not simple_mode:
        draw_ci_lines(axMTop, A_gx_m, A_glo75_m, A_ghi75_m, A_glo95_m, A_ghi95_m)
    else:
        draw_individuals(axMTop, A_m)
    if A_gx_m is not None and A_gmed_m is not None:
        axMTop.plot(A_gx_m, A_gmed_m, lw=2.2, zorder=3, color=MEDIAN_COLOR)
    axMTop.text(0.01, 0.98, "C", transform=axMTop.transAxes, ha='left', va='top',
                fontsize=12, fontweight='bold', bbox=dict(fc='white', ec='none', alpha=0.9, pad=0.4))

    # ---- Plot m: Frag (D)
    if not simple_mode:
        draw_ci_lines(axMBot, F_gx_m, F_glo75_m, F_ghi75_m, F_glo95_m, F_ghi95_m)
    else:
        draw_individuals(axMBot, F_m)
    if F_gx_m is not None and F_gmed_m is not None:
        axMBot.plot(F_gx_m, F_gmed_m, lw=2.2, zorder=3, color=MEDIAN_COLOR)
    axMBot.text(0.01, 0.98, "D", transform=axMBot.transAxes, ha='left', va='top',
                fontsize=12, fontweight='bold', bbox=dict(fc='white', ec='none', alpha=0.9, pad=0.4))

    # Legend: panel A
    handles = [
        Line2D([0, 1], [0, 0], color=MEDIAN_COLOR, linestyle='-', linewidth=2.2, label="Median"),
        Line2D([0, 1], [0, 0], color=CI95_COLOR, linestyle=CI95_LS, linewidth=CI95_LW, label=f"{CI_LABEL} 95% CI"),
        Line2D([0, 1], [0, 0], color=CI75_COLOR, linestyle=CI75_LS, linewidth=CI75_LW, label=f"{CI_LABEL} 75% CI"),
    ] if not simple_mode else [
        Line2D([0, 1], [0, 0], color=MEDIAN_COLOR, linestyle='-', linewidth=2.2, label="Median"),
        Line2D([0, 1], [0, 0], color="0.6", linestyle='-', linewidth=1.4, label="Individual runs"),
    ]
    axNeTop.legend(handles=handles, loc="upper left", bbox_to_anchor=(LEGEND_X, LEGEND_Y),
                   borderaxespad=0.0, frameon=True, framealpha=0.9, fontsize=9,
                   labelspacing=0.3, handlelength=2.2, handletextpad=0.6, borderpad=0.4)

    axNeBot.set_xlabel('Years ago (×1000)', fontweight='bold')
    axMBot.set_xlabel('Years ago (×1000)', fontweight='bold')

    fig.subplots_adjust(hspace=0.28, wspace=0.16)
    plt.tight_layout()

    out = f"timeseries_byMIS_SELECTED.{mis}.2x2"
    fig.savefig(out + ".pdf", bbox_inches='tight')
    fig.savefig(out + ".png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[write] {out}.pdf / .png")

    # ---- Write TSVs -----------------------------------------------------------
    def _write_series_tsv(out_base, mis, label, gx, gmed, glo95, ghi95, glo75, ghi75, n_runs):
        if gx is None or gmed is None:
            print(f"[{mis}] [skip TSV] {label}: no data"); return
        lo75 = glo75 if glo75 is not None else np.full_like(gx, np.nan)
        hi75 = ghi75 if ghi75 is not None else np.full_like(gx, np.nan)
        df_out = pd.DataFrame({
            "time": gx,
            "median": gmed,
            "lo95": glo95 if glo95 is not None else np.full_like(gx, np.nan),
            "hi95": ghi95 if ghi95 is not None else np.full_like(gx, np.nan),
            "lo75": lo75,
            "hi75": hi75,
            "n_runs": np.repeat(int(n_runs), len(gx))
        })
        out_path = f"{out_base}.{mis}.{label}.tsv"
        df_out.to_csv(out_path, sep="\t", index=False)
        print(f"[write] {out_path}  (rows={len(df_out)})")

    _write_series_tsv(out, mis, "Arab.N1",
                      A_gx_n, A_gmed_n, A_glo95_n, A_ghi95_n, A_glo75_n, A_ghi75_n, n_runs=len(A_n1))
    _write_series_tsv(out, mis, "Frag.N1",
                      F_gx_n, F_gmed_n, F_glo95_n, F_ghi95_n, F_glo75_n, F_ghi75_n, n_runs=len(F_n1))
    _write_series_tsv(out, mis, "Arab.m",
                      A_gx_m, A_gmed_m, A_glo95_m, A_ghi95_m, A_glo75_m, A_ghi75_m, n_runs=len(A_m))
    _write_series_tsv(out, mis, "Frag.m",
                      F_gx_m, F_gmed_m, F_glo95_m, F_ghi95_m, F_glo75_m, F_ghi75_m, n_runs=len(F_m))

    return used_row_indices, A_n1, F_n1, A_m, F_m

# ── Utility: stack hash to compare lists quickly ------------------------------

def _hash_stack(series_list):
    if not series_list:
        return "EMPTY"
    stacked = np.vstack([
        np.column_stack((s.index.values, s.values)).astype(np.float64) for s in series_list
    ])
    return hashlib.md5(stacked.tobytes()).hexdigest()

# ── Driver ────────────────────────────────────────────────────────────────────

def _suffix_selected_used(fname):
    if fname.endswith(".SELECTED_USED.tsv"):
        return fname
    return fname + ".SELECTED_USED.tsv"

def _peek(idx: dict, label: str, root: str):
    print(f"[peek] {label} (root={os.path.abspath(root)}) first 5:")
    for i, p in enumerate(sorted(idx.values())[:5], 1):
        print(f"   {i:02d}: {p}")

def main():
    print(f"[info] CI_METHOD={CI_METHOD}  (bootstrap enabled={BOOTSTRAP_ENABLED})  B={BOOTSTRAP_B}  SEED={BOOTSTRAP_SEED}  "
          f"NE_CI_SPACE={NE_CI_SPACE}  M_CI_SPACE={M_CI_SPACE}  X_MAX={X_MAX}")
    print(f"[info] USE_PASSED_ONLY={USE_PASSED_ONLY}  STRICT_IMN_GUARD={STRICT_IMN_GUARD}")
    print(f"[roots] ROOT_ARAB={ROOT_ARAB}")
    print(f"[roots] ROOT_FRAG={ROOT_FRAG}")
    print(f"[interp] N_POINTS={N_POINTS}")
    if os.path.abspath(ROOT_ARAB) == os.path.abspath(ROOT_FRAG):
        raise SystemExit("ROOT_ARAB and ROOT_FRAG are identical — aborting to avoid identical plots.")

    if FILTER_KEYS or FILTER_REGEX or FILTER_LIST_FILE:
        print(f"[select] FILTER_KEYS={FILTER_KEYS or '—'}  FILTER_REGEX={FILTER_REGEX or '—'}  "
              f"FILTER_LIST_FILE={FILTER_LIST_FILE or '—'}")

    iso = load_lr2005(CLIMATE_FILE)
    if iso is None:
        print(f"[warn] Could not load δ18O from {CLIMATE_FILE} — overlay disabled.")

    path_A, path_F = build_path_indexes()
    _assert_roots_only(path_A, ROOT_ARAB, "Arab")
    _assert_roots_only(path_F, ROOT_FRAG, "Frag")

    print(f"[paths] Arab index size={len(path_A)}  Frag index size={len(path_F)}")
    _peek(path_A, "Arab", ROOT_ARAB)
    _peek(path_F, "Frag", ROOT_FRAG)

    tsv_map = {"MIS2": MIS2_TSV, "MIS6": MIS6_TSV, "MIS8": MIS8_TSV, "MIS10": MIS10_TSV}

    for mis in MIS_LIST:
        tsv = tsv_map[mis]
        print(f"[{mis}] TSV: {tsv}")
        rows_all = load_rows_from_tsv(tsv, use_passed_only=USE_PASSED_ONLY)
        print(f"[{mis}] rows selected from TSV: {len(rows_all)}")

        used = draw_2x2_for_mis(mis, rows_all, iso, path_A, path_F)
        if used is None:
            used_rows = rows_all.iloc[[]].copy()
            A_n1 = F_n1 = A_m = F_m = []
        else:
            used_idx, A_n1, F_n1, A_m, F_m = used
            used_rows = rows_all.iloc[used_idx].copy()

        out_name = _suffix_selected_used(os.path.basename(tsv))
        used_rows.to_csv(out_name, sep="\t", index=False)
        print(f"[write] {out_name}  (rows={len(used_rows)})")

        print(f"[hash] {mis} Arab N1: {_hash_stack(A_n1)}")
        print(f"[hash] {mis} Frag N1: {_hash_stack(F_n1)}")
        print(f"[hash] {mis} Arab m : {_hash_stack(A_m)}")
        print(f"[hash] {mis} Frag m : {_hash_stack(F_m)}")

    print("[done] MIS2 / MIS6 / MIS8 / MIS10 finished.")

if __name__ == "__main__":
    main()
