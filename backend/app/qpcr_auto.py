#!/usr/bin/env python3
"""
qpcr_auto - drop in a raw PCR file, it detects what it is, analyses it,
shows the maths, and writes a publication-style draft report (HTML + Excel).

Supported inputs (auto-detected):
  * Ct/Cq result tables  : Bio-Rad CFX, Applied Biosystems / QuantStudio / StepOne,
                           Roche LightCycler, Qiagen Rotor-Gene, generic CSV/TSV/XLSX
  * Raw amplification    : matrix (Cycle x wells) or long format (Well, Cycle, Rn/RFU)
  * Melt curve / HRM     : Temperature x wells
  * Digital PCR          : Positives/Negatives (or total partitions) per well
Not supported: end-point/gel-image PCR (no numeric data to analyse).

Usage:
  python qpcr_auto.py data.csv
  python qpcr_auto.py raw_curves.csv --plate platemap.csv
  python qpcr_auto.py data.xlsx --ref GAPDH --control Control --out report.html
"""
import argparse, base64, csv, html, io, math, os, re, sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import find_peaks, savgol_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------- reading
def norm(s):
    s = str(s).strip().lower().replace("\u0442", "t").replace("\u0421", "c")  # Cyrillic look-alikes (ABI 'C\u0442')
    return re.sub(r"[^a-z0-9]", "", s)

ALIAS = {
    "well":   {"well", "wellposition", "pos", "wellno", "position", "wellid"},
    "sample": {"sample", "samplename", "name", "sampleid"},
    "target": {"target", "targetname", "detector", "detectorname", "gene", "assay", "probe"},
    "ct":     {"ct", "cq", "cp", "ctvalue", "cqvalue", "crossingpoint", "cycle threshold"},
    "task":   {"content", "task", "sampletype", "type"},
    "qty":    {"startingquantity", "sq", "quantity", "givenconc", "standardconc", "startingqty", "copies", "concentration"},
    "dye":    {"fluor", "reporter", "dye", "channel"},
    "pos":    {"positives", "positivedroplets", "positivepartitions", "positive", "npos"},
    "neg":    {"negatives", "negativedroplets", "negativepartitions", "negative", "nneg"},
    "total":  {"accepteddroplets", "totalpartitions", "partitions", "totaldroplets", "acceptedpartitions"},
    "cycle":  {"cycle", "cycles", "cyclenumber"},
    "temp":   {"temperature", "temp", "temperaturec", "tempc"},
    "fluor":  {"rn", "rfu", "fluorescence", "deltarn", "drn", "signal"},
}

def load_grid(path):
    """Return {sheet_name: list of rows} for csv/tsv/txt/xlsx."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xlsm", ".xls"):
        import openpyxl
        wb = openpyxl.load_workbook(path, data_only=True)
        return {ws.title: [["" if c is None else str(c) for c in r] for r in ws.iter_rows(values_only=True)]
                for ws in wb.worksheets}
    raw = open(path, "rb").read()
    for enc in ("utf-8-sig", "utf-16", "latin-1"):
        try:
            text = raw.decode(enc); break
        except Exception:
            continue
    lines = [l for l in text.splitlines()]
    sample = "\n".join(lines[:60])
    delim = max([",", "\t", ";"], key=lambda d: sample.count(d))
    rows = list(csv.reader(lines, delimiter=delim))
    return {os.path.basename(path): rows}

def find_header(rows, need=2):
    best, best_i = 0, None
    for i, r in enumerate(rows[:200]):
        score = sum(1 for c in r if any(norm(c) in v for v in ALIAS.values()))
        if r and norm(r[0]) in (ALIAS["cycle"] | ALIAS["temp"]) and len([c for c in r if str(c).strip()]) >= 3:
            score += 2          # Cycle/Temperature-first matrix: well names are free text
        if score > best:
            best, best_i = score, i
    return best_i if best >= need else None

def table_from(rows, hdr):
    cols = [str(c).strip() for c in rows[hdr]]
    body = []
    for r in rows[hdr + 1:]:
        if not any(str(c).strip() for c in r):
            if body: break           # first blank line after data = end of table
            continue
        body.append((list(r) + [""] * len(cols))[:len(cols)])
    return pd.DataFrame(body, columns=cols)

def pick(df, key):
    for c in df.columns:
        if norm(c) in ALIAS[key]:
            return c
    if key == "qty":
        for c in df.columns:
            if norm(c).startswith(("startingquant", "quantity", "givenconc")):
                return c
    return None

def num(s):
    return pd.to_numeric(s.astype(str).str.replace(",", ".", regex=False).str.strip(), errors="coerce")

# --------------------------------------------------------------- detection
def detect(grid, path):
    """Return dict(pcr_type, instrument, confidence, evidence, df, kind)."""
    ev = []
    alltext = " ".join(str(c) for rows in grid.values() for r in rows[:60] for c in r).lower()
    inst = "Generic / unknown instrument"
    sig = [("cfx", "Bio-Rad CFX"), ("quantstudio", "Applied Biosystems QuantStudio"), ("stepone", "Applied Biosystems StepOne"),
           ("7500", "Applied Biosystems 7500"), ("lightcycler", "Roche LightCycler"), ("rotor-gene", "Qiagen Rotor-Gene"),
           ("rotorgene", "Qiagen Rotor-Gene"), ("quantasoft", "Bio-Rad QX/QuantaSoft ddPCR"), ("quantstudio 3d", "QuantStudio 3D dPCR"),
           ("qiacuity", "QIAcuity dPCR"), ("quantagene", "QuantaGene"), ("mic qpcr", "Bio Molecular Systems Mic"),
           ("magnetic induction", "Bio Molecular Systems Mic")]
    for k, v in sig:
        if k in alltext:
            inst = v; ev.append(f"Instrument keyword '{k}' found in file header"); break

    for name, rows in grid.items():
        hdr = find_header(rows)
        if hdr is None:
            continue
        df = table_from(rows, hdr)
        cols = {k: pick(df, k) for k in ALIAS}
        ev.append(f"Header row found (sheet '{name}', row {hdr + 1}): " + ", ".join(map(str, df.columns[:10])))

        # dPCR
        if cols["pos"] and (cols["neg"] or cols["total"]):
            ev.append("Columns for positive/negative partitions present -> digital PCR")
            return dict(kind="dpcr", pcr_type="Digital PCR (ddPCR / dPCR, absolute quantification)", instrument=inst,
                        confidence="high", evidence=ev, df=df, cols=cols)
        # melt curve
        if cols["temp"] and not cols["cycle"] and df.shape[1] >= 3 and num(df.iloc[:, 1]).notna().mean() > .8:
            ev.append("First column is temperature and the other columns are numeric signals -> melt curve / HRM")
            return dict(kind="melt", pcr_type="Melt-curve / HRM analysis (post-qPCR)", instrument=inst,
                        confidence="high", evidence=ev, df=df, cols=cols)
        # long-format amplification
        if cols["well"] and cols["cycle"] and cols["fluor"]:
            ev.append("Well + Cycle + fluorescence columns (long format) -> raw amplification curves")
            return dict(kind="amp_long", pcr_type="Real-time PCR (qPCR) - raw amplification curves", instrument=inst,
                        confidence="high", evidence=ev, df=df, cols=cols)
        # wide amplification matrix
        if cols["cycle"] and df.shape[1] >= 3:
            ev.append("Cycle column followed by numeric well columns -> raw amplification curves (matrix)")
            return dict(kind="amp_wide", pcr_type="Real-time PCR (qPCR) - raw amplification curves", instrument=inst,
                        confidence="high", evidence=ev, df=df, cols=cols)
        # Ct result table
        if cols["ct"]:
            if inst.startswith("Generic"):
                if cols["dye"] and cols["task"] and norm(cols["task"]) == "content":
                    inst = "Bio-Rad CFX (inferred from Fluor/Content columns)"
                elif cols["dye"] and cols["task"] and norm(cols["task"]) == "task":
                    inst = "Applied Biosystems QuantStudio/StepOne (inferred from Reporter/Task columns)"
                elif norm(df.columns[0]) == "include" or "pos" in [norm(c) for c in df.columns]:
                    inst = "Roche LightCycler (inferred from Include/Pos/Cp columns)"
            ev.append(f"Ct/Cq column '{cols['ct']}' present -> real-time PCR result table")
            has_ref = True
            return dict(kind="ct", pcr_type="Real-time quantitative PCR (qPCR / RT-qPCR) - Ct/Cq result table", instrument=inst,
                        confidence="high" if cols["sample"] or cols["well"] else "medium", evidence=ev, df=df, cols=cols)
    return dict(kind="unknown", pcr_type="Unrecognised", instrument=inst, confidence="none",
                evidence=ev + ["No recognisable PCR header (Ct/Cq, Cycle, Temperature, Positives/Negatives) found. "
                               "End-point/gel PCR has no numeric data to analyse."], df=None, cols={})

# ------------------------------------------------------------ small helpers
def f(x, d=3):
    return "n/a" if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))) else f"{x:.{d}f}"

from cycler import cycler
OKABE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]   # colour-blind safe
GREY = "#8c8c8c"

def pub_style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 9, "axes.titlesize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 7.5, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8, "lines.linewidth": 1.2,
        "axes.prop_cycle": cycler(color=OKABE),
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
pub_style()
# Keep raster exports within small-container memory limits. Vector PDF/SVG retain
# publication-quality line art; set PCR_EXPORT_DPI=600 for high-resolution PNGs.
FIG_DPI = max(120, min(600, int(os.environ.get("PCR_EXPORT_DPI", "300"))))
EXPORT_TIFF = os.environ.get("PCR_EXPORT_TIFF", "0").lower() in ("1", "true", "yes")
PREVIEW_DPI = 90

def stars(p):
    return "" if p is None or p != p else ("***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns")

REF_NAMES = {"gapdh", "actb", "bactin", "betaactin", "18s", "18srrna", "rn18s", "hprt", "hprt1", "rplp0", "36b4", "tbp",
             "b2m", "ppia", "gusb", "ubc", "pgk1", "rpl13a", "hmbs", "tubb", "ywhaz", "u6", "5s", "rps18", "rpl4"}
CTRL_RE = re.compile(r"^(ctrl|control|ctl|untreated|wt|wildtype|vehicle|mock|normal|healthy|calibrator|nt|dmso|pbs)", re.I)

def group_of(sample):
    g = re.sub(r"[\s_\-\.]*(rep|r|n|#)?\d+$", "", str(sample), flags=re.I).strip(" _-.")
    return g if g else str(sample)

# --------------------------------------------------------- raw curve maths
def ct_from_curve(cyc, y):
    """ABI-style: baseline = mean of cycles 3..15, threshold = baseline + 10*SD (floored at 10% of plateau)."""
    y = np.asarray(y, float); n = len(y)
    b0, b1 = 2, max(5, min(15, n // 2))
    base = y[b0:b1]; mu, sd = base.mean(), base.std(ddof=1)
    yb = y - mu
    thr = max(10 * sd, 0.10 * yb.max())
    steps = dict(baseline_cycles=(int(cyc[b0]), int(cyc[b1 - 1])), baseline_mean=mu, baseline_sd=sd, threshold=thr,
                 plateau=yb.max())
    if yb.max() < 5 * sd or yb.max() < thr:
        return None, steps
    for i in range(1, n - 1):
        if yb[i] >= thr and yb[i + 1] >= thr and yb[i - 1] < thr:
            c = cyc[i - 1] + (thr - yb[i - 1]) / (yb[i] - yb[i - 1]) * (cyc[i] - cyc[i - 1])
            steps.update(x0=cyc[i - 1], x1=cyc[i], y0=yb[i - 1], y1=yb[i])
            return float(c), steps
    return None, steps

# ----------------------------------------------------------------- analysis
class Analysis:
    def __init__(self, det, args):
        self.det, self.args = det, args
        self.math, self.qc, self.figs, self.tables, self.results_text = [], [], [], {}, []
        self.meta = {}
        self.figdir = None

    def add_fig(self, key, title, legend, fig):
        """Write compact raster + scalable vector figures without retaining large buffers."""
        import gc
        key = re.sub(r"[^A-Za-z0-9_.-]+", "_", key); files = {}
        try:
            if self.figdir:
                os.makedirs(self.figdir, exist_ok=True)
                formats = [("png", dict(dpi=FIG_DPI)), ("pdf", {}), ("svg", {})]
                if EXPORT_TIFF:
                    formats.insert(1, ("tiff", dict(dpi=FIG_DPI, pil_kwargs={"compression": "tiff_lzw"})))
                for ext, kw in formats:
                    fp = os.path.join(self.figdir, f"{key}.{ext}")
                    fig.savefig(fp, format=ext, **kw)
                    files[ext] = fp
            # Reports use small previews; avoid a second high-resolution raster in RAM.
            b = io.BytesIO()
            fig.savefig(b, format="png", dpi=PREVIEW_DPI)
            preview = "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()
            b.close()
            self.figs.append(dict(key=key, title=title, legend=legend, files=files, preview=preview))
        finally:
            plt.close(fig)
            gc.collect()

    def M(self, title, lines):
        self.math.append((title, lines))

    # ---- Ct table -> tidy
    def tidy_from_ct(self, df, cols):
        t = pd.DataFrame()
        t["Well"] = df[cols["well"]] if cols["well"] else ""
        t["Sample"] = df[cols["sample"]] if cols["sample"] else t["Well"]
        t["Target"] = df[cols["target"]] if cols["target"] else (df[cols["dye"]] if cols["dye"] else "Target1")
        t["Ct"] = num(df[cols["ct"]]).where(lambda s: s > 0)
        t["Task"] = df[cols["task"]].astype(str) if cols["task"] else ""
        t["Qty"] = num(df[cols["qty"]]) if cols["qty"] else np.nan
        t = t[(t["Sample"].astype(str).str.strip() != "") | (t["Well"].astype(str).str.strip() != "")].copy()
        tk = t["Task"].str.lower() + " " + t["Sample"].astype(str).str.lower()
        t["Role"] = np.where(tk.str.contains(r"ntc|no.?template|neg(?:ative)? ctrl|negative control|water"), "NTC",
                    np.where(tk.str.contains(r"std|standard") | t["Qty"].notna() & tk.str.contains("std|standard"), "Standard", "Unknown"))
        t.loc[(t["Role"] == "Unknown") & t["Qty"].notna() & (t["Qty"] > 0) & ~tk.str.contains("unkn|unknown"), "Role"] = "Standard"
        miss = (t["Role"] == "Standard") & t["Qty"].isna()
        if miss.any():
            def q_from_name(x):
                m = re.findall(r"\d+\.?\d*(?:[eE][+-]?\d+)?", str(x))
                return float(m[-1]) if m else np.nan
            t.loc[miss, "Qty"] = t.loc[miss, "Sample"].map(q_from_name)
            self.qc.append("Standard quantities were read from the sample names (no quantity column found) - please verify.")
        return t

    def tidy_from_curves(self, d):
        cols = d["cols"]; df = d["df"].copy()
        out = []
        curves = {}
        if d["kind"] == "amp_wide":
            cyc = num(df[cols["cycle"]]).values
            for c in df.columns:
                if c == cols["cycle"]: continue
                y = num(df[c]).values
                if np.isnan(y).all(): continue
                curves[str(c)] = (cyc, y)
        else:
            df["_c"] = num(df[cols["cycle"]]); df["_f"] = num(df[cols["fluor"]])
            key = df[cols["well"]].astype(str) + (" | " + df[cols["target"]].astype(str) if cols["target"] else "")
            for k, g in df.groupby(key, sort=False):
                g = g.sort_values("_c"); curves[k] = (g["_c"].values, g["_f"].values)
        for k, (cyc, y) in curves.items():
            m = ~np.isnan(y) & ~np.isnan(cyc)
            ct, st = ct_from_curve(cyc[m], y[m])
            out.append(dict(Well=k.split(" | ")[0], Key=k, Ct=ct, **st))
        res = pd.DataFrame(out)
        self.curves = curves
        self.meta["curve_steps"] = res
        # plate map
        t = res[["Well", "Ct"]].copy()
        t["Sample"] = t["Well"]; t["Target"] = "Target1"; t["Task"] = ""; t["Qty"] = np.nan; t["Role"] = "Unknown"
        if self.args.plate:
            pm = pd.read_csv(self.args.plate, sep=None, engine="python")
            pm.columns = [c.strip() for c in pm.columns]
            low = {norm(c): c for c in pm.columns}
            t = t.drop(columns=["Sample", "Target", "Task", "Qty", "Role"]).merge(
                pm.rename(columns={low.get("well", "Well"): "Well", low.get("sample", "Sample"): "Sample",
                                   low.get("target", "Target"): "Target"}), on="Well", how="left")
            for c, dflt in (("Sample", ""), ("Target", "Target1")):
                t[c] = t[c].fillna(dflt)
            t["Qty"] = num(t[low["quantity"]]) if "quantity" in low else np.nan
            tk = (t[low["task"]].astype(str).str.lower() if "task" in low else pd.Series("", index=t.index))
            t["Role"] = np.where(tk.str.contains("ntc|neg"), "NTC", np.where(tk.str.contains("std|stand"), "Standard", "Unknown"))
            self.meta["plate"] = True
        else:
            self.qc.append("No plate map supplied: raw files do not contain sample/target names, so Ct values are reported per well. "
                           "Provide --plate map.csv (columns Well,Sample,Target,Task,Quantity) for gene-expression analysis.")
        # figure
        tmap = t.drop_duplicates("Well").set_index("Well")["Target"].to_dict()
        tgs = sorted(set(map(str, tmap.values()))); col = {g: OKABE[i % len(OKABE)] for i, g in enumerate(tgs)}
        fig, ax = plt.subplots(figsize=(3.6, 2.8)); seen = set()
        for k, (cyc, y) in list(curves.items())[:384]:
            g = str(tmap.get(k.split(" | ")[0], "Target1"))
            ax.plot(cyc, y, lw=.7, color=col.get(g, GREY), alpha=.85, label=g if g not in seen and len(tgs) > 1 else None); seen.add(g)
        thr = (res["baseline_mean"] + res["threshold"]).median()
        if thr == thr: ax.axhline(thr, ls="--", lw=.7, color="black", label="threshold")
        ax.set_xlabel("Cycle"); ax.set_ylabel("Fluorescence (a.u.)")
        if len(tgs) > 1 or thr == thr: ax.legend(loc="upper left")
        self.add_fig("amplification_curves", "Amplification curves",
                     "Raw real-time fluorescence versus cycle number. Dashed line: automatic threshold (baseline mean + 10 SD, "
                     "floored at 10 % of the baseline-corrected plateau). Colours denote targets.", fig)
        # math for the first well
        if len(res):
            r = res.dropna(subset=["Ct"]).head(1)
            if len(r):
                r = r.iloc[0]
                self.M("Ct from raw fluorescence (example: well " + str(r["Well"]) + ")", [
                    f"Baseline = mean fluorescence of cycles {r['baseline_cycles'][0]}-{r['baseline_cycles'][1]} = {f(r['baseline_mean'], 4)}",
                    f"Baseline SD = {f(r['baseline_sd'], 4)}",
                    f"Threshold = max(10 x SD, 10% of baseline-corrected plateau) = max({f(10 * r['baseline_sd'], 4)}, {f(0.1 * r['plateau'], 4)}) = {f(r['threshold'], 4)}",
                    f"Ct = first cycle where baseline-subtracted signal crosses threshold, by linear interpolation between cycle {r['x0']:.0f} "
                    f"(F={f(r['y0'], 4)}) and cycle {r['x1']:.0f} (F={f(r['y1'], 4)}): "
                    f"Ct = {r['x0']:.0f} + ({f(r['threshold'], 4)} - {f(r['y0'], 4)}) / ({f(r['y1'], 4)} - {f(r['y0'], 4)}) = {f(r['Ct'])}",
                    "The same procedure is applied to every well."])
        return t

    # ---- main entry
    def run(self):
        d = self.det; k = d["kind"]
        if k == "ct":
            t = self.tidy_from_ct(d["df"], d["cols"]); self.analyse_ct(t)
        elif k in ("amp_wide", "amp_long"):
            t = self.tidy_from_curves(d); self.analyse_ct(t)
        elif k == "melt":
            self.analyse_melt()
        elif k == "dpcr":
            self.analyse_dpcr()

    # ---- qPCR core
    def analyse_ct(self, t):
        self.tables["Raw_Ct"] = t[["Well", "Sample", "Target", "Role", "Qty", "Ct"]]
        n_wells = len(t); n_und = int(t["Ct"].isna().sum())
        self.meta.update(n_wells=n_wells, n_undetermined=n_und,
                         targets=sorted(t["Target"].astype(str).unique()), samples=sorted(t["Sample"].astype(str).unique()))
        # NTC check
        ntc = t[t["Role"] == "NTC"]
        for tg, g in ntc.groupby("Target"):
            amp = g[g["Ct"].notna() & (g["Ct"] < 35)]
            if len(amp):
                self.qc.append(f"NTC for {tg} amplified (Ct {', '.join(f(x, 2) for x in amp['Ct'])}) - possible contamination / primer-dimer.")
        if len(ntc) == 0:
            self.qc.append("No no-template control (NTC) identified in the file.")
        unk = t[t["Role"] == "Unknown"].copy()
        # replicate stats
        rep = unk.groupby(["Sample", "Target"])["Ct"].agg(["count", "mean", "std", "min", "max"]).reset_index()
        rep["CV%"] = rep["std"] / rep["mean"] * 100
        self.tables["Replicates"] = rep
        for _, r in rep[rep["std"] > 0.3].iterrows():
            self.qc.append(f"High technical-replicate SD (>0.30 cycles): {r['Sample']} / {r['Target']} SD = {f(r['std'], 2)}.")
        for _, r in unk[unk["Ct"] > 35].iterrows():
            self.qc.append(f"Late Ct (>35, low confidence): {r['Sample']} / {r['Target']} Ct = {f(r['Ct'], 2)}.")
        if n_und:
            self.qc.append(f"{n_und} well(s) have no Ct (undetermined / no amplification).")
        # example mean maths
        ex = unk.dropna(subset=["Ct"]).groupby(["Sample", "Target"]).filter(lambda g: len(g) > 1)
        if len(ex):
            s0, t0 = ex.iloc[0]["Sample"], ex.iloc[0]["Target"]
            v = ex[(ex["Sample"] == s0) & (ex["Target"] == t0)]["Ct"].values
            self.M(f"Technical-replicate mean (example: {s0} / {t0})", [
                f"Mean Ct = ({' + '.join(f(x) for x in v)}) / {len(v)} = {f(v.mean())}",
                f"SD = sqrt( sum((Ct_i - mean)^2) / (n-1) ) = {f(v.std(ddof=1))}",
                "Calculated identically for every sample/target pair (see Replicates sheet in the Excel file)."])
        # standard curves
        self.eff = {}
        std = t[(t["Role"] == "Standard") & t["Qty"].notna() & (t["Qty"] > 0) & t["Ct"].notna()]
        sc_rows = []
        for tg, g in std.groupby("Target"):
            if g["Qty"].nunique() < 3: continue
            x = np.log10(g["Qty"].values); y = g["Ct"].values
            lr = stats.linregress(x, y)
            E = 10 ** (-1 / lr.slope) - 1
            self.eff[tg] = 1 + E
            sc_rows.append(dict(Target=tg, Slope=lr.slope, Intercept=lr.intercept, R2=lr.rvalue ** 2, Efficiency_pct=E * 100,
                                Levels=g["Qty"].nunique()))
            self.M(f"Standard curve & amplification efficiency ({tg})", [
                f"Linear regression of Ct against log10(quantity), {len(g)} points / {g['Qty'].nunique()} dilution levels:",
                f"Ct = slope x log10(Q) + intercept  ->  Ct = {f(lr.slope, 4)} x log10(Q) + {f(lr.intercept, 4)}",
                f"R\u00b2 = {f(lr.rvalue ** 2, 4)}",
                f"Efficiency E = 10^(-1/slope) - 1 = 10^(-1/{f(lr.slope, 4)}) - 1 = {f(E, 4)}  ({f(E * 100, 1)} %)",
                "Acceptable range (MIQE practice): 90-110 % and R\u00b2 \u2265 0.98."])
            if not (0.9 <= 1 + E <= 1.1 + 0.0) and not (90 <= E * 100 + 0 <= 110):
                self.qc.append(f"{tg}: efficiency {f(E * 100, 1)} % outside 90-110 %.")
            if lr.rvalue ** 2 < 0.98:
                self.qc.append(f"{tg}: standard-curve R\u00b2 = {f(lr.rvalue ** 2, 4)} (<0.98).")
            fig, ax = plt.subplots(figsize=(3.2, 2.8))
            ax.scatter(x, y, s=14, color=OKABE[0], zorder=3); xs = np.array([x.min(), x.max()])
            ax.plot(xs, lr.slope * xs + lr.intercept, color=OKABE[3], zorder=2)
            ax.set_xlabel("log$_{10}$ quantity"); ax.set_ylabel("Ct"); ax.set_title(str(tg), fontstyle="italic")
            ax.text(0.97, 0.97, f"y = {lr.slope:.3f}x + {lr.intercept:.2f}\nR$^2$ = {lr.rvalue ** 2:.4f}\nE = {E * 100:.1f} %",
                    transform=ax.transAxes, ha="right", va="top", fontsize=7.5)
            self.add_fig(f"standard_curve_{tg}", f"Standard curve - {tg}",
                         f"Standard curve for {tg}: Ct versus log10 starting quantity ({len(g)} points, {g['Qty'].nunique()} dilution levels). "
                         f"Line: least-squares fit; efficiency E = 10^(-1/slope) - 1.", fig)
            # absolute quantities
            a = unk[unk["Target"] == tg].dropna(subset=["Ct"]).copy()
            if len(a):
                a["Quantity"] = 10 ** ((a["Ct"] - lr.intercept) / lr.slope)
                self.tables.setdefault("Absolute_Quantity", pd.DataFrame())
                self.tables["Absolute_Quantity"] = pd.concat([self.tables["Absolute_Quantity"],
                                                              a[["Sample", "Target", "Well", "Ct", "Quantity"]]])
                r0 = a.iloc[0]
                self.M(f"Absolute quantification ({tg}, example: {r0['Sample']})", [
                    "Q = 10^((Ct - intercept) / slope)",
                    f"Q = 10^(({f(r0['Ct'])} - {f(lr.intercept, 4)}) / {f(lr.slope, 4)}) = {r0['Quantity']:.4g}"])
        if sc_rows: self.tables["Standard_Curves"] = pd.DataFrame(sc_rows)
        # relative quantification
        self.heatmap(unk)
        self.relative(unk)

    def heatmap(self, unk):
        m = unk.groupby(["Sample", "Target"])["Ct"].mean().unstack("Target")
        if m.shape[0] < 2 or m.shape[0] > 48 or m.shape[1] < 1: return
        fig, ax = plt.subplots(figsize=(max(3.0, 0.75 * m.shape[1] + 1.9), max(2.2, 0.27 * m.shape[0] + 1.0)))
        im = ax.imshow(m.values, aspect="auto", cmap="viridis_r")
        ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=30, ha="right", fontstyle="italic")
        ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index)
        lo, hi = np.nanmin(m.values), np.nanmax(m.values)
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                v = m.values[i, j]
                if v == v: ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.5,
                                   color="black" if (v - lo) / (hi - lo + 1e-9) < .45 else "white")
        for sp in ax.spines.values(): sp.set_visible(False)
        cb = fig.colorbar(im, ax=ax, fraction=.06, pad=.03); cb.set_label("Mean Ct"); cb.outline.set_visible(False)
        self.add_fig("ct_heatmap", "Mean Ct overview",
                     "Mean Ct (technical replicates) for every sample and target. Lower Ct (yellow) = higher template abundance.", fig)

    def relative(self, unk):
        tg_all = sorted(unk["Target"].astype(str).unique())
        ref = self.args.ref
        if not ref:
            cand = [t for t in tg_all if norm(t) in REF_NAMES]
            ref = cand[0] if cand else None
        if not ref or ref not in tg_all or len(tg_all) < 2:
            self.meta["relative"] = "skipped"
            self.qc.append("Relative expression (2^-\u0394\u0394Ct) not run: no reference gene found (use --ref NAME) or only one target present.")
            return
        s = unk.groupby(["Sample", "Target"])["Ct"].mean().unstack("Target")
        if ref not in s.columns: return
        groups = pd.Series({smp: group_of(smp) for smp in s.index})
        ctrl = self.args.control
        if not ctrl:
            cg = [g for g in sorted(groups.unique()) if CTRL_RE.match(g)]
            ctrl = cg[0] if cg else sorted(groups.unique())[0]
            self.qc.append(f"Control/calibrator group assumed to be '{ctrl}' (override with --control). "
                           f"Groups inferred from sample names by removing trailing replicate numbers.")
        self.meta.update(relative="done", ref=ref, ctrl=ctrl, groups=sorted(groups.unique()))
        rows, per_sample = [], []; plot_data = {}
        use_pfaffl = ref in self.eff
        for tg in [t for t in s.columns if t != ref]:
            d = (s[tg] - s[ref]).dropna()
            dct = pd.DataFrame({"dCt": d, "Group": groups.reindex(d.index)})
            for smp, r in dct.iterrows():
                per_sample.append(dict(Target=tg, Sample=smp, Group=r["Group"], Ct_target=s.loc[smp, tg], Ct_ref=s.loc[smp, ref], dCt=r["dCt"]))
            c = dct[dct["Group"] == ctrl]["dCt"]
            if c.empty: continue
            lines = [f"Reference gene: {ref}.  Calibrator group: {ctrl} (n = {len(c)}).",
                     "\u0394Ct = Ct(target) - Ct(reference), per sample",
                     f"Mean \u0394Ct of {ctrl} = {f(c.mean())}"]
            for gname, g in dct.groupby("Group"):
                x = g["dCt"]
                ddct = x.mean() - c.mean()
                fold = 2 ** (-ddct)
                sd = x.std(ddof=1) if len(x) > 1 else np.nan
                lo, hi = (2 ** (-(ddct + sd)), 2 ** (-(ddct - sd))) if not np.isnan(sd) else (np.nan, np.nan)
                if gname == ctrl or len(x) < 2 or len(c) < 2:
                    p = np.nan
                else:
                    p = stats.ttest_ind(x, c, equal_var=False).pvalue
                pf = np.nan
                if use_pfaffl and tg in self.eff and gname != ctrl:
                    cm = s.loc[[i for i in s.index if groups[i] == ctrl]]; gm = s.loc[[i for i in s.index if groups[i] == gname]]
                    pf = self.eff[tg] ** (cm[tg].mean() - gm[tg].mean()) / self.eff[ref] ** (cm[ref].mean() - gm[ref].mean())
                plot_data.setdefault(tg, []).append(dict(group=gname, fold=fold, lo=lo, hi=hi, p=p, pts=2 ** (-(x.values - c.mean())), ctrl=(gname == ctrl)))
                rows.append(dict(Target=tg, Group=gname, n=len(x), mean_dCt=x.mean(), SD_dCt=sd, ddCt=ddct, Fold_2_ddCt=fold,
                                 Fold_low=lo, Fold_high=hi, p_Welch=p, Fold_Pfaffl=pf))
                if gname != ctrl:
                    lines += [f"-- {gname} (n = {len(x)}) --",
                              f"Mean \u0394Ct = {f(x.mean())},  SD = {f(sd)}",
                              f"\u0394\u0394Ct = mean \u0394Ct({gname}) - mean \u0394Ct({ctrl}) = {f(x.mean())} - {f(c.mean())} = {f(ddct)}",
                              f"Fold change = 2^(-\u0394\u0394Ct) = 2^({f(-ddct)}) = {f(fold)}   (range 2^-(\u0394\u0394Ct\u00b1SD): {f(lo)} - {f(hi)})",
                              f"Welch two-sample t-test on \u0394Ct values: p = {p:.3g}" if not np.isnan(p) else
                              "Statistics: need n \u2265 2 per group for a t-test."]
                    if not np.isnan(pf):
                        lines.append(f"Pfaffl (efficiency-corrected) ratio = E_target^(\u0394Ct_target) / E_ref^(\u0394Ct_ref) = "
                                     f"{f(self.eff[tg], 3)}^({f(cm[tg].mean() - gm[tg].mean())}) / {f(self.eff[ref], 3)}^({f(cm[ref].mean() - gm[ref].mean())}) = {f(pf)}")
            self.M(f"Relative expression of {tg} (Livak 2^-\u0394\u0394Ct)", lines)
        if rows:
            R = pd.DataFrame(rows); self.tables["Relative_Expression"] = R
            self.tables["dCt_per_sample"] = pd.DataFrame(per_sample)
            tg_list = list(plot_data); n = len(tg_list)
            fig, axes = plt.subplots(1, n, figsize=(min(7.2, 2.2 * n + 1.0), 2.9), squeeze=False)
            for ax, tg in zip(axes[0], tg_list):
                items = plot_data[tg]
                for i, it in enumerate(items):
                    colr = GREY if it["ctrl"] else OKABE[(i % (len(OKABE) - 1)) + 0]
                    ax.bar(i, it["fold"], width=.62, color=colr, edgecolor="black", linewidth=.6, zorder=2)
                    top = it["fold"]
                    if it["lo"] == it["lo"]:
                        ax.errorbar(i, it["fold"], yerr=[[max(it["fold"] - it["lo"], 0)], [max(it["hi"] - it["fold"], 0)]],
                                    fmt="none", ecolor="black", capsize=3, lw=.9, zorder=3); top = max(top, it["hi"])
                    pts = it["pts"]; jit = np.linspace(-.16, .16, len(pts)) if len(pts) > 1 else np.zeros(1)
                    ax.scatter(i + jit, pts, s=10, facecolor="white", edgecolor="black", linewidth=.6, zorder=4); top = max(top, pts.max())
                    ax.text(i, top * 1.05, stars(it["p"]), ha="center", va="bottom", fontsize=9)
                ax.axhline(1, ls="--", lw=.6, color="black", zorder=1)
                ax.set_xticks(range(len(items))); ax.set_xticklabels([it["group"] for it in items], rotation=30 if max(len(it["group"]) for it in items) > 6 else 0,
                                                                     ha="right" if max(len(it["group"]) for it in items) > 6 else "center")
                ax.set_title(str(tg), fontstyle="italic"); ax.margins(y=.15); ax.set_ylim(bottom=0)
            axes[0][0].set_ylabel(f"Relative expression\n(fold change vs {ctrl})")
            fig.tight_layout()
            self.add_fig("relative_expression", f"Relative expression normalised to {ref}",
                         f"Relative expression ({', '.join(map(str, tg_list))}) normalised to {ref} by the 2^-\u0394\u0394Ct method; calibrator group: {ctrl} (grey, = 1). "
                         "Bars: mean fold change; error bars: 2^-(\u0394\u0394Ct \u00b1 SD of \u0394Ct); open circles: individual biological samples. "
                         "Welch's t-test on \u0394Ct vs calibrator: * p < 0.05, ** p < 0.01, *** p < 0.001, ns = not significant.", fig)

    # ---- melt
    def analyse_melt(self):
        d = self.det["df"]; cols = self.det["cols"]
        T = num(d[cols["temp"]]).values; rows = []
        fig, ax = plt.subplots(figsize=(3.6, 2.8)); ncurve = 0
        for c in d.columns:
            if c == cols["temp"]: continue
            y = num(d[c]).values
            if np.isnan(y).all(): continue
            m = ~np.isnan(y) & ~np.isnan(T)
            Tm_, ym = T[m], y[m]
            w = min(11, (len(ym) // 2) * 2 - 1)
            ys = savgol_filter(ym, w, 3) if w >= 5 else ym
            dy = -np.gradient(ys, Tm_)
            pk, pr = find_peaks(dy, prominence=0.1 * dy.max())
            rows.append(dict(Well=c, Peaks=len(pk), Tm1=Tm_[pk[0]] if len(pk) else np.nan,
                             Tm2=Tm_[pk[1]] if len(pk) > 1 else np.nan))
            ax.plot(Tm_, dy, lw=1.0, label=str(c)); ncurve += 1
        ax.set_xlabel("Temperature (\u00b0C)"); ax.set_ylabel("-dF/dT")
        if ncurve <= 8: ax.legend()
        self.add_fig("melt_curve", "Melt-curve analysis",
                     "Negative first derivative of fluorescence with respect to temperature (-dF/dT) after Savitzky-Golay smoothing. "
                     "Peak maxima give the melting temperature (Tm) of each amplicon.", fig)
        R = pd.DataFrame(rows); self.tables["Melt_Tm"] = R
        for _, r in R[R["Peaks"] > 1].iterrows():
            self.qc.append(f"Well {r['Well']} has {int(r['Peaks'])} melt peaks (possible non-specific product / primer-dimer).")
        self.M("Melting temperature (Tm)", [
            "Fluorescence is smoothed (Savitzky-Golay, 3rd order) and the negative derivative -dF/dT is computed numerically.",
            "Tm = temperature at each peak maximum of -dF/dT (peaks must exceed 10 % of the tallest peak).",
            "A single sharp peak indicates one specific amplicon; extra peaks at lower Tm suggest primer-dimers."])

    # ---- digital PCR
    def analyse_dpcr(self):
        d = self.det["df"]; c = self.det["cols"]
        pos = num(d[c["pos"]])
        if c["neg"]: neg = num(d[c["neg"]]); tot = pos + neg
        else: tot = num(d[c["total"]]); neg = tot - pos
        vol = self.args.dpcr_vol
        out = pd.DataFrame({"Well": d[c["well"]] if c["well"] else range(1, len(d) + 1),
                            "Sample": d[c["sample"]] if c["sample"] else "",
                            "Target": d[c["target"]] if c["target"] else "Target1",
                            "Positives": pos, "Negatives": neg, "Total": tot})
        p = out["Positives"] / out["Total"]
        lam = -np.log(out["Negatives"] / out["Total"])
        out["lambda"] = lam
        out["copies_per_uL"] = lam / vol
        z = 1.96; n = out["Total"]
        centre = (p + z * z / (2 * n)) / (1 + z * z / n)
        half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
        out["CI_low"] = -np.log(1 - (centre - half)) / vol
        out["CI_high"] = -np.log(1 - (centre + half)) / vol
        self.tables["dPCR_Concentration"] = out
        lab = (out["Sample"].astype(str) + "\n" + out["Target"].astype(str)) if out["Target"].astype(str).nunique() > 1 else out["Sample"].astype(str)
        fig, ax = plt.subplots(figsize=(max(2.6, 0.6 * len(out) + 1.2), 2.8)); tg_u = list(out["Target"].astype(str).unique())
        cols_ = [OKABE[tg_u.index(str(t)) % len(OKABE)] for t in out["Target"]]
        yerr = [(out["copies_per_uL"] - out["CI_low"]).clip(lower=0).fillna(0), (out["CI_high"] - out["copies_per_uL"]).clip(lower=0).fillna(0)]
        ax.bar(range(len(out)), out["copies_per_uL"], yerr=yerr, color=cols_, edgecolor="black", linewidth=.6, capsize=3, error_kw=dict(lw=.9))
        ax.set_xticks(range(len(out))); ax.set_xticklabels(lab, rotation=30 if len(out) > 4 else 0, ha="right" if len(out) > 4 else "center")
        ax.set_ylabel("Concentration (copies/\u00b5L)")
        self.add_fig("dpcr_concentration", "Digital PCR concentrations",
                     "Absolute target concentration per reaction from Poisson-corrected positive-partition fractions. Error bars: 95 % Wilson confidence interval.", fig)
        r = out.iloc[0]
        self.M("Absolute quantification by Poisson statistics (example: first well)", [
            f"Positive fraction p = positives / total = {int(r['Positives'])} / {int(r['Total'])} = {f(r['Positives'] / r['Total'], 5)}",
            f"Mean copies per partition \u03bb = -ln(1 - p) = -ln({int(r['Negatives'])}/{int(r['Total'])}) = {f(r['lambda'], 5)}",
            f"Concentration = \u03bb / partition volume = {f(r['lambda'], 5)} / {vol} \u00b5L = {r['copies_per_uL']:.4g} copies/\u00b5L",
            "95 % CI: Wilson interval on p, transformed through the same Poisson equation.",
            f"Partition volume assumed {vol} \u00b5L (0.85 nL, Bio-Rad droplet). Change with --dpcr-vol."])
        for _, r in out[out["Total"] < 10000].iterrows():
            self.qc.append(f"Well {r['Well']}: only {int(r['Total'])} accepted partitions (<10,000 recommended).")
        self.meta["n_wells"] = len(out)

# ------------------------------------------------------------------ report
CSS = """body{font-family:Georgia,serif;max-width:900px;margin:30px auto;padding:0 18px;color:#222;line-height:1.5}
h1{font-size:1.7em;border-bottom:2px solid #333}h2{margin-top:1.6em;border-bottom:1px solid #aaa}
table{border-collapse:collapse;font-size:.85em;margin:8px 0}td,th{border:1px solid #bbb;padding:3px 8px;text-align:right}
th{background:#eee}td:first-child,th:first-child{text-align:left}.warn{background:#fff4d6;padding:4px 10px;margin:3px 0}
details{margin:6px 0;background:#f7f7f7;padding:6px 10px}code,pre{font-family:Consolas,monospace;font-size:.88em}
.box{background:#eef5ff;padding:10px 14px;border-left:4px solid #4c78a8}figure{margin:12px 0}figcaption{font-size:.85em;color:#555}"""

def df_html(df, fmt=3):
    return df.to_html(index=False, float_format=lambda x: f"{x:.{fmt}f}" if abs(x) >= 0.001 or x == 0 else f"{x:.3g}", na_rep="-", border=0)

def draft_text(a):
    d = a.det; m = a.meta; T = a.tables; out = []
    out.append(f"<b>Methods.</b> Raw data were exported from {html.escape(d['instrument'])} and automatically classified as "
               f"<i>{html.escape(d['pcr_type'])}</i>. ")
    if d["kind"] in ("amp_wide", "amp_long"):
        out.append("Cycle-threshold (Ct) values were determined from raw fluorescence using a baseline of cycles 3-15 and a threshold of "
                   "10 standard deviations of baseline fluorescence (floored at 10 % of the baseline-corrected plateau), with linear interpolation between cycles. ")
    if "Standard_Curves" in T:
        out.append("Amplification efficiency was calculated from the slope of a ten-fold dilution standard curve as E = 10<sup>(-1/slope)</sup> - 1. ")
    if m.get("relative") == "done":
        out.append(f"Relative expression was calculated with the 2<sup>-\u0394\u0394Ct</sup> method (Livak &amp; Schmittgen, 2001) using {html.escape(m['ref'])} "
                   f"as reference gene and the {html.escape(m['ctrl'])} group as calibrator; "
                   + ("efficiency-corrected ratios (Pfaffl, 2001) are also given. " if any(r == r for r in T["Relative_Expression"]["Fold_Pfaffl"]) else "")
                   + "Groups were compared on \u0394Ct values with Welch's t-test (two-sided). ")
    if d["kind"] == "dpcr":
        out.append("Target concentration was computed from the fraction of positive partitions using Poisson statistics (\u03bb = -ln(1-p)); "
                   "95 % confidence intervals used the Wilson interval (Dube et al., 2008). ")
    out.append("<br><br><b>Results.</b> ")
    if m.get("relative") == "done":
        R = T["Relative_Expression"]; R = R[R["Group"] != m["ctrl"]]
        for _, r in R.iterrows():
            p = "" if r["p_Welch"] != r["p_Welch"] else f", p = {r['p_Welch']:.3g}"
            direc = "up-regulated" if r["Fold_2_ddCt"] > 1 else "down-regulated"
            out.append(f"{html.escape(str(r['Target']))} was {direc} {r['Fold_2_ddCt']:.2f}-fold in {html.escape(str(r['Group']))} relative to "
                       f"{html.escape(m['ctrl'])} (n = {int(r['n'])}{p}). ")
    elif "Absolute_Quantity" in T:
        out.append("Absolute quantities were interpolated from the standard curve (table below). ")
    elif "dPCR_Concentration" in T:
        out.append("Absolute copy concentrations with 95 % CI are given in the table below. ")
    elif "Melt_Tm" in T:
        out.append("Melting temperatures are given in the table below. ")
    else:
        out.append("Ct values per well/sample are given in the tables below. ")
    out.append("<br><br><i>Draft text generated automatically - verify against your experimental design before submission.</i>")
    return "".join(out)

def build_report(a, path, outdir):
    d = a.det; m = a.meta
    out_html = os.path.join(outdir, "report.html"); xl = os.path.join(outdir, "tables.xlsx")
    h = [f"<html><head><meta charset='utf-8'><title>PCR analysis report</title><style>{CSS}</style></head><body>",
         "<h1>PCR Analysis Report (auto-generated draft)</h1>",
         f"<p>Source file: <code>{html.escape(os.path.basename(path))}</code></p>",
         "<h2>1. File recognition</h2><div class='box'>"
         f"<b>Detected PCR type:</b> {html.escape(d['pcr_type'])}<br><b>Instrument / software:</b> {html.escape(d['instrument'])}<br>"
         f"<b>Confidence:</b> {d['confidence']}<br><b>Evidence:</b><ul>" + "".join(f"<li>{html.escape(e)}</li>" for e in d["evidence"]) + "</ul>"
         "<small>Note: a Ct table alone cannot distinguish qPCR from RT-qPCR (the reverse-transcription step is an assay-design fact); "
         "it is described here as qPCR/RT-qPCR.</small></div>"]
    if m.get("targets"):
        h.append(f"<p>{m['n_wells']} wells, {len(m['targets'])} target(s): {html.escape(', '.join(m['targets']))}; "
                 f"{len(m['samples'])} sample name(s); {m['n_undetermined']} undetermined.</p>")
    h.append("<h2>2. Drafted Methods and Results</h2><p>" + draft_text(a) + "</p>")
    h.append("<h2>3. Quality control</h2>" + ("".join(f"<div class='warn'>\u26a0 {html.escape(q)}</div>" for q in a.qc) or "<p>No QC flags raised.</p>"))
    h.append("<h2>4. Result tables</h2>")
    for name, t in a.tables.items():
        if name == "Raw_Ct" and len(t) > 400: t = t.head(400)
        h.append(f"<h3>{name.replace('_', ' ')}</h3>" + df_html(t))
    h.append(f"<h2>5. Figures (publication files: {FIG_DPI} dpi PNG/TIFF, vector PDF/SVG)</h2>")
    for i, fg in enumerate(a.figs, 1):
        links = " ".join(f"<a href='figures/{fg['key']}.{e}'>{e.upper()}</a>" for e in fg["files"])
        h.append(f"<figure><img src='{fg['preview']}' style='max-width:100%'><figcaption><b>Figure {i}. {html.escape(fg['title'])}.</b> "
                 f"{html.escape(fg['legend'])}<br><small>Download: {links}</small></figcaption></figure>")
    h.append("<h2>6. The mathematics, step by step</h2><p>Every formula with your numbers substituted in. Worked on one representative case each; "
             "the same operation is applied to all rows (full numbers are in tables.xlsx).</p>")
    for title, lines in a.math:
        h.append(f"<details open><summary><b>{html.escape(title)}</b></summary><pre>" + html.escape("\n".join(lines)) + "</pre></details>")
    h.append("<h2>7. References</h2><ol><li>Livak KJ, Schmittgen TD. Methods 2001;25:402-8.</li><li>Pfaffl MW. Nucleic Acids Res 2001;29:e45.</li>"
             "<li>Bustin SA et al. The MIQE guidelines. Clin Chem 2009;55:611-22.</li><li>Dube S et al. PLoS ONE 2008;3:e2876 (dPCR Poisson statistics).</li></ol>")
    h.append("</body></html>")
    open(out_html, "w", encoding="utf-8").write("\n".join(h))
    with pd.ExcelWriter(xl) as w:
        for name, t in a.tables.items():
            t.to_excel(w, sheet_name=name[:31], index=False)
    with open(os.path.join(outdir, "calculations.txt"), "w", encoding="utf-8") as fh:
        for title, lines in a.math: fh.write(title + "\n" + "\n".join(lines) + "\n\n")
    return out_html, xl

def run_pipeline(file, plate=None, ref=None, control=None, dpcr_vol=0.00085, outdir=None):
    """Detect -> analyse -> write report, tables, calculations and compact raster/vector figures."""
    args = argparse.Namespace(file=file, plate=plate, ref=ref, control=control, dpcr_vol=dpcr_vol)
    det = detect(load_grid(file), file)
    if det["kind"] == "unknown":
        return det, None, None
    outdir = outdir or (os.path.splitext(file)[0] + "_pcr_output"); os.makedirs(outdir, exist_ok=True)
    a = Analysis(det, args); a.figdir = os.path.join(outdir, "figures"); a.run()
    rep, xl = build_report(a, file, outdir)
    return det, a, dict(outdir=outdir, report=rep, xlsx=xl)

def main():
    ap = argparse.ArgumentParser(description="Universal PCR raw-file analyser")
    ap.add_argument("file"); ap.add_argument("--plate"); ap.add_argument("--ref"); ap.add_argument("--control")
    ap.add_argument("--dpcr-vol", type=float, default=0.00085, help="partition volume in uL (default 0.85 nL)")
    ap.add_argument("--outdir", help="output folder (default: <file>_pcr_output)")
    g = ap.parse_args()
    det, a, paths = run_pipeline(g.file, g.plate, g.ref, g.control, g.dpcr_vol, g.outdir)
    print(f"Detected : {det['pcr_type']}\nInstrument: {det['instrument']}\nConfidence: {det['confidence']}")
    for e in det["evidence"]: print("  -", e)
    if a is None: sys.exit("Cannot analyse this file.")
    print(f"Output   : {paths['outdir']}  (report.html, tables.xlsx, calculations.txt, figures/ PNG at {FIG_DPI} dpi + PDF/SVG)")
    for q in a.qc: print("  QC:", q)

if __name__ == "__main__":
    main()