"""
Shared loader for the Normative iEEG Sleep & Wake Atlas (Pattnaik & Litt 2024,
Pennsieve DOI 10.26275/xhte-d11l). Public, CC-BY-NC-SA. No patient access required.

Facts locked from the live dataset:
  - Two sites: HUP (Penn, sub-RID*, sleep-staged) and MNI (wake-only).
  - Sleep/wake EDF clips live in files/processed_eeg_final/<pt>_<STATE>_<idx>.edf
    STATE in {W, N2, N3, R}; each clip = 30 s, 204 Hz, multi-channel bipolar montage.
  - Channel EDF names (e.g. 'LA01-LA02') match metadata column 'name'.
  - Anatomy in 'final_label' (HUP, fine e.g. 'left entorhinal') and 'reg' (AAL).
  - Line noise 60 Hz (both sites). Nyquist 102 Hz -> no ripple/HFO band.
"""
from __future__ import annotations
import json, os, urllib.request
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
EEGDIR = os.path.join(RAW, "eeg")
DERIVED = os.path.join(ROOT, "data", "derived")
META_CSV = os.path.join(RAW, "metadata", "combined_atlas_metadata.csv")
EXIST_CSV = os.path.join(RAW, "metadata", "existing_files.csv")
API = "https://api.pennsieve.io/discover/datasets/414/versions/1/files/download-manifest"
LINE_HZ = 60.0
SFREQ = 204.0
STATES = ["W", "N2", "N3", "R"]

# ROI groups by substring match on final_label (HUP) with reg (AAL) fallback.
ROI_RULES = {
    "entorhinal":     (["entorhinal"], []),
    "parahippocampal":(["parahippocampal"], ["ParaHippocampal"]),
    "hippocampal":    (["hippocamp"], ["Hippocampus"]),
    "amygdala":       (["amygdala"], ["Amygdala"]),
    # neocortical temporal comparison (near the mesiotemporal contacts)
    "temporal_neocortex": (["middle temporal", "superior temporal", "inferior temporal"],
                            ["Temporal_Mid", "Temporal_Sup", "Temporal_Inf"]),
}
MESIOTEMPORAL = ["entorhinal", "parahippocampal", "hippocampal", "amygdala"]


def _pre(paths):
    req = urllib.request.Request(API, data=json.dumps({"paths": paths}).encode(),
                                 headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=90))["data"]


def download(paths, dest=EEGDIR):
    """Presign via Pennsieve manifest and download any not-yet-cached files."""
    os.makedirs(dest, exist_ok=True)
    need = [p for p in paths if not os.path.exists(os.path.join(dest, os.path.basename(p)))]
    got = [os.path.join(dest, os.path.basename(p)) for p in paths if p not in need]
    for i in range(0, len(need), 80):  # presigned URLs are short-lived; batch
        for d in _pre(need[i:i + 80]):
            out = os.path.join(dest, d["fileName"])
            urllib.request.urlretrieve(d["url"], out)
            got.append(out)
    return got


def load_metadata(normative_only=True):
    """Channel table with an added 'roi_group'. normative_only drops SOZ/resected."""
    df = pd.read_csv(META_CSV)
    fl = df["final_label"].astype(str).str.lower()
    rg = df["reg"].astype(str)
    df["roi_group"] = pd.Series([pd.NA] * len(df), index=df.index, dtype=object)
    for group, (fl_keys, reg_keys) in ROI_RULES.items():
        m = fl.str.contains("|".join(fl_keys), na=False)
        if reg_keys:
            m = m | rg.str.contains("|".join(reg_keys), na=False)
        df.loc[m & df["roi_group"].isna(), "roi_group"] = group
    if normative_only:
        bad = (df.get("ch1_soz") == True) | (df.get("ch2_soz") == True) | \
              (df.get("ch1_resected") == True) | (df.get("ch2_resected") == True)
        df = df[~bad.fillna(False)]
    return df


def clip_index():
    """DataFrame of every available EDF clip: pt, state, idx, filename, path."""
    ex = pd.read_csv(EXIST_CSV)
    fn = ex.iloc[:, 0].astype(str).str.rsplit("/", n=1).str[-1]
    rec = fn.str.extract(r"(?P<pt>sub-RID\d+)_(?P<state>W|N2|N3|R)_(?P<idx>\d+)\.edf")
    rec["filename"] = fn
    rec["path"] = "files/processed_eeg_final/" + fn
    return rec.dropna(subset=["pt"]).astype({"idx": int})


def read_clip(path_or_pt, state=None, idx=None):
    """Return (data[n_ch, n_t], ch_names, sfreq) for a cached/remote clip."""
    import mne
    if state is not None:
        fname = f"{path_or_pt}_{state}_{idx}.edf"
        local = os.path.join(EEGDIR, fname)
        if not os.path.exists(local):
            download([f"files/processed_eeg_final/{fname}"])
        path = local
    else:
        path = path_or_pt
    r = mne.io.read_raw_edf(path, preload=True, verbose="ERROR")
    return r.get_data(), r.ch_names, r.info["sfreq"]
