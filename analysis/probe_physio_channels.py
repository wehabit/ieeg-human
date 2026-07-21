"""One-off: does HUP phaseII on iEEG.org expose an EKG/physio channel?
Decides feasibility of the LC-infraslow RR<->spindle coherence (brief 3A)."""
import json, os, re
from ieeg.auth import Session

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
c = json.load(open(os.path.join(ROOT, "data", "ieeg_secret", "credentials.json")))
s = Session(c["username"], c["password"])

PHYS = re.compile(r"ekg|ecg|emg|eog|resp|sao2|spo2|pleth|pulse|chin|heart|"
                  r"\bhr\b|abd|thor|snore|flow|pupil|cz|fz|pz|oz|c3|c4|o1|o2|fp", re.I)
DEPTH = re.compile(r"^[A-Za-z]{1,4}\d+$")

for ds_name in ["HUP165_phaseII", "HUP157_phaseII", "HUP130_phaseII"]:
    try:
        ds = s.open_dataset(ds_name)
        labels = list(ds.get_channel_labels())
        phys = [l for l in labels if PHYS.search(l)]
        nondepth = [l for l in labels if not DEPTH.match(l)]
        print(f"=== {ds_name}: {len(labels)} channels ===")
        print("  PHYSIO candidates :", phys if phys else "NONE")
        print("  non-depth labels  :", nondepth[:50])
        print("  first 8           :", labels[:8])
        print("  last 15           :", labels[-15:])
    except Exception as e:
        print(f"=== {ds_name}: ERROR {type(e).__name__}: {e}")
