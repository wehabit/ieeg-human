"""
Probe iEEG.org for datasets THIS account can open that have MTL depth channels + long duration.
Tries known public naming conventions (Penn HUP phase-II continuous, Mayo I001_P*, Study 0**),
reports the reachable ones with sample rate, hours, and mesial-temporal channel count.
"""
import json, os, re, sys
from ieeg.auth import Session

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CRED = os.path.join(ROOT, "data", "ieeg_secret", "credentials.json")
MTL = re.compile(r"^[LR]?(A|AH|AMY|H|HC|HIP|AH|PH|EC|ENT|PHG|MH|MT|TP)\d", re.I)

HUP_N = [64,65,68,70,72,73,78,82,86,87,88,89,111,112,116,117,130,133,134,136,138,139,
         140,141,142,143,150,151,152,153,157,160,165,168,171,172,173,177,178,179,181,
         182,185,187,191,199,205,208,211,212,215,216]
CANDIDATES = ([f"HUP{n}_phaseII" for n in HUP_N] +
              [f"HUP{n}_phaseII_D01" for n in HUP_N[:20]] +
              [f"Study {i:03d}" for i in range(1, 31)] +
              [f"I001_P{p:03d}_D01" for p in range(1, 40)])


def main():
    c = json.load(open(CRED))
    s = Session(c["username"], c["password"])
    hits = []
    for name in CANDIDATES:
        try:
            ds = s.open_dataset(name)
            labels = ds.get_channel_labels()
            d = ds.get_time_series_details(labels[0])
            sf = d.sample_rate
            dur = (getattr(d, "duration", 0) or d.number_of_samples / sf * 1e6) / 3.6e9
            mtl = [l for l in labels if MTL.match(l.replace(" ", ""))]
            hits.append((name, len(labels), sf, dur, len(mtl), mtl[:6]))
            print(f"OK  {name:22s} {len(labels):3d}ch {sf:6.0f}Hz {dur:6.2f}h  MTL={len(mtl):2d} {mtl[:5]}")
            s.close_dataset(ds)
        except Exception as e:
            msg = str(e).lower()
            if "authenti" in msg or "401" in msg:
                sys.exit(f"AUTH FAILED: {e}")
            # silent on not-found / no-access
    print(f"\n{len(hits)} reachable datasets. Best MTL+long candidates:")
    for name, n, sf, dur, nmtl, ex in sorted(hits, key=lambda h: (-h[4], -h[3]))[:12]:
        print(f"  {name:22s} {nmtl:2d} MTL ch, {dur:.1f} h, {sf:.0f} Hz  e.g. {ex}")


if __name__ == "__main__":
    main()
