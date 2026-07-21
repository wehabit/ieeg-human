"""
iEEG.org (IEEG Portal) downloader — continuous full-night iEEG for this project.

Credentials live in data/ieeg_secret/credentials.json (gitignored), never in code.

Commands:
  python analysis/ieeg_portal.py test
      Verify login works (tries a few known public demo datasets).
  python analysis/ieeg_portal.py info  <DATASET_NAME>
      Print channel labels, sample rate, and total duration for a dataset.
  python analysis/ieeg_portal.py pull  <DATASET_NAME> --start 0 --dur 300 \
      [--chans LAH1,LAH2 | --chan-idx 0,1,2 | --all] --out data/ieeg_portal/seg.npz
      Download a time window (seconds) for chosen channels -> .npz (data, ch_names, sfreq).

Find dataset names by browsing https://www.ieeg.org after logging in (each study/session has
an ID like 'I001_P034_D01'). Paste that ID as DATASET_NAME.
"""
import argparse, json, os, sys
import numpy as np
from ieeg.auth import Session

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CRED = os.path.join(ROOT, "data", "ieeg_secret", "credentials.json")
# a few historically public demo datasets, tried by `test`
DEMO = ["I001_P034_D01", "I521_A0001_D002", "I521_A0001_D001", "Study 005"]


def load_session():
    if not os.path.exists(CRED):
        sys.exit(f"Missing {CRED}. Create it with {{\"username\":..., \"password\":...}}")
    c = json.load(open(CRED))
    return Session(c["username"], c["password"])


def details(ds):
    """Return (labels, sfreq, duration_sec)."""
    labels = ds.get_channel_labels()
    d = ds.get_time_series_details(labels[0])
    sf = d.sample_rate
    dur_usec = getattr(d, "duration", None)
    dur_sec = (dur_usec / 1e6) if dur_usec else (d.number_of_samples / sf)
    return labels, sf, dur_sec


def cmd_test(_):
    s = load_session()
    print("Logged in as:", json.load(open(CRED))["username"])
    for name in DEMO:
        try:
            ds = s.open_dataset(name)
            labels, sf, dur = details(ds)
            print(f"  OK  '{name}': {len(labels)} ch, {sf:.0f} Hz, {dur/3600:.2f} h")
            s.close_dataset(ds)
            print("\nAUTH OK. Now browse ieeg.org, copy a dataset ID, and run `info <ID>`.")
            return
        except Exception as e:
            msg = str(e).lower()
            if "authenti" in msg or "401" in msg or "unauthor" in msg:
                sys.exit(f"AUTH FAILED — check username/password. ({e})")
            print(f"  (no access to demo '{name}': {type(e).__name__})")
    print("\nAUTH appears OK (no auth error), but none of the demo datasets were reachable.\n"
          "That's fine — browse ieeg.org, copy a real dataset ID you can see, and run `info <ID>`.")


def cmd_info(a):
    s = load_session()
    ds = s.open_dataset(a.dataset)
    labels, sf, dur = details(ds)
    print(f"Dataset: {a.dataset}")
    print(f"  sample rate: {sf:.1f} Hz | duration: {dur/3600:.2f} h ({dur:.0f} s)")
    print(f"  channels ({len(labels)}):")
    for i, l in enumerate(labels):
        print(f"    [{i:3d}] {l}")
    s.close_dataset(ds)


def cmd_pull(a):
    s = load_session()
    ds = s.open_dataset(a.dataset)
    labels, sf, dur = details(ds)
    if a.all:
        idx = list(range(len(labels)))
    elif a.chan_idx:
        idx = [int(i) for i in a.chan_idx.split(",")]
    elif a.chans:
        want = a.chans.split(",")
        idx = ds.get_channel_indices(want)
    else:
        sys.exit("choose --all, --chans, or --chan-idx")
    start_usec = int(a.start * 1e6)
    dur_usec = int(a.dur * 1e6)
    print(f"pulling {len(idx)} ch, {a.start}-{a.start+a.dur}s from '{a.dataset}' ({sf:.0f} Hz)...")
    data = ds.get_data(start_usec, dur_usec, idx)   # samples x channels
    ch_names = [labels[i] for i in idx]
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    np.savez_compressed(a.out, data=data.astype(np.float32), ch_names=ch_names,
                        sfreq=sf, start_sec=a.start)
    print(f"saved {data.shape} -> {a.out}")
    s.close_dataset(ds)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("test").set_defaults(func=cmd_test)
    p = sub.add_parser("info"); p.add_argument("dataset"); p.set_defaults(func=cmd_info)
    p = sub.add_parser("pull"); p.add_argument("dataset")
    p.add_argument("--start", type=float, default=0.0); p.add_argument("--dur", type=float, default=300.0)
    p.add_argument("--chans"); p.add_argument("--chan-idx"); p.add_argument("--all", action="store_true")
    p.add_argument("--out", default=os.path.join(ROOT, "data", "ieeg_portal", "segment.npz"))
    p.set_defaults(func=cmd_pull)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
