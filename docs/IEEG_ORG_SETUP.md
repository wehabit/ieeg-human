# iEEG.org (IEEG Portal) — download setup

Status: **working.** Login verified, `info` and `pull` tested end-to-end (downloaded a real
5000 Hz segment). This is the source for **continuous full-night** iEEG that the atlas/Falach/Zurich
sets couldn't give (they were clips).

## Credentials
Stored in `data/ieeg_secret/credentials.json` (gitignored, never committed):
```json
{"username": "<IEEG_USERNAME>", "password": "<IEEG_PASSWORD>"}
```
Never place a real account identifier, password, access token, or credential history in this file.

## Client
`ieegpy` (the real one) installed from GitHub `ieeg-portal/ieegpy` into `.venv`.
(The PyPI package named `ieeg` is a *different, unrelated* package — don't use it.)

## Commands
```bash
source .venv/bin/activate
python analysis/ieeg_portal.py test                 # verify login
python analysis/ieeg_portal.py info  <DATASET_ID>   # channels, sample rate, duration
python analysis/ieeg_portal.py pull  <DATASET_ID> --start 0 --dur 300 \
       --chans LAH1,LAH2,LA1 --out data/ieeg_portal/night1.npz
# or --chan-idx 0,1,2  or  --all
```
Downloads save as `.npz` (data [samples×channels], ch_names, sfreq, start_sec).
Times are in **seconds**.

## Verified datasets on THIS account (probed 2026-07)

`analysis/ieeg_find_datasets.py` found 80 reachable datasets. Best continuous MTL + sleep:

| Dataset | Duration | Rate | MTL depth |
|---|---|---|---|
| **HUP165_phaseII** | **459 h (~19 days)** | **1024 Hz** | LA/LB/LC (amygdala + hippocampal), 47 ch |
| HUP157_phaseII | ~190 h | 1024 Hz | 36 ch |
| HUP130_phaseII | ~126 h | 1024 Hz | 32 ch |
| HUP211_phaseII | ~210 h | 512 Hz | 48 ch |
| HUP171/177/182_phaseII | 190–286 h | 512 Hz | 31–35 ch |

These are the **same HUP patients as our Pennsieve atlas** (Slow-power–Ripple-coupling), but the full continuous
phase-II monitoring — many full nights. Use 1024 Hz ones (HUP165/157/130) for ripples. Rerun the
probe anytime with `python analysis/ieeg_find_datasets.py`.

## Finding a good sleep + MTL dataset
There is no public catalog repo, so browse the portal:
1. Log in at https://www.ieeg.org → **Tools → Data** (or the dataset browser).
2. **Most multi-day epilepsy-monitoring datasets contain sleep** (patients are recorded 24/7 for
   days). The demo `I001_P034_D01` is 32 h continuous — a full circadian span.
3. Pick one with **depth / SEEG electrodes in MTL** — channel names like `LAH*/RAH*` (hippocampus),
   `LA*/RA*` (amygdala), `LEC*/REC*` (entorhinal), `LPHG*` (parahippocampal). Grid/strip names
   (`Grid*`) are neocortical ECoG — skip for MTL work.
4. Copy the dataset ID and run `info <ID>` to check channels + duration before pulling.

## Practical notes on size
5000 Hz × many channels × hours = large. Don't pull a whole night at once. Pull **windows**
(e.g. 5–10 min blocks across the night, or the NREM periods), a handful of MTL channels at a time.
Then feed the `.npz` into an Ripple-coupling/Nesting-timecourse-style pipeline (bipolar → SO phase → spindle/ripple/fast-ripple
coupling), now with true overnight coverage and state transitions.

## Useful ieeg-portal GitHub repos
- **ieegpy** — the Python client (installed).
- **ieeg-2edf** — stream portal data straight to EDF (if you prefer EDF over npz).
- **ieeg-analytic-tools** — analysis helpers.
- EDF/Nicolet/Nervus/Persyst **readers** — for local clinical files, not the portal.
