# iEEG — Human Test of the Kipnis Sleep-Coordination Mechanism

Phase-1 analysis of open, normative human intracranial EEG to test whether the Kipnis /
Jiang-Xie slow-wave sleep-coordination mechanism appears — and is state-dependent — in human
hippocampus and entorhinal cortex.

Sister project to the Buzsáki mouse pilot (`hpaticmousebuzsakilab`). The mouse pilot showed 50 Hz
vibrotactile drive organized dorsal-hippocampal firing at **theta** (6–10 Hz) but not at the Kipnis
**0.5–4 Hz** slow-wave band — because it had no natural NREM state. This project uses open overnight
iEEG that *does* contain natural N2/N3/REM/wake to close that gap.

- Plan: [docs/PLAN.md](docs/PLAN.md)
- Data: open MNI Open iEEG Atlas (wake) + Multicenter iEEG Sleep Atlas (sleep). No new patient
  access or committee approval required.

## Layout
```
data/raw       downloaded atlas clips (gitignored)
data/derived   tidy per-channel/per-state index + band powers
analysis        one script per module M1–M5
outputs         figures + csv
docs            plan + per-module writeups (DEC4 house style)
env             requirements.txt
```
