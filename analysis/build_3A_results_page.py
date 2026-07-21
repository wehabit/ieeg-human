import base64, os, json

SCR = "/private/tmp/claude-501/-Users-paris-Documents-Buzsakli-Lab-Github/6c568d4c-5b04-4605-bece-490c8d0e3e2b/scratchpad"
OUTD = "/Users/paris/Documents/iEEG/outputs/results_3A_tutorial_style"

def b64(name):
    with open(os.path.join(SCR, name), "rb") as f:
        return base64.b64encode(f.read()).decode()

IMGS = {
    "mtl": b64("HUP165_phaseII_3A_results.jpg"),
    "ied": b64("HUP165_phaseII_iedmasked_3A_results.jpg"),
    "ctx": b64("HUP165_phaseII_cortical_3A_results.jpg"),
}

CSS = """
<title>Test 3A — spindle–heart infraslow coupling in HUP165</title>
<style>
:root{
  --teal:#0f9b96; --rose:#d6415f; --amber:#c47d15; --violet:#6d5f9c;
  --bg:#f6f9f9; --surface:#ffffff; --surface-2:#eef3f3;
  --text:#1c2226; --muted:#5b686e; --border:#d9e3e3;
  --chip-null-bg:#e8eeee; --chip-null-fg:#4d5a60;
  --chip-marg-bg:#fbf0dc; --chip-marg-fg:#8a5806;
  --chip-ok-bg:#dff0ef; --chip-ok-fg:#0b6f6b;
  --shadow:0 1px 2px rgba(20,35,35,.06),0 8px 24px rgba(20,35,35,.05);
}
@media (prefers-color-scheme:dark){
  :root{
    --bg:#121618; --surface:#191f21; --surface-2:#212a2c;
    --text:#e4ecec; --muted:#9aa9ad; --border:#2a3437;
    --teal:#2fbdb6; --rose:#e8687f; --amber:#d79a34;
    --chip-null-bg:#232c2e; --chip-null-fg:#a5b3b7;
    --chip-marg-bg:#3a2f18; --chip-marg-fg:#e0a94a;
    --chip-ok-bg:#123331; --chip-ok-fg:#4fd0c8;
    --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px rgba(0,0,0,.28);
  }
}
:root[data-theme="dark"]{
  --bg:#121618; --surface:#191f21; --surface-2:#212a2c;
  --text:#e4ecec; --muted:#9aa9ad; --border:#2a3437;
  --teal:#2fbdb6; --rose:#e8687f; --amber:#d79a34;
  --chip-null-bg:#232c2e; --chip-null-fg:#a5b3b7;
  --chip-marg-bg:#3a2f18; --chip-marg-fg:#e0a94a;
  --chip-ok-bg:#123331; --chip-ok-fg:#4fd0c8;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px rgba(0,0,0,.28);
}
:root[data-theme="light"]{
  --bg:#f6f9f9; --surface:#ffffff; --surface-2:#eef3f3;
  --text:#1c2226; --muted:#5b686e; --border:#d9e3e3;
  --teal:#0f9b96; --rose:#d6415f; --amber:#c47d15;
  --chip-null-bg:#e8eeee; --chip-null-fg:#4d5a60;
  --chip-marg-bg:#fbf0dc; --chip-marg-fg:#8a5806;
  --chip-ok-bg:#dff0ef; --chip-ok-fg:#0b6f6b;
  --shadow:0 1px 2px rgba(20,35,35,.06),0 8px 24px rgba(20,35,35,.05);
}

*{box-sizing:border-box;}
body{margin:0;background:var(--bg);color:var(--text);
  font-family:"Iowan Old Style","Palatino Linotype",Palatino,Georgia,serif;
  font-size:17px;line-height:1.62;-webkit-font-smoothing:antialiased;}
.wrap{max-width:1180px;margin:0 auto;padding:56px 24px 80px;display:flex;flex-direction:column;gap:44px;}
.prose{max-width:66ch;}

h1,h2,h3,.mono,.chip,th{font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;}
h1{font-size:clamp(1.85rem,3.6vw,2.7rem);line-height:1.12;letter-spacing:-.021em;font-weight:700;
   margin:.35em 0 .3em;text-wrap:balance;}
h2{font-size:1.32rem;letter-spacing:-.012em;font-weight:650;margin:0 0 .15em;text-wrap:balance;}
p{margin:0 0 1em;}
.mono,code,td.num,th.num{font-family:ui-monospace,SFMono-Regular,"SF Mono",Menlo,Consolas,monospace;
  font-variant-numeric:tabular-nums;}
.eyebrow{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.735rem;letter-spacing:.13em;
  text-transform:uppercase;color:var(--muted);margin:0;}

header .lede{font-size:1.12rem;color:var(--muted);max-width:64ch;margin-top:.2em;}

.banner{border-left:3px solid var(--amber);background:var(--surface);border-radius:0 10px 10px 0;
  padding:18px 22px;box-shadow:var(--shadow);}
.banner p{margin:0;}
.banner strong{color:var(--amber);}

.tablewrap{overflow-x:auto;background:var(--surface);border:1px solid var(--border);
  border-radius:12px;box-shadow:var(--shadow);}
table{border-collapse:collapse;width:100%;min-width:620px;font-size:.94rem;}
th,td{text-align:left;padding:13px 18px;border-bottom:1px solid var(--border);}
thead th{font-size:.7rem;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);font-weight:600;
  background:var(--surface-2);}
tbody tr:last-child td{border-bottom:none;}
td.num,th.num{text-align:right;}
tr.highlight td{background:color-mix(in srgb,var(--teal) 7%,transparent);}
td:first-child{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.86rem;}

.chip{display:inline-block;font-size:.7rem;font-weight:650;letter-spacing:.05em;text-transform:uppercase;
  padding:3px 9px;border-radius:999px;white-space:nowrap;}
.chip.null{background:var(--chip-null-bg);color:var(--chip-null-fg);}
.chip.marg{background:var(--chip-marg-bg);color:var(--chip-marg-fg);}
.chip.ok{background:var(--chip-ok-bg);color:var(--chip-ok-fg);}

.fig{background:var(--surface);border:1px solid var(--border);border-radius:14px;overflow:hidden;
  box-shadow:var(--shadow);}
.fig .head{padding:20px 24px 16px;display:flex;flex-wrap:wrap;gap:12px 18px;align-items:baseline;
  border-bottom:1px solid var(--border);}
.fig .head .grow{flex:1 1 260px;}
.fig figure{margin:0;background:#fff;}
.fig img{display:block;width:100%;height:auto;}
.fig .read{padding:16px 24px 20px;color:var(--muted);font-size:.95rem;margin:0;}
.fig .read b{color:var(--text);font-weight:600;}

.stats{display:flex;flex-wrap:wrap;gap:8px 22px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
  font-size:.8rem;font-variant-numeric:tabular-nums;color:var(--muted);}
.stats b{color:var(--text);font-weight:600;}

.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(290px,1fr));gap:18px;}
.card{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:20px 22px;
  box-shadow:var(--shadow);}
.card h3{margin:.1em 0 .5em;font-size:1rem;letter-spacing:-.005em;}
.card ul{margin:0;padding-left:1.1em;}
.card li{margin-bottom:.5em;}
.card li::marker{color:var(--teal);}
.card.warn li::marker{color:var(--rose);}

.legend{display:flex;flex-wrap:wrap;gap:6px 20px;font-size:.85rem;color:var(--muted);align-items:center;}
.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:6px;vertical-align:baseline;}
hr{border:none;border-top:1px solid var(--border);margin:0;}
a{color:var(--teal);}
:focus-visible{outline:2px solid var(--teal);outline-offset:3px;border-radius:4px;}
@media (prefers-reduced-motion:no-preference){.fig{transition:box-shadow .2s ease;}}
</style>
"""

def fig_block(key, eyebrow, title, chip_cls, chip_txt, stats, read):
    return f"""
<section class="fig">
  <div class="head">
    <div class="grow">
      <p class="eyebrow">{eyebrow}</p>
      <h2>{title}</h2>
    </div>
    <span class="chip {chip_cls}">{chip_txt}</span>
  </div>
  <div class="head" style="border-bottom:1px solid var(--border);padding-top:14px;padding-bottom:14px;">
    <div class="stats">{stats}</div>
  </div>
  <figure><img src="data:image/jpeg;base64,{IMGS[key]}" alt="{title} — five-panel result figure"></figure>
  <p class="read">{read}</p>
</section>"""

BODY = f"""
<div class="wrap">

<header>
  <p class="eyebrow">iEEG · HUP165_phaseII · 140 min · 93% NREM</p>
  <h1>Does spindling rise and fall with the heart every ~50 seconds?</h1>
  <p class="lede">Test 3A of the locus-coeruleus infraslow hypothesis, run on real intracranial data
  and drawn in the same visual language as the synthetic teaching figure — so the measured result can
  be read against the idealized one panel for panel.</p>
</header>

<div class="banner">
  <p><strong>Status: suggestive, not established.</strong> One subject, one night. The cortical result
  clears the per-bin chance threshold by 0.104 vs 0.101 and <em>fails</em> the band-corrected
  threshold of 0.184. Region, frequency and control all point the same way — but a margin that thin
  is exactly what a single subject cannot settle.</p>
</div>

<section>
  <p class="eyebrow" style="margin-bottom:10px;">Coherence at 0.02 Hz, by channel set</p>
  <div class="tablewrap">
    <table>
      <thead>
        <tr><th>Channel set</th><th class="num">Spindle peak</th><th class="num">Coherence @0.02 Hz</th>
        <th class="num">Chance</th><th>Verdict</th></tr>
      </thead>
      <tbody>
        <tr><td>MTL depth</td><td class="num">12.5 Hz</td><td class="num">0.080</td><td class="num">0.101</td>
            <td><span class="chip null">not significant</span></td></tr>
        <tr><td>MTL depth · IED-masked</td><td class="num">12.5 Hz</td><td class="num">0.070</td><td class="num">0.101</td>
            <td><span class="chip null">not significant</span></td></tr>
        <tr class="highlight"><td>Lateral neocortex · IED-masked</td><td class="num">12.9 Hz</td>
            <td class="num">0.104</td><td class="num">0.101</td>
            <td><span class="chip marg">marginal · p≈0.05</span></td></tr>
        <tr><td>SWA control (cortical)</td><td class="num">—</td><td class="num">0.006</td><td class="num">0.101</td>
            <td><span class="chip ok">null, as predicted</span></td></tr>
      </tbody>
    </table>
  </div>
  <p class="legend" style="margin-top:12px;">
    <span><span class="dot" style="background:var(--teal)"></span>spindle / sigma power</span>
    <span><span class="dot" style="background:var(--rose)"></span>heart rate</span>
    <span><span class="dot" style="background:var(--amber)"></span>slow-wave activity (control)</span>
    <span>Heart data throughout: 12,811 beats · 92 bpm · 0.02% dropped</span>
  </p>
</section>

{fig_block("mtl","Run 1 · mesial temporal depth","Spindles measured in the hippocampus and amygdala",
  "null","not significant",
  "<span><b>0.080</b> coherence @0.02 Hz</span><span>chance <b>0.101</b></span>"
  "<span>band peak <b>0.0195 Hz</b></span><span>spindle peak <b>12.5 Hz</b></span>",
  "Panel ③ shows the two slow lines never lock together, and panel ④ never leaves the grey chance band. "
  "Panel ① is dominated by sigma transients reaching 2500% of median — the reason the next run exists.")}

{fig_block("ied","Run 2 · same channels, epileptiform activity removed","Ruling out epileptic spikes as the cause",
  "null","not significant",
  "<span><b>0.070</b> coherence @0.02 Hz</span><span>chance <b>0.101</b></span>"
  "<span><b>12.1%</b> of samples masked</span><span>top-50 power bins cut by <b>0.4%</b></span>",
  "Masking 20–80 Hz spike activity barely touched the transients (peak 5240% → 5220%), so they are "
  "<b>not epileptiform</b> — they are large sigma-band events. The negative here is a real property of "
  "the mesial temporal signal, not contamination.")}

{fig_block("ctx","Run 3 · lateral neocortical contacts","Measuring where spindles actually live",
  "marg","marginal · p≈0.05",
  "<span><b>0.104</b> coherence @0.02 Hz</span><span>chance <b>0.101</b></span>"
  "<span>band-corrected chance <b>0.184</b></span><span>SWA control <b>0.006</b></span>",
  "Spindles are thalamo<b>cortical</b>, and Lecci's effect is parietal-maximal — so this is the region the "
  "hypothesis actually predicts. Coherence rises above chance, but by a hair, and panel ⑤ still shows "
  "<b>no isolated 0.02 Hz peak</b>, which in the published figure is unmistakable.")}

<div class="grid2">
  <div class="card">
    <h3>Three things that line up with the literature</h3>
    <ul>
      <li><b>Cortex beats mesial temporal</b> (0.104 vs 0.070) — the predicted direction, since spindles
      are thalamocortical.</li>
      <li><b>The peak sits at 0.0195 Hz</b> in every single run. Lecci reports 0.019 Hz in humans.</li>
      <li><b>Slow-wave activity stays null</b> at 0.006 — roughly 19× below sigma, reproducing the
      sigma-yes / SWA-no dissociation of Lecci Fig. 1. This is a control on real data, which is far
      stronger evidence the pipeline works than the synthetic tests were.</li>
    </ul>
  </div>
  <div class="card warn">
    <h3>Why this still isn't a finding</h3>
    <ul>
      <li>The cortical value clears chance by <b>0.003</b> and fails band correction — the same
      marginal-crossing pattern I refused to trust when the mesial result produced it.</li>
      <li><b>No isolated 0.02 Hz bump</b> in either power spectrum (panel ⑤).</li>
      <li><b>n = 1</b>, single night, single window, an epilepsy patient at 92 bpm — tachycardic for
      NREM, which may itself blunt autonomic modulation.</li>
    </ul>
  </div>
</div>

<hr>

<section class="prose">
  <p class="eyebrow">What settles it</p>
  <p>The cohort. Twenty-five HUP subjects carry both depth electrodes and EKG, with OpenNeuro
  <span class="mono">ds003848</span> as an independent replication — different site, cortical grids,
  and real EOG/EMG sleep staging. If a cortical coherence near 0.1 at 0.019 Hz is real, it appears as a
  consistent across-subject effect. If it is noise, it scatters around chance. One subject cannot tell
  those apart; twenty-five can.</p>
</section>

</div>
"""

out = os.path.join(SCR, "hup165_3A_results.html")
with open(out, "w") as f:
    f.write(CSS + BODY)
print(out, f"{os.path.getsize(out)/1024:.0f} KB")
