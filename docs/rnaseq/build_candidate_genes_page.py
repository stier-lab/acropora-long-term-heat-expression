#!/usr/bin/env python3
"""Build docs/team_summary/LTH_candidate_genes.html from candidate_genes_reference.csv.

Run from the repo root:  python3 docs/rnaseq/build_candidate_genes_page.py
The CSV is the source of truth for the full gene list; the summary sections,
prediction sketches and timeline are written here. Rebuild after editing either.
No dependencies beyond Python 3.
"""
import csv
import html
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "docs" / "rnaseq" / "candidate_genes_reference.csv"
OUT = ROOT / "docs" / "team_summary" / "LTH_candidate_genes.html"
UPDATED = "19 November 2026"

e = html.escape

# Gene sets: (name, genes, D1 prediction, D10 prediction)
SETS = [
    ("Early wound response",
     "c-Fos, c-Jun, JNK and ERK kinases, Toll-like receptors, NF-κB, TRAF6, lectins, complement",
     "Wounded differ from unwounded at both temperatures. Heat changed this little in <i>A. aspera</i>; c-Jun rose more under heat.",
     "A second rise is possible: c-Fos, ERK2 and TRAF6 rose again at day 10 in <i>A. aspera</i>."),
    ("Skeleton building",
     "Galaxin 2, Saarp2, mucin, Usomp, carbonic anhydrase 2, collagens, calicoblast markers",
     "Little change yet. These genes rose in mid-to-late regrowth in <i>A. muricata</i>.",
     "The key test. If heat blocks regrowth, wounding raises these genes less at 31 °C than at 28 °C."),
    ("Heat stress and energy",
     "Ribosomal and mitochondrial genes, antioxidant genes, heat shock proteins",
     "Differ between temperatures in wounded and unwounded fragments alike.",
     "Same. Chronic 32 °C did not raise HSP70 or HSP90 in <i>A. aspera</i>, so expect a weak heat shock signal."),
]

# Studies: (short name, url, one-line design, one-line finding)
STUDIES = [
    ("Xu et al. 2023", "https://doi.org/10.3389/fevo.2022.979278",
     "<i>A. millepora</i>, aquaria, 6 h to 10 days",
     "Wnt, FGF, ADAMTS, galaxin and c-Fos peak 6 h after wounding."),
    ("van de Water et al. 2015a", "https://doi.org/10.1111/mec.13257",
     "<i>A. aspera</i>, on the reef, 16 h to 10 days",
     "Immune genes respond in waves; c-Fos rises at days 4 and 10; c-Jun never changes."),
    ("van de Water et al. 2015b", "https://doi.org/10.1007/s10750-015-2243-z",
     "<i>A. aspera</i>, heat then wounding, 24 h",
     "Heat changes the wound response little; c-Jun rises more under heat. Two aquaria per temperature."),
    ("Bonesso et al. 2017", "https://doi.org/10.7717/peerj.3719",
     "<i>A. aspera</i>, tips removed, 26 vs 32 °C, 12 days",
     "Heat halves tip regrowth; HSP70 and HSP90 do not rise."),
    ("Han et al. 2025", "https://doi.org/10.1038/s42003-025-08089-6",
     "<i>A. muricata</i>, tips cut, to 39 days; single cells at day 21",
     "Skeleton-building genes rise in mid-to-late regrowth, which peaks at days 15–27."),
]

CAUTIONS = [
    "Heat also cut growth in unwounded fragments (10.2% vs 16.0%), so heated versus ambient wounded fragments alone cannot show that heat blocks regrowth.",
    "Heated fragments lost most of their symbionts. Fewer symbiont-hosting cells can shift host expression without any gene changing its regulation, so check key results against each sample's symbiont density.",
    "Our D1 sample (24 h) may miss the earliest wound peak, which came at 6 h in <i>A. millepora</i>.",
    "The skeleton-building genes come from one single-cell study at day 21 in a slower-regenerating species.",
    "Fluorescent protein results in the wounding studies are protein, not RNA, so they do not predict our data.",
]


def interaction_svg(kind):
    """Small sketch of skeleton-gene expression: unwounded vs wounded at 28 and 31 C."""
    # x: unwounded=70, wounded=230; y grows downward (higher expression = smaller y)
    if kind == "block":
        amb = (140, 45)
        heat = (158, 146)
    else:
        amb = (140, 45)
        heat = (165, 70)
    return f"""<svg viewBox="0 0 300 210" role="img" aria-labelledby="t-{kind}">
  <title id="t-{kind}">{'Wounding raises skeleton genes at 28 °C but barely at 31 °C' if kind == 'block' else 'Wounding raises skeleton genes by the same amount at both temperatures; 31 °C is lower overall'}</title>
  <line x1="40" y1="20" x2="40" y2="175" class="ax"/>
  <line x1="40" y1="175" x2="280" y2="175" class="ax"/>
  <text x="70" y="196" class="tick" text-anchor="middle">Unwounded</text>
  <text x="230" y="196" class="tick" text-anchor="middle">Wounded</text>
  <text x="18" y="100" class="tick" text-anchor="middle" transform="rotate(-90 18 100)">Skeleton genes</text>
  <line x1="70" y1="{amb[0]}" x2="230" y2="{amb[1]}" class="l28"/>
  <circle cx="70" cy="{amb[0]}" r="5" class="p28"/><circle cx="230" cy="{amb[1]}" r="5" class="p28"/>
  <line x1="70" y1="{heat[0]}" x2="230" y2="{heat[1]}" class="l31"/>
  <circle cx="70" cy="{heat[0]}" r="5" class="p31"/><circle cx="230" cy="{heat[1]}" r="5" class="p31"/>
  <text x="238" y="{amb[1] + 4}" class="lab28">28 °C</text>
  <text x="238" y="{heat[1] + 4}" class="lab31">31 °C</text>
</svg>"""


def timeline_svg():
    """Days 0-15 on one axis: our phenotype, our samples, and when earlier studies saw responses."""
    x0, x1 = 170, 860
    def X(d):
        return x0 + (x1 - x0) * d / 15
    rows = []
    y = 34
    def row(label, body, sub=""):
        nonlocal y
        s = f'<text x="0" y="{y + 4}" class="rl">{label}</text>'
        if sub:
            s += f'<text x="0" y="{y + 19}" class="rs">{sub}</text>'
        s += f'<line x1="{x0}" y1="{y}" x2="{x1}" y2="{y}" class="guide"/>' + body
        rows.append(s)
        y += 46
    # Our phenotype
    row("Our corals", (
        f'<circle cx="{X(1)}" cy="34" r="6" class="both"/>'
        f'<text x="{X(1)}" y="20" class="ev" text-anchor="middle">wound covered</text>'
        f'<rect x="{X(6)}" y="29" width="{X(8) - X(6)}" height="10" rx="5" class="both"/>'
        f'<text x="{X(7)}" y="20" class="ev" text-anchor="middle">axial corallite</text>'
        f'<rect x="{X(12)}" y="29" width="{X(15) - X(12)}" height="10" rx="5" class="amb"/>'
        f'<text x="{X(15)}" y="20" class="ev" text-anchor="end">radial corallites: 12 of 12 at 28 °C, 4 of 12 at 31 °C</text>'))
    # Our samples
    marks = "".join(
        f'<line x1="{X(d)}" y1="{y - 9}" x2="{X(d)}" y2="{y + 9}" class="samp"/>'
        f'<text x="{X(d) + 6}" y="{y + 4}" class="sl">D{d}</text>' for d in (1, 3, 10))
    row("Our RNA samples", marks)
    row("Xu 2023", f'<circle cx="{X(0.25)}" cy="{y}" r="6" class="prior"/><text x="{X(0.25) + 12}" y="{y + 4}" class="ev">wound genes peak at 6 h</text>', "A. millepora")
    vw = "".join(f'<rect x="{X(a)}" y="{y - 5}" width="{X(b) - X(a)}" height="10" rx="5" class="prior"/>' for a, b in ((0.7, 2), (3.6, 4.4), (9.6, 10.4)))
    row("van de Water 2015a", vw + f'<text x="{X(5.2)}" y="{y + 4}" class="ev">immune genes in waves</text>', "A. aspera")
    row("Bonesso 2017", f'<rect x="{X(0)}" y="{y - 5}" width="{X(12) - X(0)}" height="10" rx="5" class="priorlight"/><text x="{X(12) + 8}" y="{y + 4}" class="ev">regrowth done at 26 °C</text>', "A. aspera")
    row("Han 2025", f'<rect x="{X(13)}" y="{y - 5}" width="{X(15) - X(13) + 6}" height="10" rx="5" class="prior"/><text x="{X(13) - 8}" y="{y + 4}" class="ev" text-anchor="end">regrowth peaks at days 15–27, past this axis</text>', "A. muricata")
    ticks = "".join(
        f'<line x1="{X(d)}" y1="{y - 18}" x2="{X(d)}" y2="{y - 12}" class="ax"/>'
        f'<text x="{X(d)}" y="{y}" class="tick" text-anchor="middle">{d}</text>' for d in (0, 1, 3, 5, 8, 10, 12, 15))
    axis = (f'<line x1="{x0}" y1="{y - 18}" x2="{x1}" y2="{y - 18}" class="ax"/>' + ticks +
            f'<text x="{(x0 + x1) / 2}" y="{y + 20}" class="tick" text-anchor="middle">Days after wounding</text>')
    h = y + 30
    return (f'<svg viewBox="0 0 900 {h}" role="img" aria-label="Timeline from day 0 to 15 comparing our coral milestones and RNA sampling days with when earlier studies saw gene responses">'
            + "".join(rows) + axis + "</svg>")


def gene_rows(rows):
    out = []
    for r in rows:
        checked = r.get("source_check", "")
        ok = checked and not checked.startswith("Not")
        out.append(
            "<tr data-text=\"{t}\"><td><b>{g}</b><span class=\"grp\">{c}</span></td>"
            "<td>{k}</td><td>{u}</td><td>{cav}</td><td class=\"src\">{s}<span class=\"{cls}\">{chk}</span></td></tr>".format(
                t=e(" ".join(r.values()).lower()), g=e(r["gene_or_family"]), c=e(r["candidate_group"]),
                k=e(r["key_evidence"]), u=e(r["lth_use"]), cav=e(r["caveat"]), s=e(r["primary_sources"]),
                cls="chk" if ok else "nochk", chk="Checked against the paper" if ok else "Not rechecked"))
    return "\n".join(out)


def main():
    rows = list(csv.DictReader(open(CSV, newline="")))
    sets = "\n".join(
        f'<tr{" class=\"key\"" if n == "Skeleton building" else ""}><th scope="row">{n}<span class="grp">{g}</span></th><td>{d1}</td><td class="d10">{d10}</td></tr>'
        for n, g, d1, d10 in SETS)
    studies = "\n".join(
        f'<li><a href="{u}">{n}</a> <span class="dz">{d}.</span> {f}</li>' for n, u, d, f in STUDIES)
    cautions = "\n".join(f"<li>{c}</li>" for c in CAUTIONS)
    nchk = sum(1 for r in rows if r.get("source_check", "") and not r["source_check"].startswith("Not"))
    page = TEMPLATE.format(updated=UPDATED, sets=sets, studies=studies, cautions=cautions,
                           genes=gene_rows(rows), n=len(rows), nchk=nchk,
                           block=interaction_svg("block"), stress=interaction_svg("stress"),
                           timeline=timeline_svg())
    OUT.write_text(page, encoding="utf-8")
    print(f"Wrote {OUT.relative_to(ROOT)} ({len(rows)} gene families)")


TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>What earlier work predicts for our expression data</title>
<!-- Built by docs/rnaseq/build_candidate_genes_page.py from candidate_genes_reference.csv. Do not edit by hand. -->
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Instrument+Sans:wght@400;500;600;700&family=Literata:ital,opsz,wght@0,7..72,400;0,7..72,600;1,7..72,400&display=swap" rel="stylesheet">
<style>
  :root {{ --ink:#16233B; --sub:#55657A; --rule:#DCE2E8; --tint:#F3F6F8; --sand:#EFE9DD;
    --amb:#56B4E9; --amb-ink:#0B5A86; --heat:#D55E00; --heat-ink:#A84800; --prior:#7C8A9C;
    --serif:"Literata",Georgia,serif; --sans:"Instrument Sans","Helvetica Neue",Arial,sans-serif; }}
  * {{ box-sizing:border-box; }}
  html {{ font-size:17px; }}
  body {{ margin:0; background:#fff; color:var(--ink); font-family:var(--serif); line-height:1.65; -webkit-font-smoothing:antialiased; }}
  .col {{ max-width:42rem; margin:0 auto; padding:0 1.1rem; }}
  .wide {{ max-width:60rem; margin:0 auto; padding:0 1.1rem; }}
  header {{ padding-top:2.6rem; }}
  nav {{ font-family:var(--sans); font-size:0.9rem; display:flex; gap:1.2rem; margin-bottom:2rem; }}
  h1 {{ font-family:var(--sans); font-weight:700; font-size:2.25rem; line-height:1.12; letter-spacing:-0.018em; margin:0 0 0.9rem; }}
  h2 {{ font-family:var(--sans); font-size:1.2rem; font-weight:600; margin:0 0 0.4rem; }}
  section {{ margin-top:3.2rem; }}
  p {{ margin:0 0 1rem; }}
  .lead {{ font-size:1.16rem; line-height:1.55; }}
  .meta {{ font-family:var(--sans); color:var(--sub); font-size:0.88rem; }}
  .intro {{ color:var(--sub); font-size:0.98rem; }}
  a {{ color:var(--amb-ink); text-underline-offset:2px; }}
  a:hover {{ color:var(--ink); }}
  :focus-visible {{ outline:2px solid var(--amb-ink); outline-offset:3px; border-radius:2px; }}

  /* The one bold element: the D10 test drawn both ways */
  .test {{ margin-top:2.4rem; background:var(--sand); border-radius:6px; padding:1.8rem 1.6rem 1.4rem; }}
  .test h2 {{ font-size:1.35rem; }}
  .pair {{ display:grid; grid-template-columns:1fr 1fr; gap:1.6rem; margin-top:1rem; }}
  .pair figure {{ margin:0; }}
  .pair h3 {{ font-family:var(--sans); font-size:1.02rem; font-weight:600; margin:0 0 0.2rem; }}
  .pair figcaption {{ font-family:var(--sans); font-size:0.88rem; line-height:1.5; color:var(--ink); margin-top:0.2rem; }}
  .pair svg {{ width:100%; height:auto; display:block; }}
  svg .ax {{ stroke:var(--ink); stroke-width:1.2; }}
  svg .tick {{ font-family:var(--sans); font-size:12px; fill:var(--sub); }}
  svg .l28 {{ stroke:var(--amb); stroke-width:3.5; }}
  svg .l31 {{ stroke:var(--heat); stroke-width:3.5; stroke-dasharray:7 5; }}
  svg .p28 {{ fill:var(--amb); }} svg .p31 {{ fill:var(--heat); }}
  svg .lab28 {{ font-family:var(--sans); font-size:13px; font-weight:600; fill:var(--amb-ink); }}
  svg .lab31 {{ font-family:var(--sans); font-size:13px; font-weight:600; fill:var(--heat-ink); }}
  .test .why {{ font-family:var(--sans); font-size:0.9rem; color:var(--sub); margin:1.1rem 0 0; }}

  /* Timeline */
  .tl .scroll {{ overflow-x:auto; }}
  .tl svg {{ width:100%; min-width:760px; height:auto; display:block; margin-top:0.8rem; }}
  svg .rl {{ font-family:var(--sans); font-size:13.5px; font-weight:600; fill:var(--ink); }}
  svg .rs {{ font-family:var(--sans); font-size:11.5px; font-style:italic; fill:var(--sub); }}
  svg .guide {{ stroke:var(--rule); stroke-width:1; }}
  svg .ev {{ font-family:var(--sans); font-size:12px; fill:var(--sub); }}
  svg .both {{ fill:var(--ink); }} svg .amb {{ fill:var(--amb); }}
  svg .samp {{ stroke:var(--ink); stroke-width:3; }}
  svg .sl {{ font-family:var(--sans); font-size:12.5px; font-weight:600; fill:var(--ink); }}
  svg .prior {{ fill:var(--prior); }} svg .priorlight {{ fill:#C9D1DB; }}

  /* Predictions table */
  table {{ width:100%; border-collapse:collapse; font-family:var(--sans); font-size:0.9rem; line-height:1.48; }}
  th, td {{ padding:0.75rem 0.7rem 0.75rem 0; border-bottom:1px solid var(--rule); text-align:left; vertical-align:top; }}
  thead th {{ font-weight:600; font-size:0.85rem; color:var(--sub); border-bottom:1.5px solid var(--ink); }}
  tbody th {{ font-weight:600; width:30%; }}
  td.d10 {{ font-weight:500; }}
  tr.key td.d10 {{ color:var(--heat-ink); }}
  .grp {{ display:block; font-weight:400; color:var(--sub); font-size:0.8rem; margin-top:0.2rem; }}

  ul.plain {{ list-style:none; padding:0; margin:0.6rem 0 0; font-family:var(--sans); font-size:0.92rem; }}
  ul.plain li {{ padding:0.55rem 0; border-bottom:1px solid var(--rule); }}
  ul.plain a {{ font-weight:600; }}
  .dz {{ color:var(--sub); }}
  ul.cau {{ padding-left:1.1rem; margin:0.6rem 0 0; }}
  ul.cau li {{ margin-bottom:0.55rem; }}

  details {{ margin-top:3.2rem; border-top:1.5px solid var(--ink); padding-top:1rem; }}
  summary {{ font-family:var(--sans); font-weight:600; font-size:1.2rem; cursor:pointer; }}
  summary span {{ font-weight:400; font-size:0.88rem; color:var(--sub); margin-left:0.4rem; }}
  input[type=search] {{ font-family:var(--sans); font-size:0.95rem; width:100%; max-width:24rem; padding:0.5rem 0.7rem; border:1px solid var(--rule); border-radius:4px; margin:1rem 0 0.6rem; }}
  #genes {{ font-size:0.82rem; }}
  #genes thead th {{ position:sticky; top:0; background:#fff; }}
  td.src {{ color:var(--sub); font-size:0.78rem; }}
  .chk, .nochk {{ display:block; margin-top:0.3rem; font-weight:600; }}
  .chk {{ color:var(--amb-ink); }} .nochk {{ color:#8A94A3; font-weight:500; }}
  footer {{ font-family:var(--sans); font-size:0.84rem; color:var(--sub); margin:3rem auto 3rem; }}
  @media (max-width:700px) {{
    html {{ font-size:16px; }} h1 {{ font-size:1.75rem; }}
    .pair {{ grid-template-columns:1fr; }}
    .test {{ padding:1.3rem 1rem; }}
    tbody th {{ width:auto; }}
  }}
</style>
</head>
<body>
<header class="col">
  <nav><a href="../index.html">Project home</a><a href="LTH_results_summary.html">Results summary</a></nav>
  <h1>What earlier work predicts for our expression data</h1>
  <p class="lead">Heat did not slow healing, but it stopped most wounded fragments from forming radial corallites. Earlier <i>Acropora</i> studies suggest the expression data can tell us why: look at skeleton-building genes 10 days after wounding.</p>
  <p class="meta">RNA-seq from our long-term heating (LTH) experiment, sampled at days 1, 3 and 10. Updated {updated}.</p>
</header>

<div class="wide">
  <section class="test" aria-labelledby="test-h">
    <h2 id="test-h">Two ways the day 10 data could look</h2>
    <p class="intro" style="margin:0">Sketches of expected patterns, not data.</p>
    <div class="pair">
      <figure>
        <h3>Heat blocks regrowth</h3>
        {block}
        <figcaption>Wounding switches on skeleton building at 28 °C but barely at 31 °C. This would explain the missing radial corallites.</figcaption>
      </figure>
      <figure>
        <h3>Heat is general stress</h3>
        {stress}
        <figcaption>Wounding raises these genes by the same amount at both temperatures, and heat lowers everything. Regrowth would stall from general stress or low energy rather than a specific block on skeleton building.</figcaption>
      </figure>
    </div>
    <p class="why">Heat also cut growth in unwounded fragments (10.2% vs 16.0%), so the test is whether the lines diverge, not whether 31 °C sits lower.</p>
  </section>
</div>

<div class="wide">
  <section class="tl" aria-labelledby="tl-h">
    <h2 id="tl-h">Our samples against earlier timelines</h2>
    <p class="intro">Our day 10 sample falls before radial corallites form and inside the second wave of wound genes seen in <i>A. aspera</i>. Han's skeleton genes peaked later, in a slower species.</p>
    <div class="scroll" tabindex="0" aria-label="Timeline, scrolls sideways on small screens">{timeline}</div>
  </section>
</div>

<div class="wide">
  <section aria-labelledby="sets-h">
    <h2 id="sets-h">What each gene set should do</h2>
    <p class="intro">We compiled these sets in September 2026, before expression counts arrived, and corrected them against the original papers on {updated}.</p>
    <table>
      <thead><tr><th scope="col">Gene set</th><th scope="col">Day 1</th><th scope="col">Day 10</th></tr></thead>
      <tbody>
{sets}
      </tbody>
    </table>
  </section>
</div>

<div class="col">
  <section aria-labelledby="cau-h">
    <h2 id="cau-h">Before reading too much into a result</h2>
    <ul class="cau">
{cautions}
    </ul>
  </section>

  <section aria-labelledby="src-h">
    <h2 id="src-h">Sources</h2>
    <p class="intro" style="margin:0">Each checked against its full text on {updated}.</p>
    <ul class="plain">
{studies}
    </ul>
  </section>
</div>

<div class="wide">
  <details>
    <summary>Full candidate list<span>{n} gene families, {nchk} checked against the papers</span></summary>
    <label for="q" class="meta" style="display:block;margin-top:1rem">Filter by gene, study or keyword</label>
    <input id="q" type="search" placeholder="galaxin, Han, heat, D10">
    <table id="genes">
      <thead><tr><th scope="col" style="width:16%">Gene or family</th><th scope="col" style="width:30%">Evidence</th><th scope="col" style="width:20%">Use here</th><th scope="col" style="width:18%">Caveat</th><th scope="col">Sources</th></tr></thead>
      <tbody>
{genes}
      </tbody>
    </table>
    <p class="meta">Source table: <code>docs/rnaseq/candidate_genes_reference.csv</code>. Analysis plan: <code>docs/rnaseq/expression_integration_analysis_plan.md</code>.</p>
  </details>
</div>

<footer class="col">Stier Lab, UC Santa Barbara, with the Bay Lab, UC Davis.</footer>
<script>
  (function() {{
    var q = document.getElementById('q'), rows = document.querySelectorAll('#genes tbody tr');
    q.addEventListener('input', function() {{
      var v = q.value.trim().toLowerCase();
      rows.forEach(function(r) {{ r.style.display = !v || r.getAttribute('data-text').indexOf(v) > -1 ? '' : 'none'; }});
    }});
  }})();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
