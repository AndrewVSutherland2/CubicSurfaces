#!/usr/bin/env python3
"""make_seed_html.py -- build a self-contained, browsable HTML page for the
cubic-surface seed database.

Reads:
  database_seed.txt   lines  label:source_coeffs:cubic
  class_info.txt      lines  label:d:t:order:groupname     (from class_info.m)
  WE6subgroups.txt    lines  label:gens:d:t                (S27 generators)
  lmfdb_fields.csv    rows   label,"{coeffs}"              (LMFDB number-field labels)

Writes:
  seed_database.html  a single self-contained file (data embedded as JSON,
                      equations rendered with KaTeX from a CDN).

Usage:  python3 make_seed_html.py
"""
import json, re, csv, collections

WE6 = 51840

def poly_to_tex(coeffs):
    terms = []
    for i in range(len(coeffs) - 1, -1, -1):
        c = coeffs[i]
        if c == 0:
            continue
        if i == 0:
            body = str(abs(c))
        else:
            xpart = "x" if i == 1 else "x^{%d}" % i
            body = xpart if abs(c) == 1 else "%d%s" % (abs(c), xpart)
        terms.append((c < 0, body))
    if not terms:
        return "0"
    out = ("-" if terms[0][0] else "") + terms[0][1]
    for neg, body in terms[1:]:
        out += (" - " if neg else " + ") + body
    return out

def cubic_to_tex(s):
    return s.replace("*", "")            # implied multiplication

def max_coeff(cubic):
    nums = re.findall(r"\d+", re.sub(r"\^\d+", "", cubic))
    return max([int(n) for n in nums] + [1])

def name_to_tex(name):
    s = re.sub(r"([A-Za-z\)])(\d+)", r"\1_{\2}", name)
    return s.replace("*", r" \times ").replace(":", r" \rtimes ")

def norm(coeffs_str):
    return coeffs_str.strip().strip("[]{}").replace(" ", "")

# ---- class metadata --------------------------------------------------------
meta = {}
with open("class_info.txt") as f:
    for ln in f:
        ln = ln.strip()
        if ln:
            label, d, t, order, gname = ln.split(":", 4)
            meta[label] = dict(d=int(d), t=int(t), order=int(order), name=gname)

# ---- S27 generators per class ----------------------------------------------
gens = {}
with open("WE6subgroups.txt") as f:
    for ln in f:
        parts = ln.rstrip("\n").split(":")
        if len(parts) >= 2:
            gens[parts[0]] = parts[1]

# ---- LMFDB number-field labels keyed by normalised coefficient list --------
nflabel = {}
with open("lmfdb_fields.csv") as f:
    for row in csv.DictReader(f):
        nflabel[norm(row["c"])] = row["label"]

# ---- seed: group cubics by class -------------------------------------------
classes = collections.OrderedDict()
n_cubics = 0
with open("database_seed.txt") as f:
    for ln in f:
        ln = ln.rstrip("\n")
        if not ln or ln.startswith("#"):
            continue
        label, coeffs_s, cubic = ln.split(":", 2)
        coeffs = json.loads(coeffs_s)
        n_cubics += 1
        classes.setdefault(label, []).append({
            "polytex": poly_to_tex(coeffs),
            "cubictex": cubic_to_tex(cubic),
            "raw": cubic,
            "maxc": max_coeff(cubic),
            "nf": nflabel.get(norm(coeffs_s), ""),
        })

# ---- assemble data ---------------------------------------------------------
data = []
for label, cubs in classes.items():
    m = meta.get(label, dict(d=0, t=0, order=WE6 // int(label.split(".")[2]), name="?"))
    cubs.sort(key=lambda c: c["maxc"])
    data.append({
        "label": label,
        "name": m["name"],
        "nametex": name_to_tex(m["name"]),
        "order": m["order"],
        "dt": "%dT%d" % (m["d"], m["t"]),
        "gens": gens.get(label, ""),
        "cubics": cubs,
    })
data.sort(key=lambda c: (c["order"], c["dt"], c["label"]))

by_order = collections.Counter(c["order"] for c in data)
summary = {"n_cubics": n_cubics, "n_classes": len(data), "n_total": 350,
           "by_order": sorted(by_order.items())}

PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cubic Surface Seed Database</title>
<script>try{var t=localStorage.getItem('cstheme');if(t==='dark')document.documentElement.dataset.theme='dark';}catch(e){}</script>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css" crossorigin="anonymous">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js" crossorigin="anonymous"></script>
<style>
  :root{
    --bg:#f6f7f9; --panel:#ffffff; --panel2:#f1f3f6; --ink:#1b1f27; --muted:#5b6573;
    --line:#e4e8ee; --accent:#2563eb; --accent2:#15803d; --chip:#eef1f6;
    --hgrad:linear-gradient(180deg,#eef2f8,#f6f7f9); --shadow:0 1px 2px rgba(20,30,50,.06);
  }
  html[data-theme="dark"]{
    --bg:#0f1115; --panel:#171a21; --panel2:#1d212b; --ink:#e6e9ef; --muted:#9aa4b2;
    --line:#2a2f3a; --accent:#6ea8fe; --accent2:#7ee787; --chip:#222735;
    --hgrad:linear-gradient(180deg,#141822,#0f1115); --shadow:none;
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--ink);
    font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
  a{color:var(--accent);text-decoration:none}
  a:hover{text-decoration:underline}
  header{padding:28px 20px 18px;border-bottom:1px solid var(--line);background:var(--hgrad)}
  .wrap{max-width:1000px;margin:0 auto;padding:0 16px}
  h1{margin:0 0 4px;font-size:24px;letter-spacing:.2px}
  .sub{color:var(--muted);font-size:14px;margin-bottom:14px}
  .stats{display:flex;flex-wrap:wrap;gap:8px;margin-top:10px}
  .stat{background:var(--chip);border:1px solid var(--line);border-radius:10px;padding:8px 12px}
  .stat b{font-size:18px;color:var(--accent2)}
  .stat span{color:var(--muted);font-size:12px;display:block}
  .ordbar{display:flex;flex-wrap:wrap;gap:6px;margin-top:12px}
  .ord{background:var(--chip);border:1px solid var(--line);border-radius:999px;
    padding:4px 10px;font-size:12.5px;color:var(--muted);cursor:pointer;user-select:none}
  .ord:hover{border-color:var(--accent)}
  .ord.active{background:var(--accent);color:#fff;border-color:var(--accent);font-weight:600}
  .controls{position:sticky;top:0;z-index:5;background:var(--bg);border-bottom:1px solid var(--line);padding:10px 0}
  .controls .wrap{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
  input[type=search]{flex:1;min-width:200px;background:var(--panel);border:1px solid var(--line);
    color:var(--ink);border-radius:10px;padding:9px 12px;font-size:14px}
  .btn{background:var(--panel2);border:1px solid var(--line);color:var(--ink);
    border-radius:10px;padding:8px 12px;cursor:pointer;font-size:13px}
  .btn:hover{border-color:var(--accent)}
  main{padding:18px 0 60px}
  .count{color:var(--muted);font-size:13px;margin:6px 2px 14px}
  .card{background:var(--panel);border:1px solid var(--line);border-radius:14px;margin:0 0 12px;overflow:hidden;box-shadow:var(--shadow)}
  .chead{display:flex;align-items:center;gap:12px;padding:14px 16px;cursor:pointer}
  .chead:hover{background:var(--panel2)}
  .gname{font-size:18px;min-width:60px}
  .ghead-main{flex:1;min-width:0}
  .glabel{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;color:var(--muted)}
  .gmeta{color:var(--muted);font-size:12.5px}
  .pill{background:var(--chip);border:1px solid var(--line);border-radius:999px;
    padding:2px 9px;font-size:12px;color:var(--muted);margin-left:6px;white-space:nowrap}
  .caret{color:var(--muted);transition:transform .15s ease}
  .card.open .caret{transform:rotate(90deg)}
  .cbody{display:none;border-top:1px solid var(--line);padding:6px 16px 14px}
  .card.open .cbody{display:block}
  .surf{padding:12px 0;border-bottom:1px dashed var(--line)}
  .surf:last-child{border-bottom:0}
  .surf .from{color:var(--muted);font-size:12.5px;margin-bottom:6px;display:flex;align-items:baseline;gap:8px;flex-wrap:wrap}
  .eqrow{display:flex;align-items:center;gap:10px}
  .eq{overflow-x:auto;padding:2px 0;flex:1;min-width:0}
  .mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
  .mc{color:var(--muted);font-size:11.5px;white-space:nowrap}
  .copy{background:var(--panel2);border:1px solid var(--line);color:var(--muted);border-radius:8px;
    padding:3px 8px;cursor:pointer;font-size:11.5px;white-space:nowrap}
  .copy:hover{border-color:var(--accent);color:var(--accent)}
  .copy.ok{color:var(--accent2);border-color:var(--accent2)}
  footer{border-top:1px solid var(--line);color:var(--muted);font-size:12.5px;padding:18px 0 40px}
  .hidden{display:none}
</style>
</head>
<body>
<header><div class="wrap">
  <h1>Cubic Surface Seed Database</h1>
  <div class="sub">Smooth cubic surfaces over <span class="mono">&#8474;</span> with prescribed Galois action on their 27 lines,
    built by the Elsenhans&ndash;Jahnel Algorithm&nbsp;5.1 &middot;
    <a href="https://github.com/AndrewVSutherland2/CubicSurfaces">repository</a></div>
  <div class="stats" id="stats"></div>
  <div class="ordbar" id="ordbar"></div>
</div></header>
<div class="controls"><div class="wrap">
  <input type="search" id="q" placeholder="Search by group (S3, C2^2&hellip;), W(E6) label, or dTt&hellip;" autocomplete="off">
  <button class="btn" id="expand">Expand all</button>
  <button class="btn" id="collapse">Collapse all</button>
  <button class="btn" id="theme">&#9728; Light</button>
</div></div>
<main><div class="wrap">
  <div class="count" id="count"></div>
  <div id="classes"></div>
</div></main>
<footer><div class="wrap">
  Each class is one conjugacy class of subgroups of <span class="mono">W(E<sub>6</sub>)</span>
  (the Galois action on the 27 lines); cubics within a class come from different number fields and are
  sorted by largest coefficient. The number fields link to the
  <a href="https://www.lmfdb.org/NumberField/">LMFDB</a>; the copy buttons yield Magma input
  (the cubic in <span class="mono">P^3</span>, and the class as a subgroup of <span class="mono">S27</span>).
  This seed covers the classes of small intrinsic discriminant &Delta;<sub>Cl</sub>; generated by
  <span class="mono">make_seed_html.py</span>.
</div></footer>

<script>
const DATA = __DATA__;
const SUMMARY = __SUMMARY__;
let activeOrder = null, cards = [];

function tex(node, s, display){
  try{ katex.render(s, node, {throwOnError:false, displayMode:!!display}); }
  catch(e){ node.textContent = s; }
}
function cubicMagma(raw){ return 'P<x,y,z,w> := PolynomialRing(RationalField(), 4);\nf := '+raw+';'; }
function subgroupMagma(gens){ return 'S27 := SymmetricGroup(27);\nG := sub<S27 | '+gens+'>;'; }
function copyBtn(text, btn){
  const done=()=>{ const o=btn.textContent; btn.textContent='copied!'; btn.classList.add('ok');
    setTimeout(()=>{btn.textContent=o; btn.classList.remove('ok');}, 1200); };
  if(navigator.clipboard && navigator.clipboard.writeText){
    navigator.clipboard.writeText(text).then(done).catch(()=>fallback(text,done));
  } else fallback(text, done);
}
function fallback(text, done){
  const ta=document.createElement('textarea'); ta.value=text; ta.style.position='fixed'; ta.style.opacity='0';
  document.body.appendChild(ta); ta.select();
  try{ document.execCommand('copy'); done(); }catch(e){} document.body.removeChild(ta);
}

function setTheme(t){
  document.documentElement.dataset.theme = (t==='dark') ? 'dark' : '';
  try{ localStorage.setItem('cstheme', t==='dark'?'dark':'light'); }catch(e){}
  const b=document.getElementById('theme');
  if(b) b.innerHTML = (t==='dark') ? '&#9728; Light' : '&#9790; Dark';
}

function renderStats(){
  document.getElementById('stats').innerHTML =
    `<div class="stat"><b>${SUMMARY.n_cubics}</b><span>cubic surfaces</span></div>`+
    `<div class="stat"><b>${SUMMARY.n_classes} / ${SUMMARY.n_total}</b><span>W(E6) classes</span></div>`+
    `<div class="stat"><b>${SUMMARY.by_order.length}</b><span>distinct group orders</span></div>`;
  const ob = document.getElementById('ordbar');
  ob.innerHTML = '<span class="ord'+(activeOrder===null?' active':'')+'" data-o="">all orders</span>'+
    SUMMARY.by_order.map(([o,n])=>`<span class="ord" data-o="${o}">order ${o} &middot; ${n}</span>`).join('');
  ob.querySelectorAll('.ord').forEach(el=>el.onclick=()=>{
    const v=el.getAttribute('data-o'); activeOrder = v===''?null:+v; apply();
  });
}

function makeCard(c){
  const card = document.createElement('div');
  card.className = 'card';
  card.dataset.hay = (c.name+' '+c.label+' '+c.dt+' order '+c.order).toLowerCase();
  const head = document.createElement('div');
  head.className = 'chead';
  head.innerHTML =
    `<span class="caret">&#9656;</span>`+
    `<span class="gname"></span>`+
    `<span class="ghead-main"><span class="gmeta">${c.dt} &middot; order ${c.order}</span>`+
    `<span class="pill">${c.cubics.length} cubic${c.cubics.length>1?'s':''}</span><br>`+
    `<span class="glabel">${c.label}</span></span>`;
  tex(head.querySelector('.gname'), c.nametex, false);
  if(c.gens){
    const sb = document.createElement('button');
    sb.className='copy'; sb.textContent='copy subgroup';
    sb.title='Copy this class as a subgroup of S27 (Magma)';
    sb.onclick=(e)=>{ e.stopPropagation(); copyBtn(subgroupMagma(c.gens), sb); };
    head.appendChild(sb);
  }
  const body = document.createElement('div');
  body.className = 'cbody';
  let built = false;
  head.onclick = ()=>{
    card.classList.toggle('open');
    if(card.classList.contains('open') && !built){
      built = true;
      c.cubics.forEach(s=>{
        const d = document.createElement('div'); d.className='surf';
        const from = document.createElement('div'); from.className='from';
        const fspan = document.createElement('span');
        from.appendChild(document.createTextNode('from '));
        from.appendChild(fspan);
        tex(fspan, 'f(x) = ' + s.polytex, false);
        if(s.nf){
          const a=document.createElement('a'); a.href='https://www.lmfdb.org/NumberField/'+s.nf;
          a.target='_blank'; a.rel='noopener'; a.textContent='LMFDB '+s.nf;
          from.appendChild(document.createTextNode(' · ')); from.appendChild(a);
        }
        const eqrow = document.createElement('div'); eqrow.className='eqrow';
        const eq = document.createElement('div'); eq.className='eq';
        tex(eq, s.cubictex+' \\,=\\, 0', true);
        const mc = document.createElement('span'); mc.className='mc'; mc.textContent='max coeff '+s.maxc;
        const cp = document.createElement('button'); cp.className='copy'; cp.textContent='copy Magma';
        cp.title='Copy the cubic as Magma input';
        cp.onclick=()=>copyBtn(cubicMagma(s.raw), cp);
        eqrow.appendChild(eq); eqrow.appendChild(mc); eqrow.appendChild(cp);
        d.appendChild(from); d.appendChild(eqrow); body.appendChild(d);
      });
    }
  };
  card.appendChild(head); card.appendChild(body);
  return card;
}

function apply(){
  const q = document.getElementById('q').value.trim().toLowerCase();
  let shown = 0;
  cards.forEach(({c,el})=>{
    const vis = (activeOrder===null || c.order===activeOrder) && (!q || el.dataset.hay.includes(q));
    el.classList.toggle('hidden', !vis); if(vis) shown++;
  });
  document.getElementById('count').textContent = `${shown} class${shown===1?'':'es'} shown`;
  document.querySelectorAll('#ordbar .ord').forEach(el=>{
    const v=el.getAttribute('data-o');
    el.classList.toggle('active', (v===''&&activeOrder===null) || (+v===activeOrder));
  });
}

function init(){
  const wrap = document.getElementById('classes');
  cards = DATA.map(c=>{ const el=makeCard(c); wrap.appendChild(el); return {c,el}; });
  document.getElementById('q').addEventListener('input', apply);
  document.getElementById('expand').onclick = ()=>cards.forEach(({el})=>{
    if(!el.classList.contains('hidden') && !el.classList.contains('open')) el.querySelector('.chead').click();
  });
  document.getElementById('collapse').onclick = ()=>cards.forEach(({el})=>el.classList.remove('open'));
  let saved='light'; try{ saved = localStorage.getItem('cstheme') || 'light'; }catch(e){}
  setTheme(saved);
  document.getElementById('theme').onclick = ()=>
    setTheme(document.documentElement.dataset.theme==='dark' ? 'light' : 'dark');
  renderStats(); apply();
}
if(document.readyState==='loading') document.addEventListener('DOMContentLoaded', init);
else init();
</script>
</body>
</html>
"""

page = PAGE.replace("__DATA__", json.dumps(data, ensure_ascii=False)) \
           .replace("__SUMMARY__", json.dumps(summary))
with open("seed_database.html", "w") as f:
    f.write(page)
print("wrote seed_database.html: %d classes, %d cubics, %d fields linked" %
      (len(data), n_cubics, sum(1 for c in data for x in c["cubics"] if x["nf"])))
