#!/usr/bin/env python3
"""make_nosplit_html.py -- build a self-contained, browsable HTML page for the
no-splitting-field cubic-surface database (the large-|Gal| classes reached by
CubicSurfaceNoSplittingField).

Reads:
  database_seed_nosplit.txt   label:source_coeffs:orbit_sizes:frob_orders:consistent:cubic
  class_info_nosplit.txt      label:d:t:order:groupname   (from class_info_nosplit.m)
  WE6subgroups.txt            label:gens:d:t              (S27 generators)
  lmfdb_fields.csv            label,"{coeffs}"            (LMFDB number-field labels)

Writes:
  seed_database_nosplit.html  a single self-contained file (surfaces embedded as
                              JSON; equations for the SOURCE polynomial rendered
                              with KaTeX; the surfaces themselves are huge and
                              non-minimized, so they are shown as a monospace
                              preview + copy-to-Magma, not typeset).

Usage:  python3 make_nosplit_html.py
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

def max_coeff_digits(cubic):
    nums = re.findall(r"\d+", re.sub(r"\^\d+", "", cubic))
    return max((len(n) for n in nums), default=1)

def name_to_tex(name):
    # C2^4.S5, C3^3.S4.C2, S6, ... -> LaTeX with subscripts/superscripts
    s = re.sub(r"([A-Za-z])(\d+)", r"\1_{\2}", name)
    s = s.replace("^", "^")
    s = s.replace("*", r" \times ").replace(":", r" \rtimes ").replace(".", r".\,")
    return s

def norm(coeffs_str):
    return coeffs_str.strip().strip("[]{}").replace(" ", "")

def clean(s):
    return s.strip().replace("  ", " ")

# ---- class metadata --------------------------------------------------------
meta = {}
with open("class_info_nosplit.txt") as f:
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
for fn in ("lmfdb_fields.csv", "lmfdb_fields_nosplit.csv"):
    try:
        with open(fn) as f:
            for row in csv.DictReader(f):
                nflabel[norm(row["c"])] = row["label"]
    except FileNotFoundError:
        pass

# ---- read the realizations -------------------------------------------------
data = []
n_field = 0
with open("database_seed_nosplit.txt") as f:
    for ln in f:
        ln = ln.rstrip("\n")
        if not ln or ln.startswith("#"):
            continue
        label, coeffs_s, orbits, frob, consistent, cubic = ln.split(":", 5)
        coeffs = json.loads(coeffs_s)
        m = meta.get(label, dict(d=0, t=0, order=WE6 // int(label.split(".")[2]), name="?"))
        nf = nflabel.get(norm(coeffs_s), "")
        if nf:
            n_field += 1
        data.append({
            "label": label,
            "name": m["name"],
            "nametex": name_to_tex(m["name"]),
            "order": m["order"],
            "dt": "%dT%d" % (m["d"], m["t"]),
            "gens": gens.get(label, ""),
            "polytex": poly_to_tex(coeffs),
            "orbits": clean(orbits),
            "frob": clean(frob),
            "digits": max_coeff_digits(cubic),
            "nf": nf,
            "raw": cubic,
        })

data.sort(key=lambda c: (c["order"], c["dt"], c["label"]))
by_order = collections.Counter(c["order"] for c in data)
summary = {"n_surfaces": len(data), "n_classes": len(data), "n_total": 350,
           "n_seeded": 58, "by_order": sorted(by_order.items()),
           "min_digits": min(c["digits"] for c in data),
           "max_digits": max(c["digits"] for c in data)}

PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cubic Surface Database &middot; no-splitting-field classes</title>
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
  .cbody{display:none;border-top:1px solid var(--line);padding:12px 16px 16px}
  .card.open .cbody{display:block}
  .row{display:flex;gap:8px;flex-wrap:wrap;align-items:baseline;margin:4px 0;color:var(--muted);font-size:13px}
  .row .k{min-width:120px;color:var(--muted)}
  .row .v{color:var(--ink)}
  .mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
  .surfbox{margin-top:10px;border-top:1px dashed var(--line);padding-top:10px}
  .surfhead{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:6px}
  .surfhead .mc{color:var(--muted);font-size:11.5px}
  pre.surf{background:var(--panel2);border:1px solid var(--line);border-radius:8px;
    padding:10px;margin:0;overflow-x:auto;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
    font-size:11.5px;white-space:pre-wrap;word-break:break-all;max-height:220px;overflow-y:auto;color:var(--ink)}
  .copy{background:var(--panel2);border:1px solid var(--line);color:var(--muted);border-radius:8px;
    padding:3px 8px;cursor:pointer;font-size:11.5px;white-space:nowrap}
  .copy:hover{border-color:var(--accent);color:var(--accent)}
  .copy.ok{color:var(--accent2);border-color:var(--accent2)}
  footer{border-top:1px solid var(--line);color:var(--muted);font-size:12.5px;padding:18px 0 40px}
  .hidden{display:none}
  .note{background:var(--chip);border:1px solid var(--line);border-radius:10px;padding:10px 12px;
    color:var(--muted);font-size:12.5px;margin:10px 0 0}
</style>
</head>
<body>
<header><div class="wrap">
  <h1>Cubic Surface Database &middot; no-splitting-field classes</h1>
  <div class="sub">Smooth cubic surfaces over <span class="mono">&#8474;</span> with prescribed Galois action on
    their 27 lines, for the large-<span class="mono">|Gal(f)|</span> classes reached <b>without</b> building the
    splitting field (<span class="mono">CubicSurfaceNoSplittingField</span>) &middot;
    <a href="https://github.com/AndrewVSutherland2/CubicSurfaces">repository</a></div>
  <div class="stats" id="stats"></div>
  <div class="ordbar" id="ordbar"></div>
  <div class="note">Each surface is an <b>existence certificate</b>: an explicit smooth cubic whose 27 lines carry the
    prescribed action, but with very large, <b>non-minimized</b> coefficients (see max-coeff digits). They are shown as a
    monospace preview with copy-to-Magma, not typeset. Certification is by mod-<span class="mono">q</span> Frobenius orders
    on the 27 lines (a necessary condition).</div>
</div></header>
<div class="controls"><div class="wrap">
  <input type="search" id="q" placeholder="Search by group (S6, C2^4.A5&hellip;), W(E6) label, or dTt&hellip;" autocomplete="off">
  <button class="btn" id="expand">Expand all</button>
  <button class="btn" id="collapse">Collapse all</button>
  <button class="btn" id="theme">&#9728; Light</button>
</div></div>
<main><div class="wrap">
  <div class="count" id="count"></div>
  <div id="classes"></div>
</div></main>
<footer><div class="wrap">
  Each class is one conjugacy class of subgroups of <span class="mono">W(E<sub>6</sub>)</span> (the Galois action on the
  27 lines). The surface comes from the number field <span class="mono">f</span>; regenerate it with
  <span class="mono">CubicSurfaceNoSplittingField(G, f)</span>. Number fields link to the
  <a href="https://www.lmfdb.org/NumberField/">LMFDB</a>; the copy buttons yield Magma input (the cubic in
  <span class="mono">P^3</span>, and the class as a subgroup of <span class="mono">S27</span>). Classes the method cannot
  yet reach are listed in <span class="mono">nosplit_unrealized.txt</span>. Generated by
  <span class="mono">make_nosplit_html.py</span>.
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
    `<div class="stat"><b>${SUMMARY.n_classes}</b><span>large-|Gal| classes</span></div>`+
    `<div class="stat"><b>${SUMMARY.n_classes + SUMMARY.n_seeded} / ${SUMMARY.n_total}</b><span>W(E6) classes with a surface</span></div>`+
    `<div class="stat"><b>${SUMMARY.by_order.length}</b><span>distinct group orders</span></div>`+
    `<div class="stat"><b>${SUMMARY.min_digits}&ndash;${SUMMARY.max_digits}</b><span>max-coeff digits</span></div>`;
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
    `<span class="pill">orbits ${c.orbits}</span><br>`+
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
      // source polynomial + LMFDB
      const r1 = document.createElement('div'); r1.className='row';
      const k1 = document.createElement('span'); k1.className='k'; k1.textContent='number field';
      const v1 = document.createElement('span'); v1.className='v'; const fspan=document.createElement('span');
      v1.appendChild(fspan); tex(fspan, 'f(x) = ' + c.polytex, false);
      if(c.nf){ const a=document.createElement('a'); a.href='https://www.lmfdb.org/NumberField/'+c.nf;
        a.target='_blank'; a.rel='noopener'; a.textContent=' LMFDB '+c.nf;
        v1.appendChild(document.createTextNode(' ·')); v1.appendChild(a); }
      r1.appendChild(k1); r1.appendChild(v1); body.appendChild(r1);
      // orbits / frobenius / digits
      const meta2 = [['27-line orbits', c.orbits],
                     ['Frobenius orders', c.frob + '  (all lie in the target group)'],
                     ['max-coeff digits', ''+c.digits]];
      meta2.forEach(([k,v])=>{ const r=document.createElement('div'); r.className='row';
        const ks=document.createElement('span'); ks.className='k'; ks.textContent=k;
        const vs=document.createElement('span'); vs.className='v mono'; vs.textContent=v;
        r.appendChild(ks); r.appendChild(vs); body.appendChild(r); });
      // surface: preview + copy
      const sb = document.createElement('div'); sb.className='surfbox';
      const sh = document.createElement('div'); sh.className='surfhead';
      const lab = document.createElement('span'); lab.textContent='cubic surface';
      const mc = document.createElement('span'); mc.className='mc';
      mc.textContent = c.raw.length.toLocaleString()+' chars';
      const cp = document.createElement('button'); cp.className='copy'; cp.textContent='copy Magma';
      cp.title='Copy the full cubic as Magma input';
      cp.onclick=()=>copyBtn(cubicMagma(c.raw), cp);
      sh.appendChild(lab); sh.appendChild(mc); sh.appendChild(cp);
      const pre = document.createElement('pre'); pre.className='surf';
      pre.textContent = c.raw.length>1600 ? c.raw.slice(0,1600)+'…' : c.raw;
      sb.appendChild(sh); sb.appendChild(pre); body.appendChild(sb);
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
with open("seed_database_nosplit.html", "w") as f:
    f.write(page)
print("wrote seed_database_nosplit.html: %d classes, %d fields linked, digits %d-%d" %
      (len(data), n_field, summary["min_digits"], summary["max_digits"]))
