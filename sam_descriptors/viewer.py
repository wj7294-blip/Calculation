"""viewer.html 생성 (3Dmol.js 를 파일 안에 포함 → 파일 하나로 열림).

  collect()   : results/molecules/*.json + summary/database.csv 에서 뷰어용 데이터 수집
  cube_text() : cubes/*.cube 를 압축된 문자열로 변환 (HOMO/LUMO 등고면용)
  build()     : viewer.html 작성
"""
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S"}
PROP_KEYS = [("HOMO", "HOMO (eV)"), ("LUMO", "LUMO (eV)"), ("gap", "gap (eV)"), ("IP", "vertical IP (eV)"),
             ("EA", "vertical EA (eV)"), ("dipole", "|μ| (D)"), ("dipole_axial", "μ·(앵커→말단) (D)"),
             ("q_anchor_atom", "q 앵커 원자"), ("q_anchor_group", "q 앵커기"), ("q_terminal_atom", "q 말단 원자"),
             ("q_terminal_group", "q 말단기"), ("zw_dist_boltz", "N⁺–X 거리, 볼츠만 평균 (Å)"),
             ("MW", "분자량"), ("TPSA", "TPSA (Å²)"), ("logP", "logP"), ("nF", "F 원자 수"),
             ("nRotB", "회전 가능 결합"), ("chain_bonds", "앵커–말단 결합 수"), ("E_min_Eh", "E_min (Eh)"),
             ("solvent", "용매(ALPB)"), ("n_conf_ok", "conformer 수"),
             ("n_conf_connectivity_changed", "결합 변화 conformer")]


def cube_text(path):
    lines = open(path).read().split("\n")
    natoms = abs(int(lines[2].split()[0]))
    head = lines[:6 + natoms]
    vals = np.array(" ".join(lines[6 + natoms:]).split(), dtype=float)
    out = []
    for i in range(0, len(vals), 6):
        chunk = []
        for v in vals[i:i + 6]:
            if abs(v) < 5e-4:
                chunk.append("0")
            else:
                s = f"{v:.3f}".replace("0.", ".").replace("-0.", "-.")
                chunk.append(s)
        out.append(" ".join(chunk))
    return "\n".join(head + out) + "\n"


def collect():
    os.chdir(HERE)
    frames = [pd.read_csv(p) for p in ("results/summary.csv", "results/database.csv") if os.path.exists(p)]
    df = pd.concat(frames, ignore_index=True).drop_duplicates("name", keep="last")
    mols = []
    for _, r in df.iterrows():
        p = f"results/molecules/{r['name']}.json"
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        sym = [SYM.get(z, "X") for z in d["Z"]]
        confs = []
        for c in d["conformers"]:
            confs.append(dict(e=round(c["rel_kcal"], 2), w=round(c["boltz_w"], 3),
                              xyz=np.round(np.array(c["positions"]), 3).tolist(),
                              q=np.round(np.array(c["charges"]), 3).tolist(),
                              mu=np.round(np.array(c["dipole_vec"]), 2).tolist(),
                              bad=bool(c["connectivity_changed"]), changes=c["connectivity_changed"],
                              zw=round(c["zw_dist"], 2) if "zw_dist" in c else None))
        props = {}
        for k, _ in PROP_KEYS:
            v = r.get(k)
            if isinstance(v, float) and np.isnan(v):
                v = None
            elif isinstance(v, (np.integer,)):
                v = int(v)
            elif isinstance(v, (np.floating,)):
                v = round(float(v), 4)
            props[k] = v
        cubes = {}
        for lab in ("homo", "lumo"):
            cp = f"cubes/{r['name']}_{lab}.cube"
            if os.path.exists(cp):
                cubes[lab] = cube_text(cp)
        mols.append(dict(name=r["name"], series=r["series"], smiles=r["smiles"], anchor=r["anchor"],
                         terminal=r["terminal"], sym=sym, a=d["anchor_idx"], t=d["terminal_idx"],
                         npos=d["nplus_idx"], props=props, confs=confs, cubes=cubes))
    return mols


TEMPLATE = r"""<!DOCTYPE html>
<html lang="ko"><head><meta charset="utf-8">
<title>SAM descriptor viewer</title>
<style>
:root{--bg:#f7f7f5;--panel:#fff;--ink:#1a1a19;--muted:#666;--line:#e3e2dd;--accent:#2a78d6}
*{box-sizing:border-box}body{margin:0;font:13px/1.45 -apple-system,"Segoe UI",Roboto,"Noto Sans KR",sans-serif;background:var(--bg);color:var(--ink)}
.app{display:grid;grid-template-columns:250px 1fr 320px;height:100vh}
.left,.right{background:var(--panel);border-right:1px solid var(--line);overflow:auto;padding:10px}
.right{border-right:0;border-left:1px solid var(--line)}
h1{font-size:15px;margin:2px 0 8px}
#search{width:100%;padding:6px 8px;border:1px solid var(--line);border-radius:6px;margin-bottom:8px}
.series{font-size:11px;color:var(--muted);margin:10px 0 3px;text-transform:uppercase;letter-spacing:.04em}
.item{padding:4px 6px;border-radius:5px;cursor:pointer;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.item:hover{background:#eef3fb}.item.sel{background:var(--accent);color:#fff}
.center{display:flex;flex-direction:column;min-width:0}
#view{flex:1;position:relative;min-height:300px}
.bar{display:flex;flex-wrap:wrap;}.bar.top{flex-wrap:nowrap}.bar{gap:8px 14px;align-items:center;padding:8px 12px;border-bottom:1px solid var(--line);background:var(--panel)}
.bar label{display:flex;gap:4px;align-items:center;cursor:pointer}
button{border:1px solid var(--line);background:#fff;border-radius:5px;padding:3px 9px;cursor:pointer}
button:disabled{opacity:.4;cursor:default}
.title{font-weight:600;font-size:14px}
.sub{color:var(--muted);font-size:11px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;min-width:0}
table{border-collapse:collapse;width:100%}td{padding:3px 4px;border-bottom:1px solid var(--line);vertical-align:top}
td:first-child{color:var(--muted);width:58%}
.legend{font-size:11px;color:var(--muted);margin-top:8px}
.conf{font-variant-numeric:tabular-nums}
.warn{color:#b3261e;font-size:11px}
</style>
<script>__3DMOL__</script>
</head><body>
<div class="app">
 <div class="left">
  <h1>SAM 후보 분자</h1>
  <input id="search" placeholder="이름·SMILES 검색…">
  <div id="list"></div>
 </div>
 <div class="center">
  <div class="bar top">
   <span class="title" id="mname">–</span>
   <span class="sub" id="msmiles"></span>
   <span style="flex:1"></span>
   <button id="prev">◀</button><span class="conf" id="confinfo"></span><button id="next">▶</button>
  </div>
  <div class="bar">
   <label><input type="checkbox" id="cbDip" checked> 쌍극자 화살표</label>
   <label><input type="checkbox" id="cbQ"> 부분전하 색칠</label>
   <label><input type="checkbox" id="cbHomo"> HOMO</label>
   <label><input type="checkbox" id="cbLumo"> LUMO</label>
   <label>등고값 <input type="range" id="iso" min="0.005" max="0.06" step="0.005" value="0.02" style="width:90px"><span id="isov">0.02</span></label>
   <label><input type="checkbox" id="cbLabel"> 원자 번호</label>
   <span class="sub" id="cubenote"></span>
  </div>
  <div id="view"></div>
 </div>
 <div class="right">
  <h1>성질 (최저 에너지 구조)</h1>
  <table id="props"></table>
  <div class="legend">
   쌍극자 화살표: 음(−)→양(+) 방향, 길이 0.3 Å/D. 부분전하: GFN2-xTB Mulliken 형 전하, 파랑(−0.5) – 흰색(0) – 빨강(+0.5).
   HOMO/LUMO 등고면: 파랑/빨강 = 위상 ±, 최저 에너지 conformer 에서만 표시(대표 분자 10개).
   N⁺–X 점선: 쯔비터이온의 N⁺ 와 음이온 중심(S 또는 C) 거리.
  </div>
 </div>
</div>
<script>
const DATA = __DATA__;
const PROP_KEYS = __PROPKEYS__;
const SERIES_LABEL = {A:"시리즈 A (앵커×말단)",B:"시리즈 B (불소화)",C:"시리즈 C (쯔비터이온)",ref:"기준 분자",val:"검증(기체상)",user:"사용자 추가"};
let viewer = $3Dmol.createViewer(document.getElementById('view'), {backgroundColor:'#ffffff'});
let cur = null, ci = 0, model = null, shapes = [], labels = [], surf = [];

function buildList(filter){
  const box = document.getElementById('list'); box.innerHTML='';
  const f = (filter||'').toLowerCase();
  for (const s of ['A','B','C','ref','val','user']){
    const items = DATA.filter(m=>m.series===s && (!f || m.name.toLowerCase().includes(f) || m.smiles.toLowerCase().includes(f)));
    if(!items.length) continue;
    const h=document.createElement('div'); h.className='series'; h.textContent=SERIES_LABEL[s]||s; box.appendChild(h);
    for(const m of items){
      const d=document.createElement('div'); d.className='item'+(cur===m?' sel':''); d.textContent=m.name; d.title=m.smiles;
      d.onclick=()=>select(m); box.appendChild(d);
    }
  }
}
function xyzString(m,c){
  let s=m.sym.length+"\n"+m.name+"\n";
  for(let i=0;i<m.sym.length;i++){const p=c.xyz[i]; s+=m.sym[i]+" "+p[0]+" "+p[1]+" "+p[2]+"\n";}
  return s;
}
function chargeColor(q){
  const t=Math.max(-1,Math.min(1,q/0.5));
  const mix=(a,b,u)=>Math.round(a+(b-a)*u);
  let r,g,b;
  if(t<0){const u=-t; r=mix(255,42,u); g=mix(255,120,u); b=mix(255,214,u);} else {r=mix(255,227,t); g=mix(255,73,t); b=mix(255,72,t);}
  return "rgb("+r+","+g+","+b+")";
}
function clearExtras(){
  for(const s of shapes) viewer.removeShape(s); shapes=[];
  for(const l of labels) viewer.removeLabel(l); labels=[];
  for(const s of surf) viewer.removeShape(s); surf=[];
}
function draw(){
  if(!cur) return;
  const c=cur.confs[ci];
  viewer.removeAllModels(); clearExtras();
  model=viewer.addModel(xyzString(cur,c),'xyz');
  const atoms=model.selectedAtoms({});
  if(document.getElementById('cbQ').checked){
    atoms.forEach((a,i)=>{ viewer.setStyle({serial:a.serial},{stick:{radius:0.14,color:'#bbb'},sphere:{radius:0.32,color:chargeColor(c.q[i])}}); });
  } else {
    viewer.setStyle({},{stick:{radius:0.17},sphere:{scale:0.26}});
  }
  // 앵커·말단 강조
  if(cur.a!=null) viewer.addStyle({serial:atoms[cur.a].serial},{sphere:{scale:0.34}});
  if(document.getElementById('cbLabel').checked){
    atoms.forEach((a,i)=>labels.push(viewer.addLabel(cur.sym[i]+i,{position:a,fontSize:10,backgroundOpacity:0.4,inFront:true})));
  }
  // 쌍극자
  if(document.getElementById('cbDip').checked){
    let cx=0,cy=0,cz=0; for(const p of c.xyz){cx+=p[0];cy+=p[1];cz+=p[2];} const n=c.xyz.length; cx/=n;cy/=n;cz/=n;
    const k=0.3, mu=c.mu, L=Math.hypot(mu[0],mu[1],mu[2]);
    shapes.push(viewer.addArrow({start:{x:cx-mu[0]*k/2,y:cy-mu[1]*k/2,z:cz-mu[2]*k/2},end:{x:cx+mu[0]*k/2,y:cy+mu[1]*k/2,z:cz+mu[2]*k/2},radius:0.09,radiusRatio:2.2,mid:0.75,color:'#6250d6'}));
    labels.push(viewer.addLabel("μ = "+L.toFixed(2)+" D",{position:{x:cx+mu[0]*k/2,y:cy+mu[1]*k/2,z:cz+mu[2]*k/2},fontSize:12,fontColor:'#6250d6',backgroundColor:'white',backgroundOpacity:0.7,inFront:true}));
  }
  // 쯔비터이온 거리
  if(cur.npos!=null && cur.t!=null){
    const p=c.xyz[cur.npos], q=c.xyz[cur.t];
    shapes.push(viewer.addLine({start:{x:p[0],y:p[1],z:p[2]},end:{x:q[0],y:q[1],z:q[2]},dashed:true,color:'#008300',linewidth:3}));
    const d=Math.hypot(p[0]-q[0],p[1]-q[1],p[2]-q[2]);
    labels.push(viewer.addLabel("N⁺–X "+d.toFixed(2)+" Å",{position:{x:(p[0]+q[0])/2,y:(p[1]+q[1])/2,z:(p[2]+q[2])/2},fontSize:12,fontColor:'#008300',backgroundColor:'white',backgroundOpacity:0.7,inFront:true}));
  }
  // 궤도
  const iso=parseFloat(document.getElementById('iso').value);
  const haveCube=Object.keys(cur.cubes).length>0;
  document.getElementById('cbHomo').disabled=!haveCube||ci!==0; document.getElementById('cbLumo').disabled=!haveCube||ci!==0;
  document.getElementById('cubenote').textContent = haveCube ? (ci===0?'':'(궤도는 최저 E conformer 에서만)') : '(이 분자는 궤도 cube 없음)';
  if(haveCube && ci===0){
    for(const lab of ['homo','lumo']){
      if(!document.getElementById('cb'+lab[0].toUpperCase()+lab.slice(1)).checked) continue;
      const txt=cur.cubes[lab];
      surf.push(viewer.addVolumetricData(txt,'cube',{isoval:iso,color:'#2a78d6',opacity:0.75,smoothness:3}));
      surf.push(viewer.addVolumetricData(txt,'cube',{isoval:-iso,color:'#e34948',opacity:0.75,smoothness:3}));
    }
  }
  viewer.render();
  document.getElementById('confinfo').textContent="conformer "+(ci+1)+"/"+cur.confs.length+"  ΔE = "+c.e.toFixed(2)+" kcal/mol  (w="+c.w.toFixed(3)+")"+(c.bad?"  ⚠ 결합 변화: "+c.changes.join(' '):"");
  document.getElementById('prev').disabled=ci===0; document.getElementById('next').disabled=ci>=cur.confs.length-1;
}
function select(m){
  const same=(cur===m);
  cur=m; ci=0;
  document.getElementById('mname').textContent=m.name;
  document.getElementById('msmiles').textContent=m.smiles+"  ·  "+(SERIES_LABEL[m.series]||m.series);
  const t=document.getElementById('props'); t.innerHTML='';
  for(const [k,lab] of PROP_KEYS){
    const v=m.props[k]; if(v===null||v===undefined||v==='') continue;
    const tr=document.createElement('tr'); tr.innerHTML="<td>"+lab+"</td><td>"+(typeof v==='number'?(Number.isInteger(v)?v:v.toFixed(3)):v)+"</td>"; t.appendChild(tr);
  }
  buildList(document.getElementById('search').value);
  draw(); if(!same){viewer.zoomTo(); viewer.render();}
}
document.getElementById('search').oninput=e=>buildList(e.target.value);
document.getElementById('prev').onclick=()=>{if(ci>0){ci--;draw();}};
document.getElementById('next').onclick=()=>{if(cur&&ci<cur.confs.length-1){ci++;draw();}};
for(const id of ['cbDip','cbQ','cbHomo','cbLumo','cbLabel']) document.getElementById(id).onchange=draw;
document.getElementById('iso').oninput=e=>{document.getElementById('isov').textContent=e.target.value; draw();};
document.addEventListener('keydown',e=>{if(e.key==='ArrowRight')document.getElementById('next').click(); if(e.key==='ArrowLeft')document.getElementById('prev').click();});
buildList(''); if(DATA.length) select(DATA.find(m=>m.series==='A')||DATA[0]);
</script>
</body></html>
"""


def build(out="viewer.html"):
    os.chdir(HERE)
    mols = collect()
    js = open("assets/3Dmol-min.js").read()
    html = (TEMPLATE.replace("__3DMOL__", js)
            .replace("__DATA__", json.dumps(mols, ensure_ascii=False, separators=(",", ":")))
            .replace("__PROPKEYS__", json.dumps(PROP_KEYS, ensure_ascii=False)))
    open(out, "w").write(html)
    print(f"{out}: 분자 {len(mols)}개, 궤도 cube {sum(bool(m['cubes']) for m in mols)}개, {os.path.getsize(out)/1e6:.1f} MB")


if __name__ == "__main__":
    build()
