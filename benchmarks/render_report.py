"""Render an offline interactive report from executed JSON/TSV results."""
import csv
import gzip
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def table_reader(path):
    return path.open() if path.exists() else gzip.open(path.with_suffix('.tsv.gz'),'rt')


def main():
    bench = json.loads((ROOT/'results/real/real_benchmark.json').read_text())
    joint = json.loads((ROOT/'results/real_topk.json').read_text())
    labels = json.loads((ROOT/'data/gene_labels.json').read_text())['labels']
    studies = []
    for s in bench['summaries']:
        acc = s['accession']
        ids = s['practical_gene_audit']['baseline_top20_gene_ids']
        with table_reader(ROOT/f'results/real/{acc}_gene_audit.tsv') as f:
            genes = {r['feature_id']:r for r in csv.DictReader(f,delimiter='\t') if r['feature_id'] in ids}
        with table_reader(ROOT/f'results/real/{acc}_full_extrema.tsv') as f:
            intervals = [r for r in csv.DictReader(f,delimiter='\t') if r['feature_id'] in ids]
        studies.append(dict(summary=s,ids=ids,genes=genes,intervals=intervals,
                            audits=[r for r in joint['datasets'] if r['accession']==acc]))
    data = json.dumps(dict(studies=studies,labels=labels),ensure_ascii=False).replace('</','<\\/')
    page = r'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ExtremaRank · 配对组学候选稳健性</title>
<style>
*{box-sizing:border-box}body{font:16px/1.65 system-ui,sans-serif;margin:0;color:#14283b;background:#f5f8fa}main{max-width:1140px;margin:auto;padding:35px 28px 60px}h1{font-size:46px;line-height:1.15;margin:12px 0}h2{font-size:24px;margin:0 0 15px}h3{font-size:18px}.eyebrow{color:#136e83;font-size:12px;letter-spacing:2px;font-weight:750}.lead{font-size:21px;color:#52687a}.grid{display:grid;grid-template-columns:1fr 1fr;gap:20px}.card{background:white;border:1px solid #dce5ea;border-radius:13px;padding:25px;margin-top:20px}.metrics{display:flex;gap:35px;flex-wrap:wrap}.metric{font-size:30px;font-weight:750;color:#136e83}.metric small{display:block;font-size:12px;font-weight:500;color:#52687a}.note{font-size:13px;color:#52687a}a{color:#136e83}label{display:inline-block;margin:0 18px 12px 0}select{padding:7px;border:1px solid #bdced7;border-radius:5px;font-size:14px}input{accent-color:#136e83}.formula{font:16px/1.6 ui-monospace,monospace;background:#edf3f5;border-radius:6px;padding:15px;overflow:auto}.two{color:#136e83;font-weight:700}.replacement{color:#b45435;font-weight:700}.pill{background:#fdece6;color:#b45435;font-size:12px;padding:4px 8px;border-radius:5px}.bar{display:block;height:8px;background:#136e83;border-radius:3px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:9px;text-align:left;border-bottom:1px solid #e5ecef}th{font-size:12px;color:#52687a}#geneTable{overflow:auto}.foot{border-top:1px solid #dce5ea;margin-top:30px;padding-top:15px}@media(max-width:800px){main{padding:25px 16px}.grid{grid-template-columns:1fr}h1{font-size:35px}.card{padding:18px}}
</style><main><div class="eyebrow">RESEARCH SOFTWARE · v0.1.0 · 2026-10-03</div><h1>ExtremaRank</h1><p class="lead">少数供体能否改变你的候选名单？给出可复算的最坏情况与共同删除反例。</p><p>每次删除都重新计算均值与方差。双族极值规则把单基因组合枚举变为排序后的线性候选扫描；共享搜索审计整个 Top-K 集合。</p>
<div class="metrics"><div class="metric">200,000<small>独立精确极值校验 · 0 不一致</small></div><div class="metric">29,925<small>完整审计的真实数据基因</small></div><div class="metric">22× / 55×<small>固定比较 · 同目标精确穷举</small></div></div>
<div class="grid"><section class="card"><h2>为何计算量会下降</h2><p class="formula">R = Σx / √(Σx²)</p><p>固定保留数时，R 与有符号配对 t 的顺序相同。极值存在于连续块和前缀＋后缀两族中；去重后，中间基数只需 D 个候选。</p><label>供体数 <input id="n" type="range" min="6" max="64" value="18"> <b id="nv"></b></label><br><label>删除数 <input id="b" type="range" min="1" max="8" value="4"> <b id="bv"></b></label><div id="complexity"></div><p class="note">这里是数学候选计数，大组合未全部运行。完整共享 Top-K 搜索仍可能指数复杂。</p></section><section class="card"><h2>两族都有实际必要性</h2><div id="ablation"></div><p class="note">恰好删除四个供体，单独采用一族会给出错误极值。比例来自全部合格基因。</p><a href="../docs/THEOREM.md">定理与证明</a> · <a href="../docs/PRIOR_ART.md">先行方法</a></section></div>
<section class="card"><h2>方向稳定，名单仍可变化</h2><label>数据 <select id="study"><option value="0">GSE87290 · PBMC / LPS · 14 配对</option><option value="1">GSE50760 · 正常 / 原发结直肠癌 · 18 配对</option></select></label><label>目标 <select id="direction"><option value="up">上调 · 有符号 t</option><option value="down">下调 · 有符号 t</option><option value="absolute">双向 · |t|</option></select></label><div id="witness"></div><p class="note">每条反例只使用一套共同删除供体，全部基因排名已独立复算。它表示输入敏感性，不是生物学真值，也不是自动剔除样本的依据。</p></section>
<section class="card"><h2>原始上调 Top-20 的分数区间</h2><label>恰好删除 <input id="r" type="range" min="0" max="4" value="0"> <b id="rv"></b> 个供体</label><p class="note">每个基因的区间端点可能由不同删除集合达到；区间重叠本身不是共同排名反例。</p><div id="geneTable"></div></section>
<section class="card"><h2>用于你的配对效应矩阵</h2><p>CSV、TSV 和 gzip；第一列 donor_id，其余为基因。支持上调、下调、绝对值排序。核心无第三方依赖，无需 GPU。</p><pre class="formula">python -m pip install .<br>extremarank effects.csv --direction absolute --budget 2 --top-k 20 --output audit</pre><p>输出逐基因最坏分数、方向稳健性、均值与方差、共同供体反例，以及 CERTIFIED / REFUTED / UNRESOLVED 状态。</p><a href="../README.md">README</a> · <a href="../docs/USAGE.md">用法</a> · <a href="../docs/REPORT_CN.md">中文报告</a></section><p class="foot note">已执行的软件及公开数据验证。精确性针对输入 binary64 值；只审计冻结配对效应和普通配对 t 排名。未证明 FDR、疾病预测或临床效用提升。检索支持具体方法区别，尚未建立全球优先权；全文缺口保留在文档中。本报告离线可查看。</p></main>
<script>const data=__DATA__;const el=id=>document.getElementById(id);const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));const label=id=>data.labels[id]?data.labels[id].symbol+' ('+id+')':id;const fmt=x=>Number.isFinite(Number(x))?Number(x).toFixed(3):x;
function comb(n,b){let out=1n;for(let i=1;i<=b;i++)out=out*BigInt(n-i+1)/BigInt(i);return out}function complexity(){const n=Number(el('n').value),limit=Math.min(8,n-2);el('b').max=limit;if(Number(el('b').value)>limit)el('b').value=limit;const b=Number(el('b').value);el('nv').textContent=n;el('bv').textContent=b;el('complexity').innerHTML=`<p>穷举：<b>${comb(n,b).toLocaleString()}</b> 个保留子集<br>双族去重：<b class="two">${n}</b> 个候选</p>`}
function render(){const study=data.studies[Number(el('study').value)],r=Number(el('r').value);el('rv').textContent=r;const a=study.audits.find(x=>x.direction===el('direction').value);el('witness').innerHTML=`<p><span class="pill">共同删除反例已复核</span> 删除 <b>${a.deleted_donor_ids.map(esc).join(', ')}</b> 后，Top-20 更换 ${a.entered_gene_ids.length} 个成员。</p><p>退出：<span class="replacement">${a.exited_gene_ids.map(x=>esc(label(x))).join('、')}</span><br>进入：<span class="two">${a.entered_gene_ids.map(x=>esc(label(x))).join('、')}</span></p><p class="note">${a.genes.toLocaleString()} 个基因，审计本机 ${a.audit_seconds.toFixed(2)} 秒；基线与反例均以独立 Fraction 运算复算。</p>`;const rows=study.ids.map(id=>{const e=study.intervals.find(x=>x.feature_id===id&&Number(x.deletions)===r),g=study.genes[id];return `<tr><td>${esc(label(id))}</td><td>${fmt(g.baseline_t)}</td><td>${fmt(e.minimum_t)}</td><td>${fmt(e.maximum_t)}</td><td>${g.robust_positive==='True'?'是':'否'}</td></tr>`}).join('');el('geneTable').innerHTML='<table><thead><tr><th>基因</th><th>原始 t</th><th>最小 t</th><th>最大 t</th><th>≤4 删除仍为正</th></tr></thead><tbody>'+rows+'</tbody></table>'}
el('ablation').innerHTML=data.studies.map(x=>{const s=x.summary,lo=s.full_spectrum_ablations.contiguous_only['4'].minimum_wrong_rate*100,hi=s.full_spectrum_ablations.ends_only['4'].maximum_wrong_rate*100;return `<h3>${s.accession}</h3><p>仅连续块：<b>${lo.toFixed(2)}%</b> 下界错误<span class="bar" style="width:${lo*3}%"></span></p><p>仅两端：<b>${hi.toFixed(2)}%</b> 上界错误<span class="bar" style="width:${hi*3}%"></span></p>`}).join('');['n','b'].forEach(id=>el(id).addEventListener('input',complexity));['study','direction','r'].forEach(id=>el(id).addEventListener('input',render));complexity();render();</script></html>'''.replace('__DATA__',data)
    (ROOT/'results/report.html').write_text(page,encoding='utf-8')
    ds=list(range(6,65,2))
    def points(values):
        return ' '.join(f'{70+(d-6)/58*580:.2f},{310-v/7*240:.2f}' for d,v in zip(ds,values))
    svg=['<svg xmlns="http://www.w3.org/2000/svg" width="740" height="390" viewBox="0 0 740 390">','<rect width="740" height="390" fill="white"/>','<g font-family="Arial,sans-serif" fill="#14283b">','<text x="70" y="30" font-size="20">Exact subset candidates · deletion count 4</text>','<text x="70" y="52" font-size="12">Mathematical counts; shared search can still be exponential</text>']
    for v in range(8):
        y=310-v/7*240
        svg += [f'<line x1="70" y1="{y}" x2="650" y2="{y}" stroke="#e2e9ed"/>',f'<text x="22" y="{y+4}" font-size="12">10^{v}</text>']
    svg += [f'<polyline points="{points([math.log10(math.comb(d,4)) for d in ds])}" fill="none" stroke="#b45435" stroke-width="3"/>',f'<polyline points="{points([math.log10(d) for d in ds])}" fill="none" stroke="#136e83" stroke-width="3"/>','<text x="420" y="100" fill="#b45435">Exhaustive choose(D,4)</text>','<text x="420" y="258" fill="#136e83">Two-family: D candidates</text>']
    for d in (6,18,32,48,64):
        svg.append(f'<text x="{70+(d-6)/58*580-8}" y="332" font-size="12">{d}</text>')
    svg += ['<text x="300" y="363" font-size="14">Biological donor count D</text></g></svg>']
    (ROOT/'results/candidate_scaling.svg').write_text('\n'.join(svg)+'\n')
    print('Wrote report.html and candidate_scaling.svg')


if __name__=='__main__':
    main()
