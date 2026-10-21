"""Offline GO BP enrichment for the fixed manuscript resources.

Run with explicit --obo, --gaf and a new --output-dir. Uses deduplicated
gene annotations, BP filtering, NOT exclusion, is_a/part_of propagation,
foreground-complement Fisher tests and BH over all background-supported terms.
"""
import gzip
import os
import pickle
import sys
from collections import defaultdict
import numpy as np
import networkx as nx
from scipy.stats import fisher_exact
from pathlib import Path
import argparse,json,hashlib,socket
cli=argparse.ArgumentParser(description="Recompute paper GO BP enrichment from explicit offline OBO/GAF inputs.")
cli.add_argument("--obo",type=Path,required=True)
cli.add_argument("--gaf",type=Path,required=True)
cli.add_argument("--output-dir",type=Path,required=True)
cli.add_argument("--code-root",type=Path,default=Path(__file__).resolve().parents[1])
args=cli.parse_args()
ROOT=args.code_root.resolve(); OUT=args.output_dir.resolve()
if OUT.exists() or OUT==ROOT or ROOT in OUT.parents:cli.error("Use a new output directory outside the code package.")
if not args.obo.is_file() or not args.gaf.is_file():cli.error("Both offline input files must exist.")
pin=ROOT/'provenance/go_enrichment/resource_identity.json'
if pin.exists():
 expected=json.loads(pin.read_text('utf-8'))
 for f,k in [(args.obo,'ontology'),(args.gaf,'annotation')]:
  actual=hashlib.sha256(f.read_bytes()).hexdigest()
  if actual!=expected[k]['sha256']:cli.error('Input differs from the fixed paper snapshot: '+f.name)
OUT.mkdir(parents=True)
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT))
from lrot_core import LIGANDS,RECEPTORS
OBO=str(args.obo.resolve()); GAF=str(args.gaf.resolve())
def no_network(*a,**k):raise RuntimeError("Offline analysis: networking is disabled")
socket.create_connection=no_network;socket.socket.connect=no_network;socket.socket.connect_ex=no_network

def parse_obo(path):
    """返回 (name, parents)：仅 biological_process 术语。"""
    name, parents, namespaces = ({}, defaultdict(set), {})
    cur, ns = (None, None)
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if line == '[Term]':
                cur, ns = (None, None)
            elif line.startswith('id: '):
                cur = line[4:].strip()
            elif line.startswith('name: '):
                if cur:
                    name[cur] = line[6:].strip()
            elif line.startswith('namespace: '):
                ns = line[11:].strip()
                namespaces[cur] = ns
            elif line.startswith('is_a: ') and cur and (ns == 'biological_process'):
                parents[cur].add(line[6:].split(' ')[0].strip())
            elif line.startswith('relationship: part_of ') and cur and (ns == 'biological_process'):
                parents[cur].add(line.split(' ')[2].strip())
    bp = {t for t in name if namespaces.get(t) == 'biological_process'}
    return ({t:name[t] for t in bp}, {t:set(p for p in parents.get(t,set()) if p in bp) for t in bp})

def ancestors_of(parents_bp):
    """按拓扑序计算每个术语的全部祖先（true-path rule）。

    注意：边建为 t -> parent(t)，拓扑序会把子节点排在父节点之前，因此必须
    **逆拓扑序** 遍历，才能保证用到的 anc[parent] 已经算好。
    """
    G = nx.DiGraph()
    for t, ps in parents_bp.items():
        for p in ps:
            G.add_edge(t, p)
    anc = {}
    for t in reversed(list(nx.topological_sort(G))):
        s = set()
        for p in G.successors(t):
            s.add(p)
            s |= anc.get(p, set())
        anc[t] = s
    return anc

def parse_gaf(path, bp_terms):
    """返回 gene_symbol -> set(直接注释的 BP 术语)。"""
    g2t = defaultdict(set)
    n = 0
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            qual, goid, aspect, sym = (c[3], c[4], c[8], c[2])
            if aspect != 'P' or goid not in bp_terms:
                continue
            if 'NOT' in qual.split('|'):
                continue
            g2t[sym].add(goid)
            n += 1
    return (g2t, n)

def run_analysis():
    print('=' * 74)
    print('  GO BP enrichment of LROT LR genes — v2 (real background, deduped)')
    print('=' * 74)
    print('\n[1/5] 解析 GO 本体 …')
    name, parents = parse_obo(OBO)
    bp_terms = set(name)
    print('    BP 术语 %d 个' % len(bp_terms))
    anc = ancestors_of({t: parents.get(t, set()) for t in bp_terms})
    print('    祖先表就绪（平均每术语 %.1f 个祖先）' % np.mean([len(v) for v in anc.values()]))
    print('\n[2/5] 解析 MGI GAF（小鼠注释）…')
    g2t, n_ann = parse_gaf(GAF, bp_terms)
    print('    注释行 %d；基因 %d 个' % (n_ann, len(g2t)))
    print('\n[3/5] 沿本体传播（true-path rule）…')

    def propagate(terms):
        s = set(terms)
        for t in terms:
            s |= anc.get(t, set())
        return s
    bg_sets = {g: propagate(ts) for g, ts in g2t.items()}
    background = set(bg_sets)
    print('    背景基因集 %d 个（有 BP 注释的小鼠基因）' % len(background))
    bg_count = defaultdict(int)
    for g, s in bg_sets.items():
        for t in s:
            bg_count[t] += 1
    lr = sorted({g for g in set(LIGANDS) | set(RECEPTORS)})
    lr_canon = lr
    print('    LR 基因 %d 个（符号已为标准小鼠符号 Artn / Nrtn 等）' % len(lr))
    fg_sets = {g: bg_sets.get(g, set()) for g in lr_canon}
    annotated = [g for g in lr_canon if fg_sets[g]]
    missing = [g for g in lr_canon if not fg_sets[g]]
    print('    注释成功 %d/%d %s' % (len(annotated), len(lr_canon), '缺失 %s' % missing if missing else ''))
    fg_count = defaultdict(int)
    fg_genes = defaultdict(list)
    for g in annotated:
        for t in sorted(fg_sets[g]):
            fg_count[t] += 1
            fg_genes[t].append(g)
    print('\n[4/5] Fisher 精确检验 + BH 校正 …')
    n_fg = len(annotated)
    n_bg = len(background)
    tests = []
    # Predefine all BP terms represented in the frozen annotation background.
    # Terms with no foreground hit remain in the multiple-testing family (p=1).
    for t, c in bg_count.items():
        a = fg_count.get(t,0)
        if c < 1:
            continue
        b = n_fg - a
        d = n_bg - n_fg - c + a
        if d < 0:
            continue
        if c < a or d < 0:raise ValueError('Invalid foreground/background partition')
        orr, p = fisher_exact([[a, b], [c-a, d]], alternative='greater')
        tests.append(dict(go=t, term=name.get(t, ''), a=a, c=c, p=p, fold=a / n_fg / (c / n_bg)))
    tests.sort(key=lambda r: (r['p'], r['go']))
    m = len(tests)
    for i, r in enumerate(tests):
        r['padj'] = min(1.0, r['p'] * m / (i + 1))
    for i in range(m - 2, -1, -1):
        tests[i]['padj'] = min(tests[i]['padj'], tests[i + 1]['padj'])
    print('    检验术语 %d 个（分子=%d 基因，分母=背景 %d 基因）' % (m, n_fg, n_bg))
    sig = [r for r in tests if r['padj'] < 0.05 and r['a'] >= 3]
    print('    padj<0.05 且命中≥3 基因：%d 个术语' % len(sig))
    print('\n[5/5] 写出结果 …')
    lines = []
    lines.append('=' * 78)
    lines.append('  GO BP enrichment — LROT LR genes (v2: real background, deduped)')
    lines.append('=' * 78)
    lines.append('')
    lines.append('Annotation source : GO go-basic.obo + mgi.gaf.gz (MGI, mouse)')
    lines.append('Aspect            : biological_process only; NOT qualifiers excluded')
    lines.append('Foreground        : %d LR genes (%d/%d annotated)' % (n_fg, len(annotated), len(lr_canon)))
    lines.append('Background        : %d mouse genes with >=1 BP annotation' % n_bg)
    lines.append('True-path rule    : annotations propagated to ancestors (is_a/part_of)')
    lines.append('Test              : Fisher exact (one-sided greater); BH over all %d tested terms' % m)
    lines.append('Reported          : padj < 0.05 AND >=3 foreground genes')
    lines.append('')
    lines.append('%-12s %-54s %4s %7s %8s %8s %10s' % ('GO ID', 'Term', 'fg', 'bg', 'fold', 'p', 'p_adj'))
    lines.append('-' * 112)
    for r in sig:
        lines.append('%-12s %-54s %4d %7d %8.1f %8.1e %10.1e' % (r['go'], r['term'][:54], r['a'], r['c'], r['fold'], r['p'], r['padj']))
    lines.append('')
    lines.append('--- Top 25 by p (all tested terms) ---')
    for r in tests[:25]:
        lines.append('%-12s %-54s fg=%3d bg=%6d fold=%6.1f p=%.2e padj=%.2e' % (r['go'], r['term'][:54], r['a'], r['c'], r['fold'], r['p'], r['padj']))
    lines.append('')
    lines.append('--- 前景命中基因（padj<0.05, >=3 基因）---')
    for r in sig[:30]:
        lines.append('%s  %s' % (r['go'], r['term']))
        lines.append('    %s' % ', '.join(sorted(fg_genes[r['go']])))
    out_txt = os.path.join(OUT, 'lrot_go_enrichment_results.txt')
    with open(out_txt, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print('    -> %s' % out_txt)
    with open(os.path.join(OUT, 'lrot_go_enrichment_data.pkl'), 'wb') as f:
        pickle.dump(dict(sig=sig, tests=tests, n_fg=n_fg, n_bg=n_bg, fg_genes=dict(fg_genes), annotated=annotated, missing=missing), f)
    print('\n★ 关键术语（供正文引用）★')
    keys = ['axon guidance', 'angiogenesis', 'canonical Wnt signaling pathway', 'Notch signaling pathway', 'positive regulation of cell population proliferation', 'cell migration', 'positive regulation of cell migration', 'ephrin receptor signaling pathway', 'semaphorin-plexin signaling pathway', 'viral process']
    for k in keys:
        hit = [r for r in tests if r['term'].lower() == k.lower()]
        if hit:
            r = hit[0]
            print('   %-52s fg=%2d bg=%5d fold=%5.1f padj=%.1e%s' % (r['term'][:52], r['a'], r['c'], r['fold'], r['padj'], '' if r['padj'] < 0.05 else '   (不显著)'))
    return 0


def source_identity(path):
 p=Path(path);return {"file_name":p.name,"size_bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()}
if __name__=="__main__":
 run_analysis()
 import csv
 d=pickle.loads((OUT/"lrot_go_enrichment_data.pkl").read_bytes())
 with (OUT/"all_terms.csv").open("w",encoding="utf-8",newline="") as f:
  w=csv.DictWriter(f,fieldnames=["go","term","a","c","p","padj","fold"]);w.writeheader();w.writerows(d["tests"])
 obo_headers=[];gaf_headers=[]
 with open(OBO,encoding="utf-8") as f:
  for line in f:
   if line.startswith("[Term]"):break
   if line.strip():obo_headers.append(line.strip())
 with gzip.open(GAF,"rt",encoding="utf-8") as f:
  for line in f:
   if not line.startswith("!"):break
   gaf_headers.append(line.rstrip())
 identity={"ontology":source_identity(OBO),"annotation":source_identity(GAF),"ontology_headers":obo_headers,"annotation_headers":gaf_headers,"foreground_genes":sorted(set(LIGANDS)|set(RECEPTORS)),"n_foreground":d["n_fg"],"n_background":d["n_bg"],"tested_terms":len(d["tests"]),"zero_hit_terms":sum(r["a"]==0 for r in d["tests"]),"significant_hit3":len(d["sig"]),"test":"one-sided Fisher exact foreground versus background complement","multiple_testing":"BH over every background-supported BP term, including zero-foreground terms","input_policy":"Frozen local files; verification date is not acquisition date; historical download history is not inferred.","code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (OUT/"resource_identity.json").write_text(json.dumps(identity,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
 print(json.dumps({k:identity[k] for k in ["n_foreground","n_background","tested_terms","zero_hit_terms","significant_hit3"]}),flush=True)
