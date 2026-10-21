"""Prepare the exact GO/MGI paper inputs offline or by explicit public download."""
from pathlib import Path
import argparse,hashlib,json,shutil,urllib.request
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output-dir',required=True,type=Path)
g=p.add_mutually_exclusive_group(required=True)
g.add_argument('--resource-dir',type=Path,help='Offline directory containing go-basic.obo and mgi.gaf.gz')
g.add_argument('--download',action='store_true',help='Explicitly download current official files and enforce fixed checksums')
p.add_argument('--proxy',help='Optional proxy for this request only; no system settings are changed')
a=p.parse_args();OUT=a.output_dir.resolve()
if OUT.exists() or OUT==ROOT or ROOT in OUT.parents:p.error('Output must be new and outside the code package.')
pin=json.loads((ROOT/'provenance/go_enrichment/resource_identity.json').read_text('utf-8'))
OUT.mkdir(parents=True)
urls={'go-basic.obo':'https://current.geneontology.org/ontology/go-basic.obo','mgi.gaf.gz':'https://current.geneontology.org/annotations/mgi.gaf.gz'}
results=[]
for name,key in [('go-basic.obo','ontology'),('mgi.gaf.gz','annotation')]:
 dest=OUT/name
 if a.download:
  handler=urllib.request.ProxyHandler({'http':a.proxy,'https':a.proxy}) if a.proxy else urllib.request.ProxyHandler()
  op=urllib.request.build_opener(handler)
  with op.open(urls[name],timeout=60) as response:dest.write_bytes(response.read())
 else:shutil.copy2(a.resource_dir/name,dest)
 actual=hashlib.sha256(dest.read_bytes()).hexdigest()
 results.append({'file':name,'sha256':actual,'matches_fixed_input':actual==pin[key]['sha256'],'public_url':urls[name],'downloaded':a.download})
 (OUT/'preparation.json').write_text(json.dumps(results,indent=2)+'\n','utf-8')
 if actual!=pin[key]['sha256']:raise ValueError('Official current resource changed or offline resource differs: '+name+'. Output retained; no substitute accepted. Use the fixed resource companion archive.')
notice='''# Fixed GO/MGI resource companion

These exact files regenerate the paper GO results and were byte-matched to official downloads on 2026-10-16. The original acquisition dates are unknown. The OBO header reports releases/2026-07-26; the GAF header reports format 2.2 and generation date 2026-05-21. Header identifiers must not be represented as recovered acquisition records.

Attribution: Gene Ontology Consortium; Mouse Genome Informatics, The Jackson Laboratory. GO data/products and MGI data/annotations declare CC BY 4.0 in their official policies:

- https://geneontology.org/docs/go-citation-policy/
- https://www.informatics.jax.org/mgihome/other/copyright.shtml
- Licence: https://creativecommons.org/licenses/by/4.0/

The resource bytes have not been modified. This companion includes annotations, not third-party website images. Preserve attribution and the SHA-256 checksums when archiving or redistributing. Current download links are mutable; a checksum mismatch must fail. The companion has no public archive identifier until the authors publish it.
'''
(OUT/'README.md').write_text(notice,'utf-8')
(OUT/'SOURCE.json').write_text(json.dumps(pin,ensure_ascii=False,indent=2)+'\n','utf-8')
(OUT/'SHA256SUMS.txt').write_text(''.join(f"{r['sha256']}  {r['file']}\n" for r in results),'utf-8')
print(json.dumps({'success':True,'downloaded':a.download,'verified_resources':results}))
