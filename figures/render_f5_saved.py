"""Render F5 from six saved per-spot arrays; this entry never solves OT.

The arrays are from the dated prospective fixed-input rerun recorded in
figures/current_release/inputs/F5_input_binding.json. They are not
claimed to be the unarchived per-spot arrays from the historical run.
"""
import argparse
import hashlib
import importlib.util
from pathlib import Path
import sys
import numpy as np

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--input-npz",type=Path,required=True)
 parser.add_argument("--output-png",type=Path,required=True)
 args=parser.parse_args()
 if args.output_png.name!="lrot_figure8_dlpfc_score_ranking.png":
  parser.error("Output filename must match the existing F5 renderer.")
 repo=Path(__file__).resolve().parents[1]
 if not (repo/"figures/gen_figure8_dlpfc_score_ranking.py").is_file():
  repo=Path.cwd()
 spec=importlib.util.spec_from_file_location("f5_original_plot",repo/"figures/gen_figure8_dlpfc_score_ranking.py")
 module=importlib.util.module_from_spec(spec)
 sys.modules[spec.name]=module
 spec.loader.exec_module(module)
 with np.load(args.input_npz,allow_pickle=False) as z:
  required=[f"{name}_{field}" for name in ("FGW","LROT") for field in ("score","correct","row_entropy")]
  if set(z.files)!=set(required):
   raise ValueError("F5 plot input members differ from the six pinned arrays")
  pooled={name:{field:z[f"{name}_{field}"] for field in ("score","correct","row_entropy")} for name in ("FGW","LROT")}
 for name in pooled:
  if any(len(a)!=2000 or not np.isfinite(a).all() for a in pooled[name].values()):
   raise ValueError("Invalid F5 saved array length or finite values")
 metrics=[]
 for pair_index,(sid_a,sid_b) in enumerate(module.SECTION_PAIRS):
  lo=pair_index*module.N_SPOTS;hi=lo+module.N_SPOTS
  for name in ("FGW","LROT"):
   score=pooled[name]["score"][lo:hi]
   correct=pooled[name]["correct"][lo:hi]
   top20=np.argsort(score)[-int(round(len(score)*0.2)):]
   metrics.append(
    f"PAIR\t{sid_a}->{sid_b}\t{name}\tn={len(score)}\t"
    f"spearman_rho={module.spearmanr(score,correct).statistic:.3f}\t"
    f"accuracy={correct.mean():.3f}\t"
    f"top20_accuracy={correct[top20].mean():.3f}"
   )
 args.output_png.parent.mkdir(parents=True,exist_ok=True)
 fig,summary=module.make_figure(pooled,metrics,args.output_png.parent)
 if fig.resolve()!=args.output_png.resolve():
  raise RuntimeError("Renderer wrote an unexpected file")
 print("saved",fig)
 print("figure_sha256",hashlib.sha256(fig.read_bytes()).hexdigest())
 print("summary_sha256",hashlib.sha256(summary.read_bytes()).hexdigest())
if __name__=="__main__":
 main()
