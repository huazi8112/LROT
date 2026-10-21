from pathlib import Path
import subprocess,sys,json,time,os
# Historical batch orchestration is opt-in, never the default formal API.
if os.environ.get("LROT_ENABLE_HISTORICAL_BATCH") != "1":
    raise RuntimeError("Historical batch execution is disabled; use the formal paper_protocol API. Explicit opt-in is required.")
if not os.environ.get("LROT_DATA_ROOT"):
    raise RuntimeError("LROT_DATA_ROOT is required before starting any historical batch jobs.")
B=Path(__file__).parent
jobs=[('run_dlpfc_baselines.py',''),('run_official_paste2.py',''),('run_paste2_multiseed.py',''),('lr_weight_sensitivity.py',''),('preproc_sensitivity.py',''),('cancer_erbb2_ablation_truecoord.py',''),('cancer_pathway_ablation_scan_truecoord.py',''),('lr_gene_set_fairness.py',''),('lr_cellchat_comparison.py',''),('lrot_3d_multi_seed.py',''),('lrot_3d_reconstruction.py',''),('run_dlpfc_3d.py',''),('batch_effect_scan.py',''),('spatial_distortion_scan.py',''),('structural_batch_scan.py',''),('robustness_failure_scan.py',''),('robustness_cost_diagnosis.py',''),('supp/_supp_s17_checks.py','sweep'),('supp/_supp_s17_checks.py','eps'),('supp/_supp_s17_checks.py','control'),('pathway_activity_gamma_sweep.py',''),('pathway_activity_dlpfc.py',''),('run_dlpfc_full_scale.py',''),('lrot_scaling.py','')]
rows=[]
for driver,stage in jobs:
    name=Path(driver).stem+('_'+stage if stage else '')
    status=B/'runs'/name/'driver_status.json'
    if status.exists():
        old=json.loads(status.read_text(encoding='utf-8'))
        if old.get('completed'):rows.append(old);continue
    print('START',name,flush=True);t=time.time()
    with (B/(name+'.log')).open('w',encoding='utf-8') as log:
        p=subprocess.run([sys.executable,'-B','-X','utf8',str(B/'portable_run_driver.py'),driver]+([stage] if stage else []),stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
    rows.append({'driver':driver,'stage':stage,'exit_code':p.returncode,'seconds':time.time()-t})
    (B/'queue_status.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
    print('END',name,p.returncode,round(time.time()-t,1),flush=True)
print('QUEUE_FINISHED',flush=True)
