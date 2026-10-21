"""Prepare pinned external inputs in a NEW directory outside the code package.

STAligner preparation is offline. CellChat preparation only downloads when
--download is explicitly supplied; --rda provides a fully offline alternative.
No dependency installation, training, solver execution or existing-file overwrite.
"""
from __future__ import annotations
import argparse, ast, hashlib, json, pickle, re, shutil
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
RECIPE = ROOT / 'provenance/f3_staligner/source_recipe'
MANIFEST = ROOT / 'data/external_sources.json'

def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def apply_patch(original, patch):
    lines = original.decode('utf-8').splitlines(keepends=True)
    changes = patch.decode('utf-8').replace('\r\n','\n').splitlines(keepends=True)
    result, cursor, index = [], 0, 0
    while index < len(changes):
        match = re.match(r'^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@', changes[index])
        if not match:
            index += 1
            continue
        start = int(match[1]) - 1
        if start < cursor:
            raise ValueError('Overlapping patch hunks')
        result.extend(lines[cursor:start]); cursor = start
        old_count = new_count = 0; index += 1
        while index < len(changes) and not changes[index].startswith('@@ '):
            line = changes[index]
            if not line or line[0] not in ' +-':
                raise ValueError('Unsupported patch record')
            prefix, value = line[0], line[1:]
            if prefix in ' -':
                if cursor >= len(lines) or lines[cursor] != value:
                    raise ValueError('Patch context differs from pinned upstream source')
                cursor += 1; old_count += 1
            if prefix in ' +':
                result.append(value); new_count += 1
            index += 1
        if old_count != int(match[2] or 1) or new_count != int(match[4] or 1):
            raise ValueError('Patch hunk length mismatch')
    result.extend(lines[cursor:])
    return ''.join(result).encode('utf-8')

def prepare_staligner(output, manifest):
    source = manifest['staligner']
    for name, digest in source['upstream_source_sha256'].items():
        if checksum(RECIPE/'upstream/STAligner'/name) != digest:
            raise ValueError('Bundled upstream source hash mismatch: ' + name)
    patch = RECIPE/'mnn_utils.patch'
    if checksum(patch) != source['patch_sha256']:
        raise ValueError('Patch hash mismatch')
    package = output/'STAligner'; package.mkdir()
    for name, digest in source['training_source_sha256'].items():
        data = (RECIPE/'upstream/STAligner'/name).read_bytes()
        if name == 'mnn_utils.py':
            data = apply_patch(data, patch.read_bytes())
        if hashlib.sha256(data).hexdigest() != digest:
            raise ValueError('Reconstructed training-source mismatch: ' + name)
        (package/name).write_bytes(data)
    shutil.copy2(RECIPE/'LICENSE', output/'LICENSE-STAligner.txt')
    # Metadata for this newly reconstructed source tree, not a claim about the
    # historical installer. It lets importlib.metadata find the pinned version
    # when this directory is used via PYTHONPATH without another STAligner install.
    metadata = output/'STAligner-1.0.0.dist-info'; metadata.mkdir()
    (metadata/'METADATA').write_bytes(
        b'Metadata-Version: 2.1\nName: STAligner\nVersion: 1.0.0\n'
        b'Summary: Fixed upstream STAligner with disclosed LROT neighbour patch\n'
        b'Home-page: https://github.com/zhoux85/STAligner\nLicense: MIT\n')
    shutil.copy2(RECIPE/'LICENSE',metadata/'LICENSE')
    (metadata/'source_recipe.json').write_text(json.dumps(source,indent=2)+'\n',encoding='utf-8')
    return {'resource':'STAligner','offline':True,'upstream_commit':source['commit'],
            'source_sha256':{n:checksum(package/n) for n in source['training_source_sha256']},
            'all_six_training_sources_match':True,'dependencies_installed':False,
            'training_executed':False}

def prepare_cellchat(output, args, manifest):
    source = manifest['cellchat']
    rda = output/'CellChatDB.mouse.rda'
    if args.download:
        # The URL is fixed in a reviewed manifest; no local content is transmitted.
        handlers = [urllib.request.ProxyHandler({'http':args.proxy,'https':args.proxy})] if args.proxy else []
        opener = urllib.request.build_opener(*handlers)
        request = urllib.request.Request(source['rda_url'], headers={'User-Agent':'LROT-pinned-resource-fetch'})
        with opener.open(request, timeout=60) as response:
            rda.write_bytes(response.read())
    else:
        shutil.copy2(args.rda, rda)
    if checksum(rda) != source['rda_sha256']:
        raise ValueError('RDA hash mismatch; output retained for inspection, not accepted')
    builder = ROOT/'tools/build_cellchat_db.py'
    if checksum(builder) != source['builder_sha256']:
        raise ValueError('Bound legacy builder source differs')
    # Execute only the hash-bound original conversion definitions and pair loop.
    # The original script's main, synthetic generator, imports and writes do not run.
    import rdata
    converted = rdata.conversion.convert(rdata.parser.parse_file(rda))
    inter = converted['CellChatDB.mouse']['interaction']
    namespace = {'lig':inter['ligand'].astype(str).values,
                 'rec':inter['receptor'].astype(str).values}
    nodes = []
    for node in ast.parse(builder.read_text(encoding='utf-8')).body:
        if isinstance(node, ast.FunctionDef) and node.name == 'mouse_symbol':
            nodes.append(node)
        elif isinstance(node, ast.Assign):
            targets = {v.id for t in node.targets for v in ast.walk(t) if isinstance(v,ast.Name)}
            if targets & {'SUBUNIT_FIX','NON_GENE','db','n_fixed','n_dropped'}:
                nodes.append(node)
        elif (isinstance(node,ast.For) and isinstance(node.iter,ast.Call)
              and isinstance(node.iter.func,ast.Name) and node.iter.func.id=='zip'
              and isinstance(node.target,ast.Tuple)
              and all(isinstance(v,ast.Name) for v in node.target.elts)
              and [v.id for v in node.target.elts]==['l','r']):
            nodes.append(node)
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(builder),'exec'),namespace)
    genes_file = ROOT/source['derived_gene_names_path']
    if checksum(genes_file) != source['derived_gene_names_sha256']:
        raise ValueError('Derived gene-list hash mismatch')
    genes = set(json.loads(genes_file.read_text(encoding='utf-8')))
    full = namespace['db']
    valid = {k:v for k,v in full.items() if k[0] in genes and k[1] in genes}
    pkl = output/'CellChatDB_mouse_lr.pkl'
    pkl.write_bytes(pickle.dumps({'db_full':full,'db_valid_data':valid},protocol=4))
    if checksum(pkl) != source['derived_pkl_sha256']:
        raise ValueError('Derived PKL differs; output retained, not accepted')
    shutil.copy2(ROOT/'data/cellchat_LICENSE.txt',output/'LICENSE-CellChat.txt')
    return {'resource':'CellChatDB.mouse','upstream_commit':source['commit'],
            'RDA_sha256':checksum(rda),'derived_PKL_sha256':checksum(pkl),
            'full_pairs':len(full),'valid_pairs':len(valid),
            'download_performed':bool(args.download),'rdata_version':rdata.__version__,
            'solver_executed':False,'historical_acquisition_date_recovered':False}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='resource',required=True)
    st = sub.add_parser('staligner',help='Offline reconstruction of the six fixed training sources')
    st.add_argument('--output-dir',type=Path,required=True)
    cc = sub.add_parser('cellchat',help='Prepare the fixed public RDA and byte-identical derived PKL')
    cc.add_argument('--output-dir',type=Path,required=True)
    inputs = cc.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--rda',type=Path,help='Offline copy of the pinned public RDA')
    inputs.add_argument('--download',action='store_true',help='Explicitly download the pinned public RDA')
    cc.add_argument('--proxy',help='Optional HTTP(S) proxy URL; does not change system settings')
    args = parser.parse_args();output = args.output_dir.resolve()
    if output.exists() or output==ROOT or ROOT in output.parents:
        parser.error('Output must be a new directory outside the code package')
    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    output.mkdir(parents=True,exist_ok=False)
    result = prepare_staligner(output,manifest) if args.resource=='staligner' else prepare_cellchat(output,args,manifest)
    (output/'preparation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))

if __name__=='__main__':
    main()
