"""Apply approved S2D parameter labels to a locked raster panel only.

Counts and bars are never redrawn. Parameter text is read from six accepted
diagnostics. This operation is neither a solver run nor experimental replication.
"""
from pathlib import Path
import json, hashlib

TEXT_MASKS=((982,183,1283,261),(982,283,1283,361))

def _number(x):
    s=format(x,'g')
    return s[1:] if s.startswith('0.') else s

def recorded_settings(records):
    if len(records)!=6 or not all(r['accepted'] is True for r in records):
        raise ValueError('Six accepted recorded diagnostics are required; no numeric fallback.')
    pairs=list(dict.fromkeys((r['parameters']['alpha'],r['parameters']['gamma']) for r in records))
    resources=list(dict.fromkeys(r['resource'] for r in records))
    if len(pairs)!=2 or len(resources)!=3 or any(sum(r['resource']==name for r in records)!=2 for name in resources):
        raise ValueError('Recorded resource/parameter pairing is incomplete.')
    if any(r['alpha']!=r['parameters']['alpha'] or r['gamma']!=r['parameters']['gamma'] for r in records):
        raise ValueError('Parameter summary disagrees with the recorded settings.')
    return pairs

def apply(template,records_path,output,font_regular):
    from PIL import Image,ImageDraw,ImageFont,ImageChops
    records=json.loads(Path(records_path).read_text(encoding='utf-8'))
    pairs=recorded_settings(records)
    source=Image.open(template).convert('RGB');im=source.copy();d=ImageDraw.Draw(im)
    font=ImageFont.truetype(str(font_regular),64)
    labels=['α='+_number(a)+', γ='+_number(g) for a,g in pairs]
    for box,label in zip(TEXT_MASKS,labels):
        if font.getlength(label)>box[2]-box[0]-8:
            raise ValueError('Recorded parameter label does not fit the approved text strip.')
        d.rectangle((box[0],box[1],box[2]-1,box[3]-1),fill='white')
        d.text((box[0]+8,box[1]+10),label,font=font,fill='black',anchor='lt')
    difference=ImageChops.difference(source,im)
    for box in TEXT_MASKS:difference.paste((0,0,0),box)
    unchanged=difference.getbbox() is None
    if not unchanged:raise ValueError('S2D pixels outside approved parameter text strips changed.')
    im.save(output,dpi=(300,300))
    return {'labels':labels,'outside_parameter_text_masks_unchanged':unchanged,
            'counts_redrawn':False,'record_keys':[r['key'] for r in records],
            'recorded_counts':[r['outer_iterations'] for r in records],
            'records_sha256':hashlib.sha256(Path(records_path).read_bytes()).hexdigest(),
            'verification':'parameter-label patch on retained current panel; not end-to-end experiment reproduction'}
