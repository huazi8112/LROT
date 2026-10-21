"""Apply approved local figure edits only. Never import experiment or solver modules."""
from pathlib import Path
import json,math,statistics
from PIL import Image,ImageDraw,ImageFont,ImageChops
FONT_REGULAR=None
FONT_BOLD=None
BUNDLE=None
def font(size,bold=False):return ImageFont.truetype(str(FONT_BOLD if bold else FONT_REGULAR),max(8,int(size)))
def text(im,xy,t,size=40,anchor='mm',fill='black',bold=False):ImageDraw.Draw(im).text(xy,t,font=font(size,bold),anchor=anchor,fill=fill)
def rotated(im,xy,t,size=40):
 f=font(size);b=f.getbbox(t);tmp=Image.new('RGBA',(b[2]-b[0]+10,b[3]-b[1]+16));ImageDraw.Draw(tmp).text((5-b[0],8-b[1]),t,font=f,fill='black');tmp=tmp.rotate(90,expand=True);im.paste(tmp,(int(xy[0]-tmp.width/2),int(xy[1]-tmp.height/2)),tmp)
def clear(im,box):ImageDraw.Draw(im).rectangle(tuple(map(int,box)),fill='white')

def apply(number,bundle,output,font_regular,font_bold):
    global FONT_REGULAR,FONT_BOLD,BUNDLE
    FONT_REGULAR=font_regular;FONT_BOLD=font_bold;BUNDLE=bundle
    payload={}
    names={'S3':'S3_go_display_rows.json','S6':'S6_serial_observations.json','S21':'S21_verified_result_rows.json','S23_FGW':'S23_1001_2001_FGW.json','S23_LROT':'S23_1001_2001_LROT.json','S23_old':'S23_archived_second_point.json'}
    selected=['S'+str(number)] if number in (3,6,21) else []
    for key in selected:payload[key]=json.loads((bundle/'inputs'/names[key]).read_text(encoding='utf-8'))
    proofs=[]
    def picture(kind,n):
        manifest=json.loads((bundle/'figure_manifest.json').read_text(encoding='utf-8'))
        row=next(r for r in manifest['figures'] if r['id']=='S'+str(n))
        return Image.open(bundle.parents[1]/row['asset']['path']).convert('RGB')
    def save(n,im,reason,checks=None):
        im.save(output,dpi=(300,300));proofs.extend(checks or [])
    def region_equal(n,im,box,kind='baseline'):
        original=picture(kind,n)
        return {'figure':'S'+str(n),'box':box,'pixel_identical':im.crop(box).tobytes()==original.crop(box).tobytes()}
    if number == 3:
        # S3: replace only the statistically corrected A region; keep all original B pixels including legend.
        im=picture('baseline',3);dr=ImageDraw.Draw(im);cut=2484
        clear(im,(0,85,cut-1,im.height))
        go=payload['S3']['rows']
        assert len(go)==15
        x0,x1,y0,y1=918,2263,96,1750; pitch=102;centers=[211+i*pitch for i in range(15)]
        dr.rectangle((x0,y0,x1,y1),outline='#333333',width=3)
        def cmap(t):
         # Same red/yellow/green palette as the baseline, interpolated in the saved corrected significance scale.
         stops=[(0,(165,0,38)),(.25,(244,109,67)),(.5,(255,255,191)),(.75,(102,189,99)),(1,(0,104,55))]
         for (a,c),(b,d) in zip(stops,stops[1:]):
          if a<=t<=b:return tuple(round(u+(v-u)*(t-a)/(b-a)) for u,v in zip(c,d))
         return stops[-1][1]
        vmin,vmax=42,54
        for row,y in zip(go,centers):
         val=-math.log10(row['padj']);colour=cmap((val-vmin)/(vmax-vmin));end=x0+(x1-x0)*row['a']/55
         dr.rectangle((x0+2,y-40,end,y+40),fill=colour,outline='#aaaa88',width=2)
         term=row['term']
         # Match the baseline's one-line labels; shrink only labels that would otherwise clip.
         size=36
         while font(size).getlength(term)>x0-35:size-=1
         text(im,(x0-22,y),term,size,'rm');text(im,(end+7,y),f"{row['a']}* [{val:.0f}]",29,'lm')
        for t in [0,10,20,30,40,50]:
         x=x0+(x1-x0)*t/55;dr.line((x,y1,x,y1+14),fill='black',width=3);text(im,(x,y1+45),str(t),38)
        for y in range(342,1505):
         t=1-(y-342)/(1505-342);dr.line((2296,y,2351,y),fill=cmap(t),width=1)
        dr.rectangle((2296,342,2351,1505),outline='black',width=2)
        for v in [42,44,46,48,50,52,54]:
         y=1505-(v-vmin)/(vmax-vmin)*(1505-342);dr.line((2352,y,2365,y),fill='black',width=2);text(im,(2373,y),str(v),32,'lm')
        rotated(im,(2455,920),'−log10(adjusted p)',35)
        text(im,((x0+x1)/2,1860),'Number of LR Genes',42)
        save(3,im,'A仅更新已核实的GO排序、计数、校正p及色标；B完整基线区域含图例保持逐像素一致。',[region_equal(3,im,(cut,0,im.width,im.height))])
        
        return output, proofs
    if number == 6:
        # S6: a single chart, baseline method colours and logarithmic axis, saved three-call means/SD only.
        im=Image.new('RGB',picture('baseline',6).size,'white');dr=ImageDraw.Draw(im);w,h=im.size
        x0,x1,y0,y1=210,w-70,95,h-330
        stats={}
        for method in {t['method'] for t in payload['S6']['trials']}:
         values=[t['seconds'] for t in payload['S6']['trials'] if t['method']==method]
         if len(values)!=3:raise ValueError('S6 requires three recorded calls per method')
         stats[method.replace('Official ','')]={'mean':statistics.mean(values),'sample_SD':statistics.stdev(values)}
        names=['LROT','Soft FGW','PASTE','PASTE2','Expression OT'];labels=['LROT\nSoft FGW + LR','Soft FGW','PASTE (local)\nKL, α=0.1','PASTE2 (local)\nKL, α=0.1, s=1','Expression OT']
        colors=['#2e8b57','#e18d8d','#d27536','#ac89bc','#929292']
        def sy(v):return y1-math.log(v)/(math.log(20))*(y1-y0)
        for v in [1,2,5,10,20]:
         y=sy(v);dr.line((x0,y,x1,y),fill='#e5e5e5',width=2);text(im,(x0-25,y),str(v),40,'rm')
        dr.line((x0,y0,x0,y1,x1,y1),fill='black',width=3)
        step=(x1-x0)/5
        for i,(n,l,c) in enumerate(zip(names,labels,colors)):
         r=stats[n];x=x0+(i+.5)*step;bw=step*.60;m=r['mean'];sd=r['sample_SD']
         dr.rectangle((x-bw/2,sy(m),x+bw/2,y1),fill=c)
         dr.line((x,sy(m-sd),x,sy(m+sd)),fill='black',width=4)
         for y in [sy(m-sd),sy(m+sd)]:dr.line((x-18,y,x+18,y),fill='black',width=4)
         text(im,(x,sy(m+sd)-44),f'{m:.3f} s',44,bold=True)
         for j,line in enumerate(l.split('\n')):text(im,(x,y1+57+j*55),line,40)
        rotated(im,(65,(y0+y1)/2),'Recorded local time (s; log scale)',43)
        save(6,im,'仅以15条已核验本地计时更新五柱、样本SD和标签；沿用单柱图、原方法配色、log轴；不含STAligner训练及加速箭头。')
        
        return output, proofs
    if number == 12:
        # S12 retains original curve pixels, with small text masks only.
        im=picture('baseline',12);w,h=im.size
        patches=[(0,0,w//2,95),(w//2,0,w,95),(12,255,112,h-185),(1968,220,2048,h-170),(920,145,1935,270),(2100,h-112,w,h)]
        for b in patches:clear(im,b)
        approved=Image.open(BUNDLE/'assets/S12.png').convert('RGB')
        for b in patches:im.paste(approved.crop(b),b[:2])
        old=picture('baseline',12);mask=Image.new('1',im.size,0);md=ImageDraw.Draw(mask)
        for b in patches:md.rectangle(b,fill=1)
        outside=ImageChops.difference(old,im);outside.paste((0,0,0),(0,0,w,h),mask)
        save(12,im,'原曲线、配色、面板及索引保留；仅将Official与objective/loss标签改为安装版CG loss和不同的加权成本诊断。',[{'figure':'S12','outside_text_masks_unchanged':outside.getbbox() is None}])
        
        return output, proofs
    if number == 20:
        # S20: only D/F labels; original A/B/C/E and all underlying scientific pixels retained.
        im=picture('baseline',20);w,h=im.size
        patches=[(w*.52,h*.34,w*.95,h*.382),(w*.55,h*.661,w*.95,h*.706)]
        for b in patches:clear(im,b)
        text(im,(w*.75,h*.356),'D: Top-30 LR compatibility',51)
        text(im,(w*.75,h*.373),'(bright = greater compatibility)',45)
        text(im,(w*.75,h*.682),'F: Same source spots: raw / mapped',49)
        text(im,(w*.75,h*.699),'(z=0 raw, z=1 mapped display state)',43)
        old=picture('baseline',20);mask=Image.new('1',im.size,0);md=ImageDraw.Draw(mask)
        for b in patches:md.rectangle(tuple(map(int,b)),fill=1)
        diff=ImageChops.difference(old,im);diff.paste((0,0,0),(0,0,w,h),mask)
        save(20,im,'保持原3行2列与原始点云、热图和直方图；仅D兼容性标签、F同一source状态标题；原轴文字不变。',[{'figure':'S20','outside_text_masks_unchanged':diff.getbbox() is None}])
        
        return output, proofs
    if number == 21:
        # S21: same original two-by-two chart types and colours; redraw only required affected metric panels from audited result points.
        im=picture('baseline',21);w,h=im.size;dr=ImageDraw.Draw(im)
        rr=payload['S21']['rows'];xx=[r['gamma'] for r in rr];assert xx==[0,.05,.1,.2,.5,1.]
        acc=[r['accuracy'] for r in rr];cost=[r['expression_cost_in_result_record'] for r in rr];ent=[r['entropy_recomputed'] for r in rr];iters=[r['outer_iterations'] for r in rr]
        def chart(rect,title,ylabel,vals,lo,hi,ticks,color,bar=False,fmt='.3f'):
         x,y,X,Y=rect;clear(im,rect);d=ImageDraw.Draw(im);a=x+235;b=X-55;top=y+125;bot=Y-140;sz=46
         text(im,((a+b)/2,y+58),title,52)
         d.rectangle((a,top,b,bot),outline='black',width=3)
         def py(v):return bot-(v-lo)/(hi-lo)*(bot-top)
         for v in ticks:
          p=py(v);d.line((a-15,p,a,p),fill='black',width=3);text(im,(a-25,p),f'{v:g}',sz,'rm')
         rotated(im,(x+80,(top+bot)/2),ylabel,sz)
         points=[];step=(b-a)/6
         for i,v in enumerate(vals):
          px=a+(i+.5)*step;yy=py(v);points.append((px,yy));text(im,(px,bot+63),'γ='+str(xx[i]).rstrip('0').rstrip('.') if xx[i] else 'γ=0',43)
          if bar:d.rectangle((px-step*.38,yy,px+step*.38,bot-2),fill=(color[i] if isinstance(color,list) else color))
          text(im,(px,max(top+20,yy-28)),format(v,fmt),43)
         if not bar:
          d.line(points,fill=color,width=8)
          for px,yy in points:d.ellipse((px-12,yy-12,px+12,yy+12),fill=color)
        chart((0,0,w//2,h//2),'A: Region-vote accuracy vs γ','Accuracy',acc,0,1,[0,.2,.4,.6,.8,1],['#6b9abe']*5+['#ffb92e'],True)
        chart((w//2,0,w,h//2),'B: Expression transport cost vs γ','Transport cost',cost,.435,.532,[.44,.46,.48,.50,.52],'#dc143c')
        chart((0,h//2,w//2,h),'C: Mean row entropy vs γ','Mean entropy (nats)',ent,2.85,4.1,[3,3.2,3.4,3.6,3.8,4],'#008000',fmt='.2f')
        chart((w//2,h//2,w,h),'D: Accepted outer iterations vs γ','Outer iterations',iters,0,14,[0,2,4,6,8,10,12,14],'#a649aa',True,fmt='d')
        save(21,im,'保留原柱/线/线/柱、2×2与各面板配色；更新已核验点及D次数；.747在图注另列历史值，不连为同输入纠错。')
        
        return output, proofs
    if number == 23:
        raise ValueError('S23 historical displacement patch is retired; use the pinned current asset.')
    raise ValueError('Only approved local revisions are supported')
