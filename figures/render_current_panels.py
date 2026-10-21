"""Only drawing operations on saved arrays/metrics; no experiment imports."""
from pathlib import Path

def render(fig_id, data, output):
    if fig_id == 'F8':
        return composite(data, output)
    import numpy as np
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':15,'axes.labelsize':16,
                         'axes.titlesize':17,'xtick.labelsize':14,'ytick.labelsize':14,
                         'svg.fonttype':'none','pdf.fonttype':42})
    blue='#4E79A7';orange='#F28E2B'
    def primary(ax,d,s):
        ax.axhspan(*s['ci95_t'],color=blue,alpha=.15)
        ax.axhline(s['mean_difference'],color=blue,lw=2)
        ax.axhline(0,color='gray',ls='--')
        ax.scatter(range(1,31),d,c=[orange if x>0 else blue for x in d],s=40)
        ax.set(xlabel='Simulation pair index',ylabel='LROT − soft FGW accuracy',ylim=(-.25,.7))
        ax.set_xticks([1,5,10,15,20,25,30]);ax.grid(axis='y',alpha=.2)
    if fig_id in ['F3','S23']:
        d=data['paired_differences'];s=data['summary']
        assert len(d)==30 and abs(np.mean(d)-s['mean_difference'])<1e-12
        if fig_id=='S23':
            fig,ax=plt.subplots(figsize=(10,5));primary(ax,d,s)
            ax.set_title('30 controlled pairs, historical 42/123 excluded')
            fig.subplots_adjust(bottom=.2,left=.12)
        else:
            fig=plt.figure(figsize=(11.5,9));gs=fig.add_gridspec(2,2,left=.08,right=.98,top=.94,bottom=.12,hspace=.63,wspace=.37)
            names=['Expr OT','Soft FGW','LROT\nγ = 0.1','LROT\nγ = 0.2','PASTE2','STAligner\nlocal NN*']
            values=data['accuracy'];cols=['#999999',blue,orange,'#C46E13','#B07AA1','#59A14F']
            ax=fig.add_subplot(gs[0,0]);ax.bar(range(6),values,color=cols,width=.7)
            ax.set_xticks(range(6),names,fontsize=14,rotation=30,ha='right')
            ax.set(ylabel='Region-label agreement',title='A  Realistic seed 42',ylim=(0,1.06))
            for i,y in enumerate(values):ax.text(i,y+.025,f'{y:.3f}',ha='center',fontsize=14)
            ax=fig.add_subplot(gs[0,1]);hh=data['entropy'];ax.bar(range(4),hh,color=cols[:4],width=.65)
            ax.set_xticks(range(4),names[:4],fontsize=14)
            ax.set(ylabel='Mean row entropy (nats)',title='B  Transport dispersion',ylim=(0,5))
            for i,y in enumerate(hh):ax.text(i,y+.12,f'{y:.3f}',ha='center',fontsize=14)
            ax=fig.add_subplot(gs[1,:]);primary(ax,d,s)
            ax.set_title('C  Controlled benchmark: 30 simulation pairs')
            fig.text(.51,.01,f"Mean paired difference {s['mean_difference']:+.3f}; 95% CI [{s['ci95_t'][0]:+.3f}, {s['ci95_t'][1]:+.3f}]; 24/30 positive",ha='center',fontsize=15)
    elif fig_id=='S6':
        trials=data['trials'];methods=['Expression OT','Soft FGW','LROT','Official PASTE','Official PASTE2']
        vals=[[r['seconds'] for r in trials if r['method']==m] for m in methods]
        assert len(trials)==15 and all(len(v)==3 for v in vals)
        fig,ax=plt.subplots(figsize=(10,5))
        ax.bar(range(5),[np.mean(v) for v in vals],yerr=[np.std(v,ddof=1) for v in vals],capsize=5,color=['#999999',blue,orange,'#59A14F','#B07AA1'])
        ax.set_xticks(range(5),['Expr OT','Soft FGW','LROT','PASTE (local)\nKL, α=.1','PASTE2 (local)\nKL, α=.1'])
        ax.set(ylabel='Preparation + computation (s)',title='Local CPU timing: three serial repetitions',ylim=(0,max(map(max,vals))*1.3))
        fig.subplots_adjust(left=.12,bottom=.21)
        for i,v in enumerate(vals):ax.text(i,np.mean(v)+np.std(v,ddof=1)+.25,f'{np.mean(v):.2f}',ha='center')
    elif fig_id=='S12':
        fig,axes=plt.subplots(1,2,figsize=(11.5,5.3));fig.subplots_adjust(left=.09,right=.97,bottom=.25,top=.89,wspace=.35)
        for r,c in zip(data['CG'],[blue,orange]):
            axes[0].plot(range(len(r['loss'])),r['loss'],'o-',color=c,label=f"PASTE α={r['parameters']['alpha']:g}",lw=2,ms=4)
        for tag,c in [('Soft FGW',blue),('LROT',orange)]:
            loss=data['soft_loss'][tag];axes[1].plot(range(1,len(loss)+1),loss,'o-',color=c,label=tag,lw=2,ms=4)
        axes[0].set(xlabel='CG iteration',ylabel='Package-reported CG loss',title='A  Tested installed PASTE')
        axes[1].set(xlabel='Outer iteration',ylabel='Weighted cost diagnostic',title='B  Soft-plan diagnostics')
        for ax in axes:ax.legend(loc='upper center',bbox_to_anchor=(.5,-.25),ncol=2,frameon=False,fontsize=14)
    else:
        raise ValueError('No validated renderer for '+fig_id)
    dest=Path(output)/('Figure_'+fig_id+'.png')
    for ext in ['png','svg','pdf']:
        fig.savefig(dest.with_suffix('.'+ext),dpi=350,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    return dest

def composite(data,output):
    from PIL import Image, ImageDraw, ImageFont
    root=Path(__file__).resolve().parents[1]
    imgs=[]
    for p,label in zip(data['panels'],['A','B','C']):
        with Image.open(root/p) as im:
            im=im.convert('RGB');im=im.resize((round(im.width*1250/im.height),1250),Image.Resampling.LANCZOS)
            imgs.append((im,label))
    canvas=Image.new('RGB',(sum(im.width for im,_ in imgs)+52*2+26*2,1302),'white')
    x=26;xs=[]
    for im,_ in imgs:canvas.paste(im,(x,26));xs.append(x);x+=im.width+52
    draw=ImageDraw.Draw(canvas)
    try:font=ImageFont.truetype('arialbd.ttf',112)
    except OSError:font=ImageFont.load_default()
    for (im,label),px in zip(imgs,xs):
        box=draw.textbbox((0,0),label,font=font);tw,th=box[2]-box[0],box[3]-box[1];lx,ly=px+14,40
        draw.rectangle([lx-6,ly-6,lx+tw+10,ly+th+8],fill='white')
        draw.text((lx,ly),label,font=font,fill=(20,30,40))
    dest=Path(output)/'Figure_F8.png';canvas.save(dest,dpi=(300,300));return dest
