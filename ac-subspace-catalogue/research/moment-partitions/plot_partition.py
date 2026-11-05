"""Static scientific illustration of the exact N=22 prescribed-moment partition."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,
                    'svg.fonttype':'none','figure.facecolor':'#f5f4ee','axes.facecolor':'#f5f4ee'})
fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
x=np.linspace(-11,11,400)
ax.plot(x,x*x,color='#b7beb9',lw=1.5,label=r'Moment curve $(m,m^2)$')
blocks=[([11,-4],'#276356','Pair: Dicke indices 0, 15'),
        ([4,-11],'#976c36','Pair: Dicke indices 7, 22'),
        ([8,1,-7],'#78619b','Triangle: Dicke indices 3, 10, 18')]
for ms,color,label in blocks:
    xx=np.array(ms,dtype=float); yy=xx*xx
    if len(xx)==3:ax.fill(xx,yy,color=color,alpha=.12)
    ax.plot(np.r_[xx,xx[0]],np.r_[yy,yy[0]],color=color,lw=2,label=label)
    ax.scatter(xx,yy,color=color,s=45,zorder=4)
ax.scatter([0],[44],s=110,marker='*',color='#172a26',zorder=8)
ax.annotate(r'Prescribed point $(0,44)$',xy=(0,44),xytext=(-1.8,86),
            arrowprops={'arrowstyle':'->','color':'#172a26'},ha='center')
ax.set(xlabel=r'Magnetic level $m=11-r$',ylabel=r'Second moment $m^2$',
       xlim=(-12,12),ylim=(-5,130),title='Three convex blocks share the AC2 moment point')
ax.grid(alpha=.15)
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2,frameon=False,fontsize=10)
fig.savefig(HERE/'moment-partition-n22.svg')
fig.savefig(HERE/'moment-partition-n22.png',dpi=160)

svg_path=HERE/'moment-partition-n22.svg'
svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines())+'\n')
