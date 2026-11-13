"""Draws Figure 8 (with panel e) and Supplementary Figures S5-S7 from results/revision."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import os
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R=os.path.join(ROOT,'results','revision')+'/'; O=os.path.join(ROOT,'supplementary')+'/'
plt.rcParams.update({'font.family':'Liberation Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.8,'savefig.dpi':300})
def lab(ax,t,x=-0.12,y=1.06): ax.text(x,y,t,transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')

# ---------- Figure 8 panel e, aligned to panel d columns
order=[8,3,6,1,7,5,4,2,0]
c=pd.read_csv(R+'Figure8e_source_cluster_MMP13.csv').set_index('cluster').loc[order]
W,H=2078,620
fig=plt.figure(figsize=(W/300,H/300),dpi=300)
x0,x1=281/W,1906/W
ax=fig.add_axes([x0,0.33,x1-x0,0.52])
x=np.arange(len(order))
ax.bar(x,c['mean'],width=0.6,color='#B2182B',alpha=0.85)
ax.errorbar(x,c['mean'],yerr=[c['mean']-c['ci_low'],c['ci_high']-c['mean']],fmt='none',ecolor='k',elinewidth=0.8,capsize=2)
ax.set_xlim(-0.5,len(order)-0.5); ax.set_xticks(x)
ax.set_xticklabels([f"{k}\n(n={int(n)})" for k,n in zip(order,c['n_cells'])],fontsize=6.5)
ax.set_ylabel('Mean $\\it{MMP13}$\n(log-normalised)',fontsize=7); ax.tick_params(labelsize=6.5); ax.set_xlabel('Leiden cluster (ordered as in d)',fontsize=7)
fig.text(205/W,0.93,'e',fontsize=11,fontweight='bold')
fig.savefig(O+'Figure8e_panel.png',dpi=300); plt.close(fig)
top=Image.open(os.path.join(ROOT,'figures','Figure8_panels_a-d.png')).convert('RGB'); e=Image.open(O+'Figure8e_panel.png').convert('RGB').resize((W,H))
new=Image.new('RGB',(W,top.height+H),'white'); new.paste(top,(0,0)); new.paste(e,(0,top.height))
new.save(O+'Figure8_with_panel_e.png',dpi=(300,300))

# ---------- Figure S5: seed comparison and centrality
s=pd.read_csv(R+'TableS6_seed_comparison.csv'); rnd=s[s.type!='named']; nm=s[s.type=='named'].set_index('seed')
net=pd.read_csv(R+'TableS6b_network_centrality.csv',index_col=0)
fig,axs=plt.subplots(1,3,figsize=(7.6,2.8),gridspec_kw={'wspace':0.75})
ax=axs[0]; ax.hist(rnd.n_robust,bins=40,color='#BDBDBD')
for g,col in [('MMP13','#B2182B'),('FAP','#2166AC'),('THY1','#4393C3'),('OLFML2B','#92C5DE'),('POSTN','#5AAE61'),('MMP3','#7F3B08')]:
    ax.axvline(max(nm.loc[g,'n_robust'],40),color=col,lw=1.4,label=f"{g} ({int(nm.loc[g,'n_robust'])})")
ax.set_xlabel('Robust partners (n)'); ax.set_ylabel('Matched random genes (n)'); ax.set_xlim(-150,7600); ax.set_xticks([0,2000,4000,6000]); ax.legend(fontsize=6.5,frameon=False,loc='upper right',handlelength=1.0,borderaxespad=0); lab(ax,'a',-0.25)
ax=axs[1]
data=[rnd.GSE114007_pos_frac.dropna(),rnd.GSE82107_pos_frac.dropna()]
vp=ax.violinplot(data,showmedians=True)
for b in vp['bodies']: b.set_facecolor('#BDBDBD'); b.set_alpha(0.8)
for k in ['cmedians','cmins','cmaxes','cbars']: vp[k].set_color('#636363'); vp[k].set_linewidth(0.8)
ax.scatter([1,2],[nm.loc['MMP13','GSE114007_pos_frac'],nm.loc['MMP13','GSE82107_pos_frac']],color='#B2182B',zorder=3,s=20,label='MMP13')
ax.set_xticks([1,2]); ax.set_xticklabels(['Cartilage\n(GSE114007)','Synovium\n(GSE82107)'])
ax.set_ylabel('Top-150 positive partners\nreplicated (fraction)'); ax.legend(fontsize=6,frameon=False,loc='upper center',bbox_to_anchor=(0.5,1.12)); ax.set_ylim(0,1.05); lab(ax,'b',-0.3)
ax=axs[2]
ax.scatter(net.degree_rank,net.eigen_rank,s=2,color='#BDBDBD',rasterized=True)
for g,col in [('MMP13','#B2182B'),('FAP','#2166AC'),('THY1','#4393C3'),('OLFML2B','#92C5DE'),('POSTN','#5AAE61')]:
    if g in net.index:
        ax.scatter(net.loc[g,'degree_rank'],net.loc[g,'eigen_rank'],color=col,s=18,zorder=3)
        ax.annotate(g,(net.loc[g,'degree_rank'],net.loc[g,'eigen_rank']),fontsize=7,xytext=(-4,4) if g in ('OLFML2B','FAP') else (4,3),ha='right' if g in ('OLFML2B','FAP') else 'left',textcoords='offset points')
ax.set_xscale('log'); ax.set_yscale('log'); ax.invert_xaxis(); ax.invert_yaxis()
ax.set_xlabel('Rank by weighted degree'); ax.set_ylabel('Rank by eigenvector centrality'); lab(ax,'c',-0.3)
fig.subplots_adjust(left=0.08,right=0.98,bottom=0.22,top=0.88); fig.savefig(O+'FigureS5_seed_centrality.png'); plt.close(fig)

# ---------- Figure S6: hypertrophy by tissue
a=pd.read_csv(R+'FigureS6_source_hypertrophy_by_tissue.csv',index_col=0); b=pd.read_csv(R+'FigureS6_source_noSPP1_ALPL.csv',index_col=0)
def ci(r,n):
    z=np.arctanh(r); se=1.06**0.5/np.sqrt(n-3); return np.tanh(z-1.96*se),np.tanh(z+1.96*se)
def nb(ds,t):
    v=b[(b.dataset==ds)&(b.tissue==t)].rho.values; return v[0] if len(v) else np.nan
TN={'cartilage':'Cartilage','meniscus':'Meniscus','synovium':'Synovium','bone':'Subchondral bone'}
rows=[]
for t in ['cartilage','meniscus','synovium','bone']:
    d=a[(a.tissue==t)&(a.dataset!='pooled')]
    for _,r in d.iterrows():
        lo,hi=ci(r.rho,r.n); rows.append((f"{TN[t]}: {r.dataset} (n={int(r.n)})",r.rho,lo,hi,nb(r.dataset,t),False))
    if len(d)>1:
        p=a[(a.tissue==t)&(a.dataset=='pooled')].iloc[0]
        rows.append((f"{TN[t]}: pooled samples (n={int(p.n)})",p.rho,p.ci_low,p.ci_high,nb('pooled',t),True))
p=a[(a.tissue=='all')].iloc[0]
rows.append((f"All tissues: pooled samples (n={int(p.n)})",p.rho,p.ci_low,p.ci_high,nb('pooled','all'),True))
fig,ax=plt.subplots(figsize=(5.4,0.3*len(rows)+0.9))
y=np.arange(len(rows))[::-1]
for yi,(l,r,lo,hi,rb,pool) in zip(y,rows):
    ax.plot([lo,hi],[yi,yi],color='k',lw=0.8)
    ax.scatter(r,yi,marker='D' if pool else 's',color='#B2182B' if pool else 'k',s=30 if pool else 16,zorder=3)
    if not np.isnan(rb): ax.scatter(rb,yi,marker='o',facecolor='white',edgecolor='#2166AC',linewidth=1,s=18,zorder=4)
ax.axvline(0,color='grey',lw=0.6,ls='--')
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows],fontsize=7)
ax.set_xlabel('Spearman ρ, MMP13 vs hypertrophy score (95% CI)'); ax.set_xlim(-0.8,1)
ax.scatter([],[],marker='s',color='k',s=16,label='Single dataset'); ax.scatter([],[],marker='D',color='#B2182B',s=30,label='Samples pooled across datasets')
ax.scatter([],[],marker='o',facecolor='white',edgecolor='#2166AC',s=18,label='Score without SPP1 and ALPL')
fig.tight_layout(rect=(0,0.07,1,1)); fig.legend(*ax.get_legend_handles_labels(),fontsize=6.5,frameon=False,loc='lower center',ncol=3,handletextpad=0.3,columnspacing=1); fig.savefig(O+'FigureS6_hypertrophy_by_tissue.png'); plt.close(fig)

# ---------- Figure S7: single-cell sensitivity
pr=pd.read_csv(R+'Figure8e_source_cluster_MMP13.csv').set_index('cluster').loc[order]
se=pd.read_csv(R+'TableS12_sensitivity_cluster_MMP13.csv').set_index('cluster').loc[order]
fig,axs=plt.subplots(1,2,figsize=(7.2,2.6))
x=np.arange(len(order)); w=0.38
for ax,col,ylabel,L in [(axs[0],'mean','Mean MMP13 (log-normalised)','a'),(axs[1],'pct_detected','Cells with detectable MMP13 (%)','b')]:
    for d,off,c_,name in [(pr,-w/2,'#B2182B','Primary analysis'),(se,w/2,'#2166AC','Sensitivity analysis')]:
        ax.bar(x+off,d[col],width=w,color=c_,alpha=0.85,label=name)
        if col=='mean': ax.errorbar(x+off,d['mean'],yerr=[d['mean']-d['ci_low'],d['ci_high']-d['mean']],fmt='none',ecolor='k',elinewidth=0.7,capsize=1.5)
    ax.set_xticks(x); ax.set_xticklabels([str(k) for k in order]); ax.set_xlabel('Leiden cluster'); ax.set_ylabel(ylabel); lab(ax,L,-0.15)
axs[0].legend(fontsize=6.5,frameon=False)
axs[1].text(0.98,1.08,'1,597 cells in both analyses\nAdjusted Rand index = 1.00',transform=axs[1].transAxes,ha='right',va='top',fontsize=6.5)
fig.tight_layout(); fig.savefig(O+'FigureS7_sc_sensitivity.png'); plt.close(fig)
print('ok', new.size)
