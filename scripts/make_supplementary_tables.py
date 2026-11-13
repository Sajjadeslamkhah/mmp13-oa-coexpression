"""Builds supplementary/Supplementary_Tables_S5-S12.xlsx from results/revision and data/."""
import pandas as pd, numpy as np
from openpyxl import load_workbook, Workbook
from openpyxl.styles import Font, Alignment
import os
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R=os.path.join(ROOT,'results','revision')+'/'
out=os.path.join(ROOT,'supplementary','Supplementary_Tables_S5-S12.xlsx')
sheets=[]  # (name, title, df)
pv=pd.read_csv(os.path.join(ROOT,'data','assessed_series.csv'))
role={'curated_assessment (not search-derived)':'Discovery','curated_assessment (pre-specified validation)':'Validation (pre-specified, held-out)','curated_assessment (single-cell localisation)':'Single-cell localisation'}
fx=lambda v:str(v).replace('not available','Not recorded').replace('not applicable','Not applicable')
df5=pd.DataFrame({'Accession':pv.accession,'Role':pv.stage.map(role),'Platform':pv.platform.map(fx),'Tissue':pv.tissue,'OA (n)':pv.n_OA.map(fx),'Control (n)':pv.n_control.map(fx),'Decision':pv.decision.str.capitalize(),'Reason / note':pv.reason})
sheets.append(('S5','Table S5. All GEO series assessed for inclusion. Datasets were selected from a pre-specified list; no systematic database search was performed.',df5))
s6=pd.read_csv(R+'TableS6_seed_comparison.csv').rename(columns={'seed':'Seed gene','type':'Seed type','n_tested':'Genes tested','n_pos':'Robust positive partners','n_neg':'Robust negative partners','n_robust':'Robust partners (total)',
 'GSE114007_pos_frac':'GSE114007: fraction of top-150 positive partners with positive ρ','GSE114007_neg_frac':'GSE114007: fraction of top-150 negative partners with negative ρ','GSE114007_genome_r':'GSE114007: genome-wide r',
 'GSE82107_pos_frac':'GSE82107: fraction of top-150 positive partners with positive ρ','GSE82107_neg_frac':'GSE82107: fraction of top-150 negative partners with negative ρ','GSE82107_genome_r':'GSE82107: genome-wide r'})
s6=s6.dropna(axis=1,how='all')
s6['Seed type']=s6['Seed type'].map({'named':'Named gene','random':'Expression-matched random gene'}).fillna(s6['Seed type'])
sheets.append(('S6a','Table S6a. Co-expression meta-analysis repeated with other seed genes (same thresholds as the primary analysis): named genes and 1,000 random genes matched to MMP13 on mean expression and variance within each cohort. Empty cells: not computed for that seed. Fractions use all 150 partners as the denominator; Table 3 of the main text uses the partners measured in each validation dataset, so the two can differ slightly (e.g., 78.7% here vs 79.7% for MMP13 in GSE114007).',s6))
n=pd.read_csv(R+'TableS6b_network_centrality.csv',index_col=0).rename_axis('Gene').reset_index().rename(columns={'weighted_degree':'Weighted degree (β = 6)','eigenvector':'Eigenvector centrality','degree_rank':'Rank by weighted degree','eigen_rank':'Rank by eigenvector centrality'})
sheets.append(('S6b','Table S6b. Genome-wide weighted co-expression network (5,000 most variable genes, which included MMP13; cross-cohort meta-correlations; soft threshold β = 6).',n))
CT={'lining_fibroblast':'Synovial lining fibroblast','T_cell':'T cell','pericyte_SMC':'Pericyte/smooth muscle cell','B_plasma':'B/plasma cell'}
fix=lambda v: CT.get(v,v.replace('_',' ').capitalize())
m=pd.read_csv(R+'TableS7a_markers.csv',index_col=0); m=m.rename_axis('Cell type').reset_index(); m.columns=['Cell type','Marker genes']; m['Cell type']=m['Cell type'].map(fix)
sheets.append(('S7a','Table S7a. Marker genes used for the cell-type composition scores (approach of MCP-counter). Marker genes that were locked MMP13 partners were removed before scoring. The panels are approximations; the osteoblast panel consists largely of osteocyte genes.',m))
c0=pd.read_csv(R+'TableS7b_mmp13_vs_composition.csv',index_col=0)
metarows=[]
for ct,sub in c0.groupby('cell_type',sort=False):
    t=sub.copy(); t['z']=np.arctanh(t.rho.clip(-.999,.999))
    mg=t[t.dataset.isin(['GSE117999','GSE98918'])]
    t=pd.concat([t[~t.dataset.isin(['GSE117999','GSE98918'])][['z','n']],pd.DataFrame([{'z':mg.z.mean(),'n':mg.n.mean()}])]) if len(mg) else t[['z','n']]
    w=t.n-3; metarows.append({'dataset':'Meta-analysis (6 cohorts)','tissue':'All','cell_type':ct,'rho':np.tanh((t.z*w).sum()/w.sum()),'p':np.nan,'n':np.nan})
c=pd.concat([c0,pd.DataFrame(metarows)],ignore_index=True)
c['cell_type']=c['cell_type'].map(fix)
c=c.rename(columns={'dataset':'Dataset','tissue':'Tissue','cell_type':'Cell type','rho':'Spearman ρ with MMP13','p':'P','n':'Samples (n)'})
sheets.append(('S7b','Table S7b. Spearman correlation of MMP13 with each cell-type composition score within OA samples, per dataset (GSE169077: OA and control samples combined, because it has fewer than eight OA samples), and the cohort-level meta-correlation (Fisher z, weights n − 3; GSE117999 and GSE98918 averaged as one cohort).',c))
COLS={'SYMBOL':'Gene','meta_rho':'Meta ρ','meta_z':'Meta Z','n_cohorts':'Cohorts (n)','same_sign_frac':'Fraction of cohorts with same sign','p':'P','fdr':'FDR (BH)','robust':'Robust partner'}
for nm,f,t in [('S7c','TableS7c_coexpr_composition_adjusted_2PC.csv','Table S7c. MMP13 co-expression meta-analysis within OA tissue, adjusted for the first two principal components of the composition scores (partial correlations).'),
               ('S7d','TableS7c_coexpr_composition_adjusted_1PC.csv','Table S7d. As Table S7c, adjusted for the first principal component only (sensitivity analysis).'),
               ('S8a','TableS8a_coexpr_control.csv','Table S8a. MMP13 co-expression meta-analysis within the control samples of the discovery datasets.')]:
    sheets.append((nm,t,pd.read_csv(R+f).rename(columns=COLS)))
d=pd.read_csv(R+'TableS8b_OA_vs_control_difference.csv'); d=d[d.SYMBOL!='MMP13'].rename(columns={'SYMBOL':'Gene','z_diff':'Z (OA minus control, Fisher z)','p':'P','fdr':'FDR (BH)'})
sheets.append(('S8b','Table S8b. Gene-wise test of whether the correlation with MMP13 differs between OA and control tissue.',d))
x=pd.read_csv(R+'TableS9_sex_models.csv',index_col=0).rename(columns={'dataset':'Dataset','formula_age':'Model 1','beta_oa_age':'Model 1: OA coefficient (log2)','p_oa_age':'Model 1: P','std_effect_age':'Model 1: standardised effect','n_age':'Model 1: n','formula_age_sex':'Model 2','beta_oa_age_sex':'Model 2: OA coefficient (log2)','p_oa_age_sex':'Model 2: P','std_effect_age_sex':'Model 2: standardised effect','n_age_sex':'Model 2: n'})
sheets.append(('S9','Table S9. Ordinary least-squares models of MMP13 with group and, where annotated, age (Model 1), and with sex added (Model 2); the exact formula of each model is given per dataset. Sex from GEO annotation or inferred from XIST and Y-chromosome genes. Standardised effect, OA coefficient divided by the residual standard deviation.',x))
HF=lambda v: v.replace('Pperoxisome','Peroxisome').replace('heme Metabolism','Heme Metabolism').replace('  ',' ')
a=pd.read_csv(R+'TableS10a_cameraPR.csv',index_col=0); a['set']=a['set'].map(HF); a=a.rename(columns={'set':'Hallmark set','n_genes':'Genes (n)','inter_gene_cor':'Inter-gene correlation (OA samples)','direction':'Direction','t':'t','p':'P','fdr':'FDR (BH)'})
sheets.append(('S10a','Table S10a. Pre-ranked CAMERA (Python implementation of the limma::cameraPR formula) on consensus meta-Z, with inter-gene correlation estimated from OA samples.',a))
b=pd.read_csv(R+'TableS10b_sample_permutation.csv',index_col=0); b['set']=b['set'].map(HF); b=b.rename(columns={'set':'Hallmark set','mean_meta_z':'Mean meta-Z','null_mean':'Null mean','null_sd':'Null SD','p_perm':'Permutation P','fdr':'FDR (BH)'})
sheets.append(('S10b','Table S10b. Sample-label permutation test (MMP13 values permuted within cohort; 1,000 permutations; full meta-analysis rerun per permutation).',b))
summ=pd.DataFrame([
 ['Case–control (Hedges\' g, random effects, HKSJ)','Primary analysis (6 cohorts)','g = 0.69 (95% CI −0.39 to 1.77); P = 0.16; I² = 67%; τ² = 0.43'],
 ['Case–control','GSE57218 controls aged 15–18 years excluded (3 adult controls left)','GSE57218 g = 0.64 (−0.55 to 1.83); pooled g = 0.70 (−0.42 to 1.81); P = 0.17; I² = 67%; τ² = 0.48'],
 ['Case–control','GSE169077 excluded (5 cohorts)','g = 0.41 (−0.21 to 1.03); P = 0.14; I² = 28%'],
 ['Co-expression within OA','Full rerun (reference for the rows below)','1,513 positive, 407 negative robust partners (1,920 in total)'],
 ['Co-expression within OA','GSE169077 excluded (Table S11a)','1,383 positive, 680 negative robust partners; genome-wide r with full analysis = 0.98; 1,680 of the 1,920 robust partners retained'],
 ['Co-expression within OA','Spearman-specific variance 1.06/(n − 3) (Table S11b)','1,424 positive, 373 negative robust partners'],
 ['Co-expression within OA','Adjusted for age and sex','1,450 positive, 477 negative robust partners; genome-wide r with unadjusted analysis = 0.97; 150/150 top positive and 149/150 top negative partners retained (gene-level table not saved)'],
],columns=['Analysis','Sensitivity analysis','Result'])
sheets.append(('S11','Table S11. Summary of sensitivity analyses. Gene-level results are given in Tables S11a and S11b.',summ))
sheets.append(('S11a','Table S11a. MMP13 co-expression meta-analysis within OA tissue without GSE169077.',pd.read_csv(R+'TableS11a_coexpr_without_GSE169077.csv').rename(columns=COLS)))
sheets.append(('S11b','Table S11b. MMP13 co-expression meta-analysis within OA tissue with Spearman-specific variance 1.06/(n − 3).',pd.read_csv(R+'TableS11b_coexpr_fieller_variance.csv').rename(columns=COLS)))
par=pd.DataFrame([['Minimum genes per cell','200','200'],['Minimum cells per gene','3','3'],['Normalisation target sum','10,000','10,000'],['Highly variable genes','2,000 (Scanpy flavour "seurat")','2,000 (Scanpy flavour "seurat")'],['Principal components','50','50'],['Nearest neighbours','15','15'],['Leiden resolution','0.6','0.6'],['Random seed','7','7'],['Mitochondrial fraction filter','None','> 20% removed'],['Doublet removal','None','Scrublet'],['Harmony integration','None','Run on the sample label of the deposited matrix, which has a single level; donor effects could not be corrected'],['Cells retained','1,597','1,597'],['Clusters','9','9'],['Cells flagged by the mitochondrial filter or Scrublet','–','0'],['Adjusted Rand index vs primary','–','1.00']],columns=['Parameter','Primary analysis','Sensitivity analysis'])
sheets.append(('S12a','Table S12a. Single-cell analysis parameters (GSE104782).',par))
p=pd.read_csv(R+'Figure8e_source_cluster_MMP13.csv',index_col=0); s=pd.read_csv(R+'TableS12_sensitivity_cluster_MMP13.csv',index_col=0)
cl=p.merge(s,on='cluster',suffixes=(' (primary)',' (sensitivity)')).sort_values('cluster')
cl.columns=[col.replace('cluster','Leiden cluster').replace('n_cells','Cells (n)').replace('mean','Mean MMP13 (log-normalised)').replace('ci_low','95% CI low').replace('ci_high','95% CI high').replace('pct_detected','Cells with detectable MMP13 (%)') for col in cl.columns]
sheets.append(('S12b','Table S12b. Mean MMP13 expression per cluster (95% bootstrap CI over cells) in the primary and sensitivity analyses.',cl))
idx=pd.DataFrame([(n,t.split('. ',1)[1]) for n,t,_ in sheets],columns=['Sheet','Content'])
with pd.ExcelWriter(out,engine='openpyxl') as w:
    idx.to_excel(w,sheet_name='Index',index=False)
    for n,t,df in sheets: df.to_excel(w,sheet_name=n,index=False,startrow=1)
wb=load_workbook(out)
for n,t,df in sheets:
    ws=wb[n]; ws['A1']=t; ws['A1'].font=Font(bold=True)
    for cc in ws[2]: cc.font=Font(bold=True); cc.alignment=Alignment(wrap_text=True,vertical='top')
    for i,col in enumerate(ws.iter_cols(min_row=2,max_row=min(ws.max_row,200)),1):
        L=max(len(str(v.value)) if v.value is not None else 0 for v in col)
        ws.column_dimensions[col[0].column_letter].width=min(max(10,L*0.9),60)
    ws.freeze_panes='A3'
wb['Index'].column_dimensions['B'].width=120
wb.save(out); print('saved', [n for n,_,_ in sheets])
