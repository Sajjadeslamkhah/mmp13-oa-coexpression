"""Revision analyses for the MMP13 co-expression manuscript (CIMB, round 1).

Every function works on a common data interface:

    DATA[name] = {"expr": DataFrame genes x samples (log2, one row per gene symbol),
                  "meta": DataFrame indexed by sample with columns
                          group   : 'OA' | 'control' | 'preserved'
                          age     : float or NaN
                          sex     : 'F' | 'M' | NaN
                          joint   : 'knee' | 'hip' | NaN
                          patient : patient id (used to pair GSE117999/GSE98918)
                          platform: platform id (GSE114007 has two)
                          tissue  : 'cartilage' | 'meniscus' | 'synovium' | 'bone'}
    COHORT[name] = cohort label (GSE117999 and GSE98918 share one label)

Use the objects from the ORIGINAL analysis notebook if you still have them; the
reproduction gates (`check_reproduction`) tell you whether the inputs match the
published numbers before any new result is produced.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

SEED = 20261113
rng_global = np.random.default_rng(SEED)

HYPERTROPHY = ["RUNX2", "COL10A1", "IHH", "ALPL", "VEGFA", "SPP1", "PTH1R"]
CATABOLIC = ["MMP13", "MMP3", "MMP1", "ADAMTS4", "ADAMTS5", "MMP9"]
ANABOLIC = ["COL2A1", "ACAN", "SOX9", "PRG4", "COMP", "TIMP3"]

# Cell-type marker genes (Table S7). Markers that are locked MMP13 partners are
# dropped automatically by `composition_scores(..., exclude=locked_genes)`.
MARKERS = {
    "chondrocyte": ["COL2A1", "ACAN", "SOX9", "COMP", "CHAD", "CILP", "FMOD", "CNMD"],
    "fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "FBLN1", "MFAP5"],
    "lining_fibroblast": ["PRG4", "HAS1", "CD55", "CLIC5", "HTRA4"],
    "osteoblast": ["SP7", "BGLAP", "IBSP", "COL1A1", "DMP1", "PHEX", "SOST"],
    "osteoclast": ["ACP5", "CTSK", "CALCR", "OSCAR", "DCSTAMP", "MMP9"],
    "macrophage": ["CD68", "CD163", "CSF1R", "C1QA", "C1QB", "MARCO", "MS4A7"],
    "T_cell": ["CD3D", "CD3E", "CD2", "CD247", "TRAC"],
    "B_plasma": ["MS4A1", "CD79A", "MZB1", "JCHAIN", "IGHM"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "EMCN"],
    "pericyte_SMC": ["ACTA2", "MYH11", "RGS5", "NOTCH3", "CNN1"],
    "adipocyte": ["ADIPOQ", "PLIN1", "LEP", "FABP4", "CIDEA"],
}
Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "UTY"]


# --------------------------------------------------------------------- utils
def zrow(df):
    """z-score each gene across samples (within dataset)."""
    m = df.mean(axis=1); s = df.std(axis=1, ddof=1).replace(0, np.nan)
    return df.sub(m, axis=0).div(s, axis=0)


def _rank_std(X):
    """X: ndarray genes x n -> per-row ranks standardised to mean 0, norm 1."""
    R = np.apply_along_axis(stats.rankdata, 1, X)
    R = R - R.mean(axis=1, keepdims=True)
    nrm = np.linalg.norm(R, axis=1, keepdims=True)
    nrm[nrm == 0] = np.nan
    return R / nrm


def _residualise(R, C):
    """Regress covariates C (n x k) out of each row of R (genes x n)."""
    if C is None or C.shape[1] == 0:
        return R
    C1 = np.column_stack([np.ones(C.shape[0]), C])
    beta, *_ = np.linalg.lstsq(C1, R.T, rcond=None)
    E = R - (C1 @ beta).T
    E = E - E.mean(axis=1, keepdims=True)
    nrm = np.linalg.norm(E, axis=1, keepdims=True); nrm[nrm == 0] = np.nan
    return E / nrm


def select_samples(meta, which="OA", min_oa=8):
    """Replicates the published rule: OA samples, or all samples if < min_oa OA."""
    if which == "OA":
        oa = meta.index[meta.group == "OA"]
        if len(oa) < min_oa:
            return list(meta.index[meta.group.isin(["OA", "control"])])
        return list(oa)
    if which == "control":
        return list(meta.index[meta.group == "control"])
    raise ValueError(which)


# --------------------------------------------------- co-expression meta core
def cohort_z(DATA, COHORT, seeds, which="OA", min_oa=8, covariates=None,
             seed_vectors=None):
    """Fisher z of Spearman(seed, gene) per cohort, exactly as module C of MMP13-OA v5:
    genes with zero variance in the selected samples are dropped per dataset, rho clipped
    to +/-0.999, z averaged within cohort, cohort n = mean sample count of its datasets."""
    names = list(seed_vectors.keys()) if seed_vectors is not None else list(seeds)
    per_cohort = {s: {} for s in names}
    n_cohort = {}
    for ds, d in DATA.items():
        meta = d["meta"]; expr = d["expr"]
        smp = select_samples(meta, which, min_oa)
        if len(smp) < 4:
            continue
        X = expr[smp]
        X = X[X.notna().all(axis=1) & (X.std(axis=1) > 0)]
        G = _rank_std(X.values)
        C = covariates(ds, meta.loc[smp]) if covariates else None
        if C is not None:
            G = _residualise(G, C)
        if seed_vectors is not None:
            S = np.vstack([seed_vectors[s][ds] for s in names])
            S = _rank_std(S)
            if C is not None: S = _residualise(S, C)
        else:
            idx = [X.index.get_loc(s) if s in X.index else None for s in names]
            S = np.vstack([G[i] if i is not None else np.full(G.shape[1], np.nan) for i in idx])
        Rm = np.clip(S @ G.T, -0.999, 0.999)                 # seeds x genes
        Z = np.arctanh(Rm)
        c = COHORT[ds]
        for k, s in enumerate(names):
            per_cohort[s].setdefault(c, []).append(pd.Series(Z[k], index=X.index))
        n_cohort.setdefault(c, []).append(len(smp))
    out = {}
    for s in names:
        cols = {c: pd.concat(lst, axis=1).mean(axis=1) for c, lst in per_cohort[s].items()}
        out[s] = pd.DataFrame(cols)
    n = pd.Series({c: float(np.mean(v)) for c, v in n_cohort.items()})
    return out, n


def meta_from_z(Z, n, var_factor=1.0, sign_frac=0.75, **_ignored):
    """Cohort-level weighted Stouffer as in MMP13-OA v5 module C:
    meta_rho = tanh(sum(z w) / sum(w)), meta_z = sum(z w) / sqrt(sum(w)), w = n - 3;
    genes need >= max(2, ceil(K/2)) cohorts; robust = FDR < 0.05 and >= 75% sign agreement.
    var_factor = 1.06 applies the Fieller variance 1.06/(n - 3) (meta_z / sqrt(1.06))."""
    Z = Z[[c for c in Z.columns if c in n.index]]
    W = n[Z.columns] - 3
    av = Z.notna()
    num = (Z.fillna(0) * W).sum(axis=1)
    res = pd.DataFrame({"meta_rho": np.tanh(num / (av * W).sum(axis=1)),
                        "meta_z": num / np.sqrt((av * W).sum(axis=1)) / np.sqrt(var_factor),
                        "n_cohorts": av.sum(axis=1)})
    res["same_sign_frac"] = (np.sign(Z).eq(np.sign(res.meta_rho), axis=0) & av).sum(axis=1) / res.n_cohorts
    res["p"] = 2 * stats.norm.sf(res.meta_z.abs())
    res = res[res.n_cohorts >= max(2, int(np.ceil(Z.shape[1] / 2)))].copy()
    res["fdr"] = multipletests(res.p, method="fdr_bh")[1]
    res["robust"] = (res.fdr < 0.05) & (res.same_sign_frac >= sign_frac)
    return res


def coexpr(DATA, COHORT, seed="MMP13", **kw):
    var_factor = kw.pop("var_factor", 1.0)
    Zd, n = cohort_z(DATA, COHORT, [seed], **kw)
    res = meta_from_z(Zd[seed].drop(index=seed, errors="ignore"), n, var_factor=var_factor)
    return res, n


def summarise(res):
    pos = res[res.robust & (res.meta_z > 0)]; neg = res[res.robust & (res.meta_z < 0)]
    return {"n_tested": len(res), "n_pos": len(pos), "n_neg": len(neg)}


def top_partners(res, k=150):
    r = res[res.robust]
    pos = r[r.meta_z > 0].sort_values("meta_z", ascending=False).index[:k]
    neg = r[r.meta_z < 0].sort_values("meta_z").index[:k]
    return list(pos), list(neg)


def check_reproduction(res, published=None, expect=(1510, 407)):
    """Compare with the published counts and, if given, with the saved T3 table (MMP13-OA v5)."""
    sm = summarise(res)
    ok = (sm["n_pos"], sm["n_neg"]) == tuple(expect)
    msg = f"robust partners: {sm['n_pos']} positive / {sm['n_neg']} negative (published {expect[0]} / {expect[1]}); genes tested {sm['n_tested']}"
    if published is not None:
        common = res.index.intersection(published.index)
        r = np.corrcoef(res.loc[common, "meta_z"], published.loc[common, "meta_z"])[0, 1]
        msg += f"; r(meta_z) with saved T3 = {r:.4f} over {len(common)} genes"
        ok = ok and r > 0.999
    print("Reproduction gate:", msg, "->", "OK" if ok else "MISMATCH - do not use new numbers until fixed")
    return ok


# ------------------------------------------------ R1-1: seeds and centrality
def matched_random_genes(DATA, COHORT, target="MMP13", k=1000, exclude=(), seed=SEED):
    """Genes matched to `target` on mean and SD percentile, averaged over datasets."""
    pm, ps = [], []
    for ds, d in DATA.items():
        smp = select_samples(d["meta"], "OA")
        X = d["expr"][smp]
        pm.append(X.mean(axis=1).rank(pct=True)); ps.append(X.std(axis=1).rank(pct=True))
    PM = pd.concat(pm, axis=1); PS = pd.concat(ps, axis=1)
    full = PM.notna().all(axis=1)
    PM, PS = PM[full].mean(axis=1), PS[full].mean(axis=1)
    d = np.sqrt((PM - PM[target]) ** 2 + (PS - PS[target]) ** 2)
    d = d.drop(index=[target, *exclude], errors="ignore").sort_values()
    pool = list(d.index[: max(5 * k, k)])
    rng = np.random.default_rng(seed)
    return list(rng.choice(pool, size=k, replace=False))


def validation_corr(VAL, seed, center_platform=False, group="OA"):
    """Spearman(seed, gene) within one group (default OA) of each validation dataset."""
    out = {}
    for ds, d in VAL.items():
        meta = d["meta"]; smp = list(meta.index[meta.group == group])
        X = d["expr"][smp]
        X = X[X.notna().all(axis=1)]
        if center_platform and meta.loc[smp, "platform"].nunique() > 1:
            pl = meta.loc[smp, "platform"]
            X = X.T.groupby(pl.values).transform(lambda v: v - v.mean()).T
        if seed not in X.index:
            continue
        G = _rank_std(X.values); s = G[X.index.get_loc(seed)]
        out[ds] = pd.Series(G @ s, index=X.index).drop(index=seed)
    return out


def seed_comparison(DATA, COHORT, VAL, named=("MMP13", "OLFML2B", "THY1", "POSTN", "FAP", "MMP3"),
                    n_random=1000, k=150):
    """Neighbourhood size and replication for named seeds and matched random genes."""
    rand = matched_random_genes(DATA, COHORT, "MMP13", n_random, exclude=named)
    seeds = list(named) + rand
    Zd, n = cohort_z(DATA, COHORT, seeds)
    rows = []
    vstats = {}
    for ds, d in VAL.items():   # precompute validation ranks once
        meta = d["meta"]; smp = list(meta.index[meta.group == "OA"])
        X = d["expr"][smp]; X = X[X.notna().all(axis=1)]
        vstats[ds] = (X.index, _rank_std(X.values))
    for s in seeds:
        res = meta_from_z(Zd[s].drop(index=s, errors="ignore"), n)
        sm = summarise(res); pos, neg = top_partners(res, k)
        row = {"seed": s, "type": "named" if s in named else "random", **sm}
        for ds, (idx, G) in vstats.items():
            if s not in idx:
                row[f"{ds}_pos_frac"] = np.nan; row[f"{ds}_genome_r"] = np.nan; continue
            v = pd.Series(G @ G[idx.get_loc(s)], index=idx).drop(index=s)
            pp = [g for g in pos if g in v.index]; nn = [g for g in neg if g in v.index]
            row[f"{ds}_pos_frac"] = float((v[pp] > 0).mean()) if pp else np.nan
            row[f"{ds}_neg_frac"] = float((v[nn] < 0).mean()) if nn else np.nan
            common = res.index.intersection(v.index)
            row[f"{ds}_genome_r"] = stats.spearmanr(res.loc[common, "meta_rho"], v[common])[0]
        rows.append(row)
    T = pd.DataFrame(rows).set_index("seed")
    T["n_robust"] = T.n_pos + T.n_neg
    rnd = T[T.type == "random"]
    pct = {c: float((rnd[c] < T.loc["MMP13", c]).mean() * 100)
           for c in T.columns if c not in ("type",) and T[c].dtype.kind in "fi"}
    return T, pct


def network_centrality(DATA, COHORT, n_top=5000, beta=6, force=("MMP13",)):
    """WGCNA-style weighted degree and eigenvector centrality on meta-correlations."""
    var = []
    for ds, d in DATA.items():
        smp = select_samples(d["meta"], "OA")
        var.append(d["expr"][smp].var(axis=1).rank(pct=True))
    V = pd.concat(var, axis=1)
    V = V[V.notna().all(axis=1)].mean(axis=1).sort_values(ascending=False)
    genes = list(dict.fromkeys(list(V.index[:n_top]) + [g for g in force if g in V.index]))
    num = np.zeros((len(genes), len(genes))); den = 0.0
    cohort_mats = {}
    for ds, d in DATA.items():
        smp = select_samples(d["meta"], "OA")
        G = np.nan_to_num(_rank_std(d["expr"].loc[genes, smp].values))
        cohort_mats.setdefault(COHORT[ds], []).append((np.arctanh(np.clip(G @ G.T, -0.999, 0.999)), len(smp)))
    for c, lst in cohort_mats.items():
        Zc = np.mean([m for m, _ in lst], axis=0); nc = np.mean([nn for _, nn in lst])
        num += (nc - 3) * Zc; den += (nc - 3)
    Rmeta = np.tanh(num / den)
    A = np.abs(Rmeta) ** beta; np.fill_diagonal(A, 0)
    k = A.sum(axis=1)
    vals, vecs = np.linalg.eigh(A)
    ev = np.abs(vecs[:, -1])
    T = pd.DataFrame({"weighted_degree": k, "eigenvector": ev}, index=genes)
    T["degree_rank"] = T.weighted_degree.rank(ascending=False).astype(int)
    T["eigen_rank"] = T.eigenvector.rank(ascending=False).astype(int)
    return T.sort_values("degree_rank")


# ------------------------------------------------- R1-2: cell-type composition
def composition_scores(DATA, exclude=(), min_markers=2):
    out = {}
    for ds, d in DATA.items():
        Z = zrow(d["expr"])
        sc = {}
        for ct, mk in MARKERS.items():
            g = [m for m in mk if m in Z.index and m not in exclude]
            if len(g) >= min_markers:
                sc[ct] = Z.loc[g].mean(axis=0)
        out[ds] = pd.DataFrame(sc)
    return out


def composition_covariates(scores, n_pc=2):
    """Returns a covariate function for cohort_z: first n_pc PCs of composition scores."""
    def f(ds, meta_sel):
        S = scores[ds].loc[meta_sel.index].dropna(axis=1)
        S = (S - S.mean()) / S.std(ddof=1).replace(0, np.nan)
        S = S.dropna(axis=1)
        U, s, Vt = np.linalg.svd(S.values, full_matrices=False)
        k = min(n_pc, max(0, S.shape[0] - 5))
        return U[:, :k] * s[:k]
    return f


def mmp13_vs_composition(DATA, scores):
    rows = []
    for ds, d in DATA.items():
        smp = select_samples(d["meta"], "OA")
        if "MMP13" not in d["expr"].index: continue
        m = d["expr"].loc["MMP13", smp]
        for ct in scores[ds].columns:
            r, p = stats.spearmanr(m, scores[ds].loc[smp, ct])
            rows.append({"dataset": ds, "tissue": d["meta"].tissue.iloc[0], "cell_type": ct,
                         "rho": r, "p": p, "n": len(smp)})
    return pd.DataFrame(rows)


def compare_meta(res_a, res_b, locked_pos, locked_neg):
    common = res_a.index.intersection(res_b.index)
    r = stats.spearmanr(res_a.loc[common, "meta_rho"], res_b.loc[common, "meta_rho"])[0]
    keep_pos = sum((g in res_b.index) and res_b.at[g, "fdr"] < 0.05 and res_b.at[g, "meta_z"] > 0 for g in locked_pos)
    keep_neg = sum((g in res_b.index) and res_b.at[g, "fdr"] < 0.05 and res_b.at[g, "meta_z"] < 0 for g in locked_neg)
    return {"genome_r": r, "pos_kept": keep_pos, "pos_total": len(locked_pos),
            "neg_kept": keep_neg, "neg_total": len(locked_neg), **summarise(res_b)}


def nnls_deconvolution(bulk_log2, ref_linear, genes=None):
    """Proportions of single-cell reference clusters in bulk samples (cartilage)."""
    from scipy.optimize import nnls
    B = 2 ** bulk_log2
    g = [x for x in (genes or ref_linear.index) if x in B.index and x in ref_linear.index]
    Rf = ref_linear.loc[g].values; Rf = Rf / Rf.sum(axis=0, keepdims=True)
    P = {}
    for s in B.columns:
        y = B.loc[g, s].values; y = y / y.sum()
        w, _ = nnls(Rf, y); P[s] = w / w.sum() if w.sum() > 0 else w
    return pd.DataFrame(P, index=ref_linear.columns).T


# ---------------------------------------------------- R1-3: hypertrophy by tissue
def hypertrophy_by_tissue(DATA, genes=HYPERTROPHY):
    rows, pooled = [], []
    for ds, d in DATA.items():
        smp = list(d["meta"].index[d["meta"].group == "OA"])
        if len(smp) < 8: continue       # as in Figure 4b (n = 117)
        Z = zrow(d["expr"][smp])
        g = [x for x in genes if x in Z.index and x != "MMP13"]
        score = Z.loc[g].mean(axis=0); m = Z.loc["MMP13"]
        r, p = stats.spearmanr(m, score)
        t = d["meta"].tissue.iloc[0]
        rows.append({"dataset": ds, "tissue": t, "n": len(smp), "rho": r, "p": p})
        pooled.append(pd.DataFrame({"tissue": t, "mmp13": m, "score": score}))
    P = pd.concat(pooled)
    for t, sub in list(P.groupby("tissue")) + [("all", P)]:
        r, p = stats.spearmanr(sub.mmp13, sub.score)
        lo, hi = np.tanh(np.arctanh(r) + np.array([-1, 1]) * 1.96 * np.sqrt(1.06 / (len(sub) - 3)))
        rows.append({"dataset": "pooled", "tissue": t, "n": len(sub), "rho": r, "p": p,
                     "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows), P


# --------------------------------------------- R1-4: control-tissue co-expression
def oa_vs_control_difference(DATA, COHORT, seed="MMP13"):
    Zoa, noa = cohort_z(DATA, COHORT, [seed], which="OA", min_oa=0)
    Zc, nc = cohort_z(DATA, COHORT, [seed], which="control")
    A, B = Zoa[seed], Zc[seed]
    cohorts = [c for c in A.columns if c in B.columns]
    D = pd.DataFrame({c: (A[c] - B[c]) / np.sqrt(1 / (noa[c] - 3) + 1 / (nc[c] - 3)) for c in cohorts})
    w = pd.Series({c: 1 / (1 / (noa[c] - 3) + 1 / (nc[c] - 3)) for c in cohorts})
    meas = D.notna()
    Zcomb = (D.fillna(0) * meas.mul(np.sqrt(w), axis=1)).sum(axis=1) / np.sqrt(meas.mul(w, axis=1).sum(axis=1))
    p = 2 * stats.norm.sf(np.abs(Zcomb))
    keep = meas.sum(axis=1) >= np.ceil(0.5 * len(cohorts))
    out = pd.DataFrame({"z_diff": Zcomb[keep], "p": p[keep.values]})
    out["fdr"] = multipletests(out.p, method="fdr_bh")[1]
    return out


# ------------------------------------------------------------- R1-6: sex
def infer_sex(expr):
    Z = zrow(expr)
    x = Z.loc["XIST"] if "XIST" in Z.index else pd.Series(0, index=Z.columns)
    y = [g for g in Y_GENES if g in Z.index]
    ys = Z.loc[y].mean(axis=0) if y else pd.Series(0, index=Z.columns)
    s = ys - x
    return pd.Series(np.where(s > 0, "M", "F"), index=expr.columns), s


def case_control_models(DATA, covars=("age",), extra=()):
    import statsmodels.formula.api as smf
    rows = []
    for ds, d in DATA.items():
        meta = d["meta"]; smp = list(meta.index[meta.group.isin(["OA", "control"])])
        df = meta.loc[smp].copy(); df["y"] = d["expr"].loc["MMP13", smp].values
        df["oa"] = (df.group == "OA").astype(int)
        cv = [c for c in (*covars, *extra) if c in df and df[c].notna().all() and df[c].nunique() > 1]
        f = "y ~ oa" + "".join(f" + C({c})" if df[c].dtype == object else f" + {c}" for c in cv)
        fit = smf.ols(f, df).fit()
        rows.append({"dataset": ds, "formula": f, "beta_oa": fit.params["oa"], "p_oa": fit.pvalues["oa"],
                     "std_effect": fit.params["oa"] / np.sqrt(fit.scale), "n": len(df)})
    return pd.DataFrame(rows)


# ------------------------------------------------- case-control meta-analysis
def hedges_g(a, b):
    n1, n2 = len(a), len(b)
    sp = np.sqrt(((n1 - 1) * np.var(a, ddof=1) + (n2 - 1) * np.var(b, ddof=1)) / (n1 + n2 - 2))
    d = (np.mean(a) - np.mean(b)) / sp
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    g = J * d
    v = (n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2))
    return g, v


def dl_hksj(g, v):
    g, v = np.asarray(g), np.asarray(v); k = len(g)
    w = 1 / v; mu_f = np.sum(w * g) / np.sum(w)
    Q = np.sum(w * (g - mu_f) ** 2)
    tau2 = max(0, (Q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
    ws = 1 / (v + tau2); mu = np.sum(ws * g) / np.sum(ws)
    se = np.sqrt(np.sum(ws * (g - mu) ** 2) / ((k - 1) * np.sum(ws)))
    t = stats.t.ppf(0.975, k - 1)
    p = 2 * stats.t.sf(abs(mu / se), k - 1)
    I2 = max(0, (Q - (k - 1)) / Q) if Q > 0 else 0
    return {"g": mu, "ci_low": mu - t * se, "ci_high": mu + t * se, "p": p, "tau2": tau2, "I2": I2, "k": k}


def case_control_meta(DATA, COHORT, exclude=(), sample_filter=None):
    per = []
    for ds, d in DATA.items():
        if ds in exclude: continue
        meta = d["meta"]
        if sample_filter is not None: meta = sample_filter(ds, meta)
        a = d["expr"].loc["MMP13", meta.index[meta.group == "OA"]].values
        b = d["expr"].loc["MMP13", meta.index[meta.group == "control"]].values
        g, v = hedges_g(a, b)
        per.append({"dataset": ds, "cohort": COHORT[ds], "g": g, "v": v, "n_oa": len(a), "n_ctrl": len(b)})
    P = pd.DataFrame(per)
    C = P.groupby("cohort").agg(g=("g", "mean"), v=("v", lambda x: np.mean(np.sqrt(x)) ** 2))  # r = 1
    return P, dl_hksj(C.g.values, C.v.values)


# -------------------------------------------- R1-7: correlation-aware gene sets
def inter_gene_correlation(DATA, genes, which="OA"):
    vals, wts = [], []
    for ds, d in DATA.items():
        smp = select_samples(d["meta"], which)
        g = [x for x in genes if x in d["expr"].index]
        if len(g) < 3: continue
        X = d["expr"].loc[g, smp].dropna()
        # remove the linear MMP13 effect so the estimate reflects residual correlation
        Xc = X.values - X.values.mean(1, keepdims=True)
        if "MMP13" in d["expr"].index:
            m = d["expr"].loc["MMP13", smp].values.astype(float); m = (m - m.mean()) / m.std()
            Xc = Xc - np.outer(Xc @ m / len(m), m)
        C = np.corrcoef(Xc); iu = np.triu_indices_from(C, 1)
        vals.append(np.nanmean(C[iu])); wts.append(len(smp))
    return float(np.average(vals, weights=wts)) if vals else 0.0


def camera_pr(stat, sets, igc=None, DATA=None):
    """Pre-ranked CAMERA (Wu & Smyth 2012, parametric version of limma::cameraPR)."""
    stat = stat.dropna(); G = len(stat); rows = []
    for name, genes in sets.items():
        idx = stat.index.intersection(genes); m = len(idx)
        if m < 10 or m > 500: continue
        rho = igc[name] if isinstance(igc, dict) else (inter_gene_correlation(DATA, idx) if DATA is not None else (igc or 0.01))
        vif = 1 + (m - 1) * rho
        mean_set = stat[idx].mean(); mean_rest = stat.drop(idx).mean()
        s2 = stat.var(ddof=1)
        t = (mean_set - mean_rest) / np.sqrt(s2 * (vif / m + 1 / (G - m)))
        p = 2 * stats.t.sf(abs(t), G - 2)
        rows.append({"set": name, "n_genes": m, "inter_gene_cor": rho, "direction": "Up" if t > 0 else "Down", "t": t, "p": p})
    T = pd.DataFrame(rows); T["fdr"] = multipletests(T.p, method="fdr_bh")[1]
    return T.sort_values("p")


def permutation_gene_sets(DATA, COHORT, sets, n_perm=1000, seed=SEED, batch=100):
    """Sample-label permutation: permute MMP13 within each cohort (patients kept aligned)."""
    rng = np.random.default_rng(seed)
    base, n = coexpr(DATA, COHORT, "MMP13")
    obs = {k: base.meta_z.reindex(v).dropna().mean() for k, v in sets.items()}
    null = {k: [] for k in sets}
    # MMP13 vectors on the selected samples
    vec, patients = {}, {}
    for ds, d in DATA.items():
        smp = select_samples(d["meta"], "OA")
        X = d["expr"][smp]; X = X[X.notna().all(axis=1)]
        vec[ds] = X.loc["MMP13"].values
        patients[ds] = d["meta"].loc[smp, "patient"].astype(str).values
    for b0 in range(0, n_perm, batch):
        sv = {}
        for b in range(b0, min(n_perm, b0 + batch)):
            perm_by_cohort = {}
            sv[f"perm{b}"] = {}
            for ds in DATA:
                c = COHORT[ds]
                if c not in perm_by_cohort:
                    pts = np.unique(patients[ds]); perm_by_cohort[c] = dict(zip(pts, rng.permutation(pts)))
                mp = perm_by_cohort[c]
                order = {p: i for i, p in enumerate(patients[ds])}
                src = [order.get(mp.get(p, p), i) for i, p in enumerate(patients[ds])]
                sv[f"perm{b}"][ds] = vec[ds][src]
        Zd, nn = cohort_z(DATA, COHORT, None, seed_vectors=sv)
        for s in sv:
            mz = meta_from_z(Zd[s], nn).meta_z
            for k, v in sets.items():
                null[k].append(mz.reindex(v).dropna().mean())
    rows = []
    for k in sets:
        nl = np.array(null[k]); o = obs[k]
        if np.isnan(o) or len(nl) == 0: continue
        p = (1 + np.sum(np.abs(nl - nl.mean()) >= abs(o - nl.mean()))) / (1 + len(nl))
        rows.append({"set": k, "mean_meta_z": o, "null_mean": nl.mean(), "null_sd": nl.std(), "p_perm": p})
    T = pd.DataFrame(rows); T["fdr"] = multipletests(T.p_perm, method="fdr_bh")[1]
    return T.sort_values("p_perm")


# ---------------------------------------------------------- single-cell
def sc_pipeline(adata, n_hvg=2000, flavor="seurat", n_pcs=50, n_neighbors=15, resolution=0.6,
                mito_max=None, doublets=False, harmony_key=None, seed=7):
    """Defaults reproduce MMP13-OA validation module V3: min_genes 200, min_cells 3, normalize_total 1e4 + log1p
    (only if counts), 2,000 HVGs (seurat flavour), scale (max 10), PCA (50), neighbours (15),
    UMAP, Leiden (resolution 0.6, scanpy default flavour), seed 7."""
    import scanpy as sc
    ad = adata.copy()
    if mito_max is not None:
        ad.var["mt"] = ad.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        ad = ad[ad.obs.pct_counts_mt < mito_max].copy()
    X = ad.X.max()
    is_counts = float(X) > 50
    sc.pp.filter_cells(ad, min_genes=200); sc.pp.filter_genes(ad, min_cells=3)
    if doublets and is_counts:
        sc.pp.scrublet(ad, random_state=seed)
        ad = ad[~ad.obs.predicted_doublet].copy()
    if is_counts:
        sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    elif float(X) > 20:
        sc.pp.log1p(ad)
    ad.raw = ad
    sc.pp.highly_variable_genes(ad, n_top_genes=n_hvg, flavor=flavor)
    adh = ad[:, ad.var.highly_variable].copy(); sc.pp.scale(adh, max_value=10)
    sc.tl.pca(adh, n_comps=n_pcs, random_state=seed)
    rep = "X_pca"
    if harmony_key:
        import harmonypy
        ho = harmonypy.run_harmony(adh.obsm["X_pca"], adh.obs, [harmony_key], random_state=seed)
        Zc = np.asarray(ho.Z_corr)
        adh.obsm["X_pca_harmony"] = Zc if Zc.shape[0] == adh.n_obs else Zc.T
        rep = "X_pca_harmony"
    sc.pp.neighbors(adh, n_neighbors=n_neighbors, use_rep=rep, random_state=seed)
    sc.tl.umap(adh, random_state=seed); sc.tl.leiden(adh, resolution=resolution, random_state=seed)
    ad.obs["leiden"] = adh.obs["leiden"]; ad.obsm["X_umap"] = adh.obsm["X_umap"]
    ad.uns["params"] = dict(min_genes=200, min_cells=3, target_sum=1e4, n_hvg=n_hvg, flavor=flavor, n_pcs=n_pcs,
                            n_neighbors=n_neighbors, resolution=resolution, mito_max=mito_max, doublets=doublets,
                            harmony_key=harmony_key, seed=seed, input_was_counts=is_counts)
    return ad


def cluster_mmp13(ad, gene="MMP13", n_boot=1000, seed=0):
    rng = np.random.default_rng(seed)
    x = ad.raw[:, gene].X
    x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    rows = []
    for c, idx in ad.obs.groupby("leiden").indices.items():
        v = x[idx]
        boots = [rng.choice(v, len(v)).mean() for _ in range(n_boot)]
        rows.append({"cluster": c, "n_cells": len(v), "mean": v.mean(), "ci_low": np.percentile(boots, 2.5),
                     "ci_high": np.percentile(boots, 97.5), "pct_detected": 100 * (v > 0).mean()})
    return pd.DataFrame(rows).sort_values("mean", ascending=False)


# --------------------------------------------------------- PRISMA GEO search

