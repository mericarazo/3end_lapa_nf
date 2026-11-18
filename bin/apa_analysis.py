#!/usr/bin/env python3
"""
Analisis APA a partir de la salida de LAPA (un pool, dos condiciones).
Reproduce los notebooks lapa_processing.ipynb + _dPUR.ipynb:

  1. une counts + usage de LAPA con polyA_clusters.bed (Feature, gene_id)
  2. filtra sitios:  presente en todas las muestras de una condicion (opcional),
                     media de counts >= --min-mean en alguna condicion,
                     usage medio >= --min-usage en alguna condicion,
                     no intergenico, y genes con >= 2 sitios
  3. dPUR = usage medio (case) - usage medio (control)
  4. Fisher exacto por sitio (counts del sitio vs resto del gen, case vs control) + BH
  5. por gen: ancla (menor padj) + socio de signo opuesto -> shorter / longer, dDistalPUR
  6. volcano plot

Uso:
  apa_analysis.py --clusters polyA_clusters.bed --counts lapa_counts.tsv --usage lapa_usage.tsv \
      --gtf anotacion.gtf --case TDP43 --control NT --outprefix pool1
"""
import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

CLUSTER_COLS = ["Chromosome", "Start", "End", "polyA_site", "count", "Strand", "Feature", "gene_id",
                "tpm", "gene_count", "usage", "fracA", "signal", "annotated_site"]


def bh(p):
    """Benjamini-Hochberg (equivalente a statsmodels fdr_bh)."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def load_gene_names(gtf):
    ids, names = [], []
    rx_id = re.compile(r'gene_id "([^"]+)"')
    rx_nm = re.compile(r'gene_name "([^"]+)"')
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            mi, mn = rx_id.search(f[8]), rx_nm.search(f[8])
            if mi:
                ids.append(mi.group(1).split(".")[0])
                names.append(mn.group(1) if mn else None)
    return pd.DataFrame({"gene_id": ids, "gene_name": names}).drop_duplicates("gene_id")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--usage", required=True)
    ap.add_argument("--gtf", required=True, help="GTF solo para los nombres de gen")
    ap.add_argument("--case", required=True, help="condicion de interes (p. ej. TDP43)")
    ap.add_argument("--control", required=True, help="condicion control (p. ej. NT)")
    ap.add_argument("--min-mean", type=float, default=5.0, help="media de counts en alguna condicion (def. 5)")
    ap.add_argument("--min-usage", type=float, default=0.0025,
                    help="usage medio en alguna condicion, en fraccion (def. 0.0025 = 0.25 %%, lo usado en el notebook)")
    ap.add_argument("--require-detected-all", action="store_true",
                    help="exige counts > 0 en todas las muestras de una condicion (en el notebook este filtro no tenia efecto)")
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--delta", type=float, default=0.05)
    ap.add_argument("--highlight", default="STMN2", help="gen a resaltar en el volcano (vacio = ninguno)")
    ap.add_argument("--outprefix", required=True)
    a = ap.parse_args()

    # ---- 1. carga y union ------------------------------------------------------------------
    cl = pd.read_csv(a.clusters, sep="\t", header=None)
    cl.columns = CLUSTER_COLS
    cl["site"] = cl["Chromosome"].astype(str) + ":" + cl["polyA_site"].astype(str) + ":" + cl["Strand"]
    cl = cl.sort_values("usage", ascending=False).drop_duplicates(subset=["site"], keep="first")

    counts = pd.read_csv(a.counts, sep="\t").fillna(0)
    usage = pd.read_csv(a.usage, sep="\t").fillna(0)

    samples = [c for c in counts.columns if c != "site"]
    case_cols = [c for c in samples if a.case in c]
    ctrl_cols = [c for c in samples if a.control in c]
    if not case_cols or not ctrl_cols:
        sys.exit(f"No encuentro muestras *_{a.case} / *_{a.control} en las columnas: {samples}")
    print(f"case ({a.case}): {case_cols}\ncontrol ({a.control}): {ctrl_cols}", file=sys.stderr)
    case_u = [c + "_usage" for c in case_cols]
    ctrl_u = [c + "_usage" for c in ctrl_cols]

    merged = (counts
              .merge(cl[["site", "Feature", "gene_id", "polyA_site", "Strand"]], on="site", how="inner")
              .merge(usage, on="site", suffixes=("", "_usage")))

    # ---- 2. filtros ------------------------------------------------------------------------
    if a.require_detected_all:
        f_present = (merged[case_cols] > 0).all(axis=1) | (merged[ctrl_cols] > 0).all(axis=1)
    else:
        f_present = pd.Series(True, index=merged.index)
    f_mean = (merged[case_cols].mean(axis=1) >= a.min_mean) | (merged[ctrl_cols].mean(axis=1) >= a.min_mean)
    f_usage = (merged[case_u].mean(axis=1) >= a.min_usage) | (merged[ctrl_u].mean(axis=1) >= a.min_usage)
    f_feat = merged["Feature"] != "intergenic"

    steps = [("sitios LAPA (con anotacion)", len(merged)),
             ("presentes en todas las muestras de una condicion", int(f_present.sum())),
             (f"+ media >= {a.min_mean:g} en una condicion", int((f_present & f_mean).sum())),
             (f"+ usage >= {a.min_usage:g} en una condicion", int((f_present & f_mean & f_usage).sum())),
             ("+ no intergenico", int((f_present & f_mean & f_usage & f_feat).sum()))]
    filt = merged[f_present & f_mean & f_usage & f_feat]
    n_per_gene = filt["gene_id"].value_counts()
    filt = filt[filt["gene_id"].isin(n_per_gene[n_per_gene >= 2].index)].reset_index(drop=True)
    steps.append(("+ >= 2 sitios por gen", len(filt)))
    pd.DataFrame(steps, columns=["filtro", "sitios"]).to_csv(f"{a.outprefix}.filter_summary.tsv", sep="\t", index=False)
    for s, n in steps:
        print(f"{s}: {n}", file=sys.stderr)
    if filt.empty:
        sys.exit("Ningun sitio supera los filtros")

    # ---- 3. dPUR ---------------------------------------------------------------------------
    filt["dPUR"] = filt[case_u].mean(axis=1) - filt[ctrl_u].mean(axis=1)

    # ---- 4. Fisher por sitio ---------------------------------------------------------------
    sample_cols = ctrl_cols + case_cols
    gene_tot = filt.groupby("gene_id")[sample_cols].transform("sum")
    k = pd.DataFrame({"case": filt[case_cols].sum(axis=1), "ctrl": filt[ctrl_cols].sum(axis=1)})
    n = pd.DataFrame({"case": gene_tot[case_cols].sum(axis=1), "ctrl": gene_tot[ctrl_cols].sum(axis=1)})
    keep = ((n > 0).sum(axis=1) > 1) & (k != n).any(axis=1)     # igual que _k_n de LAPA
    k, n = k[keep], n[keep]

    odds, pval = [], []
    for kc, kt, nc, nt in zip(k["case"], k["ctrl"], n["case"], n["ctrl"]):
        o, p = fisher_exact([[kc, kt], [nc - kc, nt - kt]])
        odds.append(o)
        pval.append(p)
    res = pd.DataFrame({"odds_ratio": odds, "pval": pval}, index=k.index)
    res["padj"] = bh(res["pval"].values)

    final = filt.loc[keep[keep].index, ["site", "gene_id", "dPUR", "Strand"]].join(res)
    final = final.sort_values("padj")
    final["gene_id"] = final["gene_id"].str.split(".").str[0]
    final = final.merge(load_gene_names(a.gtf), on="gene_id", how="left")
    final["coord"] = final["site"].str.split(":").str[1].astype(int)
    final = final.rename(columns={"dPUR": "∆dPUR"})
    final.drop(columns="coord").to_csv(f"{a.outprefix}.fisher_events_individual.tsv", sep="\t", index=False)
    print(f"Sitios testados: {len(final)} | padj<{a.padj}: {(final['padj'] < a.padj).sum()}", file=sys.stderr)

    # ---- 5. pares proximal / distal por gen ------------------------------------------------
    anchor = final.loc[final.groupby("gene_id")["padj"].idxmin()].reset_index(drop=True)
    sg = final.merge(anchor[["gene_id", "∆dPUR"]].rename(columns={"∆dPUR": "anchor_delta"}), on="gene_id")
    opp = sg[np.sign(sg["∆dPUR"]) != np.sign(sg["anchor_delta"])]
    if opp.empty:
        sys.exit("Ningun gen tiene un socio de signo opuesto")
    partner = opp.loc[opp.groupby("gene_id")["padj"].idxmin()].reset_index(drop=True)
    pair = anchor.merge(partner, on="gene_id", suffixes=("_anchor", "_partner"))

    anchor_pos = pair["∆dPUR_anchor"] > 0
    up_site = np.where(anchor_pos, pair["site_anchor"], pair["site_partner"])
    up_delta = np.where(anchor_pos, pair["∆dPUR_anchor"], pair["∆dPUR_partner"])
    up_coord = np.where(anchor_pos, pair["coord_anchor"], pair["coord_partner"])
    up_strand = np.where(anchor_pos, pair["Strand_anchor"], pair["Strand_partner"])
    down_site = np.where(anchor_pos, pair["site_partner"], pair["site_anchor"])
    down_coord = np.where(anchor_pos, pair["coord_partner"], pair["coord_anchor"])

    # strand +: externo (distal) = coordenada mayor ; strand -: coordenada menor
    up_is_distal = np.where(up_strand == "+", up_coord > down_coord, up_coord < down_coord)
    delta_distal = np.where(up_is_distal, up_delta, -up_delta)
    direction = np.where(delta_distal > 0, "longer", "shorter")

    paired = pd.DataFrame({
        "gene_id": pair["gene_id"],
        "gene_name": pair["gene_name_anchor"],
        "APA_direction": direction,
        "site_down": down_site,
        "site_up": up_site,
        "∆distalPUR": delta_distal,
        "padj": pair["padj_anchor"],
    }).sort_values("padj").reset_index(drop=True)
    paired.to_csv(f"{a.outprefix}.fisher_events_paired.tsv", sep="\t", index=False)

    # ---- 6. volcano -------------------------------------------------------------------------
    r = paired.dropna(subset=["padj", "∆distalPUR"]).copy()
    r["nl"] = -np.log10(r["padj"].clip(lower=1e-300))
    r["up"] = (r["padj"] < a.padj) & (r["∆distalPUR"] > a.delta)
    r["down"] = (r["padj"] < a.padj) & (r["∆distalPUR"] < -a.delta)
    r["sig"] = r["up"] | r["down"]

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(r.loc[~r["sig"], "∆distalPUR"], r.loc[~r["sig"], "nl"], c="grey", s=10, alpha=0.5)
    ax.scatter(r.loc[r["up"], "∆distalPUR"], r.loc[r["up"], "nl"], c="red", s=12, alpha=0.7,
               label="Distal polyA site")
    ax.scatter(r.loc[r["down"], "∆distalPUR"], r.loc[r["down"], "nl"], c="blue", s=12, alpha=0.7,
               label="Proximal polyA site")
    sig = r[r["sig"]]
    top = pd.concat([sig.sort_values("padj").head(10), sig.sort_values("∆distalPUR").head(5)]).drop_duplicates()
    for _, row in top.iterrows():
        if pd.notna(row["gene_name"]):
            ax.annotate(row["gene_name"], (row["∆distalPUR"], row["nl"]), xytext=(3, 3),
                        textcoords="offset points", fontsize=7, alpha=0.8)
    if a.highlight:
        h = r[r["gene_name"].eq(a.highlight)]
        if len(h):
            ax.scatter(h["∆distalPUR"], h["nl"], c="black", s=45, marker="*", zorder=5, label=a.highlight)
            for _, row in h.iterrows():
                ax.annotate(a.highlight, (row["∆distalPUR"], row["nl"]), xytext=(5, 5),
                            textcoords="offset points", fontsize=9, fontweight="bold")
    ax.axhline(-np.log10(a.padj), color="black", linestyle="--", linewidth=0.8)
    ax.axvline(a.delta, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-a.delta, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel(f"∆distalPUR ({a.case} - {a.control})")
    ax.set_ylabel("-log10(padj)")
    ax.set_title("Volcano plot: differential polyA site usage")
    ax.legend(loc="upper right", frameon=False)
    plt.tight_layout()
    plt.savefig(f"{a.outprefix}.volcano.png", dpi=300)
    plt.savefig(f"{a.outprefix}.volcano.pdf")

    n_up, n_down = int(r["up"].sum()), int(r["down"].sum())
    print(f"Genes con par proximal/distal: {len(r)} | significativos: {int(r['sig'].sum())} "
          f"(distal {n_up}, proximal {n_down})", file=sys.stderr)


if __name__ == "__main__":
    main()
