import sys
import argparse
import pandas as pd
import scanpy as sc
from redeviz.enhance.index import spot2sce

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--input", type=str, dest="input", metavar="spot.tsv", required=True)
    base_group.add_argument("--mask", type=str, dest="mask", metavar="mask.npz", required=True)
    base_group.add_argument("--output", type=str, dest="output", metavar="gene_li.txt", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_mask = args.mask
    f_out = args.output

    # f_in = "/home/wangdehe/HPC/project/RedeViz_ms/MouseEmbr/data/organ/E16.5_E1S1.lung.tsv"
    # f_mask = "/home/wangdehe/HPC/project/RedeViz_ms/MouseEmbr/analysis/organ_mask_ST_data/E16.5/E16.5_E1S1/lung/mask.npz"

    sce = spot2sce(f_in, f_mask, "x", "y", "MIDCounts", "geneID")
    sc.pp.normalize_total(sce, target_sum=1e4)
    sc.pp.log1p(sce)
    sc.pp.highly_variable_genes(sce, min_mean=0.01, max_mean=3, min_disp=0.5)
    HVG = sce.var_names[sce.var.highly_variable]

    spot_df = pd.read_csv(f_in, sep="\t")
    uniq_pos = spot_df[["x", "y"]].drop_duplicates()
    uniq_gene_pos = spot_df[["x", "y", "geneID"]].drop_duplicates().groupby(["geneID"]).size().reset_index(name="SpotNum")
    uniq_gene_pos["ExprRatio"] = uniq_gene_pos["SpotNum"] / uniq_pos.shape[0]
    uniq_gene_pos = uniq_gene_pos.sort_values(by=["SpotNum"], ascending=False)
    expr_gene_li = uniq_gene_pos["geneID"].to_numpy()[uniq_gene_pos["ExprRatio"]>1e-4]

    select_genes = set(expr_gene_li) & set(HVG)
    with open(f_out, "w") as f:
        for gid in select_genes:
            f.write(f"{gid}\n")

def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()
