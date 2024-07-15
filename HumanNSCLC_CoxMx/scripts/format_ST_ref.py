import sys
import argparse
import pandas as pd
import scanpy as sc


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--expr", type=str, dest="expr", metavar="expr.csv", required=True)
    base_group.add_argument("--meta", type=str, dest="meta", metavar="meta.csv", required=True)
    base_group.add_argument("--output", type=str, dest="output", metavar="output.h5ad", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_expr = args.expr
    f_meta = args.meta
    f_out = args.output

    expr_df = pd.read_csv(f_expr)
    meta_df = pd.read_csv(f_meta)

    both_cell_ID = expr_df[["fov", "cell_ID"]].merge(meta_df[["fov", "cell_ID"]], how="inner")
    expr_df = both_cell_ID.merge(expr_df, how="inner")
    meta_df = both_cell_ID.merge(meta_df, how="inner")

    expr_arr = expr_df[expr_df.columns[2:]].to_numpy()

    sce = sc.AnnData(X=expr_arr, obs=meta_df)
    sce.var_names = expr_df.columns[2:]
    sce.obs_names = [f"{fov}.{cell_ID}" for (fov, cell_ID) in zip(
        both_cell_ID["fov"].to_numpy(), both_cell_ID["cell_ID"].to_numpy()
    )]
    sce.write_h5ad(f_out)

def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()
