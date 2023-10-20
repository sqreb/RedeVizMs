import pandas as pd
from scipy.sparse import coo_matrix, save_npz
import argparse
import sys
import os


def plot_gene_bin_expr_main(args):
    f_in = args.input
    x_label = args.x_label
    y_label = args.y_label
    f_out = args.output

    spot_df = pd.read_csv(f_in, sep="\t")
    pos_df = spot_df.groupby([x_label, y_label]).size().reset_index(name="GeneNum")
    x_range = spot_df[x_label].max()
    y_range = spot_df[y_label].max()

    sig_arr = coo_matrix(([True]*pos_df.shape[0], (pos_df[x_label].to_numpy(), pos_df[y_label].to_numpy())), (x_range+1, y_range+1))
    save_npz(f_out, sig_arr)

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=plot_gene_bin_expr_main)
    parser.add_argument("--input", type=str, dest="input",
                        metavar="spot.tsv", required=True)
    parser.add_argument("--x-label", type=str, dest="x_label", metavar="x_label", required=True)
    parser.add_argument("--y-label", type=str, dest="y_label", metavar="y_label", required=True)
    parser.add_argument("--output", type=str, dest="output",
                        metavar="output.npz", required=True, help="Output.npz")
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    args.func(args)

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()


