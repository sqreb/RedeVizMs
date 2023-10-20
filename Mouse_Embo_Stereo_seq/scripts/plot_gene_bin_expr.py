import numpy as np
import pandas as pd
import cv2 as cv
from matplotlib import pyplot as plt
from scipy.sparse import coo_matrix
import argparse
import sys
import os



def mat2image_arr(gid, spot_df, nbg_label_arr, x_range, y_range, gene_name_label, UMI_label, x_label, y_label):
    if gid is None:
        return None
    tmp_spot_df = spot_df[spot_df[gene_name_label]==gid]
    expr_arr = coo_matrix((tmp_spot_df[UMI_label].to_numpy(), (tmp_spot_df[x_label].to_numpy(), tmp_spot_df[y_label].to_numpy())), (x_range+1, y_range+1)).toarray()
    nzo_arr = expr_arr[expr_arr>0]
    cutoff = np.percentile(nzo_arr, 95)
    expr_arr[expr_arr>cutoff] = cutoff
    norm_mat = expr_arr / cutoff
    expr_arr = norm_mat * 215
    expr_arr[~nbg_label_arr] = 0
    return expr_arr

def single_panel2multi_panel(x_range, y_range, nbg_label_arr, r_arr=None, g_arr=None, b_arr=None):
    if r_arr is None:
        r_arr = np.zeros((x_range+1, y_range+1), dtype=np.float64)
    if g_arr is None:
        g_arr = np.zeros((x_range+1, y_range+1), dtype=np.float64)
    if b_arr is None:
        b_arr = np.zeros((x_range+1, y_range+1), dtype=np.float64)
    r_arr = r_arr + 40
    g_arr = g_arr + 40
    b_arr = b_arr + 40
    r_arr[~nbg_label_arr] = 0
    g_arr[~nbg_label_arr] = 0
    b_arr[~nbg_label_arr] = 0
    rgb_arr = np.stack([r_arr, g_arr, b_arr], -1)
    rgb_arr = rgb_arr.astype(np.uint8)
    return rgb_arr

def plot_gene_bin_expr_main(args):
    f_in = args.input
    f_gene_list = args.gene_list
    x_label = args.x_label
    y_label = args.y_label
    gene_name_label = args.gene_name_label
    UMI_label = args.UMI_label
    f_out = args.output

    spot_df = pd.read_csv(f_in, sep="\t")
    spot_df[x_label] = spot_df[x_label] - spot_df[x_label].min()
    spot_df[y_label] = spot_df[y_label] - spot_df[y_label].min()
    x_range = spot_df[x_label].max()
    y_range = spot_df[y_label].max()

    pos_df = spot_df[[x_label, y_label]]
    pos_df = pos_df.drop_duplicates()
    nbg_label_arr = coo_matrix(([True] * pos_df.shape[0], (pos_df[x_label].to_numpy(), pos_df[y_label].to_numpy())), (x_range+1, y_range+1))
    nbg_label_arr = nbg_label_arr.toarray()

    all_gene_name_li = list(spot_df[gene_name_label].unique())

    with open(f_gene_list, "r") as f:
        for line in f.readlines():
            gid = line.rstrip("\n")
            if not gid in all_gene_name_li:
                continue
            expr_arr = mat2image_arr(gid, spot_df, nbg_label_arr, x_range, y_range, gene_name_label, UMI_label, x_label, y_label)
            img_arr = single_panel2multi_panel(x_range, y_range, nbg_label_arr, r_arr=expr_arr, g_arr=expr_arr, b_arr=None)
            plt.imsave(os.path.join(f_out, f"{gid}.png"), cv.rotate(img_arr, cv.ROTATE_90_COUNTERCLOCKWISE))


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=plot_gene_bin_expr_main)
    parser.add_argument("--input", type=str, dest="input",
                        metavar="spot.pred.tsv", required=True, help="Enhanced file generate by RedeViz")
    parser.add_argument("--x-label", type=str, dest="x_label", metavar="x_label", required=True)
    parser.add_argument("--y-label", type=str, dest="y_label", metavar="y_label", required=True)
    parser.add_argument("--gene-name-label", type=str, dest="gene_name_label", metavar="gene_name_label", required=True)
    parser.add_argument("--UMI-label", type=str, dest="UMI_label", metavar="UMI_label", required=True)
    parser.add_argument("--output", type=str, dest="output",
                        metavar="output", required=True, help="Output image")
    parser.add_argument("--gene-list", type=str, dest="gene_list", metavar="gene_list.txt", help="Gene list with TXT format")
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    args.func(args)

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()


