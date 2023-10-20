import numpy as np
import os
import sys
import argparse
import pandas as pd
import cv2 as cv
from matplotlib import pyplot as plt
from scipy.sparse import coo_matrix

def mat2image(expr_arr, nbg_label_arr):
    nzo_arr = expr_arr[expr_arr>0]
    cutoff = np.percentile(nzo_arr, 95)
    expr_arr[expr_arr>cutoff] = cutoff
    norm_mat = expr_arr / cutoff
    
    r_arr = norm_mat * 215 + 40
    g_arr = norm_mat * 215 + 40
    b_arr = np.ones_like(norm_mat) * 40

    r_arr[~nbg_label_arr] = 0
    g_arr[~nbg_label_arr] = 0
    b_arr[~nbg_label_arr] = 0

    rgb_arr = np.stack([r_arr, g_arr, b_arr], -1)
    rgb_arr = rgb_arr.astype(np.uint8)
    return rgb_arr

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--spot", type=str, dest="spot",
                        metavar="spot.tsv", required=True, help="Spot information")
    base_group.add_argument("--gene-id-label", type=str, dest="gene_id_label",
                        metavar="gene_id_label", default="geneID", required=False)
    base_group.add_argument("--UMI-label", type=str, dest="UMI_label",
                        metavar="UMI_label", default="MIDCounts", required=False)
    base_group.add_argument("--output", type=str, dest="output",
                        metavar="output", required=True, help="Output folder")
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_spot = args.spot
    f_out = args.output
    gene_id_label = args.gene_id_label
    UMI_label = args.UMI_label

    spot_df = pd.read_csv(f_spot, sep="\t")
    spot_df["bin_x_index"] = spot_df["bin_x_index"] - spot_df["bin_x_index"].min()
    spot_df["bin_y_index"] = spot_df["bin_y_index"] - spot_df["bin_y_index"].min()
    x_range = spot_df["bin_x_index"].max()
    y_range = spot_df["bin_y_index"].max()

    pos_df = spot_df[["bin_x_index", "bin_y_index"]]
    pos_df = pos_df.drop_duplicates()
    nbg_label_arr = coo_matrix(([True] * pos_df.shape[0], (pos_df["bin_x_index"].to_numpy(), pos_df["bin_y_index"].to_numpy())), (x_range+1, y_range+1))
    nbg_label_arr = nbg_label_arr.toarray()

    all_gene_li = np.array(spot_df[gene_id_label].unique())
    for gid in all_gene_li:
        f_fig = os.path.join(f_out, f"{gid}.png")
        if os.path.exists(f_fig):
            continue
        tmp_spot_df = spot_df[spot_df[gene_id_label]==gid]
        UMI_mat = coo_matrix((tmp_spot_df[UMI_label].to_numpy(), (tmp_spot_df["bin_x_index"].to_numpy(), tmp_spot_df["bin_y_index"].to_numpy())), (x_range+1, y_range+1))
        img_mat = mat2image(UMI_mat.toarray(), nbg_label_arr)
        plt.imsave(f_fig, cv.rotate(img_mat, cv.ROTATE_90_COUNTERCLOCKWISE))

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
