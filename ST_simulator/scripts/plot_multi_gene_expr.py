import numpy as np
import os
import sys
import argparse
import pandas as pd
import cv2 as cv
from matplotlib import pyplot as plt
from scipy.sparse import coo_matrix

def mat2image_arr(gid, spot_df, nbg_label_arr, x_range, y_range):
    tmp_spot_df = spot_df[spot_df["Gid"]==gid]
    expr_arr = coo_matrix((tmp_spot_df["UMI"].to_numpy(), (tmp_spot_df["bin_x_index"].to_numpy(), tmp_spot_df["bin_y_index"].to_numpy())), (x_range+1, y_range+1)).toarray()
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

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--spot", type=str, dest="spot",
                        metavar="spot.tsv", required=True, help="Spot information")
    base_group.add_argument("--output", type=str, dest="output",
                        metavar="output", required=True, help="Output folder")
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_spot = args.spot
    f_out = args.output

    spot_df = pd.read_csv(f_spot, sep="\t")
    spot_df["bin_x_index"] = spot_df["bin_x_index"] - spot_df["bin_x_index"].min()
    spot_df["bin_y_index"] = spot_df["bin_y_index"] - spot_df["bin_y_index"].min()
    x_range = spot_df["bin_x_index"].max()
    y_range = spot_df["bin_y_index"].max()

    pos_df = spot_df[["bin_x_index", "bin_y_index"]]
    pos_df = pos_df.drop_duplicates()
    nbg_label_arr = coo_matrix(([True] * pos_df.shape[0], (pos_df["bin_x_index"].to_numpy(), pos_df["bin_y_index"].to_numpy())), (x_range+1, y_range+1))
    nbg_label_arr = nbg_label_arr.toarray()

    KRT12_arr = mat2image_arr("KRT12", spot_df, nbg_label_arr, x_range, y_range)
    KRT19_arr = mat2image_arr("KRT19", spot_df, nbg_label_arr, x_range, y_range)
    CD74_arr = mat2image_arr("CD74", spot_df, nbg_label_arr, x_range, y_range)
    img_arr = single_panel2multi_panel(x_range, y_range, nbg_label_arr, r_arr=KRT12_arr, g_arr=KRT19_arr, b_arr=CD74_arr)
    plt.imsave(os.path.join(f_out, "KRT12_KRT19_CD74.png"), cv.rotate(img_arr, cv.ROTATE_90_COUNTERCLOCKWISE))


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
