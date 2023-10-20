import sys
import argparse
import os
import pandas as pd
import cv2
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_array


def plot_each_cell_type(f_dir, spot_df, x_label, y_label, cell_type_label):
    if not os.path.exists(f_dir):
        os.makedirs(f_dir)
    x_range = int(spot_df[x_label].max() + 1)
    y_range = int(spot_df[y_label].max() + 1)
    cell_type_li = spot_df[cell_type_label].unique()
    all_cell_pos = coo_array((np.ones(len(spot_df)), (spot_df[x_label].to_numpy(), spot_df[y_label].to_numpy())), shape=(x_range, y_range)).toarray()
    all_cell_type_R = 40 * all_cell_pos
    all_cell_type_G = 40 * all_cell_pos
    all_cell_type_B = 40 * all_cell_pos
    for cell_type in cell_type_li:
        tmp_spot_df = spot_df[spot_df[cell_type_label]==cell_type]
        cell_type_pos_arr = coo_array((np.ones(len(tmp_spot_df)), (tmp_spot_df[x_label].to_numpy(), tmp_spot_df[y_label].to_numpy())), shape=(x_range, y_range)).toarray().astype(bool)
        tmp_img_R = all_cell_type_R.copy()
        tmp_img_G = all_cell_type_G.copy()
        tmp_img_B = all_cell_type_B.copy()
        tmp_img_R[cell_type_pos_arr] = 255
        tmp_img_G[cell_type_pos_arr] = 255
        tmp_img_B[cell_type_pos_arr] = 0
        tmp_img = np.stack([tmp_img_R, tmp_img_G, tmp_img_B], -1).astype(np.uint8)
        fname = cell_type.replace("-", "_").replace("/", "_")
        plt.imsave(os.path.join(f_dir, f"{fname}.label.png"), cv2.rotate(tmp_img, cv2.ROTATE_90_COUNTERCLOCKWISE))

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--spot", type=str, dest="spot", metavar="spot.predict.tsv", required=True)
    base_group.add_argument("--output", type=str, dest="output", metavar="output", required=True, help="Output folder")
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_spot = args.spot
    f_out = args.output

    if not os.path.exists(f_out):
        os.makedirs(f_out)

    spot_df = pd.read_csv(f_spot, sep="\t")

    no_bg_spot = spot_df[spot_df["LabelTransfer"]!="Background"]
    spot_df["x"] = spot_df["x"] - no_bg_spot["x"].min()
    spot_df["y"] = spot_df["y"] - no_bg_spot["y"].min()
    spot_df = spot_df[spot_df["x"]>=0]
    spot_df = spot_df[spot_df["y"]>=0]

    plot_each_cell_type(f_out, spot_df, "x", "y", "ArgMaxCellType")

def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()
