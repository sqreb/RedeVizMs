import sys
import argparse
import os
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from matplotlib import pyplot as plt
import cv2 as cv


def get_img_arr(plot_df, RGB_flag, other_pos=None):
    RGB_li = plot_df[RGB_flag].to_numpy()
    R_li = [int(f"{x[1:3]}", 16) for x in RGB_li]
    G_li = [int(f"{x[3:5]}", 16) for x in RGB_li]
    B_li = [int(f"{x[5:7]}", 16) for x in RGB_li]
    x_range = plot_df["x"].max() + 1
    y_range = plot_df["y"].max() + 1
    R_arr = coo_matrix((R_li, (plot_df["x"].to_numpy(), plot_df["y"].to_numpy())), (x_range, y_range)).toarray()
    G_arr = coo_matrix((G_li, (plot_df["x"].to_numpy(), plot_df["y"].to_numpy())), (x_range, y_range)).toarray()
    B_arr = coo_matrix((B_li, (plot_df["x"].to_numpy(), plot_df["y"].to_numpy())), (x_range, y_range)).toarray()
    if other_pos is not None:
        R_arr[other_pos] = 150
        G_arr[other_pos] = 150
        B_arr[other_pos] = 150
    RGB_arr = np.stack([R_arr, G_arr, B_arr], -1)
    RGB_arr = RGB_arr.astype(np.uint8)
    return RGB_arr

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--input", type=str, dest="input", metavar="spot.tsv", required=True)
    base_group.add_argument("--emb-info", type=str, dest="emb_info", metavar="emb_info.tsv", required=True)
    base_group.add_argument("--output", type=str, dest="output", metavar="output", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_emb_info = args.emb_info
    f_out = args.output

    if not os.path.exists(f_out):
        os.makedirs(f_out)

    pred_df = pd.read_csv(f_in, sep="\t")
    emb_info_df = pd.read_csv(f_emb_info, sep="\t")
    emb_name_li = list(emb_info_df.columns)
    emb_name_li[3] = "EmbeddingState"
    emb_info_df.columns = emb_name_li
    plot_df = pred_df.merge(emb_info_df)
    plot_df["x"] = plot_df["x"] - plot_df["x"].min()
    plot_df["y"] = plot_df["y"] - plot_df["y"].min()

    celltype_RGB_arr = get_img_arr(plot_df, "cell_type_RGB")
    level1_RGB_arr = get_img_arr(plot_df, "level1_RGB")
    level0_RGB_arr = get_img_arr(plot_df, "level0_RGB")
    plt.imsave(os.path.join(f_out, "cell_type.all.png"), cv.rotate(celltype_RGB_arr, cv.ROTATE_90_COUNTERCLOCKWISE))
    plt.imsave(os.path.join(f_out, "level1.all.png"), cv.rotate(level1_RGB_arr, cv.ROTATE_90_COUNTERCLOCKWISE))
    plt.imsave(os.path.join(f_out, "level0.all.png"), cv.rotate(level0_RGB_arr, cv.ROTATE_90_COUNTERCLOCKWISE))

    plot_HQ_df = plot_df[plot_df["LabelTransfer"]!="Background"]
    other_df = plot_HQ_df[plot_HQ_df["LabelTransfer"]=="Other"]
    other_pos = coo_matrix(([1]*other_df.shape[0], (other_df["x"].to_numpy(), other_df["y"].to_numpy())), (plot_HQ_df["x"].max()+1, plot_HQ_df["y"].max()+1)).toarray()
    other_pos = other_pos.astype(bool)
    celltype_HQ_RGB_arr = get_img_arr(plot_HQ_df, "cell_type_RGB", other_pos)
    level1_HQ_RGB_arr = get_img_arr(plot_HQ_df, "level1_RGB", other_pos)
    level0_HQ_RGB_arr = get_img_arr(plot_HQ_df, "level0_RGB", other_pos)
    plt.imsave(os.path.join(f_out, "cell_type.HQ.png"), cv.rotate(celltype_HQ_RGB_arr, cv.ROTATE_90_COUNTERCLOCKWISE))
    plt.imsave(os.path.join(f_out, "level1.HQ.png"), cv.rotate(level1_HQ_RGB_arr, cv.ROTATE_90_COUNTERCLOCKWISE))
    plt.imsave(os.path.join(f_out, "level0.HQ.png"), cv.rotate(level0_HQ_RGB_arr, cv.ROTATE_90_COUNTERCLOCKWISE))


def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()

