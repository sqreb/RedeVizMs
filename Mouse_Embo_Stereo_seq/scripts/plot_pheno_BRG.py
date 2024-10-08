import numpy as np
import cv2 as cv
from matplotlib import pyplot as plt
import pandas as pd
from scipy.sparse import coo_matrix
from redeviz.posttreatment.utils import filter_pred_df
import argparse
import sys
import os



def plot_phenotype_main(args):
    f_spot = args.input
    f_out = args.output
    keep_other = args.keep_other
    white_bg = args.white_bg

    spot_df = pd.read_csv(f_spot, sep="\t")
    spot_df = spot_df[spot_df["RefCellTypeScore"] > spot_df["BackgroundScore"]]
    spot_df["x"] = spot_df["x"] - spot_df["x"].min()
    spot_df["y"] = spot_df["y"] - spot_df["y"].min()

    x_range = spot_df["x"].max() + 1
    y_range = spot_df["y"].max() + 1
    nbg_spot_df = spot_df[spot_df["LabelTransfer"]!="Background"]
    if args.denoise:
        nbg_spot_df = filter_pred_df(nbg_spot_df, min_spot_in_region=args.min_spot_num)

    if keep_other:
        sig_df = nbg_spot_df[nbg_spot_df["ArgMaxCellType"]!="Other"]
        other_df = nbg_spot_df[nbg_spot_df["ArgMaxCellType"]=="Other"]
    else:
        sig_df = nbg_spot_df[nbg_spot_df["LabelTransfer"]!="Other"]
        other_df = nbg_spot_df[nbg_spot_df["LabelTransfer"]=="Other"]
        
    other_pos = coo_matrix(([1]*other_df.shape[0], (other_df["x"].to_numpy(), other_df["y"].to_numpy())), (x_range, y_range)).toarray()
    other_pos = other_pos.astype(bool)
    R_arr = coo_matrix((sig_df["Embedding2"].to_numpy().astype(int)+20, (sig_df["x"].to_numpy(), sig_df["y"].to_numpy())), (x_range, y_range)).toarray()
    G_arr = coo_matrix((sig_df["Embedding3"].to_numpy().astype(int)+20, (sig_df["x"].to_numpy(), sig_df["y"].to_numpy())), (x_range, y_range)).toarray()
    B_arr = coo_matrix((sig_df["Embedding1"].to_numpy().astype(int)+20, (sig_df["x"].to_numpy(), sig_df["y"].to_numpy())), (x_range, y_range)).toarray()
    R_arr[other_pos] = 150
    G_arr[other_pos] = 150
    B_arr[other_pos] = 150
    RGB_arr = np.stack([R_arr, G_arr, B_arr], -1)
    RGB_arr = RGB_arr.astype(np.uint8)
    if white_bg:
        bg_pos = RGB_arr.sum(-1) == 0
        RGB_arr[bg_pos] = 255
    plt.imsave(f_out, cv.rotate(RGB_arr, cv.ROTATE_90_COUNTERCLOCKWISE))


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=plot_phenotype_main)
    parser.add_argument("--input", type=str, dest="input",
                        metavar="spot.pred.tsv", required=True, help="Enhanced file generate by RedeViz")
    parser.add_argument("--output", type=str, dest="output",
                        metavar="output.png", required=True, help="Output image")
    parser.add_argument("--keep-other", dest="keep_other", action="store_true", help="Force to impute other cell types")
    parser.add_argument("--min-denoise-spot-num", type=int, dest="min_spot_num", default=200)
    parser.add_argument("--denoise", dest="denoise", action="store_true", help="Remove scattered pixel")
    parser.add_argument("--white-bg", dest="white_bg", action="store_true")
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    args.func(args)

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()



