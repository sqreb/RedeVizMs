import numpy as np
from scipy.sparse import load_npz
import cv2 as cv
from matplotlib import pyplot as plt
import argparse
import sys


def load_expr_arr(fname):
    if fname is None:
        return None
    expr_mat = load_npz(fname)
    expr_mat = expr_mat.toarray()
    nzo_arr = expr_mat[expr_mat>0]
    cutoff = np.percentile(nzo_arr, 95)
    expr_mat[expr_mat>cutoff] = cutoff
    norm_mat = expr_mat / cutoff * 215
    return norm_mat

def single_panel2multi_panel(x_range, y_range, r_arr=None, g_arr=None, b_arr=None, sig_arr=None):
    if r_arr is None:
        r_arr = np.zeros((x_range, y_range), dtype=np.float64)
    else:
        r_arr = r_arr[:x_range, :y_range]
    if g_arr is None:
        g_arr = np.zeros((x_range, y_range), dtype=np.float64)
    else:
        g_arr = g_arr[:x_range, :y_range]
    if b_arr is None:
        b_arr = np.zeros((x_range, y_range), dtype=np.float64)
    else:
        b_arr = b_arr[:x_range, :y_range]
    r_arr = r_arr + 40
    g_arr = g_arr + 40
    b_arr = b_arr + 40
    if sig_arr is not None:
        r_arr[~sig_arr] = 0
        g_arr[~sig_arr] = 0
        b_arr[~sig_arr] = 0     
    rgb_arr = np.stack([r_arr, g_arr, b_arr], -1)
    rgb_arr = rgb_arr.astype(np.uint8)
    return rgb_arr

def plot_imputation(args):
    f_R = args.R
    f_G = args.G
    f_B = args.B
    f_sig = args.signal
    f_out = args.output

    assert not all([x is None for x in [f_R, f_G, f_B]])

    R_expr_arr = load_expr_arr(f_R)
    G_expr_arr = load_expr_arr(f_G)
    B_expr_arr = load_expr_arr(f_B)

    sig_arr = load_npz(f_sig).toarray()
    x_index, y_index = np.where(sig_arr>0)
    min_x = min(x_index)
    min_y = min(y_index)
    x_range, y_range = sig_arr.shape
    img_arr = single_panel2multi_panel(x_range, y_range, r_arr=R_expr_arr, g_arr=G_expr_arr, b_arr=B_expr_arr, sig_arr=sig_arr)
    img_arr = img_arr[min_x:, min_y:, :]
    plt.imsave(f_out, cv.rotate(img_arr, cv.ROTATE_90_COUNTERCLOCKWISE))


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=plot_imputation)
    parser.add_argument("-R", type=str, dest="R", metavar="R.expr.npz", required=False, default=None)
    parser.add_argument("-G", type=str, dest="G", metavar="G.expr.npz", required=False, default=None)
    parser.add_argument("-B", type=str, dest="B", metavar="B.expr.npz", required=False, default=None)
    parser.add_argument("--signal", type=str, dest="signal", metavar="signal.npz", required=True)
    parser.add_argument("-o", "--output", type=str, dest="output", metavar="output.png", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    args.func(args)

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()


