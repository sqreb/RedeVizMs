import tifffile
import sys
import argparse
import re
import os
import matplotlib.pyplot as plt
import numpy as np

def load_image(f_img):
    img_dict = dict()
    shape = None
    with tifffile.TiffFile(f_img) as tif:
        for page in tif.pages:
            page_desc = page.description
            if page_desc.find("FullResolution") == -1:
                continue
            name_li = list(set(re.findall("<Biomarker>([ \S]+)</Biomarker>", page.description)))
            if len(name_li)==0:
                continue
            name = name_li[0]
            img_dict[name] = page.asarray()
            shape = img_dict[name].shape
    return img_dict, shape

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--input", type=str, dest="input", metavar="input.qptiff", required=True)
    base_group.add_argument("--output", type=str, dest="output", metavar="output.png", required=True)
    base_group.add_argument("-R", type=str, dest="R", metavar="R", required=False, default=None)
    base_group.add_argument("-G", type=str, dest="G", metavar="G", required=False, default=None)
    base_group.add_argument("-B", type=str, dest="B", metavar="B", required=False, default=None)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output

    img_dict, shape = load_image(f_in)
    if args.R is None:
        R_arr = np.zeros(shape=shape, dtype=np.uint8)
    else:
        R_arr = img_dict[args.R]

    if args.G is None:
        G_arr = np.zeros(shape=shape, dtype=np.uint8)
    else:
        G_arr = img_dict[args.G]

    if args.B is None:
        B_arr = np.zeros(shape=shape, dtype=np.uint8)
    else:
        B_arr = img_dict[args.B]
    
    RGB_arr = np.stack([R_arr, G_arr, B_arr], -1).astype(np.uint8)
    plt.imsave(f_out, RGB_arr)

def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()
