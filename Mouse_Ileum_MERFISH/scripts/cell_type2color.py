import sys
import argparse
import pandas as pd
from matplotlib.colors import hsv_to_rgb
import numpy as np

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--emb-info", type=str, dest="emb_info", metavar="Embedding.bin.info.tsv", required=True)
    base_group.add_argument("--output", type=str, dest="output", metavar="color.tsv", required=True, help="Cell type color")
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_emb_info = args.emb_info
    f_out = args.output

    emb_info_df = pd.read_csv(f_emb_info, sep="\t")
    uniq_cell_type = sorted(emb_info_df["MainCellType"].unique())
    cell_type_num = len(uniq_cell_type)
    H_li = [0.9*i/cell_type_num for i in range(cell_type_num)]
    rgb_arr = np.array([hsv_to_rgb((x, 0.9, 0.9)) for x in H_li])
    color_df = pd.DataFrame(
        {
            "CellType":uniq_cell_type, 
            "R":(rgb_arr[:,0]*255).astype(int),
            "G":(rgb_arr[:,1]*255).astype(int),
            "B":(rgb_arr[:,2]*255).astype(int)
        }
    )
    color_df.to_csv(f_out, sep="\t", index_label=False, index=False)

def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()
