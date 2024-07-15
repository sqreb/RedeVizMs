import pickle
import sys
import argparse
import os
import numpy as np
import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--index", type=str, dest="index", metavar="index.pkl", required=True)
    base_group.add_argument("--celltype", type=str, dest="cell_type", metavar="cell_type.tsv", required=True)
    base_group.add_argument("--output", type=str, dest="output", metavar="output", required=True)
    return parser.parse_args(args)

def fillnan(arr, fill):
    res = list()
    for x in arr:
        if x is np.nan:
            res.append(fill)
        else:
            res.append(x)
    return res

def main(args):
    args = parse_args(args)
    f_index = args.index
    f_cell_type = args.cell_type
    f_out = args.output

    with open(f_index, "rb") as f:
        embedding_dict = pickle.load(f)

    cell_info_df = embedding_dict["cell_info"]
    cell_type_df = pd.read_csv(f_cell_type, sep="\t")
    cell_type_df["CellName"] = ["{0}.{1}".format(fov, id) for (fov, id) in zip(cell_type_df["fov"].to_numpy(), cell_type_df["cell_ID"].to_numpy())]
    cell_info_df = cell_info_df.merge(cell_type_df, on="CellName")
    cell_info_df["cell_type"] = fillnan(cell_info_df["cell_type"].to_numpy(), "LQ")
    cell_info_df["niche"] = fillnan(cell_info_df["niche"].to_numpy(), "LQ")
    cell_info_df["level1"] = fillnan(cell_info_df["level1"].to_numpy(), "LQ")
    cell_info_df["level0"] = fillnan(cell_info_df["level0"].to_numpy(), "LQ")
    cell_info_df["cell_type_RGB"] = fillnan(cell_info_df["cell_type_RGB"].to_numpy(), "#AAAAAA")
    cell_info_df["level1_RGB"] = fillnan(cell_info_df["level1_RGB"].to_numpy(), "#AAAAAA")
    cell_info_df["level0_RGB"] = fillnan(cell_info_df["level0_RGB"].to_numpy(), "#AAAAAA")
    cell_info_df.to_csv(os.path.join(f_out, "CellInfo.WithAnn.tsv"), sep="\t", index_label=False, index=False)

    cell_info_df["BinIndex"] = np.nan_to_num(cell_info_df["BinIndex"].to_numpy(), nan=-1)
    cell_info_df = cell_info_df[cell_info_df["BinIndex"]>=0]
    cell_info_df["BinIndex"] = cell_info_df["BinIndex"].astype(int)
    cell_info_no_LQ = cell_info_df[cell_info_df["cell_type"]!="LQ"]

    major_cell_type = cell_info_no_LQ.groupby(["BinIndex", "cell_type", "cell_type_RGB"]).size().reset_index(name='cell_type_num').sort_values(["BinIndex", "cell_type_num"], ascending=False).groupby(["BinIndex"]).head(1)
    major_level0 = cell_info_no_LQ.groupby(["BinIndex", "level0", "level0_RGB"]).size().reset_index(name='level0_num').sort_values(["BinIndex", "level0_num"], ascending=False).groupby(["BinIndex"]).head(1)
    major_level1 = cell_info_no_LQ.groupby(["BinIndex", "level1", "level1_RGB"]).size().reset_index(name='level1_num').sort_values(["BinIndex", "level1_num"], ascending=False).groupby(["BinIndex"]).head(1)

    emb_info = embedding_dict["embedding_info"]
    emb_info = emb_info.merge(major_cell_type, how="left").merge(major_level1, how="left").merge(major_level0, how="left")
    emb_info = emb_info.sort_values("BinIndex")

    emb_info["cell_type"] = fillnan(emb_info["cell_type"].to_numpy(), "LQ")
    emb_info["level1"] = fillnan(emb_info["level1"].to_numpy(), "LQ")
    emb_info["level0"] = fillnan(emb_info["level0"].to_numpy(), "LQ")
    emb_info["cell_type_RGB"] = fillnan(emb_info["cell_type_RGB"].to_numpy(), "#AAAAAA")
    emb_info["level1_RGB"] = fillnan(emb_info["level1_RGB"].to_numpy(), "#AAAAAA")
    emb_info["level0_RGB"] = fillnan(emb_info["level0_RGB"].to_numpy(), "#AAAAAA")
    emb_info.to_csv(os.path.join(f_out, "EmbInfo.WithAnn.tsv"), sep="\t", index_label=False, index=False)


def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()

