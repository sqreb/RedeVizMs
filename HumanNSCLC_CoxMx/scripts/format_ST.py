import pickle
import sys
import argparse
import os
import numpy as np
import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--input", type=str, dest="input", metavar="tx_file.csv", required=True)
    base_group.add_argument("--output", type=str, dest="output", metavar="output", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output

    if not os.path.exists(f_out):
        os.makedirs(f_out)
    raw_data_df = pd.read_csv(f_in)
    raw_data_df["x_um"] = raw_data_df["x_global_px"].to_numpy() * 0.18
    raw_data_df["y_um"] = raw_data_df["y_global_px"].to_numpy() * 0.18
    raw_data_df["x_um"] = raw_data_df["x_um"] - raw_data_df["x_um"].min()
    raw_data_df["y_um"] = raw_data_df["y_um"] - raw_data_df["y_um"].min()
    raw_data_df.to_csv(os.path.join(f_out, "RawSignal.tsv"), sep="\t", index_label=False, index=False)

    raw_data_df["x_1um_index"] = raw_data_df["x_um"].astype(int)
    raw_data_df["y_1um_index"] = raw_data_df["y_um"].astype(int)

    raw_data_df["x_500nm_index"] = (raw_data_df["x_um"].to_numpy() * 2).astype(int)
    raw_data_df["y_500nm_index"] = (raw_data_df["y_um"].to_numpy() * 2).astype(int)

    no_Membrane_df = raw_data_df[raw_data_df["CellComp"]!="Membrane"]

    all_dict = {
        "All": raw_data_df,
        "NoMembrane": no_Membrane_df
    }
    for label, df in all_dict.items():
        cnt_500_nm_df = df.groupby(["x_500nm_index", "y_500nm_index", "target"]).size().reset_index()
        cnt_500_nm_df.columns = ["x", "y", "Gid", "UMI"]
        cnt_500_nm_df.to_csv(os.path.join(f_out, f"{label}.500nm.tsv"), sep="\t", index_label=False, index=False)

        cnt_1_um_df = df.groupby(["x_1um_index", "y_1um_index", "target"]).size().reset_index()
        cnt_1_um_df.columns = ["x", "y", "Gid", "UMI"]
        cnt_1_um_df.to_csv(os.path.join(f_out, f"{label}.1um.tsv"), sep="\t", index_label=False, index=False)


def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()

