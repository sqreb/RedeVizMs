import numpy as np
import sys
import argparse
import pandas as pd
from matplotlib import pyplot as plt
from scipy.sparse import coo_matrix
import pickle

def hex2rgb(hex):
    r = int(hex[1:3], 16)
    g = int(hex[3:5], 16)
    b = int(hex[5:7], 16)
    return np.array([r, g, b], dtype=np.uint8)

def rgb2hex(r, g, b):
    return ('#{:02X}{:02X}{:02X}').format(r, g, b)

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--index", type=str, dest="index",
                        metavar="pretreat.pkl", required=True)
    base_group.add_argument("--color", type=str, dest="color",
                        metavar="color.tsv", required=True)
    base_group.add_argument("--embedding-plot", type=str, dest="embedding_plot",
                        metavar="embedding_plot.png", required=True)
    base_group.add_argument("--cell-type-plot", type=str, dest="cell_type_plot",
                        metavar="cell_type_plot.png", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_index = args.index
    f_color = args.color
    f_embedding_plot = args.embedding_plot
    f_cell_type_plot = args.cell_type_plot

    color_df = pd.read_csv(f_color, sep="\t")
    if "CellType" == color_df.columns[0]:
        color_df = color_df[["CellType", "R", "G", "B"]]
    else:
        color_df = color_df[["Annotation", "R", "G", "B"]]
    color_df.columns = ["MainCellType", "CellTypeR", "CellTypeG", "CellTypeB"]
    color_df["MainCellType"] = color_df["MainCellType"].astype(str)
    color_df["CellTypeRGB"] = [rgb2hex(r, g, b) for (r, g, b) in zip(
        color_df["CellTypeR"].to_numpy(),
        color_df["CellTypeG"].to_numpy(),
        color_df["CellTypeB"].to_numpy()
    )]

    with open(f_index, "rb") as f:
        index_dict = pickle.load(f)

    embedding_info_df = index_dict["embedding_info"]
    embedding_resolution = index_dict["embedding_resolution"]
    embedding_dim = index_dict["embedding_dim"]
    embedding_info_df["Embedding1"] = embedding_info_df["Embedding1BinIndex"] * embedding_resolution
    embedding_info_df["Embedding2"] = embedding_info_df["Embedding2BinIndex"] * embedding_resolution
    embedding_info_df["Embedding3"] = embedding_info_df["Embedding3BinIndex"] * embedding_resolution
    embedding_info_df["R"] = embedding_info_df["Embedding1"] + 20
    embedding_info_df["G"] = embedding_info_df["Embedding2"] + 20
    embedding_info_df["B"] = embedding_info_df["Embedding3"] + 20
    embedding_info_df["R"] = embedding_info_df["R"].astype(int)
    embedding_info_df["G"] = embedding_info_df["G"].astype(int)
    embedding_info_df["B"] = embedding_info_df["B"].astype(int)

    embedding_info_df = embedding_info_df.merge(color_df)
    x_range = embedding_info_df["Embedding1BinIndex"].max() + 1
    y_range = embedding_info_df["Embedding2BinIndex"].max() + 1
    if embedding_dim == 2:
        R_arr = coo_matrix((embedding_info_df["R"].to_numpy(), (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        G_arr = coo_matrix((embedding_info_df["G"].to_numpy(), (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        B_arr = coo_matrix((embedding_info_df["B"].to_numpy(), (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        RGB_arr = np.stack([R_arr.toarray(), G_arr.toarray(), B_arr.toarray()], -1)
        RGB_arr = RGB_arr.astype(np.uint8)
        plt.imsave(f_embedding_plot, RGB_arr)

        cell_type_R_arr = coo_matrix((embedding_info_df["CellTypeR"].to_numpy(), (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        cell_type_G_arr = coo_matrix((embedding_info_df["CellTypeG"].to_numpy(), (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        cell_type_B_arr = coo_matrix((embedding_info_df["CellTypeB"].to_numpy(), (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        cell_type_RGB_arr = np.stack([cell_type_R_arr.toarray(), cell_type_G_arr.toarray(), cell_type_B_arr.toarray()], -1)
        cell_type_RGB_arr = cell_type_RGB_arr.astype(np.uint8)
        plt.imsave(f_cell_type_plot, cell_type_RGB_arr)
    elif embedding_dim == 3:
        x = embedding_info_df["Embedding1BinIndex"].to_numpy()
        y = embedding_info_df["Embedding2BinIndex"].to_numpy()
        z = embedding_info_df["Embedding3BinIndex"].to_numpy()
        cell_tye_RGB = embedding_info_df["CellTypeRGB"].to_numpy()
        ax = plt.figure().add_subplot(projection='3d')
        ax.scatter(x, y, zs=z, c=cell_tye_RGB)
        ax.legend()
        ax.set_xlabel('Embedding1')
        ax.set_ylabel('Embedding2')
        ax.set_zlabel('Embedding3')
        plt.savefig(f_cell_type_plot, dpi=1200)

        RGB = [f"#{R:02X}{G:02X}{B:02X}" for (R, G, B) in zip(embedding_info_df["R"].to_numpy(), embedding_info_df["G"].to_numpy(), embedding_info_df["B"].to_numpy())]
        ax = plt.figure().add_subplot(projection='3d')
        ax.scatter(x, y, zs=z, c=RGB)
        ax.legend()
        ax.set_xlabel('Embedding1')
        ax.set_ylabel('Embedding2')
        ax.set_zlabel('Embedding3')
        plt.savefig(f_embedding_plot, dpi=1200)

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
