import numpy as np
import sys
import argparse
import pandas as pd
from matplotlib import pyplot as plt
from scipy.sparse import coo_matrix
import pickle
import os

def hex2rgb(hex):
    r = int(hex[1:3], 16)
    g = int(hex[3:5], 16)
    b = int(hex[5:7], 16)
    return np.array([r, g, b], dtype=np.uint8)


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--emb-color", type=str, dest="emb_color",
                        metavar="emb_color.tsv", required=True)
    base_group.add_argument("--embedding-resolution", type=float, dest="embedding_resolution",
                        metavar="embedding_resolution", required=True)
    base_group.add_argument("--embedding-dim", type=int, dest="embedding_dim",
                        metavar="embedding_dim", required=True)
    base_group.add_argument("--output", type=str, dest="output",
                        metavar="output", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_emb_color = args.emb_color
    embedding_resolution = args.embedding_resolution
    embedding_dim = args.embedding_dim
    f_out = args.output
    embedding_info_df = pd.read_csv(f_emb_color, sep="\t")

    embedding_info_df["Embedding1"] = embedding_info_df["Embedding1BinIndex"] * embedding_resolution
    embedding_info_df["Embedding2"] = embedding_info_df["Embedding2BinIndex"] * embedding_resolution
    embedding_info_df["Embedding3"] = embedding_info_df["Embedding3BinIndex"] * embedding_resolution
    embedding_info_df["R"] = embedding_info_df["Embedding1"] + 20
    embedding_info_df["G"] = embedding_info_df["Embedding2"] + 20
    embedding_info_df["B"] = embedding_info_df["Embedding3"] + 20
    embedding_info_df["R"] = embedding_info_df["R"].astype(int)
    embedding_info_df["G"] = embedding_info_df["G"].astype(int)
    embedding_info_df["B"] = embedding_info_df["B"].astype(int)

    x_range = embedding_info_df["Embedding1BinIndex"].max() + 1
    y_range = embedding_info_df["Embedding2BinIndex"].max() + 1
    if embedding_dim == 2:
        R_arr = coo_matrix((embedding_info_df["R"].to_numpy(), (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        G_arr = coo_matrix((embedding_info_df["G"].to_numpy(), (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        B_arr = coo_matrix((embedding_info_df["B"].to_numpy(), (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        RGB_arr = np.stack([R_arr.toarray(), G_arr.toarray(), B_arr.toarray()], -1)
        RGB_arr = RGB_arr.astype(np.uint8)
        plt.imsave(os.path.join(f_out, "Emb.Color.png"), RGB_arr)

        cell_type_RGB = np.array([hex2rgb(x) for x in embedding_info_df["cell_type_RGB"].to_numpy()])
        cell_type_R_arr = coo_matrix((cell_type_RGB[:,0], (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        cell_type_G_arr = coo_matrix((cell_type_RGB[:,1], (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        cell_type_B_arr = coo_matrix((cell_type_RGB[:,2], (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        cell_type_RGB_arr = np.stack([cell_type_R_arr.toarray(), cell_type_G_arr.toarray(), cell_type_B_arr.toarray()], -1)
        cell_type_RGB_arr = cell_type_RGB_arr.astype(np.uint8)
        plt.imsave(os.path.join(f_out, "CellType.Color.png"), cell_type_RGB_arr)

        level0_RGB = np.array([hex2rgb(x) for x in embedding_info_df["level0_RGB"].to_numpy()])
        level0_R_arr = coo_matrix((level0_RGB[:,0], (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        level0_G_arr = coo_matrix((level0_RGB[:,1], (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        level0_B_arr = coo_matrix((level0_RGB[:,2], (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        level0_RGB_arr = np.stack([level0_R_arr.toarray(), level0_G_arr.toarray(), level0_B_arr.toarray()], -1)
        level0_RGB_arr = level0_RGB_arr.astype(np.uint8)
        plt.imsave(os.path.join(f_out, "level0.Color.png"), level0_RGB_arr)

        level1_RGB = np.array([hex2rgb(x) for x in embedding_info_df["level1_RGB"].to_numpy()])
        level1_R_arr = coo_matrix((level1_RGB[:,0], (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        level1_G_arr = coo_matrix((level1_RGB[:,1], (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        level1_B_arr = coo_matrix((level1_RGB[:,2], (embedding_info_df["Embedding1BinIndex"].to_numpy(), embedding_info_df["Embedding2BinIndex"].to_numpy())), [x_range, y_range])
        level1_RGB_arr = np.stack([level1_R_arr.toarray(), level1_G_arr.toarray(), level1_B_arr.toarray()], -1)
        level1_RGB_arr = level1_RGB_arr.astype(np.uint8)
        plt.imsave(os.path.join(f_out, "level1.Color.png"), level1_RGB_arr)

    elif embedding_dim == 3:
        x = embedding_info_df["Embedding1BinIndex"].to_numpy()
        y = embedding_info_df["Embedding2BinIndex"].to_numpy()
        z = embedding_info_df["Embedding3BinIndex"].to_numpy()
        ax = plt.figure().add_subplot(projection='3d')
        ax.scatter(y, x, zs=z, c=embedding_info_df["cell_type_RGB"].to_numpy())
        ax.legend()
        ax.set_xlabel('Embedding2')
        ax.set_ylabel('Embedding1')
        ax.set_zlabel('Embedding3')
        plt.savefig(os.path.join(f_out, "CellType.Color.png"), dpi=1200)

        ax = plt.figure().add_subplot(projection='3d')
        ax.scatter(y, x, zs=z, c=embedding_info_df["level0_RGB"].to_numpy())
        ax.legend()
        ax.set_xlabel('Embedding2')
        ax.set_ylabel('Embedding1')
        ax.set_zlabel('Embedding3')
        plt.savefig(os.path.join(f_out, "level0.Color.png"), dpi=1200)

        ax = plt.figure().add_subplot(projection='3d')
        ax.scatter(y, x, zs=z, c=embedding_info_df["level1_RGB"].to_numpy())
        ax.legend()
        ax.set_xlabel('Embedding2')
        ax.set_ylabel('Embedding1')
        ax.set_zlabel('Embedding3')
        plt.savefig(os.path.join(f_out, "level1.Color.png"), dpi=1200)

        RGB = [f"#{R:02X}{G:02X}{B:02X}" for (R, G, B) in zip(embedding_info_df["R"].to_numpy(), embedding_info_df["G"].to_numpy(), embedding_info_df["B"].to_numpy())]
        ax = plt.figure().add_subplot(projection='3d')
        ax.scatter(y, x, zs=z, c=RGB)
        ax.legend()
        ax.set_xlabel('Embedding2')
        ax.set_ylabel('Embedding1')
        ax.set_zlabel('Embedding3')
        plt.savefig(os.path.join(f_out, "Emb.Color.png"), dpi=1200)

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
