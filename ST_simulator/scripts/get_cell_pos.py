import sys
import pandas as pd
import argparse
import os
import numpy as np
from shapely.prepared import prep
from shapely import wkt
from shapely.geometry import Point

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--input", type=str, dest="ref_cell_shape",
                        metavar="CellShapeInfo.tsv", required=True)
    base_group.add_argument("--output", type=str, dest="output",
                        metavar="output.npy", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_ref_cell_shape = args.ref_cell_shape
    f_out = args.output

    if not os.path.exists(f_out):
        os.makedirs(f_out)

    x_range = 2000
    y_range = 2000
    reselution = 500

    ref_cell_df = pd.read_csv(f_ref_cell_shape, sep="\t")

    cell_name_arr = -1 * np.ones([x_range, y_range], dtype=np.int64)
    cell_type_arr = np.zeros([x_range, y_range], dtype="U64")
    bin_index_arr = -1 * np.ones([x_range, y_range], dtype=np.int64)
    dist_to_boundary_arr = 1e6 * np.ones([x_range, y_range], dtype=np.float64)
    for cell_name, cell_type, bin_index, cell_shape in ref_cell_df[["CellName", "CellType", "BinIndex", "CellShape"]].to_numpy():
        poly = wkt.loads(cell_shape)
        poly_boundary = poly.boundary
        xmin, ymin, xmax, ymax = poly.bounds
        prep_poly = prep(poly)
        for x_index in range(max(0, int(xmin/reselution) - 10), int(xmax/reselution) + 10):
            if x_index >= x_range:
                continue
            for y_index in range(max(0, int(ymin/reselution)-10), int(ymax/reselution) + 10):
                if y_index >= y_range:
                    continue
                point = Point(x_index*reselution+reselution/2, y_index*reselution+reselution/2)
                dist2boundary = poly_boundary.distance(point)
                dist_to_boundary_arr[x_index, y_index] = min(dist2boundary / reselution, dist_to_boundary_arr[x_index, y_index])
                if prep_poly.contains(point):
                    cell_name_arr[x_index, y_index] = cell_name
                    cell_type_arr[x_index, y_index] = cell_type
                    bin_index_arr[x_index, y_index] = bin_index
    dist_to_boundary_arr[dist_to_boundary_arr>2000] = -1

    np.save(os.path.join(f_out, "Cell.Index.pos.npy"), cell_name_arr)
    np.save(os.path.join(f_out, "Cell.CellType.pos.npy"), cell_type_arr)
    np.save(os.path.join(f_out, "Cell.BinIndex.pos.npy"), bin_index_arr)
    np.save(os.path.join(f_out, "Cell.Dist2Bound.pos.npy"), dist_to_boundary_arr)

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
