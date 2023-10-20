import os
import anndata
import sys
import argparse
import numpy as np


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--input", type=str, dest="input",
                        metavar="input.h5ad", required=True, help="Reference scRNA-seq UMI matrix with H5AD format")
    base_group.add_argument("--gene-list", type=str, dest="gene_li",
                        metavar="gene_li.txt", required=True, 
                        help="Select gene list with TSV format.")
    base_group.add_argument("--gene-name-label", type=str, dest="gene_name_label",
                        metavar="gene_name_label", required=False, default="GeneSymbol", 
                        help="Gene label (default: GeneSymbol)")
    base_group.add_argument("--output", type=str, dest="output",
                        metavar="output.h5ad", required=True, help="Output file with H5AD format")
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    sce = anndata.read_h5ad(args.input)
    select_gene_set = set()
    with open(args.gene_li, "r") as f:
        for line in f.readlines():
            select_gene_set.add(line.rstrip())
    sce_gene_li = sce.var[args.gene_name_label].to_numpy()
    is_in_sce_gene = np.array([x in select_gene_set for x in sce_gene_li])
    select_sce = sce[:, is_in_sce_gene]
    select_sce.write_h5ad(args.output)


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
