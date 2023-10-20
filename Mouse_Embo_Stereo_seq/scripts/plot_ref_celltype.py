import sys
import argparse
import scanpy as sc
import pandas as pd
from plotnine import ggplot, geom_point, scale_color_manual, element_rect, geom_bar, geom_rect, geom_text, annotate, scale_y_log10, scale_fill_manual, aes, labs, theme_bw, theme, element_text, element_blank, ggsave
from plotnine.scales import scale_x_continuous, scale_y_continuous

def RGB2hex(R, G, B):
    return f"#{R:02X}{G:02X}{B:02X}"

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--input", type=str, dest="input", metavar="input.h5ad", required=True)
    base_group.add_argument("--color", type=str, dest="color", metavar="color.tsv", required=True)
    base_group.add_argument("--output", type=str, dest="output", metavar="gene_li.png", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_color = args.color
    f_out = args.output

    sce = sc.read_h5ad(f_in)
    sce.obs["sp_x"] = sce.obsm["spatial"][:, 0]
    sce.obs["sp_y"] = sce.obsm["spatial"][:, 1]
    sce.obs["sp_x"] = sce.obs["sp_x"] - sce.obs["sp_x"].min()
    sce.obs["sp_y"] = sce.obs["sp_y"] - sce.obs["sp_y"].min()

    st_df = sce.obs[["sp_x", "sp_y", "annotation"]]
    st_df = st_df[st_df["annotation"]!="Cavity"]

    color_df = pd.read_csv(f_color, sep="\t")
    color_df["hex"] = [RGB2hex(R, G, B) for (R, G, B) in zip(
        color_df["R"].to_numpy(),
        color_df["G"].to_numpy(),
        color_df["B"].to_numpy()
    )]
    cell_color_dict = {key: val for (key, val) in zip(color_df["CellType"].to_numpy(), color_df["hex"].to_numpy())}

    x_range = st_df["sp_x"].max() - st_df["sp_x"].min()
    y_range = st_df["sp_y"].max() - st_df["sp_y"].min()

    p = ggplot(st_df, aes(x="sp_x", y="sp_y", color="annotation")) + \
        geom_point(size=0.3) + \
        scale_color_manual(values=cell_color_dict) + \
        scale_x_continuous(expand=(0, 0)) + \
        scale_y_continuous(expand=(0, 0)) + \
        labs(x="x", y="y", color="Cell type") + \
        theme_bw() + \
        theme(
            text = element_text(family="Arial", size=5),
            title = element_text(family="Arial", size=6),
            axis_text = element_text(color = "black"),
            panel_background = element_rect(fill="black"),
            legend_title = element_text(family = "Arial", size=7, color="black"),
            legend_text = element_text(family = "Arial", size=6, color="black"),
            legend_key_size = 3,
            panel_grid = element_blank()
        )
    ggsave(p, filename=f_out, width=x_range/10, height=y_range/10, limitsize=False, units="cm")


def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()
