import sys
import pandas as pd
import argparse
import os
from scipy.sparse import coo_matrix
from collections import defaultdict
import numpy as np
from shapely.prepared import prep
from shapely import wkt
from shapely.geometry import Point
from matplotlib import pyplot as plt
from plotnine import ggplot, scale_x_log10, scale_x_continuous, scale_y_continuous, scale_fill_gradient, geom_point, geom_bar, geom_rect, geom_text, annotate, scale_y_log10, scale_fill_manual, aes, labs, theme_bw, theme, element_text, element_blank, ggsave


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--pred", type=str, dest="pred",
                        metavar="pred.tsv", required=True)
    base_group.add_argument("--simi", type=str, dest="simi",
                        metavar="simi.tsv", required=True)
    base_group.add_argument("--output", type=str, dest="output",
                        metavar="output", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_simi = args.simi
    f_pred = args.pred
    f_out = args.output

    if not os.path.exists(f_out):
        os.makedirs(f_out)

    x_range = 2000
    y_range = 2000
    reselution = 500

    simi_df = pd.read_csv(f_simi, sep="\t")
    simi_df["x_index"] = (simi_df["x"]/reselution).astype(int)
    simi_df["y_index"] = (simi_df["y"]/reselution).astype(int)

    pred_df = pd.read_csv(f_pred, sep="\t")

    all_cell_type_li = sorted(set(list(simi_df["CellType"].unique()) + list(pred_df["LabelTransfer"].unique())))
    all_cell_type_dict = {ct: index for (index, ct) in enumerate(all_cell_type_li)}
    all_cell_type_dict["Background"] = -1

    pred_df["CellTypeIndex"] = [all_cell_type_dict[ct] for ct in pred_df["LabelTransfer"].to_numpy()]

    cell_type_index_arr = coo_matrix((pred_df["CellTypeIndex"].to_numpy()+1, (pred_df["x"].to_numpy(), pred_df["y"].to_numpy())), [x_range, y_range]).toarray()
    cell_type_index_arr = cell_type_index_arr - 1

    cell_center_cnt_dict = defaultdict(int)
    for simi_ct, x, y in simi_df[["CellType", "x_index", "y_index"]].to_numpy():
        if x >= x_range:
            continue
        if y >= y_range:
            continue
        pred_ct = all_cell_type_li[cell_type_index_arr[x, y]]
        cell_center_cnt_dict[(simi_ct, pred_ct)] += 1
    
    cell_center_df = pd.DataFrame(
        [(simi_ct, pred_ct, cnt) for (simi_ct, pred_ct), cnt in cell_center_cnt_dict.items()],
        columns=["SimiCellType", "PredCellType", "CellNum"]
        )
    SimiSumNum = cell_center_df.groupby(["SimiCellType"])["CellNum"].sum().reset_index(name="AllSimuCellNum")
    cell_center_df = cell_center_df.merge(SimiSumNum)
    cell_center_df["Ratio"] = cell_center_df["CellNum"] / cell_center_df["AllSimuCellNum"]
    cell_center_df.to_csv(os.path.join(f_out, "CellType.Center.cnt.tsv"), sep="\t", index_label=False, index=False)

    cell_center_df["SimiCellTypeIndex"] = [all_cell_type_dict[ct] for ct in cell_center_df["SimiCellType"].to_numpy()]
    cell_center_df["PredCellTypeIndex"] = [all_cell_type_dict[ct] for ct in cell_center_df["PredCellType"].to_numpy()]

    cell_center_df["xmin"] = cell_center_df["SimiCellTypeIndex"] - 0.5
    cell_center_df["xmax"] = cell_center_df["SimiCellTypeIndex"] + 0.5
    cell_center_df["ymin"] = cell_center_df["PredCellTypeIndex"] - 0.5
    cell_center_df["ymax"] = cell_center_df["PredCellTypeIndex"] + 0.5
    cell_center_df["Label"] = [f"{r*100:.0f}" for r in cell_center_df["Ratio"].to_numpy()]

    p = ggplot(cell_center_df, aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="Ratio", x="SimiCellTypeIndex", y="PredCellTypeIndex", label="Label")) + \
        geom_rect() + \
        geom_text(size=3) + \
        scale_x_continuous(breaks=list(range(len(all_cell_type_li))), labels=all_cell_type_li) + \
        scale_y_continuous(breaks=list(range(len(all_cell_type_li))), labels=all_cell_type_li) + \
        scale_fill_gradient(limits=(0, 1), low="white", high="red") + \
        labs(x="Simulted cell types", y="Predicted cell types", fill="%Simulted") + \
        theme_bw() + \
        theme(
            text = element_text(family="Arial", size=5),
            title = element_text(family="Arial", size=6),
            axis_text_x= element_text(angle=45, hjust=1),
            axis_text = element_text(color = "black"),
            legend_title = element_text(family = "Arial", size=7, color="black"),
            legend_text = element_text(family = "Arial", size=6, color="black"),
            legend_key_size = 3,
            panel_grid = element_blank()
        )
    ggsave(p, filename=os.path.join(f_out, "CellType.SimiPct.center.pdf"), width=10, height=10, limitsize=False, units="cm")

    cell_center_same_celltype_df = cell_center_df[cell_center_df["SimiCellType"]==cell_center_df["PredCellType"]]
    p = ggplot(cell_center_same_celltype_df, aes(x="SimiCellTypeIndex", y="Ratio", label="Label")) + \
        geom_bar(stat="identity", fill="black") + \
        geom_text(size=3) + \
        scale_x_continuous(breaks=list(range(len(all_cell_type_li))), labels=all_cell_type_li) + \
        scale_y_continuous(limits=(0, 1), breaks=(0, 0.5, 1), labels=("0%", "50%", "100%")) + \
        labs(x="Simulted cell types", y="Precision") + \
        theme_bw() + \
        theme(
            text = element_text(family="Arial", size=5),
            title = element_text(family="Arial", size=6),
            axis_text_x= element_text(angle=45, hjust=1),
            axis_text = element_text(color = "black"),
            legend_title = element_text(family = "Arial", size=7, color="black"),
            legend_text = element_text(family = "Arial", size=6, color="black"),
            legend_key_size = 3,
            panel_grid = element_blank()
        )
    ggsave(p, filename=os.path.join(f_out, "CellType.Precision.center.pdf"), width=10, height=3, limitsize=False, units="cm")

    p = ggplot(cell_center_same_celltype_df, aes(x="AllSimuCellNum", y="Ratio", label="Label")) + \
        geom_point(size=0.5) + \
        scale_y_continuous(limits=(0, 1), breaks=(0, 0.5, 1), labels=("0%", "50%", "100%")) + \
        scale_x_log10() + \
        labs(x="#Cell", y="Precision") + \
        theme_bw() + \
        theme(
            text = element_text(family="Arial", size=5),
            title = element_text(family="Arial", size=6),
            axis_text_x= element_text(angle=45, hjust=1),
            axis_text = element_text(color = "black"),
            legend_title = element_text(family = "Arial", size=7, color="black"),
            legend_text = element_text(family = "Arial", size=6, color="black"),
            legend_key_size = 3,
            panel_grid = element_blank()
        )
    ggsave(p, filename=os.path.join(f_out, "CellType.Precision.CellNum.center.pdf"), width=5, height=5, limitsize=False, units="cm")

    
    simi_cell_type_index_arr = -1 * np.ones_like(cell_type_index_arr)
    for simi_ct, x, y, poly_str in simi_df[["CellType", "x_index", "y_index", "CellShape"]].to_numpy():
        poly = wkt.loads(poly_str)
        xmin, ymin, xmax, ymax = poly.bounds
        prep_poly = prep(poly)
        simi_ct_index = all_cell_type_dict[simi_ct]
        for x_index in range(int(xmin/reselution), int(xmax/reselution)):
            if x_index >= x_range:
                continue
            for y_index in range(int(ymin/reselution), int(ymax/reselution)):
                if y_index >= y_range:
                    continue
                point = Point(x_index*reselution+reselution/2, y_index*reselution+reselution/2)
                if prep_poly.contains(point):
                    simi_cell_type_index_arr[x_index, y_index] = simi_ct_index

    correct_pos_arr = cell_type_index_arr == simi_cell_type_index_arr
    correct_R_arr = np.zeros_like(correct_pos_arr, dtype=np.uint8)
    correct_G_arr = np.zeros_like(correct_pos_arr, dtype=np.uint8)
    correct_B_arr = 255 * np.ones_like(correct_pos_arr, dtype=np.uint8)
    correct_R_arr[correct_pos_arr] = 255
    correct_B_arr[correct_pos_arr] = 0
    correct_RGB_arr = np.stack([correct_R_arr, correct_G_arr, correct_B_arr], -1)

    plt.imsave(os.path.join(f_out, "CellType.AccPos.png"), correct_RGB_arr)

    flat_simi_cell_type_index = np.reshape(simi_cell_type_index_arr, -1)
    flat_pred_cell_type_index = np.reshape(cell_type_index_arr, -1)
    spot_cell_type_df = pd.DataFrame({"SimiCellTypeIndex": flat_simi_cell_type_index, "PredCellTypeIndex": flat_pred_cell_type_index})
    spot_cell_type_df["SimiCellType"] = np.array(all_cell_type_li)[spot_cell_type_df["SimiCellTypeIndex"].to_numpy()]
    spot_cell_type_df["PredCellType"] = np.array(all_cell_type_li)[spot_cell_type_df["PredCellTypeIndex"].to_numpy()]
    spot_cell_type_df["SimiCellType"][spot_cell_type_df["SimiCellTypeIndex"]==-1] = "Background"
    spot_cell_type_df["PredCellType"][spot_cell_type_df["PredCellTypeIndex"]==-1] = "Background"
    spot_cell_type_cnt_df = spot_cell_type_df.groupby(["SimiCellType", "PredCellType"]).size().reset_index(name="SpotNum")
    SimiSpotSumNum = spot_cell_type_cnt_df.groupby(["SimiCellType"])["SpotNum"].sum().reset_index(name="AllSimuSpotNum")
    spot_cell_type_cnt_df = spot_cell_type_cnt_df.merge(SimiSpotSumNum)
    spot_cell_type_cnt_df["Ratio"] = spot_cell_type_cnt_df["SpotNum"] / spot_cell_type_cnt_df["AllSimuSpotNum"]
    spot_cell_type_cnt_df.to_csv(os.path.join(f_out, "CellType.Spot.cnt.tsv"), sep="\t", index_label=False, index=False)

    spot_cell_type_cnt_df["SimiCellTypeIndex"] = [all_cell_type_dict[ct] for ct in spot_cell_type_cnt_df["SimiCellType"].to_numpy()]
    spot_cell_type_cnt_df["PredCellTypeIndex"] = [all_cell_type_dict[ct] for ct in spot_cell_type_cnt_df["PredCellType"].to_numpy()]

    spot_cell_type_cnt_df["xmin"] = spot_cell_type_cnt_df["SimiCellTypeIndex"] - 0.5
    spot_cell_type_cnt_df["xmax"] = spot_cell_type_cnt_df["SimiCellTypeIndex"] + 0.5
    spot_cell_type_cnt_df["ymin"] = spot_cell_type_cnt_df["PredCellTypeIndex"] - 0.5
    spot_cell_type_cnt_df["ymax"] = spot_cell_type_cnt_df["PredCellTypeIndex"] + 0.5
    spot_cell_type_cnt_df["Label"] = [f"{r*100:.0f}" for r in spot_cell_type_cnt_df["Ratio"].to_numpy()]

    p = ggplot(spot_cell_type_cnt_df, aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="Ratio", x="SimiCellTypeIndex", y="PredCellTypeIndex", label="Label")) + \
        geom_rect() + \
        geom_text(size=3) + \
        scale_x_continuous(breaks=list(range(len(all_cell_type_li))), labels=all_cell_type_li) + \
        scale_y_continuous(breaks=list(range(len(all_cell_type_li))), labels=all_cell_type_li) + \
        scale_fill_gradient(limits=(0, 1), low="white", high="red") + \
        labs(x="Simulted cell types", y="Predicted cell types", fill="%Simulted") + \
        theme_bw() + \
        theme(
            text = element_text(family="Arial", size=5),
            title = element_text(family="Arial", size=6),
            axis_text_x= element_text(angle=45, hjust=1),
            axis_text = element_text(color = "black"),
            legend_title = element_text(family = "Arial", size=7, color="black"),
            legend_text = element_text(family = "Arial", size=6, color="black"),
            legend_key_size = 3,
            panel_grid = element_blank()
        )
    ggsave(p, filename=os.path.join(f_out, "CellType.SimiPct.spot.pdf"), width=10, height=10, limitsize=False, units="cm")

    spot_cell_type_same_celltype_df = spot_cell_type_cnt_df[spot_cell_type_cnt_df["SimiCellType"]==spot_cell_type_cnt_df["PredCellType"]]
    p = ggplot(spot_cell_type_same_celltype_df, aes(x="SimiCellTypeIndex", y="Ratio", label="Label")) + \
        geom_bar(stat="identity", fill="black") + \
        geom_text(size=3) + \
        scale_x_continuous(breaks=list(range(len(all_cell_type_li))), labels=all_cell_type_li) + \
        scale_y_continuous(limits=(0, 1), breaks=(0, 0.5, 1), labels=("0%", "50%", "100%")) + \
        labs(x="Simulted cell types", y="Precision") + \
        theme_bw() + \
        theme(
            text = element_text(family="Arial", size=5),
            title = element_text(family="Arial", size=6),
            axis_text_x= element_text(angle=45, hjust=1),
            axis_text = element_text(color = "black"),
            legend_title = element_text(family = "Arial", size=7, color="black"),
            legend_text = element_text(family = "Arial", size=6, color="black"),
            legend_key_size = 3,
            panel_grid = element_blank()
        )
    ggsave(p, filename=os.path.join(f_out, "CellType.Precision.spot.pdf"), width=10, height=3, limitsize=False, units="cm")

    nbg_spot_cell_type_same_celltype_df = spot_cell_type_same_celltype_df[spot_cell_type_same_celltype_df["SimiCellType"]!="Background"]
    p = ggplot(nbg_spot_cell_type_same_celltype_df, aes(x="AllSimuSpotNum", y="Ratio", label="Label")) + \
        geom_point(size=0.5) + \
        scale_y_continuous(limits=(0, 1), breaks=(0, 0.5, 1), labels=("0%", "50%", "100%")) + \
        scale_x_log10() + \
        labs(x="#Spot", y="Precision") + \
        theme_bw() + \
        theme(
            text = element_text(family="Arial", size=5),
            title = element_text(family="Arial", size=6),
            axis_text_x= element_text(angle=45, hjust=1),
            axis_text = element_text(color = "black"),
            legend_title = element_text(family = "Arial", size=7, color="black"),
            legend_text = element_text(family = "Arial", size=6, color="black"),
            legend_key_size = 3,
            panel_grid = element_blank()
        )
    ggsave(p, filename=os.path.join(f_out, "CellType.Precision.CellNum.spot.pdf"), width=5, height=5, limitsize=False, units="cm")

    cell_same_ct_num = cell_center_same_celltype_df["CellNum"].sum()
    cell_total_num = cell_center_df["CellNum"].sum()
    cell_bg_num = cell_center_df["CellNum"][cell_center_df["PredCellType"]=="Background"].sum()
    cell_diff_ct_num = cell_total_num - cell_same_ct_num - cell_bg_num

    nbg_spot_cell_type_cnt_df = spot_cell_type_cnt_df[spot_cell_type_cnt_df["SimiCellType"]!="Background"]
    spot_same_ct_num = nbg_spot_cell_type_same_celltype_df["SpotNum"].sum()
    spot_total_num = nbg_spot_cell_type_cnt_df["SpotNum"].sum()
    spot_bg_num = nbg_spot_cell_type_cnt_df["SpotNum"][nbg_spot_cell_type_cnt_df["PredCellType"]=="Background"].sum()
    spot_diff_ct_num = spot_total_num - spot_same_ct_num - spot_bg_num

    stat_df = pd.DataFrame({
        "Type": ["Same cell type", "Diff cell type", "Background", "Same cell type", "Diff cell type", "Background"],
        "Stat": ["Cell Center", "Cell Center", "Cell Center", "Spot", "Spot", "Spot"],
        "Number": [cell_same_ct_num, cell_diff_ct_num, cell_bg_num, spot_same_ct_num, spot_bg_num, spot_diff_ct_num],
        "Ratio": [cell_same_ct_num/cell_total_num, cell_diff_ct_num/cell_total_num, cell_bg_num/cell_total_num, spot_same_ct_num/spot_total_num, spot_bg_num/spot_total_num, spot_diff_ct_num/spot_total_num]
    })
    stat_df.to_csv(os.path.join(f_out, "Stat.cnt.tsv"), sep="\t", index_label=False, index=False)

    p = ggplot(stat_df, aes(x="Stat", y="Ratio", fill="Type")) + \
        geom_bar(stat="identity") + \
        scale_y_continuous(limits=(0, 1), breaks=(0, 0.5, 1), labels=("0%", "50%", "100%")) + \
        labs(y="Ratio") + \
        theme_bw() + \
        theme(
            text = element_text(family="Arial", size=5),
            title = element_text(family="Arial", size=6),
            axis_text = element_text(color = "black"),
            legend_title = element_text(family = "Arial", size=7, color="black"),
            legend_text = element_text(family = "Arial", size=6, color="black"),
            legend_key_size = 3,
            panel_grid = element_blank()
        )
    ggsave(p, filename=os.path.join(f_out, "Stat.cnt.pdf"), width=5, height=5, limitsize=False, units="cm")

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
