library(ggplot2)
library(dplyr)

load("~/HPC/project/RedeViz_ms/HumanNSCLC_CoxMx/data/ST/SMI_Giotto_Object.RData")
cell_info_df <- gem@cell_metadata$rna
color_info <- read_delim("~/HPC/project/RedeViz_ms/HumanNSCLC_CoxMx/data/ST/color.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
color_li <- color_info$cell_type_RGB
names(color_li) <- color_info$cell_type

level1_color_df <- color_info[!duplicated(color_info$level1),]
level1_color_li <- level1_color_df$level1_RGB
names(level1_color_li) <- level1_color_df$level1

level0_color_df <- color_info[!duplicated(color_info$level0),]
level0_color_li <- level0_color_df$level0_RGB
names(level0_color_li) <- level0_color_df$level0

f_out <- "/home/wangdehe/HPC/project/RedeViz_ms/HumanNSCLC_CoxMx/analysis/ref_anno"

for(tissue_name in unique(cell_info_df$Run_Tissue_name)){
  tmp_cell_info_df <- cell_info_df[cell_info_df$Run_Tissue_name==tissue_name,]
  tmp_cell_info_df$cell_ID <- sapply(strsplit(tmp_cell_info_df$cell_ID, "_"), function(x){
    return(as.integer(x[4]))
  }) 
  tmp_cell_meta_df <- read_csv(sprintf("~/HPC/project/RedeViz_ms/HumanNSCLC_CoxMx/data/ST/%s/%s-Flat_files_and_images/%s_metadata_file.csv", tissue_name, tissue_name, tissue_name))
  tmp_cell_meta_df <- left_join(tmp_cell_meta_df, tmp_cell_info_df[, c("cell_ID", "fov", "cell_type", "niche")])
  tmp_cell_meta_df <- left_join(tmp_cell_meta_df, color_info)
  write_tsv(tmp_cell_meta_df, file.path(f_out, sprintf("%s.cellinfo.tsv", tissue_name)))
  
  tmp_plot_df <- tmp_cell_meta_df[!is.na(tmp_cell_meta_df$cell_type),]
  x_max <- max(tmp_plot_df$CenterX_global_px)
  x_min <- min(tmp_plot_df$CenterX_global_px)
  y_max <- max(tmp_plot_df$CenterY_global_px)
  y_min <- min(tmp_plot_df$CenterY_global_px)
  p <- ggplot(tmp_plot_df, aes(x=CenterX_global_px, y=CenterY_global_px, color=cell_type)) +
    geom_point(size=0.1) +
    scale_color_manual(values = color_li) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(
      text = element_text(family="ArialMT", size=5),
      title = element_text(family="ArialMT", size=5),
      axis.text = element_text(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill="black"),
      legend.title = element_text(family="ArialMT", size=5),
      legend.text = element_text(family="ArialMT", size=5),
      legend.key.size = unit(2, "mm")
    )
  ggsave(p, filename=file.path(f_out, sprintf("%s.celltype.pdf", tissue_name)), width=3+(x_max-x_min)/1000, height=1.5+(y_max-y_min)/1000, limitsize=FALSE, units="cm")
  
  p <- ggplot(tmp_plot_df, aes(x=CenterX_global_px, y=CenterY_global_px, color=level1)) +
    geom_point(size=0.1) +
    scale_color_manual(values = level1_color_li) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(
      text = element_text(family="ArialMT", size=5),
      title = element_text(family="ArialMT", size=5),
      axis.text = element_text(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill="black"),
      legend.title = element_text(family="ArialMT", size=5),
      legend.text = element_text(family="ArialMT", size=5),
      legend.key.size = unit(2, "mm")
    )
  ggsave(p, filename=file.path(f_out, sprintf("%s.level1.pdf", tissue_name)), width=3+(x_max-x_min)/1000, height=1.5+(y_max-y_min)/1000, limitsize=FALSE, units="cm")
  
  p <- ggplot(tmp_plot_df, aes(x=CenterX_global_px, y=CenterY_global_px, color=level0)) +
    geom_point(size=0.1) +
    scale_color_manual(values = level0_color_li) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(
      text = element_text(family="ArialMT", size=5),
      title = element_text(family="ArialMT", size=5),
      axis.text = element_text(color = "black"),
      panel.background = element_rect(fill="black"),
      panel.grid = element_blank(),
      legend.title = element_text(family="ArialMT", size=5),
      legend.text = element_text(family="ArialMT", size=5),
      legend.key.size = unit(2, "mm")
    )
  ggsave(p, filename=file.path(f_out, sprintf("%s.level0.pdf", tissue_name)), width=3+(x_max-x_min)/1000, height=1.5+(y_max-y_min)/1000, limitsize=FALSE, units="cm")
}
