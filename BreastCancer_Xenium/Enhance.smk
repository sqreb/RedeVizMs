import os

embedding_dict = {
    "tSNE": {"embedding": "tSNE", "embedding_dim": 2},
    "UMAP2": {"embedding": "UMAP", "embedding_dim": 2},
    "UMAP3": {"embedding": "UMAP", "embedding_dim": 3}
    }

REF_SAMPLE_li = ["3_scRNA", "5_scRNA", "Xenium"]

domain_params = {
    "small": {"smooth_radius": 0, "receptive_radius": 24, "emb_radius": 5, "merge_emb_dist": 5, "min_spot_num": 300},
    "large": {"smooth_radius": 50, "receptive_radius": 100, "emb_radius": 10, "merge_emb_dist": 20, "min_spot_num": 4000}
}

RGB_gene_dict = {
    "CD3E_None_VWF": ("CD3E", None, "VWF"),
    "CD40LG_None_FLT4": ("CD40LG", None, "FLT4") 
    }

rule scRNA_select_gene:
    input:
        h5ad = os.path.join("data", "ref_scRNA", "{ref_sample}.cnt.h5ad"),
        gene_li = os.path.join("data", "Xenium.gene.txt"),
    output:
        h5ad = os.path.join("analysis", "ref_scRNA", "{ref_sample}.selected.cnt.h5ad"),
    params:
        scRNA_select_gene = "python scripts/scRNA_select_gene.py",
        gene_name_label = "GeneSymbol",
    shell:
        """
{params.scRNA_select_gene} \
    --input {input.h5ad} \
    --gene-list {input.gene_li} \
    --gene-name-label {params.gene_name_label} \
    --output {output.h5ad}
        """

rule filter_scRNA:
    input:
        h5ad = rules.scRNA_select_gene.output.h5ad
    output:
        h5ad = os.path.join("analysis", "ref_scRNA", "{ref_sample}.filtered.cnt.h5ad"),
    params:
        total_UMI = 100,
        expr_gene_num = 20,
    threads: 1
    shell:
        """
RedeViz pretreatment filter_scRNA \
    --input {input.h5ad} \
    --output {output.h5ad} \
    --min-expr-gene-per-cell {params.expr_gene_num} \
    --min-UMI-per-cell {params.total_UMI}
        """

rule RedeViz_pretreat:
    input:
        h5ad = rules.filter_scRNA.output.h5ad,
    output:
        index = os.path.join("analysis", "RedeViz_embedding_info", "{embedding}", "{ref_sample}", "pretreat.pkl"),
        emb_info = os.path.join("analysis", "RedeViz_embedding_info", "{embedding}", "{ref_sample}", "Embedding.bin.info.tsv"),
    params:
        min_cell_num = 2,
        max_expr_ratio = 0.2,
        cell_type_label = "SpaSegClu",
        gene_id_label = "GeneSymbol",
        out_dir = os.path.join("analysis", "RedeViz_embedding_info", "{embedding}", "{ref_sample}"),
        embedding = lambda wildcards: embedding_dict[wildcards.embedding]["embedding"],
        embedding_dim = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_dim"],
    threads: 28
    shell:
        """
RedeViz enhance pretreat \
    -i {input.h5ad} \
    --cell-type-label {params.cell_type_label} \
    --gene-id-label {params.gene_id_label} \
    --embedding {params.embedding} \
    --embedding-dim {params.embedding_dim} \
    --min-cell-num {params.min_cell_num} \
    --max-expr-ratio {params.max_expr_ratio} \
    -o {params.out_dir}
        """

rule plot_embedding_color:
    input:
        index = rules.RedeViz_pretreat.output.index,
        color = os.path.join("data", "celltype.color.tsv"),
    output:
        embedding_plot = os.path.join("analysis", "embedding_color", "{embedding}", "{ref_sample}.EmbeddingColor.png"),
        cell_type_plot = os.path.join("analysis", "embedding_color", "{embedding}", "{ref_sample}.CellTypeColor.png"),
    params:
        plot_embedding_color = "python scripts/plot_embedding_color.py",
    shell:
        """
{params.plot_embedding_color} --index {input.index} --color {input.color} --embedding-plot {output.embedding_plot} --cell-type-plot {output.cell_type_plot}
        """

rule point2tsv:
    input:
        point = os.path.join("data", "ST", "Xenium_FFPE_Human_Breast_Cancer_Rep1_transcripts.csv"),
    output:
        point = os.path.join("analysis", "ST", "ST.points.tsv"),
    run:
        import pandas as pd
        df = pd.read_csv(input.point)
        df.to_csv(output.point, sep="\t", index_label=False, index=False)

rule point2bin:
    input:
        point = rules.point2tsv.output.point,
    output:
        res = os.path.join("analysis", "ST", "ST.bin.tsv"),
    params:
        x = "x_location",
        y = "y_location",
        gene_name_label = "feature_name",
        bin_size = 0.5,
    shell:
        """
RedeViz pretreatment point2bin \
    --point {input.point} \
    --bin-size {params.bin_size} \
    --x-label {params.x} \
    --y-label {params.y} \
    --gene-name-label {params.gene_name_label} \
    --output {output.res}
        """

rule point2bin50:
    input:
        point = rules.point2tsv.output.point,
    output:
        res = os.path.join("analysis", "ST", "ST.bin50.tsv"),
    params:
        x = "x_location",
        y = "y_location",
        gene_name_label = "feature_name",
        bin_size = 25,
    shell:
        """
RedeViz pretreatment point2bin \
    --point {input.point} \
    --bin-size {params.bin_size} \
    --x-label {params.x} \
    --y-label {params.y} \
    --gene-name-label {params.gene_name_label} \
    --output {output.res}
        """

rule infer_spot_expr:
    input:
        spot = rules.point2bin.output.res,
    output:
        cov = os.path.join("analysis", "ST", "ST.cov.txt"),
    params:
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        UMI_label = "UMI",
    shell:
        """
RedeViz pretreatment infer_spot_expr \
    --spot {input.spot} \
    --x-index-label {params.x_index_label} \
    --y-index-label {params.y_index_label} \
    --UMI-label {params.UMI_label} \
    --output {output.cov}
        """

rule RedeViz_build_enhance_index:
    input:
        index = rules.RedeViz_pretreat.output.index,
        cov = rules.infer_spot_expr.output.cov
    output:
        index = os.path.join("analysis", "RedeViz_enhance_index", "{embedding}", "{ref_sample}.index.pkl"),
    threads: 28
    shell:
        """
RedeViz enhance build \
    -i {input.index} \
    --UMI-per-spot $(cat {input.cov} | awk '{{print $1*0.75}}') \
    -o {output.index}
        """

rule RedeViz_enhance:
    input:
        index = rules.RedeViz_build_enhance_index.output.index,
        spot = rules.infer_spot_expr.input.spot,
    output:
        enhance = os.path.join("analysis", "RedeViz_enhance", "{embedding}", "{ref_sample}.enhance.tsv"),
    params:
        ave_bin_dist_cutoff = 4,
        max_expr_ratio = 0.2,
        neighbor_close_label_fct = 0.2,
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        gene_name_label = "feature_name",
        UMI_label = "UMI",
    threads: 28
    shell:
        """
RedeViz enhance run \
    -i {input.index} \
    -s {input.spot} \
    -o {output.enhance} \
    --x-index-label {params.x_index_label} \
    --y-index-label {params.y_index_label} \
    --gene-name-label {params.gene_name_label} \
    --UMI-label {params.UMI_label} \
    --ave-bin-dist-cutoff {params.ave_bin_dist_cutoff} \
    --max-expr-ratio {params.max_expr_ratio} \
    --neighbor-close-label-fct {params.neighbor_close_label_fct} \
    --window-size 100"""


rule RedeViz_phenotype_visualization:
    input:
        enhance = rules.RedeViz_enhance.output.enhance,
    output:
        img = os.path.join("analysis", "RedeViz_enhance", "{embedding}", "{ref_sample}.enhance.phenotype.{QC}.png"),
    params:
        other_params = lambda wildcards: {"HQ": "", "all": "--keep-other"}[wildcards.QC]
    shell:
        """
RedeViz posttreatment plot_phenotype \
    --input {input.enhance} \
    --output {output.img} {params.other_params} --denoise
        """

rule RedeViz_celltype_visualization:
    input:
        enhance = rules.RedeViz_enhance.output.enhance,
        color = os.path.join("data", "celltype.color.tsv")
    output:
        img = os.path.join("analysis", "RedeViz_enhance", "{embedding}", "{ref_sample}.enhance.celltype.{QC}.png"),
    params:
        other_params = lambda wildcards: {"HQ": "", "all": "--keep-other"}[wildcards.QC]
    shell:
        """
RedeViz posttreatment plot_cell_type \
    --color {input.color} \
    --input {input.enhance} \
    --output {output.img} {params.other_params} --denoise
        """

rule RedeViz_segment:
    input:
        enhance = rules.RedeViz_enhance.output.enhance,
    output:
        res = os.path.join("analysis", "RedeViz_segment", "{embedding}", "{ref_sample}", "{domain_param}_{QC}_domain", "domain.label.tsv"),
    params:
        smooth_radius = lambda wildcards: domain_params[wildcards.domain_param]["smooth_radius"],
        receptive_radius = lambda wildcards: domain_params[wildcards.domain_param]["receptive_radius"],
        emb_radius = lambda wildcards: domain_params[wildcards.domain_param]["emb_radius"],
        merge_emb_dist = lambda wildcards: domain_params[wildcards.domain_param]["merge_emb_dist"],
        min_spot_num = lambda wildcards: domain_params[wildcards.domain_param]["min_spot_num"],
        out_dir = os.path.join("analysis", "RedeViz_segment", "{embedding}", "{ref_sample}", "{domain_param}_{QC}_domain"),
        other_params = lambda wildcards: {"HQ": "", "all": "--keep-other"}[wildcards.QC]
    threads: 10
    shell:
        """
RedeViz posttreatment segment \
    --input {input.enhance} \
    --smooth-radius {params.smooth_radius} \
    --receptive-radius {params.receptive_radius} \
    --embedding-radius {params.emb_radius} \
    --merge-embedding-dist {params.merge_emb_dist} \
    --min-spot-per-domain {params.min_spot_num} \
    --output {params.out_dir} {params.other_params} --denoise
        """


rule build_RedeViz_imputation_index:
    input:
        index = rules.RedeViz_pretreat.output.index,
        sce = rules.scRNA_select_gene.input.h5ad,
    output:
        index = os.path.join("analysis", "RedeViz_imputation_index", "{embedding}", "{ref_sample}.pkl"),
    params:
        gene_id_label = "GeneSymbol",
    shell:
        """
RedeViz posttreatment impute build \
    --index {input.index} \
    --sce {input.sce} \
    --gene-name-label {params.gene_id_label} \
    --output {output.index}
        """

rule RedeViz_imputation:
    input:
        enhance = rules.RedeViz_enhance.output.enhance,
        spot = rules.infer_spot_expr.input.spot,
        index = rules.build_RedeViz_imputation_index.output,
        gene_li = os.path.join("data", "gene_li.txt")
    output:
        touch(os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{ref_sample}", "done.flag")),
    params:
        x_label = "bin_x_index",
        y_label = "bin_y_index",
        UMI_label = "UMI",
        embedding_smooth_sigma = 3.0,
        UMI_smooth_sigma = 10.0,
        out_dir = os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{ref_sample}"),
    threads: 10
    shell:
        """
RedeViz posttreatment impute run \
    --input {input.enhance} \
    --index {input.index} \
    --gene-list {input.gene_li} \
    --spot {input.spot} \
    --spot-pos-pos-label {params.x_label} {params.y_label} \
    --spot-UMI-label {params.UMI_label} \
    --embedding-smooth-sigma {params.embedding_smooth_sigma} \
    --UMI-smooth-sigma {params.UMI_smooth_sigma} \
    --output {params.out_dir} --denoise
        """

rule plot_RedeViz_imputation_gene_expr_RB:
    input:
        flag = rules.RedeViz_imputation.output,
    output:
        img = os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{ref_sample}", "{RGB}.png")
    params:
        data_dir = rules.RedeViz_imputation.params.out_dir,
        R = lambda wildcards: RGB_gene_dict[wildcards.RGB][0],
        G = lambda wildcards: RGB_gene_dict[wildcards.RGB][1],
        B = lambda wildcards: RGB_gene_dict[wildcards.RGB][2]
    shell:
        """
RedeViz posttreatment plot_gene_expr \
    -R {params.data_dir}/{params.R}.imputate.npz \
    -B {params.data_dir}/{params.B}.imputate.npz \
    -o {output.img}
        """

rule plot_bin_expr_RB:
    input:
        signal = os.path.join("analysis", "ST", "ST.{binN}.tsv"),
    output:
        img = os.path.join("analysis", "bin_expr", "{binN}", "{RGB}.png"),
    params:
        x = "bin_x_index",
        y = "bin_y_index",
        gene_name_label = "feature_name",
        UMI = "UMI",
        R = lambda wildcards: RGB_gene_dict[wildcards.RGB][0],
        G = lambda wildcards: RGB_gene_dict[wildcards.RGB][1],
        B = lambda wildcards: RGB_gene_dict[wildcards.RGB][2]
    shell:
        """
RedeViz posttreatment plot_gene_bin_expr \
    --input {input.signal} \
    --gene-name-label {params.gene_name_label} \
    --UMI-label {params.UMI} \
    --x-label {params.x} \
    --y-label {params.y} \
    -R {params.R} \
    -B {params.B} \
    --output {output.img}
        """

rule all:
    input:
        expand(rules.plot_embedding_color.output, embedding=embedding_dict.keys(), ref_sample=REF_SAMPLE_li),
        expand(rules.RedeViz_phenotype_visualization.output, embedding=embedding_dict.keys(), ref_sample=REF_SAMPLE_li, QC=["HQ", "all"]),
        expand(rules.RedeViz_celltype_visualization.output, embedding=embedding_dict.keys(), ref_sample=REF_SAMPLE_li, QC=["HQ", "all"]),
        expand(rules.RedeViz_segment.output, embedding=embedding_dict.keys(), ref_sample=REF_SAMPLE_li, QC=["HQ", "all"], domain_param=domain_params.keys()),
        expand(rules.plot_RedeViz_imputation_gene_expr_RB.output, embedding=embedding_dict.keys(), ref_sample=REF_SAMPLE_li, RGB=RGB_gene_dict.keys()),
        expand(rules.plot_bin_expr_RB.output, embedding=embedding_dict.keys(), ref_sample=REF_SAMPLE_li, binN=["bin50", "bin"], RGB=["CD3E_None_VWF"]),
