import os

embedding_dict = {
    "UMAP2": {"embedding": "UMAP", "embedding_dim": 2},
    "UMAP3": {"embedding": "UMAP", "embedding_dim": 3},
    "tSNE": {"embedding": "tSNE", "embedding_dim": 2}
    }

rule point2bin:
    input:
        point = os.path.join("data", "ST", "points.tsv"),
    output:
        res = os.path.join("analysis", "ST_bin", "ST.bin0.5.tsv"),
    params:
        x = "x_um",
        y = "y_um",
        gene_name_label = "gene",
        bin_size = 0.5,
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
    threads: 1
    shell:
        """
RedeViz pretreatment point2bin\
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
        cov = os.path.join("analysis", "ST_bin", "ST.bin0.5.cov.txt"),
    params:
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        UMI_label = "UMI",
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
    threads: 1
    shell:
        """
RedeViz pretreatment infer_spot_expr \
    --spot {input.spot} \
    --x-index-label {params.x_index_label} \
    --y-index-label {params.y_index_label} \
    --UMI-label {params.UMI_label} \
    --output {output.cov}
        """


rule RedeViz_embedding_pretreat:
    input:
        h5ad = os.path.join("data", "ref_scRNA", "sc_data.SCP1038.new_label.cnt.h5ad"),
    output:
        index = os.path.join("analysis", "embedding_info", "{embedding}", "SCP1038", "pretreat.pkl"),
    params:
        min_cell_num = 2,
        cell_type_label = "annotation",
        embedding = lambda wildcards: embedding_dict[wildcards.embedding]["embedding"],
        embedding_dim = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_dim"],
        out_dir = os.path.join("analysis", "embedding_info", "{embedding}", "SCP1038"),
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
    threads: 1
    shell:
        """
RedeViz enhance pretreat -i {input.h5ad} \
    --cell-type-label {params.cell_type_label} \
    --n-neighbors 30 \
    --embedding {params.embedding} \
    --embedding-dim {params.embedding_dim} \
    --min-cell-num {params.min_cell_num} \
    -o {params.out_dir}
        """


rule emb_info_select_gene:
    input:
        index = os.path.join("analysis", "embedding_info", "{embedding}", "SCP1038", "pretreat.pkl"),
        gene_li = os.path.join("data", "ST", "gene_list.txt"),
    output:
        index = os.path.join("analysis", "RedeViz_NGS_embedding_info", "{embedding}", "SCP1038", "pretreat.pkl"),
    params:
        min_gene_num = 40,
        min_UMI_num = 200,
    shell:
        """
RedeViz pretreatment emb_info_select_gene \
    --input {input.index} \
    --gene-list {input.gene_li} \
    --min-gene-per-bin {params.min_gene_num} \
    --min-UMI-per-bin {params.min_UMI_num} \
    --output {output.index}
        """

    
rule RedeViz_embedding_build_NGS_index:
    input:
        index = rules.emb_info_select_gene.output.index,
        cov = os.path.join("analysis", "ST_bin", "ST.bin0.5.cov.txt"),
    output:
        index = os.path.join("analysis", "RedeViz_embedding_NGS_index", "{embedding}", "SCP1038.index.pkl"),
    params:
        bin_size = [12, 21],
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
    threads: 28
    shell:
        """
RedeViz enhance build -i {input.index} \
    --bin-size {params.bin_size} \
    --UMI-per-spot $(cat {input.cov} | awk '{{print $1 * 0.75}}') \
    -o {output.index}
        """

rule RedeViz_NGS_enhance:
    input:
        index = rules.RedeViz_embedding_build_NGS_index.output.index,
        spot = os.path.join("analysis", "ST_bin", "ST.bin0.5.tsv"),
    output:
        res = os.path.join("analysis", "RedeViz_NGS_enhance", "{embedding}", "SCP1038.predict.tsv"),
    params:
        cell_radius = 11,
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        UMI_label = "UMI",
        cov = 0.35,
        gene_name_label = "gene",
        queue = "gpu1,gpu2",
        gpu_gres = "--gres=gpu:1"
    threads: 28
    shell:
        """
RedeViz enhance run \
    -i {input.index} -s {input.spot} -o {output.res} \
    --x-index-label {params.x_index_label} --y-index-label {params.y_index_label} \
    --gene-name-label {params.gene_name_label} --UMI-label {params.UMI_label} \
    --cell-radius {params.cell_radius} --mid-signal-cutoff {params.cov} \
    --window-size 150"""

rule RedeViz_phenotype_NGS_visualization:
    input:
        enhance = rules.RedeViz_NGS_enhance.output.res,
    output:
        all_img = os.path.join("analysis", "RedeViz_NGS_enhance_pheno_plot", "{embedding}", "SCP1038.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_NGS_enhance_pheno_plot", "{embedding}", "SCP1038.enhance.phenotype.HQ.png"),
    shell:
        """
RedeViz posttreatment plot_phenotype \
    --input {input.enhance} \
    --output {output.all_img} \
    --keep-other --denoise

RedeViz posttreatment plot_phenotype \
    --input {input.enhance} \
    --output {output.HQ_img} --denoise
        """

rule RedeViz_celltype_NGS_visualization:
    input:
        enhance = rules.RedeViz_NGS_enhance.output.res,
        color = os.path.join("data", "ref_scRNA", "color.tsv")
    output:
        all_img = os.path.join("analysis", "RedeViz_NGS_enhance_celltype_plot", "{embedding}", "SCP1038.enhance.celltype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_NGS_enhance_celltype_plot", "{embedding}", "SCP1038.enhance.celltype.HQ.png"),
    shell:
        """
RedeViz posttreatment plot_cell_type \
    --color {input.color} \
    --input {input.enhance} \
    --output {output.all_img} \
    --keep-other --denoise

RedeViz posttreatment plot_cell_type \
    --color {input.color} \
    --input {input.enhance} \
    --output {output.HQ_img} --denoise
        """



rule build_RedeViz_imputation_index:
    input:
        index = os.path.join("analysis", "embedding_info", "{embedding}", "SCP1038", "pretreat.pkl"),
        sce = os.path.join("data", "ref_scRNA", "sc_data.SCP1038.new_label.cnt.h5ad"),
    output:
        index = os.path.join("analysis", "RedeViz_Ngs_imputation_index", "{embedding}", "SCP1038.pkl"),
    shell:
        """
RedeViz posttreatment impute build \
    --index {input.index} \
    --sce {input.sce} \
    --embedding-smooth-sigma 3 \
    --output {output.index}
        """

rule RedeViz_imputation:
    input:
        enhance = rules.RedeViz_NGS_enhance.output.res,
        spot = os.path.join("analysis", "ST_bin", "ST.bin0.5.tsv"),
        index = rules.build_RedeViz_imputation_index.output,
        gene_li = os.path.join("data", "impute_gene_li.txt"),
    output:
        flag = touch(os.path.join("analysis", "RedeViz_Ngs_imputation", "{embedding}", "SCP1038", "done.flag")),
    params:
        x_label = "bin_x_index",
        y_label = "bin_y_index",
        UMI_label = "UMI",
        out_dir = os.path.join("analysis", "RedeViz_Ngs_imputation", "{embedding}", "SCP1038"),
    shell:
        """
RedeViz posttreatment impute run \
    --input {input.enhance} \
    --index {input.index} \
    --spot {input.spot} \
    --gene-list {input.gene_li} \
    --spot-pos-pos-label {params.x_label} {params.y_label} \
    --spot-UMI-label {params.UMI_label} \
    --output {params.out_dir} --keep-other --denoise
        """

rule RedeViz_plot_imputation:
    input:
        flag = rules.RedeViz_imputation.output,
    output:
        flag = touch(os.path.join("analysis", "RedeViz_Ngs_imputation", "{embedding}", "SCP1038", "plot.flag")),
    params:
        out_dir = os.path.join("analysis", "RedeViz_Ngs_imputation", "{embedding}", "SCP1038"),
    shell:
        """
for f in $(ls {params.out_dir}/*.imputate.npz); do
    RedeViz posttreatment plot_gene_expr -R ${{f}} -G ${{f}} -o $(echo ${{f}} | sed 's/npz/png/g')
done
        """

rule plot_embedding_color:
    input:
        index = rules.emb_info_select_gene.output.index,
        color = os.path.join("data", "ref_scRNA", "color.tsv"),
    output:
        embedding_plot = os.path.join("analysis", "embedding_color", "{embedding}", "SCP1038.EmbeddingColor.png"),
        cell_type_plot = os.path.join("analysis", "embedding_color", "{embedding}", "SCP1038.CellTypeColor.png"),
    params:
        plot_embedding_color = "python scripts/plot_embedding_color.py",
    shell:
        """
{params.plot_embedding_color} --index {input.index} --color {input.color} --embedding-plot {output.embedding_plot} --cell-type-plot {output.cell_type_plot}
        """

rule all:
    input:
        [rules.RedeViz_phenotype_NGS_visualization.output.all_img.format(embedding=embedding) for embedding in embedding_dict.keys()],
        [rules.RedeViz_celltype_NGS_visualization.output.all_img.format(embedding=embedding) for embedding in embedding_dict.keys()],
        [rules.RedeViz_plot_imputation.output.flag.format(embedding=embedding) for embedding in embedding_dict.keys()],
        [rules.plot_embedding_color.output.embedding_plot.format(embedding=embedding) for embedding in embedding_dict.keys()],