import os

embedding_dict = {
    "UMAP2": {"embedding": "UMAP", "embedding_dim": 2},
    "UMAP3": {"embedding": "UMAP", "embedding_dim": 3},
    "tSNE": {"embedding": "tSNE", "embedding_dim": 2}
    }


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


rule mask_ST_data:
    input:
        spot = os.path.join("analysis", "ST", "ST.bin.tsv"),
    output:
        mask = os.path.join("analysis", "mask_ST_data", "mask.npz"),
    params:
        mask_model = "bin",
        cell_diameter = 25,
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        UMI_label = "UMI",
        smooth_sd = 4,
        out_dir = os.path.join("analysis", "mask_ST_data")
    shell:
        """
RedeViz pretreatment mask_ST_data \
    --spot {input.spot} \
    --x-index-label {params.x_index_label} \
    --y-index-label {params.y_index_label} \
    --UMI-label {params.UMI_label} \
    --mask-model {params.mask_model} \
    --smooth-sd {params.smooth_sd} \
    --cell-diameter {params.cell_diameter} \
    --output {params.out_dir}
        """


rule RedeViz_embedding_pretreat:
    input:
        spot = rules.mask_ST_data.input.spot,
        mask = rules.mask_ST_data.output.mask,
        black_li = os.path.join("data", "black_list.txt"),
    output:
        index = os.path.join("analysis", "self_embedding_info", "{embedding}", "pretreat.pkl"),
        emb_info = os.path.join("analysis", "self_embedding_info", "{embedding}", "Embedding.bin.info.tsv"),
    params:
        min_cell_num = 2,
        gene_id_label = "feature_name",
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        UMI_label = "UMI",
        embedding = lambda wildcards: embedding_dict[wildcards.embedding]["embedding"],
        embedding_dim = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_dim"],
        out_dir = os.path.join("analysis", "self_embedding_info", "{embedding}"),
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
    threads: 1
    shell:
        """
RedeViz enhance pretreat_by_ST \
    --spot {input.spot} \
    --mask {input.mask} \
    --gene-blacklist {input.black_li} \
    --x-index-label {params.x_index_label} \
    --y-index-label {params.y_index_label} \
    --UMI-label {params.UMI_label} \
    --gene-id-label {params.gene_id_label} \
    --embedding {params.embedding} \
    --embedding-dim {params.embedding_dim} \
    --min-cell-num {params.min_cell_num} \
    --leiden-resolution 0.35 \
    -o {params.out_dir} \
    --no-HVG
        """


rule cell_type2color:
    input:
        emb_info = rules.RedeViz_embedding_pretreat.output.emb_info,
    output:
        color = os.path.join("analysis", "self_embedding_info", "{embedding}", "color.tsv"),
    params:
        cell_type2color = "python scripts/cell_type2color.py",
    shell:
        """
{params.cell_type2color} \
    --emb-info {input.emb_info} \
    --output {output.color}
        """

rule plot_embedding_color:
    input:
        index = rules.RedeViz_embedding_pretreat.output.index,
        color = rules.cell_type2color.output.color,
    output:
        embedding_plot = os.path.join("analysis", "self_embedding_color", "{embedding}.EmbeddingColor.png"),
        cell_type_plot = os.path.join("analysis", "self_embedding_color", "{embedding}.CellTypeColor.png"),
    params:
        plot_embedding_color = "python scripts/plot_embedding_color.py",
    shell:
        """
{params.plot_embedding_color} --index {input.index} --color {input.color} --embedding-plot {output.embedding_plot} --cell-type-plot {output.cell_type_plot}
        """

rule run_plot_embedding_color:
    input:
        expand(rules.plot_embedding_color.output, embedding=["tSNE", "UMAP2", "UMAP3"])

rule RedeViz_embedding_build_index:
    input:
        index = rules.RedeViz_embedding_pretreat.output.index,
        cov = os.path.join("analysis", "ST", "ST.cov.txt"),
    output:
        index = os.path.join("analysis", "RedeViz_self_embedding_index", "{embedding}.index.pkl"),
    params:
        bin_size = [12, 18, 21],
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


rule RedeViz_enhance:
    input:
        index = rules.RedeViz_embedding_build_index.output.index,
        spot = rules.mask_ST_data.input.spot,
        cov = os.path.join("analysis", "ST", "ST.cov.txt"),
    output:
        res = os.path.join("analysis", "RedeViz_self_enhance", "{embedding}.predict.tsv"),
    params:
        cell_radius = 11,
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        UMI_label = "UMI",
        gene_name_label = "feature_name",
        queue = "gpu1,gpu2",
        gpu_gres = "--gres=gpu:1"
    threads: 28
    shell:
        """
RedeViz enhance run \
    -i {input.index} -s {input.spot} -o {output.res} \
    --x-index-label {params.x_index_label} --y-index-label {params.y_index_label} \
    --gene-name-label {params.gene_name_label} --UMI-label {params.UMI_label} \
    --cell-radius {params.cell_radius} \
    --window-size 75 --mid-signal-cutoff $(cat {input.cov} | awk '{{print $1 * 0.5}}')"""


rule RedeViz_phenotype_visualization:
    input:
        enhance = rules.RedeViz_enhance.output.res,
    output:
        all_img = os.path.join("analysis", "RedeViz_self_enhance_pheno_plog", "{embedding}.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_self_enhance_pheno_plog", "{embedding}.enhance.phenotype.HQ.png"),
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

rule RedeViz_celltype_visualization:
    input:
        enhance = rules.RedeViz_enhance.output.res,
        color = rules.cell_type2color.output.color,
    output:
        all_img = os.path.join("analysis", "RedeViz_self_enhance_celltype_plog", "{embedding}.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_self_enhance_celltype_plog", "{embedding}.enhance.phenotype.HQ.png"),
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


rule all:
    input:
        [rules.RedeViz_phenotype_visualization.output.all_img.format(embedding=embedding) for embedding in embedding_dict.keys()],
        [rules.RedeViz_celltype_visualization.output.all_img.format(embedding=embedding) for embedding in embedding_dict.keys()],
