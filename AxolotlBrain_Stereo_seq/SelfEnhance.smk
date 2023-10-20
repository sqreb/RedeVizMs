import os


ST_sample_li = ["30DPI", "60DPI"]

embedding_dict = {
    "tSNE": {"embedding": "tSNE", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP2": {"embedding": "UMAP", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP3": {"embedding": "UMAP", "embedding_dim": 3, "embedding_resolution": 4}
    }


rule infer_spot_expr:
    input:
        spot = os.path.join("data", "ST", "{sample}.gem"),
    output:
        cov = os.path.join("analysis", "cov", "{sample}.cov.txt"),
    params:
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = "MIDCounts",
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
        spot = os.path.join("data", "ST", "{sample}.gem"),
    output:
        mask = os.path.join("analysis", "mask_ST_data", "{sample}", "mask.npz"),
    params:
        mask_model = "bin",
        cell_diameter = 30,
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = "MIDCounts",
        smooth_sd = 6,
        out_dir = os.path.join("analysis", "mask_ST_data", "{sample}"),
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
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
    output:
        index = os.path.join("analysis", "self_embedding_info", "{embedding}", "{sample}", "pretreat.pkl"),
        emb_info = os.path.join("analysis", "self_embedding_info", "{embedding}", "{sample}", "Embedding.bin.info.tsv"),
    params:
        min_cell_num = 2,
        gene_id_label = "geneID",
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = "MIDCounts",
        embedding = lambda wildcards: embedding_dict[wildcards.embedding]["embedding"],
        embedding_dim = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_dim"],
        out_dir = os.path.join("analysis", "self_embedding_info", "{embedding}", "{sample}"),
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
    threads: 1
    shell:
        """
RedeViz enhance pretreat_by_ST \
    --spot {input.spot} \
    --mask {input.mask} \
    --x-index-label {params.x_index_label} \
    --y-index-label {params.y_index_label} \
    --UMI-label {params.UMI_label} \
    --gene-id-label {params.gene_id_label} \
    --embedding {params.embedding} \
    --embedding-dim {params.embedding_dim} \
    --min-cell-num {params.min_cell_num} \
    -o {params.out_dir}
        """

rule cell_type2color:
    input:
        emb_info = rules.RedeViz_embedding_pretreat.output.emb_info,
    output:
        color = os.path.join("analysis", "self_embedding_info", "{embedding}", "{sample}", "color.tsv"),
    params:
        cell_type2color = "python scripts/cell_type2color.py",
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
    shell:
        """
{params.cell_type2color} \
    --emb-info {input.emb_info} \
    --output {output.color}
        """

rule RedeViz_embedding_build_index:
    input:
        index = rules.RedeViz_embedding_pretreat.output.index,
        cov = os.path.join("analysis", "cov", "{sample}.cov.txt"),
    output:
        index = os.path.join("analysis", "RedeViz_self_embedding_index", "{embedding}", "{sample}.index.pkl"),
    params:
        bin_size = [9, 21],
        queue = "gpu1,gpu2",
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
        cov = os.path.join("analysis", "cov", "{sample}.cov.txt"),
    output:
        res = os.path.join("analysis", "RedeViz_self_enhance", "{embedding}", "{sample}.predict.tsv"),
    params:
        cell_radius = 11,
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = "MIDCounts",
        gene_name_label = "geneID",
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
    --mid-signal-cutoff $(cat {input.cov} | awk '{{print $1 * 0.75}}') \
    --window-size 100"""


rule RedeViz_phenotype_visualization:
    input:
        enhance = rules.RedeViz_enhance.output.res,
    output:
        all_img = os.path.join("analysis", "RedeViz_self_enhance_pheno_plot", "{embedding}", "{sample}.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_self_enhance_pheno_plot", "{embedding}", "{sample}.enhance.phenotype.HQ.png"),
    params:
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
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
        all_img = os.path.join("analysis", "RedeViz_self_enhance_celltype_plot", "{embedding}", "{sample}.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_self_enhance_celltype_plot", "{embedding}", "{sample}.enhance.phenotype.HQ.png"),
    params:
        queue = "cpu1,cpu2,fat,amd1",
        gpu_gres = ""
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

rule self_enhance_all:
    input:
        [
            rules.RedeViz_phenotype_visualization.output.all_img.format(
                sample=sample, embedding=embedding
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for sample in ST_sample_li
        ],
        [
            rules.RedeViz_celltype_visualization.output.all_img.format(
                sample=sample, embedding=embedding
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for sample in ST_sample_li
        ],
