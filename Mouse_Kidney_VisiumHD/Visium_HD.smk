import os
from collections import defaultdict

SAMPLE_LI = ["Kidney"]
embedding_dict = {
    "tSNE": {"embedding": "tSNE", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP2": {"embedding": "UMAP", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP3": {"embedding": "UMAP", "embedding_dim": 3, "embedding_resolution": 6}
    }


RGB_gene_dict = {
    "SmoothMuscleL5": {"R": "Acta2", "G": "Upk1b", "B": "Scd1"},
    "L1L2": {"R": "Slc34a1", "G": "S100g", "B": "Atp11a"},
    "L3L4": {"R": "Slc12a1", "G": "Akr1b3", "B": "Aqp4"},
    "L1": {"R": "Nphs2", "G": "Car15", "B": "Hbb-bs"},
    "L3": {"R": "Bst1", "G": "Slc14a2", "B": "Aqp6"},
}


BASE_DIR = os.path.join("analysis", "Visium_HD")


rule infer_spot_expr:
    input:
        spot = os.path.join("data", "ST.data.tsv"),
    output:
        cov = os.path.join(BASE_DIR, "cov", "{sample}.cov.txt"),
    params:
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = "UMI",
    threads: 1
    shell:
        """
RedeViz pretreatment infer_spot_expr\
    --spot {input.spot} \
    --x-index-label {params.x_index_label} \
    --y-index-label {params.y_index_label} \
    --UMI-label {params.UMI_label} \
    --output {output.cov}
        """

rule RedeViz_embedding_pretreat:
    input:
        h5ad = os.path.join("data", "SCE.filter.WithAnn.h5ad")
    output:
        index = os.path.join(BASE_DIR, "embedding_info", "{embedding}", "pretreat.pkl"),
        emb_info = os.path.join(BASE_DIR, "embedding_info", "{embedding}", "Embedding.bin.info.tsv"),
        emb_cell_info = os.path.join(BASE_DIR, "embedding_info", "{embedding}", "CellInfo.tsv"),
    params:
        cell_type_label = "Annotation",
        embedding = lambda wildcards: embedding_dict[wildcards.embedding]["embedding"],
        embedding_dim = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_dim"],
        embedding_resolution = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_resolution"],
        out_dir = os.path.join(BASE_DIR, "embedding_info", "{embedding}"),
    threads: 16
    shell:
        """
RedeViz enhance pretreat -i {input.h5ad} \
    --cell-type-label {params.cell_type_label} \
    --embedding {params.embedding} \
    --embedding-dim {params.embedding_dim} \
    --embedding-resolution {params.embedding_resolution} \
    -o {params.out_dir}
        """

rule RedeViz_embedding_build_index:
    input:
        index = rules.RedeViz_embedding_pretreat.output.index,
        cov = rules.infer_spot_expr.output.cov,
    output:
        index = os.path.join(BASE_DIR, "RedeViz_embedding_index", "{embedding}", "{sample}.index.pkl"),
    params:
        bin_size = [3, 9],
    threads: 16
    shell:
        """
RedeViz enhance build -i {input.index} \
    --bin-size {params.bin_size} \
    --UMI-per-spot $(cat {input.cov} | awk '{{print $1 * 0.75}}') \
    -o {output.index}  --device-name cuda
        """

rule RedeViz_enhance:
    input:
        index = rules.RedeViz_embedding_build_index.output.index,
        spot = rules.infer_spot_expr.input.spot,
        cov = rules.infer_spot_expr.output.cov,
    output:
        res = protected(os.path.join(BASE_DIR, "RedeViz_enhance", "{embedding}", "{sample}.predict.tsv")),
    params:
        cell_radius = 3,
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = "UMI",
        gene_name_label = "GeneName",
    threads: 16
    shell:
        """
RedeViz enhance run \
    -i {input.index} -s {input.spot} -o {output.res} \
    --x-index-label {params.x_index_label} --y-index-label {params.y_index_label} \
    --gene-name-label {params.gene_name_label} --UMI-label {params.UMI_label} \
    --cell-radius {params.cell_radius} \
    --mid-signal-cutoff $(cat {input.cov} | awk '{{print $1 * 0.75}}') \
    --window-size 60  --device-name cuda"""


rule RedeViz_phenotype_visualization:
    input:
        enhance = rules.RedeViz_enhance.output.res,
    output:
        all_img = os.path.join(BASE_DIR, "RedeViz_enhance_pheno_plot", "{embedding}", "{sample}.enhance.phenotype.all.png"),
        HQ_img = os.path.join(BASE_DIR, "RedeViz_enhance_pheno_plot", "{embedding}", "{sample}.enhance.phenotype.HQ.png"),
    params:
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
        color = os.path.join("data", "color.tsv")
    output:
        all_img = os.path.join(BASE_DIR, "RedeViz_enhance_celltype_plot", "{embedding}", "{sample}.enhance.celltype.all.png"),
        HQ_img = os.path.join(BASE_DIR, "RedeViz_enhance_celltype_plot", "{embedding}", "{sample}.enhance.celltype.HQ.png"),
    params:
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


rule spot2bin:
    input:
        spot = rules.infer_spot_expr.input.spot,
    output:
        res = os.path.join(BASE_DIR, "bin", "{sample}.bin{binN}.txt"),
    params:
        x_index_label = "x",
        y_index_label = "y",
        gene_name_label = "GeneName",
        UMI_label = "UMI",
        gpu_gres = ""
    threads: 1
    shell:
        """
RedeViz pretreatment spot2bin \
    --spot {input.spot} \
    --bin-size {wildcards.binN} \
    --x-index-label {params.x_index_label} \
    --y-index-label {params.y_index_label} \
    --gene-name-label {params.gene_name_label} \
    --UMI-label {params.UMI_label} \
    --output {output.res}
        """


rule plot_gene_bin_expr:
    input:
        spot = rules.spot2bin.output.res,
        gene_li = os.path.join("data", "marker.gene.txt"),
    output:
        flag = touch(os.path.join(BASE_DIR, "gene_bin_expr", "{sample}.bin{binN}", "done.flag"))
    params:
        plot_gene_bin_expr = "python scripts/plot_gene_bin_expr.py",
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        gene_name_label = "GeneName",
        UMI_label = "UMI",
        out_dir = os.path.join(BASE_DIR, "gene_bin_expr", "{sample}.bin{binN}"),
        gpu_gres = ""
    threads: 10
    shell:
        """
{params.plot_gene_bin_expr} \
    --input {input.spot} \
    --x-label {params.x_index_label} \
    --y-label {params.y_index_label} \
    --gene-name-label {params.gene_name_label} \
    --UMI-label {params.UMI_label} \
    --output {params.out_dir} \
    --gene-list {input.gene_li}
        """


def get_gene_bin_params(gene_li):
    tmp_dict = RGB_gene_dict[gene_li]
    R = tmp_dict["R"]
    G = tmp_dict["G"]
    B = tmp_dict["B"]
    params = ""
    if R:
        params += f" -R {R}"
    if G:
        params += f" -G {G}"
    if B:
        params += f" -B {B}"
    return params

rule plot_bin_expr_RGB:
    input:
        spot = rules.spot2bin.output.res,
    output:
        plot = os.path.join(BASE_DIR, "gene_bin_expr_RGB", "{sample}.bin{binN}", "{gene_li}.png")
    params:
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        gene_name_label = "GeneName",
        UMI_label = "UMI",
        RGB_params = lambda wildcards: get_gene_bin_params(wildcards.gene_li),
    shell:
        """
RedeViz posttreatment plot_gene_bin_expr \
    -i {input.spot} \
    --x-label {params.x_index_label} \
    --y-label {params.y_index_label} \
    --gene-name-label {params.gene_name_label} \
    --UMI-label {params.UMI_label} \
    {params.RGB_params} -o {output.plot}
        """

rule all:
    input:
        expand(rules.RedeViz_phenotype_visualization.output, embedding=["UMAP2", "tSNE", "UMAP3"], sample=SAMPLE_LI),
        expand(rules.RedeViz_celltype_visualization.output, embedding=["UMAP2", "tSNE", "UMAP3"], sample=SAMPLE_LI),
        expand(rules.plot_gene_bin_expr.output, sample=SAMPLE_LI, binN=[1, 5, 10]),
        expand(rules.plot_bin_expr_RGB.output, sample=SAMPLE_LI, binN=[1, 5, 10], gene_li=RGB_gene_dict),
