import os


ST_sample_dict = {
    "E16.5": ["E16.5_E1S1"],
}

all_tissue_li = ['E16.5']

organ_li = ["kidney", "lung", "ear"]

embedding_dict = {
    "tSNE": {"embedding": "tSNE", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP2": {"embedding": "UMAP", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP3": {"embedding": "UMAP", "embedding_dim": 3, "embedding_resolution": 4}
    }

def get_ST_UMI_label(sample):
    if sample in ["E10.5_E2S1", "E12.5_E2S1"]:
        return "UMICounts"
    else:
        return "MIDCounts"

rule infer_spot_expr:
    input:
        spot = os.path.join("data", "organ", "{sample}.{organ}.tsv"),
    output:
        cov = os.path.join("analysis", "organ_cov", "{tissue}", "{sample}.{organ}.cov.txt"),
    params:
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = lambda wildcards: get_ST_UMI_label(wildcards.sample),
        queue = "cpu1,cpu2,amd1,fat",
        gpu_gres = ""
    threads: 10
    shell:
        """
RedeViz pretreatment infer_spot_expr\
    --spot {input.spot} \
    --x-index-label {params.x_index_label} \
    --y-index-label {params.y_index_label} \
    --UMI-label {params.UMI_label} \
    --output {output.cov}
        """

rule mask_ST_data:
    input:
        spot = os.path.join("data", "organ", "{sample}.{organ}.tsv"),
    output:
        mask = os.path.join("analysis", "organ_mask_ST_data", "{tissue}", "{sample}", "{organ}", "mask.npz"),
    params:
        mask_model = "bin",
        cell_diameter = 20,
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = lambda wildcards: get_ST_UMI_label(wildcards.sample),
        smooth_sd = 6,
        shift_cutoff = 255,
        out_dir = os.path.join("analysis", "organ_mask_ST_data", "{tissue}", "{sample}", "{organ}"),
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
    threads: 10
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
    --shift-cutoff {params.shift_cutoff} \
    --output {params.out_dir}
        """

rule spot2bin:
    input:
        spot = rules.mask_ST_data.input.spot,
    output:
        res = os.path.join("data", "organ_self_spot2bin", "{sample}.{organ}.bin{binN}.tsv"),
    params:
        x_index_label = "x",
        y_index_label = "y",
        gene_name_label = "geneID",
        UMI_label = "MIDCounts",
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
    threads: 10
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
        gene_li = os.path.join("data", "organ", "{organ}.gene.txt")
    output:
        flag = touch(os.path.join("analysis", "organ_self_expr_gene", "{tissue}", "{sample}", "{organ}", "plot_bin{binN}", "done.flag"))
    params:
        plot_gene_bin_expr = "python scripts/plot_gene_bin_expr.py",
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        gene_name_label = "geneID",
        UMI_label = "MIDCounts",
        out_dir = os.path.join("analysis", "organ_self_expr_gene", "{tissue}", "{sample}", "{organ}", "plot_bin{binN}"),
        queue = "cpu1,cpu2,fat",
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

def get_black_list_params(tissue):
    if tissue == "E16.5":
        params = "--gene-blacklist data/blacklist.E16.txt"
    else:
        params = ""
    return params


rule RedeViz_embedding_pretreat:
    input:
        spot = rules.mask_ST_data.input.spot,
        mask = rules.mask_ST_data.output.mask,
    output:
        index = os.path.join("analysis", "organ_self_embedding_info", "{embedding}", "{tissue}", "{sample}", "{organ}", "pretreat.pkl"),
        emb_info = os.path.join("analysis", "organ_self_embedding_info", "{embedding}", "{tissue}", "{sample}", "{organ}", "Embedding.bin.info.tsv"),
    params:
        min_cell_num = 2,
        gene_id_label = "geneID",
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = lambda wildcards: get_ST_UMI_label(wildcards.sample),
        embedding = lambda wildcards: embedding_dict[wildcards.embedding]["embedding"],
        embedding_dim = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_dim"],
        out_dir = os.path.join("analysis", "organ_self_embedding_info", "{embedding}", "{tissue}", "{sample}", "{organ}"),
        black_list_params = lambda wildcards: get_black_list_params(wildcards.tissue),
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
    threads: 20
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
    -o {params.out_dir} {params.black_list_params}
        """

rule RedeViz_embedding_build_index:
    input:
        index = os.path.join("analysis", "organ_self_embedding_info", "{embedding}", "{tissue}", "{sample}", "{ref_organ}", "pretreat.pkl"),
        cov = os.path.join("analysis", "organ_cov", "{tissue}", "{sample}.{organ}.cov.txt"),
    output:
        index = os.path.join("analysis", "RedeViz_organ_self_embedding_index", "{embedding}", "{tissue}", "{sample}.{organ}.Ref{ref_organ}.index.pkl"),
    params:
        bin_size = [9, 21],
        queue = "gpu1,gpu2",
        gpu_gres = ""
    threads: 20
    shell:
        """
RedeViz enhance build -i {input.index} \
    --bin-size {params.bin_size} \
    --UMI-per-spot $(cat {input.cov} | awk '{{print $1 * 0.75 }}' ) \
    -o {output.index}
        """


rule RedeViz_enhance:
    input:
        index = os.path.join("analysis", "RedeViz_organ_self_embedding_index", "{embedding}", "{tissue}", "{sample}.{organ}.Ref{ref_organ}.index.pkl"),
        spot = rules.mask_ST_data.input.spot,
        cov = os.path.join("analysis", "organ_cov", "{tissue}", "{sample}.{organ}.cov.txt"),
    output:
        res = os.path.join("analysis", "RedeViz_organ_self_enhance", "{embedding}", "{tissue}", "{sample}.{organ}.Ref{ref_organ}.predict.tsv"),
    params:
        cell_radius = 11,
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = lambda wildcards: get_ST_UMI_label(wildcards.sample),
        gene_name_label = "geneID",
        queue = "gpu1",
        gpu_gres = "--gres=gpu:1"
    threads: 20
    shell:
        """
RedeViz enhance run \
    -i {input.index} -s {input.spot} -o {output.res} \
    --x-index-label {params.x_index_label} --y-index-label {params.y_index_label} \
    --gene-name-label {params.gene_name_label} --UMI-label {params.UMI_label} \
    --cell-radius {params.cell_radius} \
    --mid-signal-cutoff $(cat {input.cov} | awk '{{print $1 * 0.75}}' ) --batch-effect-fct 1 \
    --window-size 65  --device-name cuda"""


rule RedeViz_phenotype_visualization:
    input:
        enhance = rules.RedeViz_enhance.output.res,
    output:
        all_img = os.path.join("analysis", "RedeViz_organ_self_enhance_pheno_plot", "{embedding}", "{tissue}", "{sample}.{organ}.Ref{ref_organ}.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_organ_self_enhance_pheno_plot", "{embedding}", "{tissue}", "{sample}.{organ}.Ref{ref_organ}.enhance.phenotype.HQ.png"),
    params:
        queue = "cpu1,cpu2,fat,amd1,hygon",
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

rule all:
    input:
        [
            rules.RedeViz_phenotype_visualization.output.all_img.format(
                tissue=tissue, sample=sample, embedding=embedding, organ=organ, ref_organ=ref_organ
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for tissue in ["E16.5"]
                for sample in ST_sample_dict[tissue] 
                for organ in organ_li
                for ref_organ in organ_li
        ]
