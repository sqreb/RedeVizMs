import os


ST_sample_dict = {
    "E16.5": ["E16.5_E1S1"]
}

REF_sample_dict = {
    "E16.5": ["E16.5"],
}
all_tissue_li = ['E16.5']


REF_tissue_dict = {
    "E16.5": ["E16.5"],
}

embedding_dict = {
    "tSNE": {"embedding": "tSNE", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP2": {"embedding": "UMAP", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP3": {"embedding": "UMAP", "embedding_dim": 3, "embedding_resolution": 4}
    }

domain_params = {
    "large1": {"smooth_radius": 0, "receptive_radius": 100, "emb_radius": 10, "merge_emb_dist": 10, "min_spot_num": 4000},
    "large2": {"smooth_radius": 0, "receptive_radius": 100, "emb_radius": 15, "merge_emb_dist": 15, "min_spot_num": 4000},
    "large3": {"smooth_radius": 0, "receptive_radius": 100, "emb_radius": 20, "merge_emb_dist": 20, "min_spot_num": 4000}
}

RGB_gene_dict = {
    "eye": {"R": "Col1a1", "G": "Cnmd", "B": "Elavl4"},
    "lung": {"R": "Epcam", "G": "Col1a2", "B": "Hbb-bs"}
}

rule infer_spot_expr:
    input:
        spot = os.path.join("data", "Bin1_matrix", "{sample}_GEM_bin1.tsv"),
    output:
        cov = os.path.join("analysis", "cov", "{tissue}", "{sample}.cov.txt"),
    params:
        x_index_label = "x",
        y_index_label = "y",
        UMI_label = "MIDCounts",
        queue = "cpu1,cpu2,amd1,fat",
        gpu_gres = ""
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

def get_black_list_params(ref_sample):
    if ref_sample == "E16.5":
        params = "--gene-blacklist data/blacklist.E16.txt"
    else:
        params = ""
    return params

rule RedeViz_embedding_pretreat:
    input:
        h5ad = os.path.join("data", "stomics", "{ref_sample}.cnt.filter.h5ad")
    output:
        index = os.path.join("analysis", "embedding_info", "{embedding}", "{ref_sample}", "pretreat.pkl"),
        tSNE_info = os.path.join("analysis", "embedding_info", "{embedding}", "{ref_sample}", "Embedding.bin.info.tsv"),
        tSNE_cell_info = os.path.join("analysis", "embedding_info", "{embedding}", "{ref_sample}", "CellInfo.tsv"),
    params:
        cell_type_label = "annotation",
        gene_id_label = "Gene",
        embedding = lambda wildcards: embedding_dict[wildcards.embedding]["embedding"],
        embedding_dim = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_dim"],
        embedding_resolution = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_resolution"],
        out_dir = os.path.join("analysis", "embedding_info", "{embedding}", "{ref_sample}"),
        black_list_params = lambda wildcards: get_black_list_params(wildcards.ref_sample),
        queue = "cpu1,cpu2,amd1,fat",
        gpu_gres = ""
    threads: 10
    shell:
        """
RedeViz enhance pretreat -i {input.h5ad} \
    --cell-type-label {params.cell_type_label} \
    --gene-id-label {params.gene_id_label} \
    --embedding {params.embedding} \
    --embedding-dim {params.embedding_dim} \
    --embedding-resolution {params.embedding_resolution} \
    -o {params.out_dir} {params.black_list_params}
        """

rule RedeViz_embedding_build_index:
    input:
        index = rules.RedeViz_embedding_pretreat.output.index,
        cov = rules.infer_spot_expr.output.cov,
    output:
        index = os.path.join("analysis", "RedeViz_embedding_index", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}.index.pkl"),
    params:
        bin_size = [9, 21],
        queue = "gpu1",
        gpu_gres = "--gres=gpu:1"
    threads: 20
    shell:
        """
RedeViz enhance build -i {input.index} \
    --bin-size {params.bin_size} \
    --UMI-per-spot $(cat {input.cov} | awk '{{print $1 * 0.75}}' ) \
    --norm-method None \
    -o {output.index} --device-name cuda
        """

rule RedeViz_enhance:
    input:
        index = rules.RedeViz_embedding_build_index.output.index,
        spot = rules.infer_spot_expr.input.spot,
        cov = rules.infer_spot_expr.output.cov,
    output:
        res = os.path.join("analysis", "RedeViz_enhance", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}.predict.tsv"),
    params:
        cell_radius = 11,
        x_index_label = "x",
        y_index_label = "y",
        gene_name_label = "geneID",
        UMI_label = "MIDCounts",
        queue = "gpu1",
        gpu_gres = "--gres=gpu:1"
    threads: 20
    shell:
        """
RedeViz enhance run -i {input.index} -s {input.spot} -o {output.res} \
    --x-index-label {params.x_index_label} --y-index-label {params.y_index_label} \
    --gene-name-label {params.gene_name_label} --UMI-label {params.UMI_label} \
    --cell-radius {params.cell_radius} \
    --mid-signal-cutoff $(cat {input.cov}) \
    --window-size 70  --device-name cuda"""


rule RedeViz_phenotype_visualization:
    input:
        enhance = rules.RedeViz_enhance.output.res,
    output:
        all_img = os.path.join("analysis", "RedeViz_enhance_pheno_plot", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_enhance_pheno_plot", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}.enhance.phenotype.HQ.png"),
    params:
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
    threads: 10
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

rule RedeViz_phenotype_visualization_BRG:
    input:
        enhance = os.path.join("analysis", "RedeViz_enhance", "UMAP3", "{tissue}", "{sample}.ref_sample_{ref_sample}.predict.tsv"),
    output:
        all_img = os.path.join("analysis", "RedeViz_enhance_pheno_plot_BRG", "UMAP3", "{tissue}", "{sample}.ref_sample_{ref_sample}.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_enhance_pheno_plot_BRG", "UMAP3", "{tissue}", "{sample}.ref_sample_{ref_sample}.enhance.phenotype.HQ.png"),
    params:
        queue = "cpu1,cpu2,fat",
        gpu_gres = "",
        plot_pheno_BRG = "python scripts/plot_pheno_BRG.py"
    threads: 10
    shell:
        """
{params.plot_pheno_BRG} \
    --input {input.enhance} \
    --output {output.all_img} \
    --keep-other --denoise

{params.plot_pheno_BRG} \
    --input {input.enhance} \
    --output {output.HQ_img} --denoise
        """

rule RedeViz_celltype_visualization:
    input:
        enhance = rules.RedeViz_enhance.output.res,
        color = os.path.join("data", "color.tsv"),
    output:
        all_img = os.path.join("analysis", "RedeViz_enhance_celltype_plot", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}.enhance.celltype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_enhance_celltype_plot", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}.enhance.celltype.HQ.png"),
    params:
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
    threads: 10
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

rule RedeViz_spot2domain:
    input:
        enhance = rules.RedeViz_enhance.output.res,
    output:
        res = os.path.join("analysis", "RedeViz_domain", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}", "domain_param_{domain_param}", "domain.label.tsv"),
    params:
        smooth_radius = lambda wildcards: domain_params[wildcards.domain_param]["smooth_radius"],
        receptive_radius = lambda wildcards: domain_params[wildcards.domain_param]["receptive_radius"],
        emb_radius = lambda wildcards: domain_params[wildcards.domain_param]["emb_radius"],
        merge_emb_dist = lambda wildcards: domain_params[wildcards.domain_param]["merge_emb_dist"],
        min_spot_num = lambda wildcards: domain_params[wildcards.domain_param]["min_spot_num"],
        out_dir = os.path.join("analysis", "RedeViz_domain", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}", "domain_param_{domain_param}"),
        queue = "gpu1,gpu2",
        gpu_gres = "--gres=gpu:1"
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
    --output {params.out_dir} --keep-other --denoise
        """

rule build_RedeViz_imputation_index:
    input:
        index = rules.RedeViz_embedding_pretreat.output.index,
        sce = os.path.join("data", "stomics", "{ref_sample}.cnt.filter.h5ad")
    output:
        index = os.path.join("analysis", "RedeViz_imputation_index", "{embedding}", "{ref_sample}.pkl"),
    params:
        gene_id_label = "Gene",
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
    threads: 10
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
        enhance = rules.RedeViz_enhance.output.res,
        spot = rules.infer_spot_expr.input.spot,
        index = rules.build_RedeViz_imputation_index.output,
        gene_li = os.path.join("data", "impute.gene.txt"),
    output:
        flag = touch(os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}", "done.flag")),
    params:
        x_label = "x",
        y_label = "y",
        UMI_label = "MIDCounts",
        embedding_smooth_sigma = 3.0,
        UMI_smooth_sigma = 10.0,
        out_dir = os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}"),
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
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

rule RedeViz_plot_imputation:
    input:
        flag = rules.RedeViz_imputation.output,
    output:
        flag = touch(os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}", "plot.flag")),
    params:
        out_dir = os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}"),
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
    shell:
        """
for f in $(ls {params.out_dir}/*.imputate.npz); do
    RedeViz posttreatment plot_gene_expr -R ${{f}} -G ${{f}} -o $(echo ${{f}} | sed 's/npz/png/g')
done
        """


rule spot2bin:
    input:
        spot = rules.infer_spot_expr.input.spot,
    output:
        res = os.path.join("analysis", "bin", "{tissue}", "{sample}.bin{binN}.txt"),
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

rule plot_gene_bin_expr_interest:
    input:
        spot = rules.spot2bin.output.res,
        gene_li = os.path.join("data", "impute.gene.txt"),
    output:
        flag = touch(os.path.join("analysis", "gene_bin_expr_interest", "{tissue}", "{sample}.bin{binN}", "done.flag"))
    params:
        plot_gene_bin_expr = "python scripts/plot_gene_bin_expr.py",
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        gene_name_label = "geneID",
        UMI_label = "MIDCounts",
        out_dir = os.path.join("analysis", "gene_bin_expr_interest", "{tissue}", "{sample}.bin{binN}"),
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

def get_impute_gene_params(data_dir, gene_li):
    tmp_dict = RGB_gene_dict[gene_li]
    R = tmp_dict["R"]
    G = tmp_dict["G"]
    B = tmp_dict["B"]
    params = ""
    if R:
        params += f" -R {data_dir}/{R}.imputate.npz"
    if G:
        params += f" -G {data_dir}/{G}.imputate.npz"
    if B:
        params += f" -B {data_dir}/{B}.imputate.npz"
    return params


rule get_signal_pos:
    input:
        spot = rules.RedeViz_enhance.output.res,
    output:
        sig = os.path.join("analysis", "signal_pos", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}.npz"),
    params:
        x_index_label = "x",
        y_index_label = "y",
        get_signal_pos = "python scripts/get_signal_pos.py",
        queue = "cpu1,cpu2,amd1,fat",
        gpu_gres = ""
    threads: 1
    shell:
        """
{params.get_signal_pos} \
    --input {input.spot} \
    --x-label {params.x_index_label} \
    --y-label {params.y_index_label} \
    --output {output.sig}
        """

rule RedeViz_plot_imputation_final:
    input:
        flag = rules.RedeViz_imputation.output,
        signal = rules.get_signal_pos.output.sig,
    output:
        plot = os.path.join("analysis", "RedeViz_imputation_plot_final", "{embedding}", "{tissue}", "{sample}.ref_sample_{ref_sample}", "{gene_li}.png"),
    params:
        rgb_params = lambda wildcards: get_impute_gene_params(os.path.join("analysis", "RedeViz_imputation_final", wildcards.embedding, wildcards.tissue, f"{wildcards.sample}.ref_sample_{wildcards.ref_sample}"), wildcards.gene_li),
        plot_imputation = "python scripts/plot_imputation.py",
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
    shell:
        """
{params.plot_imputation} \
    {params.rgb_params} \
    --signal {input.signal} \
    -o {output.plot}
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

rule RedeViz_plot_bin_expr_final:
    input:
        spot = rules.spot2bin.output.res,
    output:
        plot = os.path.join("analysis", "gene_bin_expr_final", "{tissue}", "{sample}.bin{binN}", "{gene_li}.png")
    params:
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        gene_name_label = "geneID",
        UMI_label = "MIDCounts",
        out_dir = os.path.join("analysis", "gene_bin_expr_final", "{tissue}", "{sample}.bin{binN}", "{gene_li}.png"),
        RGB_params = lambda wildcards: get_gene_bin_params(wildcards.gene_li),
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
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
        [
            rules.RedeViz_phenotype_visualization.output.all_img.format(
                tissue=tissue, sample=sample, ref_sample=ref_sample, embedding=embedding
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for tissue in ST_sample_dict.keys() 
                for sample in ST_sample_dict[tissue] 
                for ref_tissue in REF_tissue_dict[tissue] 
                for ref_sample in REF_sample_dict[ref_tissue]
        ],[
            rules.RedeViz_celltype_visualization.output.all_img.format(
                tissue=tissue, sample=sample, ref_sample=ref_sample, embedding=embedding
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for tissue in ST_sample_dict.keys() 
                for sample in ST_sample_dict[tissue] 
                for ref_tissue in REF_tissue_dict[tissue] 
                for ref_sample in REF_sample_dict[ref_tissue]
        ],[
            rules.RedeViz_spot2domain.output.res.format(
                tissue=tissue, sample=sample, ref_sample=ref_sample, embedding=embedding, domain_param=domain_param
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for tissue in ST_sample_dict.keys() 
                for sample in ST_sample_dict[tissue] 
                for ref_tissue in REF_tissue_dict[tissue] 
                for ref_sample in REF_sample_dict[ref_tissue]
                for domain_param in domain_params.keys()
        ],[
            rules.RedeViz_plot_imputation.output.flag.format(
                tissue=tissue, sample=sample, ref_sample=ref_sample, embedding=embedding
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for tissue in ST_sample_dict.keys() 
                for sample in ST_sample_dict[tissue] 
                for ref_tissue in REF_tissue_dict[tissue] 
                for ref_sample in REF_sample_dict[ref_tissue]
        ],[
            rules.RedeViz_plot_bin_expr_final.output.plot.format(
                tissue=tissue, sample=sample, ref_sample=ref_sample, gene_li=gene_li, binN=binN
                ) for tissue in ST_sample_dict.keys() 
                for sample in ST_sample_dict[tissue] 
                for ref_tissue in REF_tissue_dict[tissue] 
                for ref_sample in REF_sample_dict[ref_tissue]
                for gene_li in RGB_gene_dict.keys()
                for binN in [1, 50]
        ],[
            rules.RedeViz_plot_imputation_final.output.plot.format(
                tissue=tissue, sample=sample, ref_sample=ref_sample, embedding=embedding, gene_li=gene_li
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for tissue in ST_sample_dict.keys() 
                for sample in ST_sample_dict[tissue] 
                for ref_tissue in REF_tissue_dict[tissue] 
                for ref_sample in REF_sample_dict[ref_tissue]
                for gene_li in RGB_gene_dict.keys()
        ],[
            rules.plot_gene_bin_expr_interest.output.flag.format(
                tissue=tissue, sample=sample, binN=binN
                ) 
                for tissue in ST_sample_dict.keys() 
                for sample in ST_sample_dict[tissue] 
                for binN in [1, 50]
        ]
