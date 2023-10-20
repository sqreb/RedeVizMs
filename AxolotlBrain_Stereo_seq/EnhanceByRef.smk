import os


ST_sample_li = ["30DPI", "60DPI"]

REF_time_dict = {
   "30DPI": ["30DPI"],
   "60DPI": ["30DPI", "60DPI"]
}

embedding_dict = {
    "tSNE": {"embedding": "tSNE", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP2": {"embedding": "UMAP", "embedding_dim": 2, "embedding_resolution": 2.5},
    "UMAP3": {"embedding": "UMAP", "embedding_dim": 3, "embedding_resolution": 4}
    }

domain_params = {
    "large1": {"smooth_radius": 50, "receptive_radius": 100, "emb_radius": 10, "merge_emb_dist": 20, "min_spot_num": 4000},
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


rule RedeViz_embedding_pretreat:
    input:
        h5ad = os.path.join("data", "ref_scRNA", "{ref_sample}.cnt.h5ad"),
    output:
        index = os.path.join("analysis", "embedding_info", "{embedding}", "{ref_sample}", "pretreat.pkl"),
        emb_info = os.path.join("analysis", "embedding_info", "{embedding}", "{ref_sample}", "Embedding.bin.info.tsv"),
        emb_cell_info = os.path.join("analysis", "embedding_info", "{embedding}", "{ref_sample}", "CellInfo.tsv"),
    params:
        cell_type_label = "Annotation",
        gene_id_label = "Gene",
        embedding = lambda wildcards: embedding_dict[wildcards.embedding]["embedding"],
        embedding_dim = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_dim"],
        embedding_resolution = lambda wildcards: embedding_dict[wildcards.embedding]["embedding_resolution"],
        out_dir = os.path.join("analysis", "embedding_info", "{embedding}", "{ref_sample}"),
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
    -o {params.out_dir}
        """


rule RedeViz_embedding_build_index:
    input:
        index = rules.RedeViz_embedding_pretreat.output.index,
        cov = os.path.join("analysis", "cov", "{sample}.cov.txt"),
    output:
        index = os.path.join("analysis", "RedeViz_embedding_index", "{embedding}", "{sample}_ref_{ref_sample}.index.pkl"),
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
        spot = os.path.join("data", "ST", "{sample}.gem"),
        cov = os.path.join("analysis", "cov", "{sample}.cov.txt"),
    output:
        res = os.path.join("analysis", "RedeViz_enhance", "{embedding}", "{sample}_ref_{ref_sample}.predict.tsv"),
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
        all_img = os.path.join("analysis", "RedeViz_enhance_pheno_plot", "{embedding}", "{sample}_ref_{ref_sample}.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_enhance_pheno_plot", "{embedding}", "{sample}_ref_{ref_sample}.enhance.phenotype.HQ.png"),
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
        color = os.path.join("data", "celltype.color.tsv"),
    output:
        all_img = os.path.join("analysis", "RedeViz_enhance_celltype_plot", "{embedding}", "{sample}_ref_{ref_sample}.enhance.phenotype.all.png"),
        HQ_img = os.path.join("analysis", "RedeViz_enhance_celltype_plot", "{embedding}", "{sample}_ref_{ref_sample}.enhance.phenotype.HQ.png"),
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

rule RedeViz_spot2domain:
    input:
        enhance = rules.RedeViz_enhance.output.res,
    output:
        res = os.path.join("analysis", "RedeViz_domain", "{embedding}", "{sample}.ref_sample_{ref_sample}", "domain_param_{domain_param}", "domain.label.tsv"),
    params:
        smooth_radius = lambda wildcards: domain_params[wildcards.domain_param]["smooth_radius"],
        receptive_radius = lambda wildcards: domain_params[wildcards.domain_param]["receptive_radius"],
        emb_radius = lambda wildcards: domain_params[wildcards.domain_param]["emb_radius"],
        merge_emb_dist = lambda wildcards: domain_params[wildcards.domain_param]["merge_emb_dist"],
        min_spot_num = lambda wildcards: domain_params[wildcards.domain_param]["min_spot_num"],
        out_dir = os.path.join("analysis", "RedeViz_domain", "{embedding}", "{sample}.ref_sample_{ref_sample}", "domain_param_{domain_param}"),
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

rule RedeViz_compare:
    input:
        pred1 = lambda wildcards: os.path.join("analysis", "RedeViz_enhance", "{embedding}", "60DPI_ref_30DPI.predict.tsv"),
        pred2 = lambda wildcards: os.path.join("analysis", "RedeViz_enhance", "{embedding}", "60DPI_ref_60DPI.predict.tsv"),
    output:
        res = os.path.join("analysis", "RedeViz_compare", "{embedding}", "60DPI.compare.tsv"),
    params:
        queue = "cpu1,cpu2,fat",
        gpu_gres = ""
    threads: 10
    shell:
        """
RedeViz posttreatment compare \
    --input1 {input.pred1} \
    --input2 {input.pred2} \
    --SOR-cutoff 1.0 \
    --SBR-cutoff 1.0 \
    --output {output.res} \
    --keep-other \
    --min-denoise-spot-num 200 \
    --denoise
        """

rule build_RedeViz_imputation_index:
    input:
        index = rules.RedeViz_embedding_pretreat.output.index,
        sce = os.path.join("data", "ref_scRNA", "{ref_sample}.cnt.h5ad"),
    output:
        index = os.path.join("analysis", "RedeViz_imputation_index", "{embedding}", "{ref_sample}.pkl"),
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
        enhance = rules.RedeViz_enhance.output.res,
        spot = os.path.join("data", "ST", "{sample}.gem"),
        index = rules.build_RedeViz_imputation_index.output,
        gene_li = os.path.join("data", "impute.gene.txt"),
    output:
        flag = touch(os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{sample}_ref_{ref_sample}", "done.flag")),
    params:
        x_label = "bin_x_index",
        y_label = "bin_y_index",
        UMI_label = "UMI",
        out_dir = os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{sample}_ref_{ref_sample}"),
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
        flag = touch(os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{sample}_ref_{ref_sample}", "plot.flag")),
    params:
        out_dir = os.path.join("analysis", "RedeViz_imputation", "{embedding}", "{sample}_ref_{ref_sample}"),
    shell:
        """
for f in $(ls {params.out_dir}/*.imputate.npz); do
    RedeViz posttreatment plot_gene_expr -R ${{f}} -G ${{f}} -o $(echo ${{f}} | sed 's/npz/png/g')
done
        """

rule spot2bin:
    input:
        spot = os.path.join("data", "ST", "{sample}.gem"),
    output:
        res = os.path.join("analysis", "bin", "{sample}.bin{binN}.txt"),
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
        gene_li = os.path.join("data", "impute.gene.txt"),
    output:
        flag = touch(os.path.join("analysis", "gene_bin_expr", "{sample}.bin{binN}", "done.flag"))
    params:
        plot_gene_bin_expr = "python scripts/plot_gene_bin_expr.py",
        x_index_label = "bin_x_index",
        y_index_label = "bin_y_index",
        gene_name_label = "geneID",
        UMI_label = "MIDCounts",
        out_dir = os.path.join("analysis", "gene_bin_expr", "{sample}.bin{binN}"),
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

rule enhance_all:
    input:
        [
            rules.RedeViz_phenotype_visualization.output.all_img.format(
                sample=sample, embedding=embedding, ref_sample=ref_sample
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for sample in ST_sample_li
                for ref_sample in REF_time_dict[sample]
        ],[
            rules.RedeViz_celltype_visualization.output.all_img.format(
                sample=sample, embedding=embedding, ref_sample=ref_sample
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for sample in ST_sample_li
                for ref_sample in REF_time_dict[sample]
        ],[
            rules.RedeViz_spot2domain.output.res.format(
                sample=sample, ref_sample=ref_sample, embedding=embedding, domain_param=domain_param
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for sample in ST_sample_li
                for ref_sample in REF_time_dict[sample]
                for domain_param in domain_params.keys()
        ],[
            rules.RedeViz_compare.output.res.format(
                embedding=embedding
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
        ],[
            rules.RedeViz_plot_imputation.output.flag.format(
                sample=sample, ref_sample=ref_sample, embedding=embedding
                ) for embedding in ["tSNE", "UMAP2", "UMAP3"] 
                for sample in ST_sample_li
                for ref_sample in REF_time_dict[sample]
        ],[
            rules.plot_gene_bin_expr.output.flag.format(
                sample=sample, binN=binN
                ) for sample in ST_sample_li
                for binN in ["1", "50"]
        ]



