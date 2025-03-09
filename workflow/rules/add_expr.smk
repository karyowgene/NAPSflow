import os
import glob

configfile: "../config/add_expr_config.yaml"

# Get list of VCF samples from folder if not specified in config
def get_vcf_samples():
    if config.get("samples_vcf"):
        return config["samples_vcf"]
    vcf_files = glob.glob(os.path.join(config["vcf_folder"], "*.vcf.gz"))
    return [os.path.basename(f).split("_")[0] for f in vcf_files]

VCF_SAMPLES = get_vcf_samples()

rule add_expr_all:
    input:
        config["expression_file"],
        expand("{vcf_folder}/{sample}_annotated.gx.vcf", 
               vcf_folder=config["vcf_folder"], sample=VCF_SAMPLES)

rule create_tx2gene:
    input:
        fa = config["reference_fa"]
    output:
        csv = config["tx2gene_csv"]
    shell:
        "zcat {input.fa} | grep '>' | awk 'BEGIN{{FS=\" \"}};"
        "{{print substr($1,2) \",\" substr($4,6)}};' | "
        "awk 'BEGIN{{print \"TXNAME,GENEID\"}}1' > {output.csv}"

rule tximport_expression:
    input:
        tx2gene = config["tx2gene_csv"],
        abundance_files = expand("{output_dir}/kallisto_quants/{sample}/abundance.tsv",
                               output_dir=config["output_dir"], sample=VCF_SAMPLES)
    output:
        expression = os.path.join(config["output_dir"], "gene_expression.tsv")
    script:
        "../scripts/tximport.R"

rule modify_expression:
    input:
        expression = os.path.join(config["output_dir"], "gene_expression.tsv")
    output:
        modified = config["expression_file"]
    run:
        samples = "\t".join(str(s) for s in VCF_SAMPLES)
        cmd = (
            f"awk 'BEGIN {{OFS=\"\\t\"; print \"gene_name\t{samples}\"}} "
            "NR>1 {{gsub(\"abundance.\", \"\", $1); print $0}}' {input} > {output}"
        )
        shell(cmd)

rule decompress_vcf:
    input:
        os.path.join(config["vcf_folder"], "{sample}_annotated.vcf.gz")
    output:
        temp(os.path.join(config["vcf_folder"], "{sample}_annotated.vcf"))
    shell:
        "zcat {input} > {output}"

rule annotate_vcf:
    input:
        vcf = os.path.join(config["vcf_folder"], "{sample}_annotated.vcf"),
        expr = config["expression_file"]
    output:
        os.path.join(config["vcf_folder"], "{sample}_annotated.gx.vcf")
    params:
        sample = "{sample}"
    shell:
        "vcf-expression-annotator {input.vcf} {input.expr} -s {params.sample} "
        "custom gene --id-column gene_name --expression-column {params.sample} "
        "-o {output} --ignore-ensembl-id-version"
