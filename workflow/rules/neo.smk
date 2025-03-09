import os
import glob

configfile: "../config/neo_config.yaml"

# Function to get sample names from annotated VCF files
def get_samples():
    vcf_files = glob.glob(os.path.join(config["anno_dir"], "*.gx.vcf"))
    return [os.path.basename(f).split("_annotated.gx.vcf")[0] for f in vcf_files]

SAMPLES = get_samples()

rule neo_all:
    input:
        expand(os.path.join(config["output_base_dir"], "{sample}"), sample=SAMPLES)

rule pvacseq_prediction:
    input:
        vcf = os.path.join(config["anno_dir"], "{sample}_annotated.gx.vcf")
    output:
        dir_out = directory(os.path.join(config["output_base_dir"], "{sample}")),
        final_tsv = os.path.join(config["output_base_dir"], "{sample}", "all_epitopes.tsv")
    params:
        hla_types = config["hla_types"],
        normal_sample_name = lambda wildcards: str(int(wildcards.sample) + 1),
        docker_volumes = lambda wildcards, output: (
            f"-v {config['anno_dir']}:/pvacseq_mydataall_data "
            f"-v {os.path.dirname(output[0])}:/pvacseq_outputall_mydata"
        ),
        output_dir = lambda wildcards: os.path.join(config["output_base_dir"], wildcards.sample)
    threads: config["threads"]
    shell:
        """
        mkdir -p {params.output_dir} && \
        docker run {params.docker_volumes} {config['docker_image']} \
        pvacseq run \
            /pvacseq_mydataall_data/{wildcards.sample}_annotated.gx.vcf \
            {wildcards.sample} \
            {params.hla_types} \
            {config['prediction_algorithms']} \
            /pvacseq_outputall_mydata \
            {config['other_params'].format(threads=threads)} \
            -s 500 -d 500 \
            --normal-sample-name {params.normal_sample_name} \
            {config['more_params']}
        """
