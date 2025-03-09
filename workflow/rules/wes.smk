# Snakefile
configfile: "../config/wes_config.yaml"
import pandas as pd

# Load samplesheet and pairs
samples_df = pd.read_csv("../data/samplesheet.tsv", sep="\t")
pairs_df = pd.read_csv("../data/pairs.tsv", sep="\t")

# Get sample lists
ALL_SAMPLES = samples_df["sample"].tolist()
TUMORS = pairs_df["tumor_sample"].tolist()
NORMALS = pairs_df["normal_sample"].tolist()

# Validate pairs
def get_normal(tumor):
    normal = pairs_df[pairs_df["tumor_sample"] == tumor]["normal_sample"].values
    if len(normal) == 0:
        raise ValueError(f"No normal found for tumor {tumor}")
    return normal[0]

# Final targets
rule all:
    input:
        expand(f"{config['annotated_dir']}/{{sample}}_annotated.vcf.gz", sample=TUMORS),
        f"{config['pon_dir']}/pon.vcf.gz"

# Core processing rules
rule trim_reads:
    input:
        r1 = lambda wildcards: f"{config['input_dir']}/{samples_df.loc[samples_df['sample'] == wildcards.sample, 'fastq_prefix'].values[0]}_R1.fastq.gz",
        r2 = lambda wildcards: f"{config['input_dir']}/{samples_df.loc[samples_df['sample'] == wildcards.sample, 'fastq_prefix'].values[0]}_R2.fastq.gz"
    output:
        r1_paired = temp(f"{config['trimmed_dir']}/{{sample}}_R1_paired.fastq.gz"),
        r1_unpaired = temp(f"{config['trimmed_dir']}/{{sample}}_R1_unpaired.fastq.gz"),
        r2_paired = temp(f"{config['trimmed_dir']}/{{sample}}_R2_paired.fastq.gz"),
        r2_unpaired = temp(f"{config['trimmed_dir']}/{{sample}}_R2_unpaired.fastq.gz")
    log:
        f"{config['trimmed_dir']}../logs/{{sample}}_trim.log"
    threads: config["threads"]
    shell:
        """
        trimmomatic PE -threads {threads} {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} \
        {output.r2_paired} {output.r2_unpaired} \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:151 2> {log}
        """

rule align_sort:
    input:
        r1 = rules.trim_reads.output.r1_paired,
        r2 = rules.trim_reads.output.r2_paired,
        ref = config["reference"]
    output:
        temp(f"{config['sorted_dir']}/{{sample}}_sorted.bam")
    log:
        f"{config['sorted_dir']}../logs/{{sample}}_align_sort.log"
    threads: config["threads"]
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
        samtools view -@ {threads} -b | \
        samtools sort -@ {threads} -o {output} - 2> {log}
        samtools index {output} 2>> {log}
        """

# Variant calling rules
rule mutect2_tumor_normal:
    input:
        tumor = lambda wildcards: f"{config['sorted_dir']}/{{wildcards.tumor}}_sorted.bam",
        normal = lambda wildcards: f"{config['sorted_dir']}/{{get_normal(wildcards.tumor)}}_sorted.bam",
        pon = f"{config['pon_dir']}/pon.vcf.gz"
    output:
        f"{config['variants_dir']}/{{tumor}}_somatic.vcf.gz"
    log:
        f"{config['variants_dir']}../logs/{{tumor}}_mutect2.log"
    threads: 8
    shell:
        """
        gatk --java-options '-Xmx16G' Mutect2 \
        -R {config['reference']} \
        -I {input.tumor} \
        -I {input.normal} \
        --panel-of-normals {input.pon} \
        -O {output} 2> {log}
        """

# PoN creation rules
rule mutect2_normal:
    input:
        bam = lambda wildcards: f"{config['sorted_dir']}/{{wildcards.normal}}_sorted.bam"
    output:
        f"{config['pon_dir']}/{{normal}}_normal.vcf.gz"
    log:
        f"{config['pon_dir']}../logs/{{normal}}_mutect2.log"
    threads: 8
    shell:
        """
        gatk --java-options '-Xmx16G' Mutect2 \
        -R {config['reference']} \
        -I {input.bam} \
        -max-mnp-distance 0 \
        -O {output} 2> {log}
        """

rule genomicsdb_import:
    input:
        expand(f"{config['pon_dir']}/{{normal}}_normal.vcf.gz", normal=NORMALS)
    output:
        directory(f"{config['pon_dir']}/genomicsdb_workspace")
    shell:
        """
        gatk GenomicsDBImport \
        -R {config['reference']} \
        --genomicsdb-workspace-path {output} \
        --merge-input-intervals true \
        -V {' -V '.join(input)}
        """

rule create_pon:
    input:
        rules.genomicsdb_import.output
    output:
        f"{config['pon_dir']}/pon.vcf.gz"
    shell:
        """
        gatk CreateSomaticPanelOfNormals \
        -R {config['reference']} \
        -V gendb://{input} \
        -O {output}
        """

# Annotation rule
rule vep_annotation:
    input:
        vcf = f"{config['variants_dir']}/{{tumor}}_somatic.vcf.gz",
        gff = config["gff"]
    output:
        f"{config['annotated_dir']}/{{tumor}}_annotated.vcf.gz"
    params:
        plugin_dir = config["vep_plugin_dir"]
    log:
        f"{config['annotated_dir']}../logs/{{tumor}}_vep.log"
    shell:
        """
        vep --hgvs --fasta {config['reference']} \
        --gff {input.gff} \
        --plugin Downstream --plugin Frameshift --plugin Wildtype \
        -i {input.vcf} \
        -o {output} \
        --vcf --coding_only --no_intergenic \
        --dir_plugins {params.plugin_dir} 2> {log}
        """
