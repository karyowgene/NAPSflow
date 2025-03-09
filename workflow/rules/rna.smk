configfile: "../config/rna_config.yaml"

import os
import glob

def get_samples():
    source_dir = config["source_dir"]
    samples = {}
    for r1_path in glob.glob(os.path.join(source_dir, "*_R1_001.fastq.gz")):
        filename = os.path.basename(r1_path)
        parts = filename.split('_')
        if len(parts) < 3:
            continue
        sample_id = parts[2]
        r2_path = os.path.join(source_dir, filename.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
        if os.path.exists(r2_path):
            samples[sample_id] = (r1_path, r2_path)
    return samples

samples_dict = get_samples()
SAMPLES = list(samples_dict.keys())

rule rnaall:
    input:
        # MultiQC reports for raw and trimmed FASTQC
        os.path.join(config["raw_dir"], "multiqc_report.html"),
        os.path.join(config["trimmed_dir"], "multiqc_report.html"),
        # Kallisto quantification output for each sample
        expand(os.path.join(config["kallisto_quants"], "{sample}"), sample=SAMPLES)

rule rename:
    input:
        r1 = lambda wildcards: samples_dict[wildcards.sample][0],
        r2 = lambda wildcards: samples_dict[wildcards.sample][1]
    output:
        r1 = os.path.join(config["raw_dir"], "{sample}_R1.fastq.gz"),
        r2 = os.path.join(config["raw_dir"], "{sample}_R2.fastq.gz")
    shell:
        "mv {input.r1} {output.r1}; mv {input.r2} {output.r2}"

rule fastqc_raw:
    input:
        r1 = os.path.join(config["raw_dir"], "{sample}_R1.fastq.gz"),
        r2 = os.path.join(config["raw_dir"], "{sample}_R2.fastq.gz")
    output:
        html = [os.path.join(config["report_dir"], "{sample}_R1_fastqc.html"),
                os.path.join(config["report_dir"], "{sample}_R2_fastqc.html")]
    threads: config["threads"]
    shell:
        "fastqc {input.r1} {input.r2} -t {threads} --outdir {config[report_dir]}"

rule multiqc_raw:
    input:
        expand(os.path.join(config["report_dir"], "{sample}_R{read}_fastqc.html"), sample=SAMPLES, read=[1,2])
    output:
        os.path.join(config["report_dir"], "multiqc_report.html")
    shell:
        "multiqc {config[report_dir]} -o {config[report_dir]}"

rule trimrna_reads:
    input:
        r1 = os.path.join(config["raw_dir"], "{sample}_R1.fastq.gz"),
        r2 = os.path.join(config["raw_dir"], "{sample}_R2.fastq.gz")
    output:
        r1_paired   = os.path.join(config["trimmed_dir"], "{sample}_R1_paired.fastq.gz"),
        r1_unpaired = os.path.join(config["trimmed_dir"], "{sample}_R1_unpaired.fastq.gz"),
        r2_paired   = os.path.join(config["trimmed_dir"], "{sample}_R2_paired.fastq.gz"),
        r2_unpaired = os.path.join(config["trimmed_dir"], "{sample}_R2_unpaired.fastq.gz")
    params:
        adapters = config["trimmomatic_adapters"]
    threads: config["threads"]
    shell:
        "mkdir -p $(dirname {output.r1_paired}) && "
        "trimmomatic PE -threads {threads} {input.r1} {input.r2} "
        "{output.r1_paired} {output.r1_unpaired} "
        "{output.r2_paired} {output.r2_unpaired} "
        "ILLUMINACLIP:{params.adapters} LEADING:3 TRAILING:3 "
        "SLIDINGWINDOW:4:15 MINLEN:151"

rule fastqc_trimmed:
    input:
        r1 = os.path.join(config["trimreport_dir"], "{sample}_R1_paired.fastq.gz"),
        r2 = os.path.join(config["trimreport_dir"], "{sample}_R2_paired.fastq.gz")
    output:
        html = [os.path.join(config["trimreport_dir"], "{sample}_R1_paired_fastqc.html"),
                os.path.join(config["trimreport_dir"], "{sample}_R2_paired_fastqc.html")]
    threads: config["threads"]
    shell:
        "fastqc {input.r1} {input.r2} -t {threads} --outdir {config[trimreport_dir]}"

rule multiqc_trimmed:
    input:
        expand(os.path.join(config["trimreport_dir"], "{sample}_R{read}_paired_fastqc.html"), sample=SAMPLES, read=[1,2])
    output:
        os.path.join(config["trimreport_dir"], "multiqc_report.html")
    shell:
        "multiqc {config[trimreport_dir]} -o {config[trimreport_dir]}"

rule download_cdna:
    output:
        os.path.join(config["rnaref_dir"], config["cdna_fa"])
    shell:
        "curl -o {output} ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/{config[cdna_fa]}"

rule build_index:
    input:
        os.path.join(config["rnaref_dir"], config["cdna_fa"])
    output:
        os.path.join(config["rnaref_dir"], config["cdna_idx"])
    threads: config["threads"]
    shell:
        "kallisto index -i {output} -t {threads} {input}"

rule kallisto_quant:
    input:
        r1 = os.path.join(config["trimmed_dir"], "{sample}_R1_paired.fastq.gz"),
        r2 = os.path.join(config["trimmed_dir"], "{sample}_R2_paired.fastq.gz"),
        index = os.path.join(config["rnaref_dir"], config["cdna_idx"])
    output:
        directory(os.path.join(config["kallisto_quants"], "{sample}"))
    threads: config["threads"]
    shell:
        "mkdir -p {output} && "
        "kallisto quant -i {input.index} -o {output} -t {threads} {input.r1} {input.r2}"
