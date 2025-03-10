# Snakefile
configfile: "config.yaml"

rule all:
    input:
        expand(f"{config['annotated_dir']}/{{sample}}_annotated.vcf.gz", sample=config['samples']['tumor']),
        expand(f"{config['hla_dir']}/estimation/{{sample}}.result", sample=config['samples']['normal'])

# --- Step 1/2/3: Quality Control ---
rule fastqc_raw:
    input:
        expand(f"{config['input_dir']}/{{sample}}_{read}.fastq.gz", read=["R1", "R2"])
    output:
        html = expand(f"{config['input_dir']}/{{sample}}_{read}_fastqc.html", read=["R1", "R2"]),
        zip = expand(f"{config['input_dir']}/{{sample}}_{read}_fastqc.zip", read=["R1", "R2"])
    threads: config['threads']
    shell:
        "fastqc {input} -t {threads} -o {config['input_dir']}"

rule multiqc_raw:
    input:
        expand(f"{config['input_dir']}/{{sample}}_{read}_fastqc.zip", read=["R1", "R2"])
    output:
        f"{config['input_dir']}/multiqc_report.html"
    shell:
        "multiqc {config['input_dir']} -o {config['input_dir']}"

# --- Step 2: Trimming ---
rule trimmomatic:
    input:
        r1 = f"{config['input_dir']}/{{sample}}_R1.fastq.gz",
        r2 = f"{config['input_dir']}/{{sample}}_R2.fastq.gz"
    output:
        r1_paired = temp(f"{config['trimmed_dir']}/{{sample}}_R1_paired.fastq.gz"),
        r1_unpaired = temp(f"{config['trimmed_dir']}/{{sample}}_R1_unpaired.fastq.gz"),
        r2_paired = temp(f"{config['trimmed_dir']}/{{sample}}_R2_paired.fastq.gz"),
        r2_unpaired = temp(f"{config['trimmed_dir']}/{{sample}}_R2_unpaired.fastq.gz")
    log: f"{config['trimmed_dir']}/logs/{{sample}}_trim.log"
    threads: config['threads']
    shell:
        """
        trimmomatic PE -threads {threads} {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} \
        {output.r2_paired} {output.r2_unpaired} \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:151 2> {log}
        """

# --- StepTrim Quality Control ---

rule fastqc_trimmed:
    input:
        expand(f"{config['trimmed_dir']}/{{sample}}_{read}._paired.fastq.gz", read=["R1", "R2"])
    output:
        html = expand(f"{config['trimmed_dir']}/{{sample}}_{read}_paired.fastqc.html", read=["R1", "R2"]),
        zip = expand(f"{config['trimmed_dir']}/{{sample}}_{read}_paired.fastqc.zip", read=["R1", "R2"])
    threads: config['threads']
    shell:
        "fastqc {input} -t {threads} -o {config['trimmed_dir']}"

rule multiqc_trimmed:
    input:
        expand(f"{config['trimmed_dir']}/{{sample}}_{read}_paired.fastqc.zip", read=["R1", "R2"])
    output:
        f"{config['trimmed_dir']}/multiqc_report.html"
    shell:
        "multiqc {config['trimmed_dir']} -o {config['trimmed_dir']}"

# --- Step 4: Alignment ---
rule bwa_index:
    input: config['reference']
    output: touch(f"{config['reference']}.bwt")
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        r1 = rules.trimmomatic.output.r1_paired,
        r2 = rules.trimmomatic.output.r2_paired,
        idx = rules.bwa_index.output
    output:
        temp(f"{config['mapped_dir']}/{{sample}}.sam")
    log: f"{config['mapped_dir']}/logs/{{sample}}_align.log"
    threads: config['threads']
    shell:
        "bwa mem -t {threads} {config['reference']} {input.r1} {input.r2} > {output} 2> {log}"

rule sam_to_bam:
    input: f"{config['mapped_dir']}/{{sample}}.sam"
    output: temp(f"{config['mapped_dir']}/{{sample}}.bam")
    threads: config['threads']
    shell:
        "samtools view -@ {threads} -b {input} > {output}"

rule sort_bam:
    input: rules.sam_to_bam.output
    output: f"{config['sorted_dir']}/{{sample}}_sorted.bam"
    threads: config['threads']
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule index_bam:
    input: rules.sort_bam.output
    output: f"{config['sorted_dir']}/{{sample}}_sorted.bam.bai"
    shell:
        "samtools index {input}"

# --- Step 5: Mark Duplicates ---
rule mark_duplicates:
    input: rules.sort_bam.output
    output:
        bam = f"{config['sorted_dir']}/{{sample}}_dupmarked.bam",
        metrics = f"{config['sorted_dir']}/{{sample}}_metrics.txt"
    log: f"{config['sorted_dir']}/logs/{{sample}}_dups.log"
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT 2> {log}"

# --- Step 6: Add Read Groups ---
rule add_read_groups:
    input: rules.mark_duplicates.output.bam
    output: f"{config['rgbam_dir']}/{{sample}}_rg.bam"
    params:
        RGID = "HV3HWDSXY.4",
        RGLB = "library1",
        RGPL = "illumina",
        RGPU = "unit1",
        RGSM = "{sample}"
    shell:
        """
        picard AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        RGID={params.RGID} \
        RGLB={params.RGLB} \
        RGPL={params.RGPL} \
        RGPU={params.RGPU} \
        RGSM={params.RGSM}
        """

# --- Step 7: Base Quality Recalibration ---
rule base_recalibrator:
    input:
        bam = rules.add_read_groups.output,
        known = config['known_sites']
    output:
        table = temp(f"{config['rgbam_dir']}/{{sample}}_recal_data.table")
    log: f"{config['rgbam_dir']}/logs/{{sample}}_bqsr.log"
    shell:
        """
        gatk BaseRecalibrator \
        -R {config['reference']} \
        -I {input.bam} \
        --known-sites {input.known} \
        -O {output.table} 2> {log}
        """

rule apply_bqsr:
    input:
        bam = rules.add_read_groups.output,
        table = rules.base_recalibrator.output.table
    output:
        bam = f"{config['rgbam_dir']}/{{sample}}_recal.bam"
    shell:
        """
        gatk ApplyBQSR \
        -R {config['reference']} \
        -I {input.bam} \
        --bqsr-recal-file {input.table} \
        -O {output.bam}
        """

# --- Step 8: Metrics Collection ---
rule collect_metrics:
    input: rules.apply_bqsr.output.bam
    output:
        alignment = f"{config['rgbam_dir']}/{{sample}}_alignmetrics.txt",
        insert = f"{config['rgbam_dir']}/{{sample}}_insertmetrics.txt"
    shell:
        """
        gatk CollectAlignmentSummaryMetrics -I {input} -R {config['reference']} -O {output.alignment}
        gatk CollectInsertSizeMetrics -I {input} -O {output.insert}
        """

# --- Step 9: Variant Calling ---
rule mutect2_normal:
    input: 
        bam = expand(f"{config['rgbam_dir']}/{{normal}}_recal.bam", normal=config['samples']['normal'])
    output:
        vcf = f"{config['variants_dir']}/normal/{{normal}}_normal.vcf.gz"
    shell:
        """
        gatk Mutect2 \
        -R {config['reference']} \
        -I {input.bam} \
        -max-mnp-distance 0 \
        -O {output.vcf}
        """

rule genomics_db_import:
    input:
        expand(f"{config['variants_dir']}/normal/{{normal}}_normal.vcf.gz", normal=config['samples']['normal'])
    output:
        directory(f"{config['pon_dir']}")
    shell:
        """
        gatk GenomicsDBImport \
        -R {config['reference']} \
        -L {config['exome_bed']} \
        --genomicsdb-workspace-path {output} \
        {input}
        """

rule create_pon:
    input:
        db = rules.genomics_db_import.output
    output:
        vcf = f"{config['variants_dir']}/pon/pon.vcf.gz"
    shell:
        """
        gatk CreateSomaticPanelOfNormals \
        -R {config['reference']} \
        -V gendb://{input.db} \
        -O {output.vcf}
        """

rule mutect2_tumor:
    input:
        tumor = f"{config['rgbam_dir']}/{{tumor}}_recal.bam",
        normal = f"{config['rgbam_dir']}/{{normal}}_recal.bam",
        pon = rules.create_pon.output.vcf
    output:
        vcf = f"{config['variants_dir']}/tumor/{{tumor}}_somatic.vcf.gz"
    params:
        normal_sample = lambda wildcards: str(int(wildcards.tumor) + 1)
    shell:
        """
        gatk Mutect2 \
        -R {config['reference']} \
        -I {input.tumor} \
        -I {input.normal} \
        --panel-of-normals {input.pon} \
        -O {output.vcf}
        """

# --- Step 10: Variant Normalization ---
rule normalize_vcf:
    input: rules.mutect2_tumor.output.vcf
    output: f"{config['variants_dir']}/normalized/{{tumor}}_normalized.vcf.gz"
    shell:
        "vt normalize {input} -r {config['reference']} -n -o {output}"

# --- Step 11: Variant Filtering ---
rule filter_variants:
    input: rules.normalize_vcf.output
    output: 
        vcf = f"{config['variants_dir']}/filtered/{{tumor}}_filtered.vcf.gz",
        idx = f"{config['variants_dir']}/filtered/{{tumor}}_filtered.vcf.gz.tbi"
    shell:
        """
        gatk SelectVariants \
        -R {config['reference']} \
        -V {input} \
        -select-type SNP \
        -select-type INDEL \
        -O {output.vcf}
        """

# --- Step 12: Annotation ---
rule vep_annotation:
    input: rules.filter_variants.output.vcf
    output: f"{config['annotated_dir']}/{{tumor}}_annotated.vcf.gz"
    params:
        gff = "/user/home/bcancer/ref/Homo_sapiens.GRCh38.110.gff3.gz",
        plugins = "/user/home/bcancer/software/VEP_plugins"
    shell:
        """
        vep --hgvs --fasta {config['reference']}.gz \
        --gff {params.gff} \
        -i {input} \
        -o {output} \
        --vcf \
        --plugin Downstream \
        --plugin Frameshift \
        --plugin Wildtype \
        --coding_only \
        --no_intergenic \
        --dir_plugins {params.plugins}
        """

# --- Step 13: HLA Typing ---
rule hlahd_typing:
    input:
        r1 = f"{config['trimmed_dir']}/{{sample}}_R1_paired.fastq.gz",
        r2 = f"{config['trimmed_dir']}/{{sample}}_R2_paired.fastq.gz"
    output:
        f"{config['hla_dir']}/estimation/{{sample}}.result"
    params:
        freq = f"{config['hla_dir']}/freq_data",
        gene_split = f"{config['hla_dir']}/HLA_gene.split.txt",
        dict_dir = f"{config['hla_dir']}/dictionary"
    threads: config['threads']
    shell:
        """
        zcat {input.r1} > {config['hla_dir']}/data/{wildcards.sample}_R1.fastq
        zcat {input.r2} > {config['hla_dir']}/data/{wildcards.sample}_R2.fastq
        hlahd.sh -t {threads} \
        {config['hla_dir']}/data/{wildcards.sample}_R1.fastq \
        {config['hla_dir']}/data/{wildcards.sample}_R2.fastq \
        {params.gene_split} \
        {params.dict_dir} \
        {wildcards.sample} \
        {config['hla_dir']}/estimation
        """