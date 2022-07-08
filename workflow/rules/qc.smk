# first rule: fastq short reads qc
rule fastqc_untrimmed:
    input:
        fq = lambda wildcards: samples.at[wildcards.sample, 'fq'+ wildcards.number]
    output:
        html="results/qc/fastq/fastqc/{sample}_{number}_fastqc.html",
        zip="results/qc/fastq/fastqc/{sample}_{number}_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    log:
        "results/logs/qc/untrimmed/fastqc/{sample}_{number}.log"
    shell:
        "fastqc {input.fq} -o results/qc/fastq &> {log}"

rule fastqc_trimmed:
    input:
        "results/fastq/trimmed/{sample}_{number}_{paired}.fastq.gz"
    output:
        html="results/qc/trimmed_short/{sample}_{number}_{paired}_fastqc.html",
        zip="results/qc/trimmed_short/{sample}_{number}_{paired}_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    log:
        "results/logs/qc/trimmed/{sample}_{number}_{paired}.log"
    shell:
        "fastqc {input} -o results/qc/trimmed &> {log}"




rule longqc_untrimmed:
    input:
        fq = lambda wildcards: samples.at[wildcards.sample, 'ONT']
    output:
        out_dir = directory("results/qc/longqc/{sample}")
    conda:
        "../envs/longqc.yaml"
    log: 
        "results/logs/qc/untrimmed/longqc/{sample}.log"
    params:
        preset = config["longqc"]["preset"],
        extra = config["longqc"]["extra"]
    shell:
        "python longQC.py sampleqc -x {params.preset} {params.extra} -o {output.out_dir}  -p {threads} {input.fq} &> {log}"

rule longqc_trimmed:
    input:
        "results/fastq_long/trimmed_ONT/{sample}.fastq.gz"  # make sure trimmed long reads are there 
    output:
        out_dir = directory("results/qc/trimmed_ONT/{sample}")
    conda:
         "../envs/longqc.yaml"
    log: 
        "results/logs/qc/trimmed_ONT/longqc/{sample}_trimmed.log"
    params:
        preset = config["longqc"]["preset"],
        extra = config["longqc"]["extra"]
    threads: 4  # TODO: add to the command
    shell:
        "python longQC.py sampleqc -x {params.preset} {params.extra} -o {output.out_dir} -p {threads} {input.fq} &> {log}"
    
    


rule trimmomatic:
    input:
        r1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
    output:
        r1_p= "results/fastq/trimmed_short/{sample}_1_P.fastq.gz", r1_u = "results/fastq/trimmed_short/{sample}_1_UP.fastq.gz",
        r2_p= "results/fastq/trimmed_short/{sample}_2_P.fastq.gz", r2_u = "results/fastq/trimmed_short/{sample}_2_UP.fastq.gz"
    params:
        trailing = config["trimmomatic"]['trailing'],
        illuminaclip = ':'.join(config["trimmomatic"]["illuminaclip"].values()),
        extra = config["trimmomatic"]['extra']
    conda:
        "../envs/qc.yaml"
    log:
        "results/logs/trimmomatic/{sample}.log"
    threads: 4
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1_p} {output.r1_u} {output.r2_p} {output.r2_u} TRAILING:{params.trailing} ILLUMINACLIP:{params.illuminaclip} -threads {threads} {params.extra} -trimlog {log}"


rule porechop:  # trimming of ONT
    input:
        lambda wildcards: samples.at[wildcards.sample, 'ONT']
    output:
        "results/fastq/trimmed_ONT/{sample}.fastq.gz"
    params:
        extra = config["porechop"]["extra"]
    conda:
        "../envs/porechop.yaml"
    log: 
        "results/logs/porechop/{sample}.log"
    threads: 8
    shell:
        "porechop -i {input} -o {output} --threads {threads} --format fastq.gz {params.extra} &> {log}"


rule multiqc:
    input:
        fqc=expand("results/qc/fastq/{sample}_{number}_fastqc.html",sample=IDS,number=['1', '2']),
        tqc=expand("results/qc/trimmed_short/{sample}_{number}_{paired}_fastqc.html",sample=IDS,number=['1', '2'],paired=['P', 'UP'])
            if config['qc']['short'] != 'False' and config['trimming']['short'] != 'False' else [],
        longqc=expand("results/qc/trimmed_ONT/{sample}.html") if config["qc"]["long"] != 'False' and config['trimming']['long'] != 'False' else []# TODO: check if correct format
    output:
        "results/qc/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "results/logs/multiqc/multiqc.log"
    threads: 1
    params:
        extra=config["multiqc"]["extra"]
    shell:
        "multiqc results/qc {params.extra} -o results/qc &> {log}" # fixed 1 mistake that we did not submit previously, which probably caused the output to the shell.
