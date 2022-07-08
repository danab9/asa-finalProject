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
    
    


rule trimmomatic_short_reads:
    input:
        r1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
    output:
        r1_p= "results/fastq/trimmed_short/{sample}_1_P.fastq.gz", r1_u = "results/fastq/trimmed_short/{sample}_1_UP.fastq.gz",
        r2_p= "results/fastq/trimmed_short/{sample}_2_P.fastq.gz", r2_u = "results/fastq/trimmed_short/{sample}_2_UP.fastq.gz"
    params:
<<<<<<< HEAD
        trailing = config["trimmomatic_short"]['trailing'],
        illuminaclip = ':'.join(config["trimmomatic_short"]["illuminaclip"].values()),
        extra = config["trimmomatic_short"]['extra']
=======
        trailing = config["trimmomatic"]['trailing'], 
        illuminaclip = ':'.join(config["trimmomatic"]["illuminaclip"].values())
>>>>>>> 09e086efa5f5a1d0c78af1fac0707fad99997332
    conda:
        "../envs/qc.yaml"
    log:
        "results/logs/trimmomatic_short_reads/{sample}.log"
    threads: 4
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1_p} {output.r1_u} {output.r2_p} {output.r2_u} TRAILING:{params.trailing} ILLUMINACLIP:{params.illuminaclip} -threads {threads} {params.extra} -trimlog {log}"


rule trimmomatic_long_reads:
    input: 
        lambda wildcards: samples.at[wildcards.sample, 'ONT']
    output: 
        "results/fastq/trimmed_ONT/{sample}.fastq.gz"
    log: 
        "results/logs/trimmomatic_ONT/{sample}.log"
    params:
        extra = config["trimmomatic_ONT"]['extra'],
        trailing = config["trimmomatic_ONT"]['trailing']
    conda:
        "../envs/qc.yaml"
    threads: 4
    shell:
        "trimmomatic SE {inpug} {output} -threads {threads} TRAILING:{params.trailing} -trimlog {log} {params.extra}"



rule multiqc:
    input:
        fqc=expand("results/qc/fastq/{sample}_{number}_fastqc.html",sample=IDS,number=['1', '2']) if config["skip_trimming"] in['False',''] else [],
        tqc=expand("results/qc/trimmed_short/{sample}_{number}_{paired}_fastqc.html",sample=IDS,number=['1', '2'],paired=['P', 'UP'])
            if config['skip_fastQC'] in ['False',''] else [],
        longqc=expand("results/qc/trimmed_ONT/{sample}.html") # TODO: check if correct format
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
