rule fastqc_untrimmed:
    """
    Quality Control using FastQC for short reads
    Optional: if qc:short is set to 'True' in the configuration file. 
    """
    input:
        fq = lambda wildcards: samples.at[wildcards.sample, 'fq'+ wildcards.number]
    output:
        html= "results/qc/fastq/fastqc/{sample}_{number}_fastqc.html",
        zip= "results/qc/fastq/fastqc/{sample}_{number}_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    log:
        "results/logs/qc/untrimmed/fastqc/{sample}_{number}.log"
    shell:
        "fastqc {input.fq} -o results/qc/fastq &> {log}"

rule fastqc_trimmed:
    """
    Quality Control using FastQC for trimmed short reads
    Optional: if qc:short and trimming:short are set to 'True' in the configuration file. 
    """
    input:
        trimmed="results/fastq/trimmed/short/{sample}_{number}_{paired}.fastq.gz"
    output:
        html="results/qc/trimmed/short/{sample}_{number}_{paired}_fastqc.html",
        zip="results/qc/trimmed/short/{sample}_{number}_{paired}_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    log:
        "results/logs/qc/trimmed/{sample}_{number}_{paired}.log"
    shell:
        "fastqc {input.trimmed} -o results/qc/trimmed &> {log}"

# rule fastqc_long:
#     """
#     Quality Control using FastQC for ONTreads
#     Optional: if qc:long is set to 'True' in the configuration file. 
#     """
#     input:
#         fq = lambda wildcards: samples.at[wildcards.sample, 'ONT']
#     output:
#         html= "results/qc/fastq/fastqc/{sample}_fastqc.html",
#         zip= "results/qc/fastq/fastqc/{sample}_fastqc.zip"
#     conda:
#         "../envs/qc.yaml"
#     log:
#         "results/logs/qc/untrimmed/fastqc/{sample}_ONT.log"
#     shell:
#         "fastqc {input.fq} -o results/qc/fastq &> {log}"
    
rule longqc_untrimmed:
    """
    Quality Control using LongQC for long reads
    Optional: if qc:long is set to 'True' in the configuration file. 
    """
    input:
        fq = lambda wildcards: samples.at[wildcards.sample, 'ONT'],
        touch_file="results/make_longqc.done"
    output:
        out_dir = directory("results/qc/longqc/{sample}"),
        html = "results/qc/longqc/{sample}/web_summary.html"
    conda:
        "../envs/longqc.yaml"
    log: 
        "results/logs/qc/untrimmed/longqc/{sample}.log"
    params:
        preset = config["longqc"]["preset"],
        extra = config["longqc"]["extra"]
    threads: 8
    shell:
        "python results/installations/LongQC/longQC.py sampleqc -x {params.preset} {params.extra} -o {output.out_dir}  -p {threads} {input.fq} &> {log}"

rule longqc_trimmed:
    """
    Quality Control using LongQC for trimmed long reads
    Optional: if qc:long and trimming:short are set to 'True' in the configuration file. 
    """
    input:
        "results/fastq/trimmed/long/{sample}.fastq.gz"  # make sure trimmed long reads are there 
    output:
        out_dir = directory("results/qc/longqc/{sample}"),
        html = "results/qc/longqc/{sample}/web_summary.html"
    conda:
         "../envs/longqc.yaml"
    log: 
        "results/logs/qc/trimmed/long/longqc/{sample}_trimmed.log"
    params:
        preset = config["longqc"]["preset"],
        extra = config["longqc"]["extra"],
        out_dir = directory("results/qc/trimmed/long/{sample}"),
    threads: 8
    shell:
        "python results/installations/LongQC/longQC.py sampleqc -x {params.preset} {params.extra} -o {params.out_dir} -p {threads} {input.fq} &> {log}"
    
    
rule trimmomatic:
    """
    Trimming short reads using trimmomatic.
    Optional: if trimming:short is set to 'True' in the configuration file. 
    """
    input:
        r1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
    output:
        r1_p= "results/fastq/trimmed/short/{sample}_1_P.fastq.gz", r1_u = "results/fastq/trimmed/short/{sample}_1_UP.fastq.gz",
        r2_p= "results/fastq/trimmed/short/{sample}_2_P.fastq.gz", r2_u = "results/fastq/trimmed/short/{sample}_2_UP.fastq.gz"
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
        "trimmomatic PE {input.r1} {input.r2} {output.r1_p} {output.r1_u} {output.r2_p} {output.r2_u} TRAILING:{params.trailing} ILLUMINACLIP:{params.illuminaclip} -threads {threads} {params.extra} -trimlog {log} &> {log}" 
        #TODO: check if it doesnt produce shell output anymore. 


rule porechop:
    """
    Trimming long reads using trimmomatic.
    Optional: if trimming:long is set to 'True' in the configuration file. 
    """
    input:
        lambda wildcards: samples.at[wildcards.sample, 'ONT']
    output:
        "results/fastq/trimmed/long/{sample}.fastq.gz"
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
    """
    Combining quality control reports of FastQC and LongQC of untrimmed and trimmed short and long reads, depending on the set workflow. 
    Optional: if qc:short or qc:long is set to 'True' in the configuration file. 
    """
    input:
        fqc=expand("results/qc/fastq/fastqc/{sample}_{number}_fastqc.html",sample=IDS,number=['1', '2']),
        tqc=expand("results/qc/trimmed/short/{sample}_{number}_{paired}_fastqc.html",sample=IDS,number=['1', '2'],paired=['P', 'UP'])
            if config['qc']['short'] != 'False' and config['trimming']['short'] != 'False' else [],
        #longqc=expand("results/qc/fastq/fastqc/{sample}_fastqc.html", sample=IDS) # TODO: sample ONT wildcard?
        #    if config["qc"]["long"] != 'False' and config['trimming']['long'] != 'False' else [] # TODO: check if correct format
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
        "multiqc results/qc {params.extra} -o results/qc &> {log}" 