rule fastqc_untrimmed:
    """
    Quality Control using FastQC for short reads
    Optional: if qc:short is set to 'True' in the configuration file. 
    """
    input:
        fq = lambda wildcards: samples.at[wildcards.sample, 'fq' + wildcards.number]
    output:
        html= "results/qc/untrimmed/short/{sample}_{number}_short_untrimmed_fastqc.html",
        zip= "results/qc/untrimmed/short/{sample}_{number}_short_untrimmed_fastqc.zip", 
    conda:
        "../envs/qc.yaml"
    params:
        output_filename_html = lambda wildcards, input: "results/qc/untrimmed/short/" + os.path.basename(input.fq).split(".")[0] + "_fastqc.html", #Rename the samples, such that they are pretty and comprehensive when viewing the MultiQC
        output_filename_zip = lambda wildcards, input: "results/qc/untrimmed/short/" + os.path.basename(input.fq).split(".")[0] + "_fastqc.zip", 
    log:
        "results/logs/qc/untrimmed/short/{sample}_{number}.log" 
    shell:
        """ 
        mkdir -p results/qc/untrimmed/short/
        fastqc {input.fq} -o results/qc/untrimmed/short &> {log}
        mv {params.output_filename_html} {output.html} \\
            && mv  {params.output_filename_zip} {output.zip}
        """ 

rule fastqc_trimmed:
    """
    Quality Control using FastQC for trimmed short reads
    Optional: if qc:short and trimming:short are set to 'True' in the configuration file. 
    """
    input:
        trimmed="results/fastq/trimmed/short/{sample}_{number}_{paired}.fastq.gz",
    output:
        html="results/qc/trimmed/short/{sample}_{number}_short_trimmed_{paired}_fastqc.html", 
        zip="results/qc/trimmed/short/{sample}_{number}_short_trimmed_{paired}_fastqc.zip", 
    conda:
        "../envs/qc.yaml"
    log:
        "results/logs/qc/trimmed/short/{sample}_{number}_{paired}.log"
    shell:
        """
        fastqc {input.trimmed} -o results/qc/trimmed/short &> {log}
        mv results/qc/trimmed/short/{wildcards.sample}_{wildcards.number}_{wildcards.paired}_fastqc.html results/qc/trimmed/short/{wildcards.sample}_{wildcards.number}_{wildcards.paired}_short_trimmed_fastqc.html
        mv results/qc/trimmed/short/{wildcards.sample}_{wildcards.number}_{wildcards.paired}_fastqc.zip results/qc/trimmed/short/{wildcards.sample}_{wildcards.number}_{wildcards.paired}_short_trimmed_fastqc.zip
        """
    
rule longqc_untrimmed:
    """
    Quality Control using LongQC for long reads
    Optional: if qc:long is set to 'True' in the configuration file. 
    """
    input:
        fq = lambda wildcards: samples.at[wildcards.sample, 'ONT'],
        touch_file="results/installations/make_longqc.done", # Triggers the installation of longQC. 
    output:
        html = "results/qc/untrimmed/long/{sample}/web_summary.html"
    conda:
        "../envs/longqc.yaml"
    log: 
        "results/logs/qc/untrimmed/long/{sample}.log"
    threads: 8
    params:
        preset = config["longqc"]["preset"],
        extra = config["longqc"]["extra"],
        out_dir = "results/qc/untrimmed/long/{sample}",
    shell:
        """ 
        rm -f -r {params.out_dir}
        python results/installations/LongQC/longQC.py sampleqc -x {params.preset} {params.extra} -o {params.out_dir} -p {threads} {input.fq} &> {log}
        """ # Existing folder should be removed to prevent LongQC from throwing errors. 

rule longqc_trimmed:
    """
    Quality Control using LongQC for trimmed long reads
    Optional: if qc:long and trimming:short are set to 'True' in the configuration file. 
    """
    input:
        trimmed_read = "results/fastq/trimmed/long/{sample}.fastq.gz",
        touch_file="results/installations/make_longqc.done", # Triggers the installation of longQC.   
    output:
        html = "results/qc/trimmed/long/{sample}/web_summary.html"
    conda:
         "../envs/longqc.yaml"
    log: 
        "results/logs/qc/trimmed/long/{sample}_trimmed.log"
    params:
        preset = config["longqc"]["preset"],
        extra = config["longqc"]["extra"],
        out_dir = "results/qc/trimmed/long/{sample}",
    threads: 8
    shell:
        """
        rm -f -r {params.out_dir}
        python results/installations/LongQC/longQC.py sampleqc -x {params.preset} {params.extra} -o {params.out_dir} -p {threads} {input.trimmed_read} &> {log}
        """
    
    
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
        untrimmed_short = expand("results/qc/untrimmed/short/{sample}_{number}_short_untrimmed_fastqc.html", sample=IDS,number=['1', '2']) 
            if config["qc"]["short"] == "True" else [],
        trimmed_short=expand("results/qc/trimmed/short/{sample}_{number}_short_trimmed_{paired}_fastqc.html",sample=IDS,number=['1', '2'],paired=['P'])
            if config['qc']['short'] == 'True' and config['trimming']['short'] == 'True' else [], # Unpaired samples for read 2 result don't exist. Since we are also not involved in the later analysis. 
        untrimmed_long = expand("results/qc/untrimmed/long/{sample}/web_summary.html", sample=IDS) 
            if config["qc"]["long"] =="True" else [],
        trimmed_long = expand("results/qc/trimmed/long/{sample}/web_summary.html", sample=IDS) 
            if config["qc"]["long"] =="True" and config['trimming']['long'] == 'True' else [],
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
        "multiqc results/qc {params.extra} -o results/qc -f --fn_as_s_name &> {log}" #-f is used to overwrite existing reports, to prevent Snakemake from resulting in an error. 