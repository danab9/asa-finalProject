rule assembly:
    """
    Performs the (hybrid) assembly using Unicycler. 
    As input it takes short reads, raw, decontaminated and/or trimmed, depending on the steps the user set in the configuration file. 
    Similarly, it takes raw or preprocessed long reads, depending on whether the user wants to perform a hybrid assembly. 
    If a short read-only de novo assembly is performed, then Spades is used. 
    """
    input:
        short_r1 = (lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config["trimming"]["short"]=='False' and config["decontamination"] != "True" else "results/fastq/decontaminated/short/{sample}_1.fq" if config["decontamination"] == "True" else "results/fastq/trimmed/short/{sample}_1_P.fastq.gz",
        short_r2 =  (lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config["trimming"]["short"]=='False' and config["decontamination"] != "True" else "results/fastq/decontaminated/short/{sample}_2.fq" if config["decontamination"] == "True" else "results/fastq/trimmed/short/{sample}_2_P.fastq.gz",
        long = ((lambda wildcards: samples.at[wildcards.sample, 'ONT']) if config["trimming"]["long"]=='False' and config["decontamination"] != "True" else"results/fastq/decontaminated/long/{sample}.fq" if config["decontamination"] == "True" else "results/fastq/trimmed/long/{sample}.fastq.gz") if config["hybrid_assembly"] == "True" else [],
    output:
        genome = "results/assembly/{sample}/assembly.fasta",
        dir = directory("results/assembly/{sample}"),
    conda:
        "../envs/assembly.yaml"
    log:
        "results/logs/assembly/{sample}.log"
    threads: 10
    params:
        extra = config["unicycler"]["extra"],
        hybrid_assembly = "-l" if config["hybrid_assembly"] == "True" else ""
    shell:
        "unicycler -1 {input.short_r1} -2 {input.short_r2} {params.hybrid_assembly} {input.long} -o {output.dir} -t {threads} {params.extra} &> {log}" 


rule combine_files:  
    """
    Puts each genomes in a seperate folder, so that they are easily accessible for the user. 
    """
    input:
        genome = "results/assembly/{sample}/assembly.fasta"
    output:
        genome = "results/genomes/{sample}.fasta"
    shell:
        "cp {input.genome} {output.genome}"
