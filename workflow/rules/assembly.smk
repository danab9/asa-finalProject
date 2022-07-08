rule assembly:
    input:
        short_r1 = (lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config["trimming"]["short"]=='False' and config["decontamination"] != "True" else "results/fastq/decontaminated/{sample}_1.fq" if config["decontamination"] == "True" else "results/fastq/trimmed/{sample}_1_P.fastq.gz",
        short_r2 = (lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config["trimming"]["short"]=='False' and config["decontamination"] != "True" else "results/fastq/decontaminated/{sample}_2.fq" if config["decontamination"] == "True" else "results/fastq/trimmed/{sample}_2_P.fastq.gz",
        long = ((lambda wildcards: samples.at[wildcards.sample, 'long']) if config["trimming"]["long"]=='False' and config["decontamination"] != "True" else"results/fastq/decontaminated/{sample}_2.fq" if config["decontamination"] == "True" else "results/fastq/trimmed/{sample}_2_P.fastq.gz") if config["hybrid_assembly"] == "True" else ""
    output:
        genome = "results/assembly/{sample}/assembly.fasta",
        dir = "results/assembly/{sample}",
    conda:
        "../envs/assembly.yaml"
    log:
        "results/logs/assembly_{sample}.log"
    threads: 10
    params:
        extra = config["unicycler"]["extra"],
        hybrid_assembly = "-l" if config["hybrid_assembly"] == "True" else ""
    shell:
        "unicycler -1 {input.short_r1} -2 {input.short_r2} {params.hybrid_assembly} {input.long} -o {output.dir} -t {threads} {params.extra} &> {log}" 