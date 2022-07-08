rule mlst:
    input:
        reference = config["contamination_reference"]["long"],
        long = "results/fastq/trimmed/{sample}_1_P.fastq.gz" if config["trimming"]["long"]=='True' else (lambda wildcards: samples.at[wildcards.sample, 'ONT']),
    output:
        bam = "results/bam_decontaminated/long/{sample}.bam"
    log:
        "results/logs/artificialreference/{sample}.log"
    threads: 4
    conda:
        "../envs/decontamination.yaml"
    shell:
        "minimap2 -t {threads} -a {input.reference} {input.long} 2> {log} | samtools view -b - > {output.bam} 2> {log}"
