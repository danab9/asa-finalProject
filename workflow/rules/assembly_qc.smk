rule quast:
    input:
        "results/assembly/{sample}/assembly.fasta"
    output:
        "results/qc/assembly/{sample}/report.txt"
    log:
        "results/logs/quast/{sample}.log"
    conda:
        "../envs/quast.yaml"
    params:
        extra = config['quast']['extra']
    threads:
        4
    shell:
        "python quast.py -o results/qc/assembly/{sample} -t {threads} {params.extra} {input} &> {log}"


rule busco:
    input:
        "results/assembly/{sample}/assembly.fasta"
    output:
        directory("results/qc/assembly/{sample}"), # TODO: add file according to output
    log:
        "results/logs/busco/{sample}.log"
    conda:
        "../envs/busco.yaml"
    params:
        lineag = config["busco"]["lineage"],
        extra = config["busco"]["extra"]
    threads: 8
    shell:
        "busco -i {input} -o {output} -l {params.lineage} -m genome {params.extra} --cpu {threads} &> {log}" 