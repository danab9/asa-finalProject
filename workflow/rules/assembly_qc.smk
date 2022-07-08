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

