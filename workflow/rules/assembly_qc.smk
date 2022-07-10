rule quast:
    input:
        "results/assembly/{sample}/assembly.fasta"
    output:
        "results/qc/assembly/{sample}/quast/report.txt"
    log:
        "results/logs/quast/{sample}.log"
    conda:
        "../envs/quast.yaml"
    params:
        extra = config['quast']['extra']
    threads:
        4
    shell:
        "quast -o results/qc/assembly/{wildcards.sample} -t {threads} {params.extra} {input} &> {log}"


rule busco:
    input:
        "results/assembly/{sample}/assembly.fasta"
    output:
        dir=directory("results/qc/assembly/{sample}/BUSCO")
    log:
        "results/logs/busco/{sample}.log"
    conda:
        "../envs/busco.yaml"
    params:
        lineag = config["busco_params"]["lineage"],
        offline = "--offline" if config["busco_params"]["offline"]=='True' else "",
        extra = config["busco_params"]["extra"]
    threads: 8
    shell:
        "busco -i {input} -o {output.dir} -l {params.lineage} -m genome {params.offline} {params.extra} --cpu {threads} &> {log}" 