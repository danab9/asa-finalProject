rule quast:
    input:
        "results/assembly/{sample}/assembly.fasta"
    output:
        "results/qc/assembly/quast/{sample}/report.txt"
    log:
        "results/logs/quast/{sample}.log"
    conda:
        "../envs/quast.yaml"
    params:
        extra = config['quast']['extra']
    threads:
        4
    shell:
        "quast -o results/qc/assembly/quast/{wildcards.sample} -t {threads} {params.extra} {input} &> {log}"

busco_lin_name = os.path.basename(config["busco_params"]["lineage"])
rule busco:
    input:
        "results/assembly/{sample}/assembly.fasta"
    output:
        #dir=directory("results/qc/assembly/BUSCO/"),
        file = "results/qc/assembly/BUSCO/{sample}/short_summary.specific." + busco_lin_name + ".{sample}.json"
    log:
        "results/logs/busco/{sample}.log"
    conda:
        "../envs/busco.yaml"
    params:
        lineage = config["busco_params"]["lineage"],
        offline = "--offline" if config["busco_params"]["offline"]=='True' else "",
        extra = config["busco_params"]["extra"]
    threads: 8
    shell:
        "busco -i {input} -o {wildcards.sample} --out_path results/qc/assembly/BUSCO/ -l {params.lineage} -m genome {params.offline} {params.extra} --cpu {threads} -f &> {log}" 