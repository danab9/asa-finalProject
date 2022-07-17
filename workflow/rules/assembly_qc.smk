rule quast:
    """
    Quality Control of the assembled genomes using QUAST. 
    """
    input:
        assembly="results/genomes/{sample}.fasta"
    output:
        report="results/qc/assembly/quast/{sample}/report.txt"
    log:
        "results/logs/quast/{sample}.log"
    conda:
        "../envs/quast.yaml"
    params:
        extra = config['quast']['extra']
    threads:
        8
    shell:
        "quast -o results/qc/assembly/quast/{wildcards.sample} -t {threads} {params.extra} {input.assembly} &> {log}"

rule busco:
    """
    Quality Control of the assembled genomes using BUSCO.
    Optional: If busco is set to 'True' in the configuration file.  
    """ #TODO: explain the params 
    input:
        assembly="results/genomes/{sample}.fasta"
    output:
        file = "results/qc/assembly/BUSCO/{sample}/short_summary_{sample}.txt"
    log:
        "results/logs/busco/{sample}.log"
    conda:
        "../envs/busco.yaml"
    params:
        lineage = config["busco_params"]["lineage"], # path to lineage / name of lineage dataset of BUSCO
        offline = "--offline" if config["busco_params"]["offline"]=='True' else "",  # offline=True prevents BUSCO from searching online
        extra = config["busco_params"]["extra"]
    threads: 8
    shell:
        """
        busco -i {input.assembly} -o {wildcards.sample} --out_path results/qc/assembly/BUSCO/ -l {params.lineage} -m genome {params.offline} {params.extra} --cpu {threads} -f &> {log}
        mv results/qc/assembly/BUSCO/{wildcards.sample}/short_summary.specific.{params.lineage}.{wildcards.sample}.txt results/qc/assembly/BUSCO/{wildcards.sample}/{wildcards.sample}.txt
        mv results/qc/assembly/BUSCO/{wildcards.sample}/run_{params.lineage}/short_summary.txt results/qc/assembly/BUSCO/{wildcards.sample}/short_summary_{wildcards.sample}.txt
        """ # Moves output files, such that they appear by sample name and not double in the MultiQC report


rule multiqc_assembly:
    """
    Accumulation of Quality Control BUSCO and Quast reports of the assembly
    """
    input:
        quast = expand("results/qc/assembly/quast/{sample}/report.txt", sample=IDS),
        busco = expand("results/qc/assembly/BUSCO/{sample}/short_summary_{sample}.txt", sample=IDS)
    output:
        report = "results/qc/assembly/multiqc_report.html"
    log:
        "results/logs/assembly_multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    params:
        extra=config["multiqc"]["extra"]
    shell:
        "multiqc results/qc/assembly {params.extra} -o results/qc/assembly -f &> {log}" # -f replaces existing multiQC
