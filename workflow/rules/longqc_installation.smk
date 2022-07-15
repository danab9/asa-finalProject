"""
Installation of longQC to allow quality control of long ONT reads.
"""

rule clone_longqc:
    """
    Cloning longQC from github: https://github.com/yfukasawa/LongQC.git
    """
    # no input needed 
    output:
        python_file = "results/installations/LongQC/longQC.py"
    log:
        "results/logs/longqc_installation/git_clone.log"
    threads: 1
    shell:
        "git clone https://github.com/yfukasawa/LongQC.git results/installations/LongQC &> {log}" #TODO: make sure it doesn't produce output. 

rule make_longqc:
    """
    Run make commmand in LongQC/minimap2-coverage/ as part of LongQC installation.
    """
    input:
        python_file = "results/installations/LongQC/longQC.py"
    output:
        touch_file=temp(touch("results/installations/make_longqc.done"))
    log:
        "results/logs/longqc_installation/make.log"
    threads: 1
    shell:
        "make -C results/installations/LongQC/minimap2-coverage &> {log}" 