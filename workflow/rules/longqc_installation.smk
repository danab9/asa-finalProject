"""
Installation of longQC to allow quality control of long ONT reads.
"""

rule clone_longqc:
    """
    Cloning longQC from github: https://github.com/yfukasawa/LongQC.git
    """
    # no input needed 
    output:
        "results/installations/LongQC/longQC.py"
    log:
        "results/logs/longqc_installation/git_clone.log"
    threads: 1
    shell:
        """
        git clone https://github.com/yfukasawa/LongQC.git results/installations/LongQC 2> {log}
        """

rule make_longqc:
    """
    run make commmand in LongQC/minimap2-coverage/ as part of LongQC installation.
    """
    input:
        "results/installations/LongQC/longQC.py"
    output:
        touch_file=temp(touch("results/make_longqc.done"))
    log:
        "results/logs/longqc_installation/make.log"
    threads: 1
    shell:
        "make -C results/installations/LongQC/minimap2-coverage 2> {log}"