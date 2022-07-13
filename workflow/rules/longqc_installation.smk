"""
Installation of longQC to allow quality control of long ONT reads.
"""

rule clone_longqc:
    """
    Cloning longQC from github: https://github.com/yfukasawa/LongQC.git
    """
    # no input needed 
    output:
        "results/installations/LongQC/longqc.py"
    log:
        "results/logs/longqc_installation.log"
    threads: 1
    shell:
        """
        git clone https://github.com/yfukasawa/LongQC.git results/installations
        """

rule make_longqc:
    """
    run make commmand in LongQC/minimap2-coverage/ as part of LongQC installation.
    """
    input:
        longqc= "results/installations/LongQC/longqc.py",
        minimap="results/installations/LongQC/minimap2-coverage/minimap2-coverage.c"
    output:
        touch_file=touch("make_longqc.done"),
        make = "results/installations/LongQC/minimap2-coverage/minimap2-coverage.o"
    log:
    threads: 1
    shell:
        "make -C results/installations/LongQC/minimap2-coverage"