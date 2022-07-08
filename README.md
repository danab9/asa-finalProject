# asa-finalProject

# QC
## Short reads trimming with `trimmomatic`
To cut adapter from short Illumina reads, add the path to the adapter fasta file in the configuration file, under trimmomatic, illuminaclip, file. Complete the other illuminaclip values as suits your research. 

See full manual by `trimmomatic` [here](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).


## ONT trimming with `porechop`
To add specific adapter to be removed from the start of the read, add `--start_adapt` followed by your adapter sequence to the configuration file, as the extra configuration. 

For example: `extra = "--start_adapt ACGCTAGCATACGT"`

For more options, see the full manual by `porechop` [here](https://github.com/rrwick/Porechop).
