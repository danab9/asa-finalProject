# asa-finalProject

# QC
## ONT trimming with "porechop" 
To add specific adapter to be removed from the start of the read, add `--start_adapt` followed by your adapter sequence to the configuration file, as the extra configuration. 

For example: `extra = "--start_adapt ACGCTAGCATACGT"`

for more options see the full manual by `porechop` [here](https://github.com/rrwick/Porechop)