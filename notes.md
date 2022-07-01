**Strategy 1:** assemble long reads first, then use short reads to polish  (error correction) https://www.melbournebioinformatics.org.au/tutorials/tutorials/hybrid_assembly/nanopore_assembly/

**Strategy 2:** assemble short reads first, then use long reads to scaffold (i.e.  combine) the contigs 

**Strategy 3:** create contigs for both separately, and combine the contigs, possibly using a *reference*  (https://en.wikipedia.org/wiki/Hybrid_genome_assembly image)

**Strategy 4:** hybrid de novo algorithm (i.e. HybridSpades, HASLR or Unicycler, https://www.nature.com/articles/s41597-019-0311-3)
- HASLR https://www.sciencedirect.com/science/article/pii/S2589004220305770 conda 
- HybridSpades https://academic.oup.com/bioinformatics/article/32/7/1009/1743807?login=true conda
- Unicycler https://github.com/rrwick/Unicycler *only for bacterial genomes* conda
- MuSuRCA  https://academic.oup.com/bioinformatics/article/29/21/2669/195975 conda

You can use *reference(s)* for quality check

------------------------------------------------------------------------

**Quality metrics**
- BUSCO assembly and annotation completeness https://academic.oup.com/bioinformatics/article/31/19/3210/211866?login=true
- QUAST quality assessment tool for genome assemblies https://pubmed.ncbi.nlm.nih.gov/23422339/ 
    - completeness, several metrics, e.g. No. of indels per 100 kb
    - contiguitiy, metrics based on contigs, e.g. no. of contigs. largest contigs. 

**Approaches comparison papers**
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8306402/
"for the hybrid assembly strategy, more continuous assemblies can be achieved when using long reads in conjunction with Illumina reads. This strategy initiated the hybrid assembly with high-quality Illumina short reads and filled the gaps with ONT or PacBio long reads."

**Tools comparison papers**
- https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07041-8 - **Benchmarking hybrid assembly approaches for genomic analyses of bacterial pathogens using Illumina and Oxford Nanopore sequencing**
"Unicycler performed the best for achieving contiguous genomes, closely followed by MaSuRCA, while all SPAdes assemblies were incomplete."

**Gene annotation** 

Gene prediction tools: 
- https://guides.library.yale.edu/bioinformatics/gene-prediction
- https://en.wikipedia.org/wiki/List_of_gene_prediction_software
- comparison paper: https://www.semanticscholar.org/paper/Comparative-Analysis-of-Gene-Prediction-Tools%3A-hmm-Jyoti-Saini/d2db652027e9862b9deffb558131311cadca526a
- GENEID
- AUGUSTUS: needs reference species
- BRAKER: might be possible, trains AUGUSTUS (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6635606/) 
- Prodigal - https://www.biocode.ltd/catalog-tooldata/prodigal (conda, unsupervised, gff)

-------------------------------------------------
**MLST**
https://en.wikipedia.org/wiki/Multilocus_sequence_typing
typing of multiple loci, using DNA sequences of internal fragments of multiple housekeeping genes to characterize isolates of microbial species.
- mlst https://github.com/tseemann/mlst conda

**Antibiotic resistance (CARD)**
https://anaconda.org/bioconda/rgi

**Virulence factors (VFDB)


**plasmids (PLSDB)**


