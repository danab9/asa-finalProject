**Strategy 1:** assemble long reads first, then use short reads to polish  (error correction) https://www.melbournebioinformatics.org.au/tutorials/tutorials/hybrid_assembly/nanopore_assembly/

**Strategy 2:** assemble short reads first, then use long reads to scaffold (i.e.  combine) the contigs 

**Strategy 3:** create contigs for both separately, and combine the contigs, possibly using a *reference*  (https://en.wikipedia.org/wiki/Hybrid_genome_assembly image)

**Strategy 4:** hybrid de novo algorithm (i.e. HybridSpades, HASLR or Unicycler, https://www.nature.com/articles/s41597-019-0311-3)
- HASLR https://www.sciencedirect.com/science/article/pii/S2589004220305770 conda 
- HybridSpades https://academic.oup.com/bioinformatics/article/32/7/1009/1743807?login=true conda
- Unicycler https://github.com/rrwick/Unicycler *only for bacterial genomes* conda
- MuSaCrA

You can use *reference(s)* for quality check

------------------------------------------------------------------------

**Quality metrics**
- BUSCO assembly and annotation completeness https://academic.oup.com/bioinformatics/article/31/19/3210/211866?login=true
- QUAST

**Comparison papers**
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8306402/
"for the hybrid assembly strategy, more continuous assemblies can be achieved when using long reads in conjunction with Illumina reads. This strategy initiated the hybrid assembly with high-quality Illumina short reads and filled the gaps with ONT or PacBio long reads."

**Gene annotation** 

Gene prediction tools: https://guides.library.yale.edu/bioinformatics/gene-prediction
GENEID
