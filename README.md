# coronavirus_variant_caller
in-house variant caller for the Sung Lab, output used for upstream analyses. basic pipeline:

- fasta file with coronavirus sequences aligned to reference using bwa
- processed with samtools, converted to pileup format
- pileup parsed and organized via python
- output file used for upstream analyses 
