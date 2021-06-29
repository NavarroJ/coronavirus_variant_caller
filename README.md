# Variant Caller for Coronavirus strains
in-house variant caller for the Sung Lab, output used for upstream analyses. basic pipeline:

- run alignSeq.sh, use output as the input for gen_tables.py
- input fasta file with coronavirus sequences is aligned to reference using bwa
- processed with samtools, converted to pileup format
- pileup parsed and organized via python
- output files used for upstream analyses 
