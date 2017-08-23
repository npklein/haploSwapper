haploSwapper takes the haplotypes from chunked VCF files and compares them to the haplotypes in a reference VCF file. 
If the haplotype is the same but the 1|0 / 0|1 swapped between the chunk VCF and reference VCF, the direction of effect
of the phASER files gets swapped.



***EXAMPLE***

    python haploSwapper.py --chunkDir testData/rnaseqGenotypes/chr22/ --vcfReference testData/WGSgenotypes/vcf_per_sample/BIOS_LLDeep_Diagnostics_phASER.AD1NW8ACXX-3-2.chr22.vcf.gz
