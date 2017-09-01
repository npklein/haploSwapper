import argparse;
import os
import vcf
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--chunkDir", help="Directory with all the chunked VCF files", required = True)
parser.add_argument("--vcfReference", help="Reference VCF to check haplotypes with.", required = True)
parser.add_argument("--linking_file", help="File that links DNA IDs of chunk VCF and refVCF (first column chunkVCF, second columns refVCF)", required=True)
parser.add_argument("--samples_file", help="File with IDs from chunkVCF that should be included",required = True)
parser.add_argument("--out_file", help="FIle to write output to. Will write switching and genotype error out",required = True)
parser.add_argument("--chr", help="Chromosome to calculate over", required = True)
parser.add_argument("--debug", help="Print some extra info", required = False, action='store_true')
args = parser.parse_args()

def flush_print(message):
    print(message)
    sys.stdout.flush()

flush_print('NOTE!!! This gets sample name from VCF file names, e.g. only works if you have the SAME naming scheme')
flush_print('<anything except phASER.>phASER.<sampleName>.chr<anything except .chr>.vcf.gz')
flush_print('chunkDir: '+args.chunkDir)
flush_print('vcfReference: '+args.vcfReference)
flush_print('chr: '+args.chr)

sample_link = {}
# contains the link between the sample name in chunk VCF and sample name in ref VCF
with open(args.linking_file) as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        sample_link[line[0]] = line[1]

# contains the samples that need to be included
with open(args.samples_file) as input_file:
    samples_to_include = set(input_file.read().split('\n'))

def construct_haplotype(vcf_file, chr, start=None, stop=None, sample_to_use = None, silent = True):
    '''Construct haplotype from vcf_file chr:start-stop
    
    vcf_file (str): Loctation of VCF file
    chr (str): Chromosome to use
    start (int): Start of block (def: None)
    stop (int): End of block (def: None)
    
    Returns:
    haplotypeA(list): Haplotype of A allele
    haplotypeB(list): Haplotype of B allele
    start(int): Start position of input VCF
    end(int): End position of input VCF
    sample_to_use(str): If not None, use this sample. Otherwise, use first sample (chunkVCF from phASER should only have one)
    '''
    vcf_reader = vcf.Reader(filename=vcf_file)
    records = vcf_reader.fetch('22',start,stop)
    # use ordered dict so that later when taking .values() they are in the same order, can use index for xomparing
    haplotypeA = {}
    haplotypeB = {}
    genotypes = {}
    samples = set([])
    start = None
    recordInfo = {}
    for record in records:
        # saving start and end so that for the ref VCF tabix indexing can be used to grab the same chunk as for the chunk VCF
        if not start:
            start = record.POS
        end = record.POS
        ref = str(record.REF[0])
        alt = str(record.ALT[0])
        if len(ref) > 1 or len(alt) > 1:
           if args.debug:
              flush_print('Skipping '+str(record.CHROM)+':'+str(record.POS)+' because it is not a bi-allelic SNP')
              continue
        recordInfo[record.POS] = ref+'_'+alt
        for call in record.samples:
            sample = call.sample
            samples.add(sample)
            # if sample_to_use is set it means that we are looking for one particular sample
            if sample_to_use:
                if sample != sample_to_use:
                    continue
            elif len(samples) > 1:
                # when sample_to_use is not set we don't know which sample to take, so for safety if there is more than one sample in the file raise error
                raise RuntimeError('Sample to use was not given but more than one sample in VCF, goes wrong')
            genotype = record.genotype(sample)
            genotypes[record.POS] = genotype
            if genotype['GT'] == '0|0':
                haplotypeA[record.POS] = ref
                haplotypeB[record.POS] = ref
            elif genotype['GT'] == '1|1':
                haplotypeA[record.POS] = alt
                haplotypeB[record.POS] = alt
            elif genotype['GT'] == '1|0':
                haplotypeA[record.POS] = ref
                haplotypeB[record.POS] = alt
            elif genotype['GT'] == '0|1':
                haplotypeA[record.POS] = alt
                haplotypeB[record.POS] = ref
            else:
                raise RuntimeError('genotype does not have correct format, probably not phased. Should be 0|0, 1|1, 1|0, or 0|1, was: '+str(genotype['GT']))
            # break cause we already got the right sample on this line
            break
    return (haplotypeA, haplotypeB, genotypes, start, end, sample, recordInfo)


# chunk files are in args.chunkDir/<CHUNK>/vcf_per_sample
flush_print('start looking for directory called vcf_per_sample')
found_dir = False
out =  open(args.out_file,'w')
out.write('chunk\tsample\tswitch\tno switch\tswitchSnps\tnoSwitchSnps\ttotalHets\tchunkVCF hom alt - refVCF hom ref\tchunkVCF hom ref - refVCF hom alt\tchunkVCF hom - refVCF het\tchunkVCF het - refCF hom\ttotalError\thaploSize\n')

# walk through all <CHUNK> dirs and get the VCF files
for root, dirs, files in os.walk(args.chunkDir):
    # add / because also have vcf_per_sample_<some other name>, this way can split on /
    root = root+'/'
    if not '/vcf_per_sample/' in root:
        continue
    found_dir = True
    flush_print('Starting processing '+root)
    total_files = len(files)
    x = 0 
    for file in files:
        x += 1
        if file.endswith(".vcf.gz"):
            chunk = '.'.join(file.split('.chr')[1].split('.')[1:3])
            switch_count = 0
            not_switch_count = 0
            genotype_error_count = {'total':0,
                                    'chunkVCF hom alt - refVCF_hom_ref':0,  # AA TT
                                    'chunkVCF hom ref - refVCF hom alt':0,  # TT AA
                                    'chunkVCF hom - refVCF het':0,          # AA AT or AA TA or TT AT or TT AT
                                    'chunkVCF het - refCF hom':0}           # AT AA or TA AA or TA TT or AT TT
            total_count = 0
            
            # this is specific for the filename
            file_sample_name = file.split('phASER.')[1].split('.chr')[0]
            if not file_sample_name in samples_to_include:
                continue
            if x % 10 == 0:
                flush_print(str(x)+'/'+str(total_files))
            if args.debug:
                flush_print('starting phaser vcf')
            haplotypeA_chunkVCF, haplotypeB_chunkVCF, genotypes_chunkVCF, start_chunkVCF, end_chunkVCF, sample, recordInfoChunk  = construct_haplotype(os.path.join(root, file), args.chr, silent=True)
            if sample != file_sample_name:
                raise RuntimeError('Returned results from wrong sample, should have been '+file_sample_name+', was '+sample)
            if args.debug:
                flush_print('starting ref VCF')
            # do start_chunkVCF-1 because pyVCF indexing starts at 0 while VCF indexing start at 1
            # and do end_chunkVCF + 1 because it is closed indexed
            haplotypeA_refVCF, haplotypeB_refVCF, genotypes_refVCF, start_refVCF, end_refVCF, sample, recordInfoRef = construct_haplotype(os.path.join(root, args.vcfReference),
                                                                                                                            args.chr, start=start_chunkVCF-1, stop=end_chunkVCF,
                                                                                                                            sample_to_use = sample_link[sample], silent = False)
            if args.debug:
                flush_print('chunkVCF: '+str(start_chunkVCF)+'-'+str(end_chunkVCF))
                flush_print('refVCF:   '+str(start_refVCF)+'-'+str(end_refVCF))
            # sanity check, shouldnt happen
            if start_refVCF < start_chunkVCF or end_refVCF > end_chunkVCF:
                raise RuntimeError('start of refVCF lower than start of chunkVCF or end of refVCF higher than end of chunkVCF,, something wrong with SNP indexing')
            keys_a = set(haplotypeA_chunkVCF.keys())
            keys_b = set(haplotypeA_refVCF.keys())
            # only get those SNPs that have been genotypes in both the chunk VCF and the ref VCF
            overlapping_snp_positions = keys_a & keys_b # '&' operator is used for set intersection
            switch_snps = []
            no_switch_snps = []
            totalHets = 0
#            flush_print(str(len(overlapping_snp_positions))+' overlapping SNP positions')
            for snp in overlapping_snp_positions:
                genotypeCall_chunkVCF = genotypes_chunkVCF[snp]
                genotypeCall_refVCF = genotypes_refVCF[snp]
                genotype_refVCF = haplotypeA_refVCF[snp]+haplotypeB_refVCF[snp]
                genotype_chunkVCF = haplotypeA_chunkVCF[snp]+haplotypeB_chunkVCF[snp]
                if genotype_chunkVCF != genotype_refVCF and genotype_chunkVCF[1]+genotype_chunkVCF[0] != genotype_refVCF:
                    if genotypeCall_chunkVCF['GT'] == '1|1' and genotypeCall_refVCF['GT'] == '0|0':
                        genotype_error_count['chunkVCF hom alt - refVCF_hom_ref'] += 1
                    elif genotypeCall_chunkVCF['GT'] == '0|0' and genotypeCall_refVCF['GT'] == '1|1':
                        genotype_error_count['chunkVCF hom ref - refVCF hom alt'] += 1
                    elif (genotypeCall_chunkVCF['GT'] == '1|1' or genotypeCall_chunkVCF['GT'] == '0|0'
                         ) and (genotypeCall_refVCF['GT'] == '0|1' or genotypeCall_refVCF['GT'] == '1|0'):
                        genotype_error_count['chunkVCF hom - refVCF het'] += 1
                    elif (genotypeCall_refVCF['GT'] == '1|1' or genotypeCall_refVCF['GT'] == '0|0'
                         ) and (genotypeCall_chunkVCF['GT'] == '0|1' or genotypeCall_chunkVCF['GT'] == '1|0'):
                        genotype_error_count['chunkVCF het - refCF hom'] += 1
                        totalHets += 1
                    else:
                        message = 'Combinations of genotypes that I did not account for:\n'
                        message += 'chunk: '+genotypeCall_chunkVCF['GT']
                        message += 'ref: '+genotypeCall_refVCF['GT']
                        raise RuntimeError(message)
                    genotype_error_count['total'] += 1
                elif genotypeCall_chunkVCF['GT'] == '1|0' or genotypeCall_chunkVCF['GT'] == '0|1':
                    totalHets += 1
                    snp_name = args.chr+'_'+str(snp)+'_'+str(recordInfoChunk[snp])
                    if genotype_chunkVCF != genotype_refVCF:
                        switch_snps.append(snp_name)
                        switch_count += 1
                    else:
                        no_switch_snps.append(snp_name)
                        not_switch_count += 1
                total_count += 1
            out.write(chunk+'\t'+file_sample_name+'\t'+str(switch_count)+'\t'+str(not_switch_count)+'\t'+','.join(switch_snps)+'\t'+','.join(no_switch_snps)+'\t'+
                      str(totalHets)+'\t'+
                      str(genotype_error_count['chunkVCF hom alt - refVCF_hom_ref'])+'\t'+str(genotype_error_count['chunkVCF hom ref - refVCF hom alt'])+'\t'+
                      str( genotype_error_count['chunkVCF hom - refVCF het'])+'\t'+str(genotype_error_count['chunkVCF het - refCF hom'])+'\t'+
                      str(genotype_error_count['total'])+'\t'+str(total_count)+'\n')


out.close()
print('written to '+args.out_file)

if not found_dir:
    flush_print('vcf_per_sample dir not found in input directory tree')
