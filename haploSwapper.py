import argparse;
import os
import vcf
import sys
import collections

parser = argparse.ArgumentParser()
parser.add_argument("--chunkDir", help="Directory with all the chunked VCF files", required = True)
parser.add_argument("--vcfReference", help="Reference VCF to check haplotypes with.", required = True)
parser.add_argument("--chr", help="Chromosome to calculate over", required = True)
parser.add_argument("--debug", help="Print some extra info", required = False, action='store_true')
args = parser.parse_args()

def flush_print(message):
    print(message)
    sys.stdout.flush()

flush_print('chunkDir: '+args.chunkDir)
flush_print('vcfReference: '+args.vcfReference)
flush_print('chr: '+args.chr)


def construct_haplotype(vcf_file, chr, start=None, stop=None):
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
    '''
    vcf_reader = vcf.Reader(filename=vcf_file)
    records = vcf_reader.fetch('22',start,stop)
    sample = None
    haplotypeA = collections.OrderedDict()
    haplotypeB = collections.OrderedDict()
    genotypes = collections.OrderedDict()
    for record in records:
        if not sample:
            call = record.samples[0]
            sample = call.sample
            start = record.POS
        genotype = record.genotype(sample)
        if len(record.ALT[0]) > 1 or len(record.REF[0]) > 1:
            if args.debug:
                flush_print('Skipping '+str(record.CHROM)+':'+str(record.POS)+' because it is not a bi-allelic SNP')
            continue
        ref = str(record.REF[0])
        alt = str(record.ALT[0])
        genotypes[record.POS] = ref+alt
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
        end = record.POS
    return (haplotypeA, haplotypeB, genotypes, start, end)


# chunk files are in args.chunkDir/<CHUNK>/vcf_per_sample
# walk through all <CHUNK> dirs and get the VCF files
flush_print('start looking for directory called vcf_per_sample')
found_dir = False
for root, dirs, files in os.walk(args.chunkDir):
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
            if x % 100 == 0:
                flush_print(str(x)+'/'+str(total_files))
            if args.debug:
                flush_print('starting phaser vcf')
            haplotypeA_chunkVCF, haplotypeB_chunkVCF, genotypes_chunkVCF, start_chunkVCF, end_chunkVCF = construct_haplotype(os.path.join(root, file), args.chr)
            if args.debug:
                flush_print('starting ref VCF')
            # do start_chunkVCF-1 because pyVCF indexing starts at 0 while VCF indexing start at 1
            # and do end_chunkVCF + 1 because it is closed indexed
            haplotypeA_refVCF, haplotypeB_refVCF, genotypes_refVCF, start_refVCF, end_refVCF = construct_haplotype(os.path.join(root, args.vcfReference),
                                                                                                     args.chr, start=start_chunkVCF-1, stop=end_chunkVCF)
            flush_print('chunkVCF: '+str(start_chunkVCF)+'-'+str(end_chunkVCF))
            flush_print('refVCF:   '+str(start_refVCF)+'-'+str(end_refVCF))
            # sanity check, shouldnt happen
            if start_refVCF < start_chunkVCF or end_refVCF > end_chunkVCF:
                raise RuntimeError('start of refVCF lower than start of chunkVCF or end of refVCF higher than end of chunkVCF,, something wrong with SNP indexing')
            keys_a = set(haplotypeA_chunkVCF.keys())
            keys_b = set(haplotypeA_refVCF.keys())
            overlapping_snp_positions = keys_a & keys_b # '&' operator is used for set intersection
            overlapping_haplotypeA_chunkVCF = {snp: haplotypeA_chunkVCF[snp] for snp in overlapping_snp_positions}
            overlapping_haplotypeA_refVCF = {snp: haplotypeA_refVCF[snp] for snp in overlapping_snp_positions}
            overlapping_haplotypeB_chunkVCF = {snp: haplotypeB_chunkVCF[snp] for snp in overlapping_snp_positions}
            overlapping_haplotypeB_refVCF = {snp: haplotypeB_refVCF[snp] for snp in overlapping_snp_positions}
            overlapping_genotypes_chunkVCF = {snp: genotypes_chunkVCF[snp] for snp in overlapping_snp_positions}
            overlapping_genotypes_refVCF = {snp: genotypes_refVCF[snp] for snp in overlapping_snp_positions}
            

            print(''.join(list(overlapping_haplotypeA_chunkVCF.values())))
            print(''.join(list(overlapping_haplotypeB_chunkVCF.values())))
            print(''.join(list(overlapping_haplotypeA_refVCF.values())))
            print(''.join(list(overlapping_haplotypeB_refVCF.values())))
            print(''.join(list(overlapping_genotypes_chunkVCF.values())))
            print(''.join(list(overlapping_genotypes_refVCF.values())))
            
            if haplotypeA_chunkVCF != haplotypeA_refVCF:
                flush_print('switch')
            exit()

if not found_dir:
    flush_print('vcf_per_sample dir not found in input directory tree')
