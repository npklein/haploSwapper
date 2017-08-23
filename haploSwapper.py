import argparse;
import os
import vcf


parser = argparse.ArgumentParser()
parser.add_argument("--chunkDir", help="Directory with all the chunked VCF files", required = True)
parser.add_argument("--vcfReference", help="Reference VCF to check haplotypes with.", required = True)
parser.add_argument("--chr", help="Chromosome to calculate over", required = True)
parser.add_argument("--debug", help="Print some extra info", required = False, action='store_true')
args = parser.parse_args()


def construct_haplotype(vcf_file, chr, start=None, stop=None, include_homs=False):
    '''Construct haplotype from vcf_file chr:start-stop
    
    vcf_file (str): Loctation of VCF file
    chr (str): Chromosome to use
    start (int): Start of block (def: None)
    stop (int): End of block (def: None)
    include_homs (bool): If true, include homozygotes in haplotype (def: False)
    
    Returns:
    haplotypeA(list): Haplotype of A allele
    haplotypeB(list): Haplotype of B allele
    start(int): Start position of input VCF
    end(int): End position of input VCF
    '''
    vcf_reference = vcf.Reader(filename=vcf_file)
    records = vcf_reference.fetch('22',start,stop)
    sample = None
    haplotypeA = []
    haplotypeB = []
    for record in records:
        if not sample:
            call = record.samples[0]
            sample = call.sample
            start = record.POS
        genotype = record.genotype(sample)
        if len(record.ALT[0]) > 1:
            if args.debug:
                print('Skipping '+str(record.CHROM)+':'+str(record.POS)+' because it is not a bi-allelic SNP')
            continue
        if include_homs:
            if genotype['GT'] == '0|0':
                haplotypeA.append(str(record.REF))
                haplotypeB.append(str(record.REF))
            elif genotype['GT'] == '1|1':
                haplotypeA.append(str(record.ALT[0]))
                haplotypeB.append(str(record.ALT[0]))
        elif genotype['GT'] == '1|0':
            haplotypeA.append(str(record.REF))
            haplotypeB.append(str(record.ALT[0]))
        elif genotype['GT'] == '0|1':
            haplotypeA.append(str(record.ALT[0]))
            haplotypeB.append(str(record.REF))
        end = record.POS
    return (haplotypeA, haplotypeB, start, end)


# chunk files are in args.chunkDir/<CHUNK>/vcf_per_sample
# walk through all <CHUNK> dirs and get the VCF files
for root, dirs, files in os.walk(args.chunkDir):
    if not 'vcf_per_sample' in root:
        continue
    print('Starting processing '+root)
    total_files = len(files)
    x = 0 
    for file in files:
        if file.endswith(".vcf.gz"):
            x += 1
            if x % 100 == 0:
                print(str(x)+'/'+str(total_files))
            haplotypeA_chunkVCF, haplotypeB_chunkVCF, start_chunkVCF, end_chunkVCF = construct_haplotype((os.path.join(root, file)), args.chr)
            # do start_chunkVCF-1 because pyVCF indexing starts at 0 while VCF indexing start at 1
            # and do end_chunkVCF + 1 because it is closed indexed
            haplotypeA_refVCF, haplotypeB_refVCF, start_refVCF, end_refVCF = construct_haplotype((os.path.join(root, file)), args.chr, start_chunkVCF-1, end_chunkVCF+1)
            if start_chunkVCF != start_refVCF or end_chunkVCF != end_refVCF:
                message = ('VCF region not the same:\nchunkVCF: '+str(start_chunkVCF)+'-'+str(end_chunkVCF)+
                                   '\nrefVCF:   '+str(start_refVCF)+'-'+str(end_refVCF))
                raise RuntimeError(message)
            if haplotypeA_chunkVCF != haplotypeA_refVCF:
                print('switch')
            


