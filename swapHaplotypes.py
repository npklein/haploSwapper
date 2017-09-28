import argparse;
import os
import vcf
import sys
import collections
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--switchErrorDir", help="Directory with the switch errors per chunk (output from switchAndErrorCounter.py)", required = True)
parser.add_argument("--haplotypeCountsDir", help="Dir with the haplotype counts from phASER.", required = True)
parser.add_argument("--out_dir", help="Dir to write haplotype counts with switch counts to",required = True)
parser.add_argument("--chunk", help="Chunk currently being processed", required = True)
args = parser.parse_args()

def flush_print(message):
    print(message)
    sys.stdout.flush()

flush_print('switchErrorDir: '+args.switchErrorDir)
flush_print('haplotypeCountsDir: '+args.haplotypeCountsDir)
flush_print('out_dir: '+args.out_dir)


#chunk    sample    switch    no switch    chunkVCF hom alt - refVCF hom ref    chunkVCF hom ref - refVCF hom alt    chunkVCF hom - refVCF het    chunkVCF het - refCF hom
# 18535735.18703682    BC48DMACXX-1-9    0    0    0    0    0    1    1
switches = {}
for switch_error_file in glob.glob(args.switchErrorDir+'/*.txt'):
    if args.chunk not in switch_error_file:
        continue
    with open(switch_error_file) as input_file:
        input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            chunk = line[0]
            sample = line[1]
            switch = int(line[2])
            no_switch = int(line[3])
            switch_snp = line[4].split(',')
            no_switch_snp = line[5].split(',')
            switches[sample] = {'switch':switch, 'no_switch':no_switch,
                                'switch_snp':switch_snp,'no_switch_snp':no_switch}

#contig    start    stop    variants    variantCount    variantsBlacklisted    variantCountBlacklisted    haplotypeA    haplotypeB    aCount    bCount    totalCount    blockGWPhase    gwStat    max_haplo_maf    bam    aReads    bReads
#22        170950   17095   22_17095    1               0                                                 G             A             3         7         10            0|1             1         0.0364909        AC1JL5ACXX-2-19_BC52YAACXX-7-19.mdup.sorted.readGroupsAdded
switch_haplotypes = {}
header = None
for haplotype_count_file in glob.glob(args.haplotypeCountsDir+'/*'):
    with open(haplotype_count_file) as input_file:
        header = input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            sample = line[15].split('.')[0]
            if sample not in switches:
                continue
            if switches[sample]['switch'] > switches[sample]['no_switch']:
                flush_print('switched '+line[0]+':'+line[1]+'-'+line[2])
                # the geneAE looks at the 1|0 / 0|1 for constructing a/b haplotype in geneAE, so should only swap genotypes not counts
#                tmp = line[9]
#                line[9] = line[10]
 #               line[10] = tmp
                if line[12] == '0|1':
                    line[12] = '1|0'
                elif line[12] == '1|0':
                    line[12] = '0|1'
                else:
                    raise RuntimeError('Thought only het sites in this file')
            outfile = haplotype_count_file.split('/')[-1].replace('.txt','.switched.txt')
            if outfile not in switch_geneAE:
                switch_geneAE[outfile] = []
            switch_haplotypes[outfile].append('\t'.join(line)+'\n')


for file in switch_geneAE:
    with open(args.out_dir+'/'+file,'w') as out:
       out.write(header)
       for line in switch_geneAE[file]:
           out.write(line)

#    print('written to '+args.out_dir+'/'+haplotype_count_file.replace('.txt','.switched.txt'))
