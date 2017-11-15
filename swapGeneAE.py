import argparse;
import os
import vcf
import sys
import collections
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--switchErrorDir", help="Directory with the switch errors per chunk (output from switchAndErrorCounter.py)", required = True)
parser.add_argument("--geneAeCountsDir", help="Dir with the gene AE counts from phASER.", required = True)
parser.add_argument("--out_dir", help="Dir to write haplotype counts with switch counts to",required = True)
parser.add_argument("--chunk", help="Chunk currently being processed", required = True)
args = parser.parse_args()

def flush_print(message):
    print(message)
    sys.stdout.flush()

flush_print('switchErrorDir: '+args.switchErrorDir)
flush_print('geneAeCountsDir: '+args.geneAeCountsDir)
flush_print('out_dir: '+args.out_dir)

def header_index(header):
    head_index = {}
    for index, c in enumerate(header):
        head_index[c] = index
    return head_index

#chunk   sample  overlapping_snps        snps_only_chunkVCF      snps_only_refVCF        switch  no switch       switchSnps      noSwitchSnps    totalHetsChunk  totalHetsRef    chunkVCF hom alt TO refVCF hom ref      
# chunkVCF hom ref - refVCF hom alt       chunkVCF hom - refVCF het       chunkVCF het - refVCF hom       chunkVCF het - refVCF het       chunkVCF hom alt TO refVCF hom ref SNPs chunkVCF hom ref - refVCF hom alt SNPs  
# chunkVCF hom - refVCF het SNPs  chunkVCF het - refVCF hom SNPs  chunkVCF het - refVCF het SNPs  totalError      haploSize
switches = {}
flush_print('counting switches')
switchErrorFiles = glob.glob(args.switchErrorDir+'/*switchAndError.txt')
flush_print(str(len(switchErrorFiles))+ ' switch error files in '+args.switchErrorDir+'/*switchAndError.txt')
for switch_error_file in switchErrorFiles:
    if args.chunk not in switch_error_file:
        continue
    with open(switch_error_file) as input_file:
        head_index = header_index(input_file.readline().strip().split('\t'))

        for line in input_file:
            line = line.strip().split('\t')
            sample = line[head_index['sample']].split('.')[1]
            switch = int(line[head_index['switch']])
            no_switch = int(line[head_index['no switch']])
            switch_snp = set(line[head_index['switchSnps']].split(','))
            no_switch_snp = set(line[head_index['noSwitchSnps']].split(','))
            switches[sample] = {'switch':switch, 'no_switch':no_switch,
                                'switch_snp':switch_snp,'no_switch_snp':no_switch_snp}
flush_print('samples: '+str(len(switches.keys())))
# contig                        start                         stop                          name                          aCount                        bCount                        totalCount                    log2_aFC                      n_variants
#['22', '31608224', '31676065', 'ENSG00000182541', '215', '220', '435', '-0.0331668639352', '3', '22_31675185_A_C,22_31675572_C_T,22_31663842_C_G', '1', 'AD1NNNACXX-5-16.mdup.sorted.readGroupsAdded']
switch_geneAE = {}
header = None
flush_print('going through AE files')
geneAE_files = glob.glob(args.geneAeCountsDir+'/*')
flush_print('number of geneAE files: '+str(len(geneAE_files)))
for geneAE_count_file in geneAE_files:
    with open(geneAE_count_file) as input_file:
        header = input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            sample = line[11].split('.')[0]
            if sample not in switches or 'inf' in line:
                continue
            snps_in_gene = set(line[9].split(','))
            switched_snps = switches[sample]['switch_snp']
            not_switched_snps = switches[sample]['no_switch_snp']
            switch_overlap = len(snps_in_gene  & switched_snps)
            no_switch_overlap = len(snps_in_gene  & not_switched_snps)
            if switch_overlap > no_switch_overlap:
                tmp = line[4]
                line[4] = line[5]
                line[5] = tmp
                line[7] = str(float(line[7])*-1)
            outfile = geneAE_count_file.split('/')[-1].replace('.txt','.switched.txt')
            if not outfile in switch_geneAE:
                switch_geneAE[outfile] = []
            switch_geneAE[outfile].append('\t'.join(line)+'\n')
flush_print('done')

for file in switch_geneAE:
    with open(args.out_dir+'/'+file,'w') as out:
       out.write(header)
       for line in switch_geneAE[file]:
           out.write(line)

#    print('written to '+args.out_dir+'/'+haplotype_count_file.replace('.txt','.switched.txt'))
