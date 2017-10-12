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

def construct_haplotype(vcf_file, chr, start=None, stop=None, sample_to_use = None, silent = True, filter_pass_only = True):
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
    filter_pass_only(bool): If true, filter out SNPs that do not have PASS in FILTER columns (def: True)
    '''
    vcf_reader = vcf.Reader(filename=vcf_file)
    records = vcf_reader.fetch(chr,start,stop)
    # use ordered dict so that later when taking .values() they are in the same order, can use index for xomparing
    haplotypeA = {}
    haplotypeB = {}
    genotypes = {}
    samples = set([])
    start = None
    recordInfo = {}
    end = None
    sample = None
    x = 0
    for record in records:
        # saving start and end so that for the ref VCF tabix indexing can be used to grab the same chunk as for the chunk VCF
        if not start:
            start = record.POS
        # if filter = PASS: record.FILTER = [], otherwise e.g. record.FILTER = ['INACCESSIBLE']
        # so continue if record.FILTER contains value
        if filter_pass_only and record.FILTER:
            continue
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
out.write('chunk'+'\t')
out.write('sample'+'\t')
out.write('overlapping_snps'+'\t')
out.write('snps_only_chunkVCF'+'\t')
out.write('snps_only_refVCF'+'\t')
out.write('switch'+'\t') 
out.write('no switch'+'\t') 
out.write('switchSnps'+'\t')
out.write('noSwitchSnps'+'\t')
out.write('totalHetsChunk'+'\t')
out.write('totalHetsRef'+'\t')
out.write('chunkVCF hom alt TO refVCF hom ref'+'\t')
out.write('chunkVCF hom ref - refVCF hom alt'+'\t')
out.write('chunkVCF hom - refVCF het'+'\t')
out.write('chunkVCF het - refVCF hom'+'\t')
out.write('chunkVCF het - refVCF het'+'\t')
out.write('chunkVCF hom alt TO refVCF hom ref SNPs'+'\t')
out.write('chunkVCF hom ref - refVCF hom alt SNPs'+'\t')
out.write('chunkVCF hom - refVCF het SNPs'+'\t')
out.write('chunkVCF het - refVCF hom SNPs'+'\t')
out.write('chunkVCF het - refVCF het SNPs'+'\t')
out.write('totalError'+'\t')
out.write('haploSize\n')

sample_snp_errors = {}
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
                                    'chunkVCF hom alt - refVCF hom ref':0,  # AA TT
                                    'chunkVCF hom ref - refVCF hom alt':0,  # TT AA
                                    'chunkVCF hom - refVCF het':0,          # AA AT or AA TA or TT AT or TT AT
                                    'chunkVCF het - refVCF hom':0,          # AT AA or TA AA or TA TT or AT TT
                                    'chunkVCF het - refVCF het':0}           # AT CG (both het both different)
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
            if not end_chunkVCF:
                flush_print('no records in chunk')
                continue
            if sample != file_sample_name:
                raise RuntimeError('Returned results from wrong sample, should have been '+file_sample_name+', was '+sample)
            if args.debug:
                flush_print('starting ref VCF')
            # do start_chunkVCF-1 because pyVCF indexing starts at 0 while VCF indexing start at 1
            # and do end_chunkVCF + 1 because it is closed indexed
            haplotypeA_refVCF, haplotypeB_refVCF, genotypes_refVCF, start_refVCF, end_refVCF, sample, recordInfoRef = construct_haplotype(os.path.join(root, args.vcfReference),
                                                                                                                            args.chr, start=start_chunkVCF-1, stop=end_chunkVCF,
                                                                                                                            sample_to_use = sample_link[sample], silent = False)
            sample_snp_errors[sample] = {'switch_snp':[],'no_switch_snp':[],
                                         'chunkVCF_homAlt__refVCF_homRef':[],
                                         'chunkVCF_homRef__refVCF_homAlt':[],
                                         'chunkVCF_hom__refVCF_het':[],
                                         'chunkVCF_het__refVCF_hom':[],
                                         'chunkVCF_het__refVCF_het':[]}
            if not end_refVCF:
                flush_print('no records in '+args.chr+':'+str(start_chunkVCF-1)+'-'+str(end_chunkVCF))
                continue
            if args.debug:
                flush_print('chunkVCF: '+str(start_chunkVCF)+'-'+str(end_chunkVCF))
                flush_print('refVCF:   '+str(start_refVCF)+'-'+str(end_refVCF))
            # sanity check, shouldnt happen
            if start_refVCF < start_chunkVCF or end_refVCF > end_chunkVCF:
                raise RuntimeError('start of refVCF lower than start of chunkVCF or end of refVCF higher than end of chunkVCF,, something wrong with SNP indexing')
            snps_a = set(haplotypeA_chunkVCF.keys())
            snps_b = set(haplotypeA_refVCF.keys())
            # only get those SNPs that have been genotypes in both the chunk VCF and the ref VCF
            overlapping_snp_positions = snps_a & snps_b # '&' operator is used for set intersection
            snps_only_chunkVCF = snps_a - snps_b
            snps_only_refVCF = snps_b - snps_a

            switch_snps = []
            no_switch_snps = []
            chunkVCF_homAlt__refVCF_homRef = []
            chunkVCF_homRef__refVCF_homAlt = []
            chunkVCF_hom__refVCF_het = []
            chunkVCF_het__refVCF_hom = []
            chunkVCF_het__refVCF_het = []
            totalHetsChunk = 0
            totalHetsRef = 0
            
            snps_only_chunkVCF_names = []
            snps_only_refVCF_names = []
            # loop through SNPs that are not called in refVCF and check if they are het in chunkVCF. If so, add them for the ouput file
            for snp in snps_only_chunkVCF:
                genotypeCall_chunkVCF = genotypes_chunkVCF[snp]
                if genotypeCall_chunkVCF['GT'] == '0|1' or genotypeCall_chunkVCF['GT'] == '1|0':
                    snps_only_chunkVCF_names.append(args.chr+'_'+str(snp)+'_'+str(recordInfoChunk[snp]))
            for snp in snps_only_refVCF:
                genotypeCall_refVCF = genotypes_refVCF[snp]
                if genotypeCall_refVCF['GT'] == '0|1' or genotypeCall_refVCF['GT'] == '1|0':
                    snps_only_refVCF_names.append(args.chr+'_'+str(snp)+'_'+str(recordInfoRef[snp]))
                    
            # the snp keys are the numeric positions, need to change to snp name for later comparison
            overlapping_snp_positions_names = []

#            flush_print(str(len(overlapping_snp_positions))+' overlapping SNP positions')
            for snp in overlapping_snp_positions:
                snp_name = args.chr+'_'+str(snp)+'_'+str(recordInfoChunk[snp])
                genotypeCall_chunkVCF = genotypes_chunkVCF[snp]
                genotypeCall_refVCF = genotypes_refVCF[snp]
                genotype_refVCF = haplotypeA_refVCF[snp]+haplotypeB_refVCF[snp]
                if len(genotype_refVCF) > 2:
                    flush_print('Deletion or insertion? refVCF = '+genotype_refVCF)
                    continue
                genotype_chunkVCF = haplotypeA_chunkVCF[snp]+haplotypeB_chunkVCF[snp]
                # First check if genotype of chunk VCF is not same as ref VCF (or switched, so not (1|1 and 1|1), (0|0 and 0|0),
                # (0|1 and 0|1) or (1|0 and 0|1)
                if genotype_chunkVCF != genotype_refVCF and genotype_chunkVCF[1]+genotype_chunkVCF[0] != genotype_refVCF:
                    # genotype error hom ref -> hom alt
                    if genotypeCall_chunkVCF['GT'] == '1|1' and genotypeCall_refVCF['GT'] == '0|0':
                        genotype_error_count['chunkVCF hom alt - refVCF hom ref'] += 1
                        chunkVCF_homAlt__refVCF_homRef.append(snp_name)
                        sample_snp_errors[sample]['chunkVCF_homAlt__refVCF_homRef'].append(snp_name)
                    # genotype error hom alt -> hom ref
                    elif genotypeCall_chunkVCF['GT'] == '0|0' and genotypeCall_refVCF['GT'] == '1|1':
                        genotype_error_count['chunkVCF hom ref - refVCF hom alt'] += 1
                        chunkVCF_homRef__refVCF_homAlt.append(snp_name)
                        sample_snp_errors[sample]['chunkVCF_homRef__refVCF_homAlt'].append(snp_name)
                    # genotype error het -> hom
                    elif (genotypeCall_chunkVCF['GT'] == '1|1' or genotypeCall_chunkVCF['GT'] == '0|0'
                         ) and (genotypeCall_refVCF['GT'] == '0|1' or genotypeCall_refVCF['GT'] == '1|0'):
                        genotype_error_count['chunkVCF hom - refVCF het'] += 1
                        chunkVCF_hom__refVCF_het.append(snp_name)
                        totalHetsRef += 1
                        sample_snp_errors[sample]['chunkVCF_hom__refVCF_het'].append(snp_name)
                    # genotype error hom -> het
                    elif (genotypeCall_refVCF['GT'] == '1|1' or genotypeCall_refVCF['GT'] == '0|0'
                         ) and (genotypeCall_chunkVCF['GT'] == '0|1' or genotypeCall_chunkVCF['GT'] == '1|0'):
                        genotype_error_count['chunkVCF het - refVCF hom'] += 1
                        chunkVCF_het__refVCF_hom.append(snp_name)
                        totalHetsChunk += 1
                        sample_snp_errors[sample]['chunkVCF_het__refVCF_hom'].append(snp_name)
                    # genotype het complete different (e.g. chunkVCF AT, refVCF GC)
                    elif (genotypeCall_refVCF['GT'] == '1|0' or genotypeCall_refVCF['GT'] == '0|1'
                            ) and (genotypeCall_chunkVCF['GT'] == '1|0' or genotypeCall_chunkVCF['GT'] == '0|1'):
                        genotype_error_count['chunkVCF_het__refVCF_het'] += 1
                        chunkVCF_het__refVCF_het.append(snp_name)
                        sample_snp_errors[sample]['chunkVCF_het__refVCF_het'].append(snp_name)
                        totalHetsChunk += 1
                        overlapping_snp_positions_names.append(args.chr+'_'+str(snp)+'_'+str(recordInfoRef[snp]))
                    else:
                        message = 'Combinations of genotypes that I did not account for:\n'
                        message += genotype_chunkVCF + '!=' +genotype_refVCF+' and '+genotype_chunkVCF[1]+genotype_chunkVCF[0] + '!=' + genotype_refVCF +'\n'
                        message += 'chunk: '+genotypeCall_chunkVCF['GT']+'\n'
                        message += 'ref: '+genotypeCall_refVCF['GT']
                        raise RuntimeError(message)
                    genotype_error_count['total'] += 1
                # know they are the same, now check if they are hets
                elif genotypeCall_chunkVCF['GT'] == '1|0' or genotypeCall_chunkVCF['GT'] == '0|1':
                    totalHetsChunk += 1
                    totalHetsRef += 1
                    # if they are not the same but are both het, they must be swapped
                    if genotype_chunkVCF != genotype_refVCF:
                        switch_snps.append(snp_name)
                        switch_count += 1
                        sample_snp_errors[sample]['switch_snp'].append(snp_name)
                        overlapping_snp_positions_names.append(args.chr+'_'+str(snp)+'_'+str(recordInfoRef[snp]))
                    else:
                        no_switch_snps.append(snp_name)
                        not_switch_count += 1
                        sample_snp_errors[sample]['no_switch_snp'].append(snp_name)
                        overlapping_snp_positions_names.append(args.chr+'_'+str(snp)+'_'+str(recordInfoRef[snp]))

                # if it is not one of the first then they are both hom ref or both hom alt
                total_count += 1
            out.write(chunk+'\t'+'chr'+args.chr+'.'+file_sample_name+'\t'+','.join(overlapping_snp_positions_names)+'\t'+','.join(snps_only_chunkVCF_names)+
                      '\t'+','.join(snps_only_refVCF_names)+'\t'+
                        str(switch_count)+'\t'+str(not_switch_count)+'\t'+
                      ','.join(switch_snps)+'\t'+','.join(no_switch_snps)+'\t'+str(totalHetsChunk)+'\t'+str(totalHetsRef)+'\t'+
                      str(genotype_error_count['chunkVCF hom alt - refVCF hom ref'])+'\t'+str(genotype_error_count['chunkVCF hom ref - refVCF hom alt'])+'\t'+
                      str( genotype_error_count['chunkVCF hom - refVCF het'])+'\t'+str(genotype_error_count['chunkVCF het - refVCF hom'])+'\t'+
                      str(genotype_error_count['chunkVCF het - refVCF het'])+'\t'+
                      ','.join(chunkVCF_homAlt__refVCF_homRef)+'\t'+','.join(chunkVCF_homRef__refVCF_homAlt)+'\t'+','.join(chunkVCF_hom__refVCF_het)+
                      '\t'+','.join(chunkVCF_het__refVCF_hom)+'\t'+','.join(chunkVCF_het__refVCF_het)+
                      '\t'+str(genotype_error_count['total'])+'\t'+str(total_count)+'\n')


out.close()
print('written to '+args.out_file)

snp_per_sample_file = args.out_file.replace('.txt','snpErrorPerSample.txt')
with open(snp_per_sample_file,'w') as out:
    out.write('sample_name\tswitch_snps\tno_switch_snps\t')
    out.write('chunkVCF homAlt to refVCF homRef\t')
    out.write('chunkVCF homRef to refVCF homAlt\t')
    out.write('chunkVCF hom to refVCF het\t')
    out.write('chunkVCF het to refVCF hom\t')
    out.write('chunkVCF het to refVCF het')
    for sample in sample_snp_errors:
        out.write('\n'+sample+'\t')
        d = sample_snp_errors[sample]
        out.write(','.join(d['switch_snp'])+'\t')
        out.write(','.join(d['no_switch_snp'])+'\t')
        out.write(','.join(d['chunkVCF_homAlt__refVCF_homRef'])+'\t')
        out.write(','.join(d['chunkVCF_homRef__refVCF_homAlt'])+'\t')
        out.write(','.join(d['chunkVCF_hom__refVCF_het'])+'\t')
        out.write(','.join(d['chunkVCF_het__refVCF_hom'])+'\t')
        out.write(','.join(d['chunkVCF_het__refVCF_het']))
print('written to '+snp_per_sample_file)


if not found_dir:
    raise RuntimeError('vcf_per_sample dir not found in input directory tree')
