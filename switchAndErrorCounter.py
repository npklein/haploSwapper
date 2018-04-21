import argparse;
import os
import vcf
import sys
from bisect import bisect
import pandas as pd 

parser = argparse.ArgumentParser()
parser.add_argument("--testVCF", help="VCF file to test against reference", required = True)
parser.add_argument("--vcfReference", help="Reference VCF to check haplotypes with.", required = True)
parser.add_argument("--linking_file", help="File that links DNA IDs of test VCF and refVCF (first column testVCF, second columns refVCF)", required=True)
parser.add_argument("--samples_file", help="File with IDs from testVCF that should be included",required = True)
parser.add_argument("--out_file", help="File to write output to. Will write switching and genotype error out",required = True)
parser.add_argument("--chr", help="Chromosome to calculate over", required = True)
parser.add_argument("--gtf", help="Gtf file with gene coordinate information", required = True)
args = parser.parse_args()


def flush_print(message):
    print(message)
    sys.stdout.flush()

flush_print('testVCF: '+args.testVCF)
flush_print('vcfReference: '+args.vcfReference)
flush_print('chr: '+args.chr)

sample_link = {}
# contains the link between the sample name in test VCF and sample name in ref VCF
with open(args.linking_file) as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        sample_link[line[0]] = line[1]

# contains the samples that need to be included
samples_to_include_testVCF = set([])
samples_to_include_refVCF = set([])
with open(args.samples_file) as input_file:
    for line in input_file:
        sample = line.strip()
        samples_to_include_testVCF.add(sample)
        samples_to_include_refVCF.add(sample_link[sample])

gtf_info_genes = []
genes_start_stop = {}
print('parse gtf')
prev_start = 0
prev_stop = 0
prev_gene = ''
prev_type = ''
prev_feature = ''
with open(args.gtf) as input_file:
    for line in input_file:
        if line.startswith('#') or line.split('\t')[0] != args.chr:
            continue
        line = line.strip().split('\t')
        type  = [line[1]]
        start = int(line[3])
        stop = int(line[4])
        gene = line[8].split('gene_id "')[1].split('"')[0]
        feature = [line[2]]
        if feature[0] == 'gene':
            if start < prev_start:
                raise RuntimeError('should be sorted on start position')
            if start < prev_stop:
                # they should be merged
                gene = prev_gene+';'+gene
                if prev_stop > stop:
                    stop = prev_stop
                start = prev_start
                gtf_info_genes = gtf_info_genes[:-1]
                type.extend(prev_type)
                feature.extend(prev_feature)
            gtf_info_genes.append([start, stop, gene, type, feature])
            prev_start = start
            prev_stop = stop
            prev_gene = gene
            prev_type = type
            prev_feature = feature
            if gene in genes_start_stop:
                raise RuntimeError("Shouldn't have multiple times same gene when filtering feature on gene") 
            genes_start_stop[gene] = [start,stop]

print(len(gtf_info_genes))

def retrieve_haplotype_info(vcf_file, chr, samples_to_use, start=None, stop=None):
    '''Retrieve which haplotype each SNP is on
    
    vcf_file (str): Loctation of VCF file
    chr (str): Chromosome to use
    samples_to_use(list): Samples to use

    
    Returns:
    haplotypeAdict(dict): Haplotype per sample, per sample a list with index 0 = hapA, index 1 = hapB, index 2 = genotype, index 3 = ref-alt
    '''
    vcf_reader = vcf.Reader(filename=vcf_file)
    records = vcf_reader.fetch(chr,start,stop)
    # use ordered dict so that later when taking .values() they are in the same order, can use index for comparing
    hap_info, genotypes,recordInfo = {},{},{}
    x = 0
    ref_and_alt_allele_per_postion = {}
    for record in records:
        if x%100 == 0:
            flush_print(str(x)+' records processed');
        x +=1
        # if filter = PASS: record.FILTER = [], otherwise e.g. record.FILTER = ['INACCESSIBLE']
        # so continue if record.FILTER contains value
        if  record.FILTER:
            continue
        ref = str(record.REF[0])
        alt = str(record.ALT[0])
        if len(ref) > 1 or len(alt) > 1:
            #flush_print('Skipping '+str(record.CHROM)+':'+str(record.POS)+' because it is not a bi-allelic SNP')
            continue
        sample_included_count= 0
        for call in record.samples:
            sample = call.sample
            # check if current sample is one of interest
            if sample not in samples_to_use:
                continue

            if record.POS not in hap_info:
                hap_info[record.POS] = {}
                genotypes[record.POS] = {}
                ref_and_alt_allele_per_postion[record.POS] = {'refAllele':ref, 'altAllele':alt}
            sample_included_count += 1
            if sample_included_count > len(samples_to_use):
                # the rest of the samples in this line are not interesting, so break
                break

            genotype = record.genotype(sample)
            genotypes[record.POS][sample] = genotype
            if record.POS == 211633 and sample == 'gonl-96b':
                print(sample)
                # gonl-96b
            if genotype['GT'] == '0|0':
                hap_info[record.POS][sample] = {'hapA':ref, 'hapB':ref}
            elif genotype['GT'] == '1|1':
                hap_info[record.POS][sample] = {'hapA':alt, 'hapB':alt}
            elif genotype['GT'] == '1|0':
                hap_info[record.POS][sample] = {'hapA':ref, 'hapB':alt}
            elif genotype['GT'] == '0|1':
                hap_info[record.POS][sample] = {'hapA':alt, 'hapB':ref}
            else:
                raise RuntimeError('genotype does not have correct format, probably not phased. Should be 0|0, 1|1, 1|0, or 0|1, was: '+str(genotype['GT']))
    return hap_info, ref_and_alt_allele_per_postion



def make_haplotypes_per_gene(gtf_info_genes, hap_info):
    '''Per gene build the haplotype per sample
    
    gtf_info_genes(list): Gtf feature data
    hap_info(dict): Which haplotype each genotype is on per sample
    
    '''
    feature_df = pd.DataFrame(gtf_info_genes)
    position_df = pd.DataFrame(hap_info, index=[0])
    hits = position_df.apply(lambda col: (feature_df [0] <= col.name) & (col.name < feature_df [1])).values.nonzero()
    haplotype_per_gene_per_sample = {}
    for feature_index, position_index in zip(*hits):
        gene = gtf_info_genes[feature_index][2]
        if gene not in haplotype_per_gene_per_sample:
            haplotype_per_gene_per_sample[gene] = {}
        position = position_df.columns[position_index]
        sample_info = haplotype_testVCF[position]
        for sample in sample_info.keys():
            if sample not in haplotype_per_gene_per_sample[gene]:
                haplotype_per_gene_per_sample[gene][sample] = {}
            haplotype_per_gene_per_sample[gene][sample][position] = sample_info[sample]
    return(haplotype_per_gene_per_sample) 
            
flush_print('starting test VCF')
haplotype_testVCF, ref_and_alt_allele_per_postion_testVCF = retrieve_haplotype_info(args.testVCF, 
                                                                        args.chr, 
                                                                        samples_to_include_testVCF)

haplotype_per_gene_testVCF = make_haplotypes_per_gene(gtf_info_genes, haplotype_testVCF)

switch_and_error_per_gene = {}
for gene in haplotype_per_gene_testVCF:
    if gene not in switch_and_error_per_gene:
        switch_and_error_per_gene[gene] = {}
    flush_print('starting ref VCF for gene '+gene)
    # run per gene chunk so we don't have to loop through the whole reference VCF
    haplotype_refVCF, ref_and_alt_allele_per_postion_refVCF = retrieve_haplotype_info(args.vcfReference,
                                                                               args.chr,
                                                                               samples_to_include_refVCF, 
                                                                               start = genes_start_stop[gene][0],
                                                                               stop = genes_start_stop[gene][1])
    for sample in haplotype_per_gene_testVCF[gene]:
        if sample not in switch_and_error_per_gene[gene]:
            # wrongHom -> AA TT or TT AA
            # hom_shouldBe_het -> AA AT or AA TA or TT AT or TT AT
            # het_shouldBe_hom -> AT AA or TA AA or TA TT or AT TT
            # wrongHet -> AT CG (both het both different)
            switch_and_error_per_gene[gene][sample] = {'overlapping_snps':0,'snps_only_testVCF':0,
                                                       'switch':0,'no_switch':0,'wrongHom':0,
                                                       'het_shouldBe_hom':0,'wrongHet':0,'totalError':0, 
                                                       'overlapping_snp_positions':[],
                                                       'hom_shouldBe_het':0}
        refSample = sample_link[sample]
        for pos in haplotype_per_gene_testVCF[gene][sample]:           
            if pos not in ref_and_alt_allele_per_postion_refVCF:
                switch_and_error_per_gene[gene][sample]['snps_only_testVCF'] += 1
                continue
            if ref_and_alt_allele_per_postion_testVCF[pos] != ref_and_alt_allele_per_postion_refVCF[pos]:
                # ref and alt not same between WGS and RNAseq genotypes, switch haps
                haplotype_refVCF[pos][refSample]['hapA'], haplotype_refVCF[pos][refSample]['hapB'] = haplotype_refVCF[pos][refSample]['hapB'], haplotype_refVCF[pos][refSample]['hapA']
                
            switch_and_error_per_gene[gene][sample]['overlapping_snps'] += 1
            
            hapA_testVCF = haplotype_per_gene_testVCF[gene][sample][pos]['hapA']
            hapB_testVCF = haplotype_per_gene_testVCF[gene][sample][pos]['hapB']

            hapA_refVCF = haplotype_refVCF[pos][refSample]['hapA']
            hapB_refVCF = haplotype_refVCF[pos][refSample]['hapB']
            switch_and_error_per_gene[gene][sample]['overlapping_snp_positions'].append(str(pos))
            
            if hapA_testVCF != hapB_testVCF:
                # testVCF is HET
                if hapA_refVCF != hapB_refVCF:
                    # refVCF is also HET
                    if hapA_testVCF == hapA_refVCF:
                        # hapA of both is the same
                        if hapB_testVCF == hapB_refVCF:
                            # hapB is also the same, so no switch error!
                            switch_and_error_per_gene[gene][sample]['no_switch'] += 1
                        else:
                            # hapA is the same but hapB is different, so different het is called
                            switch_and_error_per_gene[gene][sample]['wrongHet'] += 1
                    elif hapA_testVCF != hapB_refVCF:
                        # hapA of testVCF is different than hapA and hapB of refVCF, wrong het called
                        switch_and_error_per_gene[gene][sample]['wrongHet'] += 1
                    else:
                        assert hapA_testVCF == hapB_refVCF and hapB_testVCF == hapA_refVCF
                        # switch error
                        switch_and_error_per_gene[gene][sample]['switch'] += 1
                else:
                    # refVCF is hom
                    switch_and_error_per_gene[gene][sample]['het_shouldBe_hom'] += 1
                    switch_and_error_per_gene[gene][sample]['totalError'] += 1
            else:
                #testVCF is HOM
                if hapA_refVCF != hapB_refVCF:
                    # refVCF = het
                    switch_and_error_per_gene[gene][sample]['hom_shouldBe_het'] += 1
                    switch_and_error_per_gene[gene][sample]['totalError'] += 1
                else:
                    # refVCF is HOM 
                    if hapA_refVCF == hapB_refVCF:
                        # same homs, but don't count for phasing
                        pass 
                    else:
                        switch_and_error_per_gene[gene][sample]['wrongHom'] += 1
                        switch_and_error_per_gene[gene][sample]['totalError'] += 1
            
            # because hapA and hapB is randomly assigned, when more than half is switched it means that they assigned different haplotypes
            if switch_and_error_per_gene[gene][sample]['switch'] > switch_and_error_per_gene[gene][sample]['no_switch']:
                switch_and_error_per_gene[gene][sample]['switch'], switch_and_error_per_gene[gene][sample]['no_switch'] = switch_and_error_per_gene[gene][sample]['no_switch'], switch_and_error_per_gene[gene][sample]['switch']
            

with open(args.out_file,'w') as out:
    out.write('gene\tsample\toverlapping_snps\tsnps_only_testVCF\t')
    out.write('switch\tno_switch\ttotalError\twrongHom\thet_shouldBe_hom\twrongHet\t')
    out.write('hom_shouldBe_het\toverlapping_snp_positions\n')
    for gene in switch_and_error_per_gene:
        for sample in switch_and_error_per_gene[gene]:
            data = switch_and_error_per_gene[gene][sample]
            out.write(gene+'\t'+sample+'\t'+str(data['overlapping_snps'])+'\t')
            out.write(str(data['snps_only_testVCF'])+'\t'+str(data['switch'])+'\t'+str(data['no_switch'])+'\t')
            out.write(str(data['totalError'])+'\t'+str(data['wrongHom'])+'\t'+str(data['het_shouldBe_hom'])+'\t')
            out.write(str(data['wrongHet'])+'\t'+str(data['hom_shouldBe_het']))
            out.write('\t'+','.join(data['overlapping_snp_positions'])+'\n')
            
            
            
            
            
