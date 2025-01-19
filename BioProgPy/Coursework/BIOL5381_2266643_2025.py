import argparse
import logging
import csv
import bisect
import traceback
import sys
import math

from os.path import isfile

import vcf
import gffutils as gff
from Bio.Seq import Seq
from pyfaidx import Fasta

COMP_DICT = {'A': 'T',
             'T': 'A',
             'G': 'C',
             'C': 'G'}

variant_results_headers = ['CHROM','POS','REF','ALT','Type','Transcript','Protein Location','Ref AA','Alt AA']

# Loops through variants in vcf file
# Returns a dict of all variants which meet the minimum quatity threshold for each chromosome, and a sorted list of their positions which are used as keys
# Multiple variations can be stored with the same position
def get_variants(vcf_filename, logger, min_quality = 20):
    try:
        print(vcf_filename)
        if not isfile(vcf_filename):
            raise Exception('file does not exist')
        
        if not (vcf_filename.endswith('.vcf') or vcf_filename.endswith('.vcf.gz')):
            raise Exception('Incorrect file format, should be .vcf file')
        
        vcf_reader = vcf.Reader(filename=vcf_filename)

        q_count = 0
        non_q_count = 0
        variants_chromosome_dict = {}

        for record in vcf_reader:
            if record.QUAL > min_quality:
                q_count += 1
                if record.CHROM in variants_chromosome_dict:
                    if int(record.POS) in variants_chromosome_dict[record.CHROM]['variants']:
                        print('eh')
                        variants_chromosome_dict[record.CHROM]['variants'][int(record.POS)]['alts'].append(record.ALT)
                    else:
                        variants_chromosome_dict[record.CHROM]['variants'][int(record.POS)] = {'ref':record.REF, 'alts':[record.ALT]}
                        bisect.insort_right(variants_chromosome_dict[record.CHROM]['sorted_pos_keys'], record.POS)
                else:
                    variants_chromosome_dict[record.CHROM] = {}
                    variants_chromosome_dict[record.CHROM]['sorted_pos_keys'] = [int(record.POS)]
                    variants_chromosome_dict[record.CHROM]['found_pos'] = []
                    variants_chromosome_dict[record.CHROM]['variants'] = {}
                    variants_chromosome_dict[record.CHROM]['variants'][int(record.POS)] = {'ref':record.REF, 'alts':[record.ALT]}

                # if q_count < 10:
                #     print(record.ID)
            else:
                non_q_count += 1
        
        # for chromosome in variants_chromosome_dict:
        #     # Later variants will be removed from this list as they are found in genes, leaving only the non-coding genes
        #     variants_chromosome_dict[chromosome]['non-coding_pos'] = variants_chromosome_dict[chromosome]['sorted_pos_keys']

        logger.info(f'{non_q_count} variants with quality score below {min_quality}\n')
        return variants_chromosome_dict
    
    except Exception as err:
        logger.error(f'Error loading variants file: {err}')
        return None

# Create feature db, or connect to it if it already exists
def get_feature_db(gff_filename):
    try:
        if not isfile(gff_filename):
            raise Exception('file does not exist')
        
        if not gff_filename.endswith('.gff'):
            raise Exception('Incorrect file format, should be .gff file')
        
        db_filepath = f'{gff_filename.split('.gff')[0]}.db'

        #print(db_filepath)
        if not isfile(db_filepath):
            db = gff.create_db(gff_filename, dbfn=db_filepath, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        else:
            db = gff.FeatureDB(db_filepath, keep_order=True)
        
        return db
    except Exception as err:
        logger.error(f'Error loading gene file: {err}')
        return None

def find_variants_for_feature(start, end, pos_list):
    assert start <= end

    lower_bound_variants = bisect.bisect_left(pos_list, start)
    upper_bound_variants = bisect.bisect_right(pos_list, end, lo=lower_bound_variants)

    return lower_bound_variants, upper_bound_variants


# Returns the codon, and alternative codon of a sequence with SNP
# Input sequence should already be adjusted to start at the first codon and continue until at least the last codon the snp could be located at
def get_alt_codon(seq, pos, alt):
    # Get the sequence up until the SNP
    alt_codon = seq[:pos]

    # Get rid of everything before the codon
    seq_mod_3 = len(alt_codon) % 3
    if seq_mod_3 > 0:
        alt_codon = alt_codon[-seq_mod_3:]  
    else:
        alt_codon = ''
    
    # Add the base at the site of the SNP, we need 2 codons for the reference and the alternative
    codon = alt_codon
    codon += seq[pos]
    alt_codon += alt

    # Add the rest of the codon if missing
    seq_len = len(alt_codon)
    rest_of_codon_seq = seq[pos+1:pos+4-seq_len]
    codon += rest_of_codon_seq
    alt_codon += rest_of_codon_seq

    return codon, alt_codon

# Returns a sequence generated from a coding region, with a buffer either side to line up with codons
# Eg. If the region starts mid codon include the bases to make up that codon from the previous coding region
def get_adjusted_cds(mRNA, seq, cds):
    #print('a')

    if mRNA == '+':
        mRNA_order_by = 'start'
    else:
        mRNA_order_by = 'end'
    
    seq_segments = ['','']
    pre_seq = True
    for cds_child in db.children(mRNA.id, featuretype='CDS', order_by=mRNA_order_by):
        if cds_child.id == cds.id:
            pre_seq = False
        else:
            if pre_seq:
                seq_segments[0] += cds_child.sequence(genome_fasta,use_strand=True)
            else:
                if len(seq_segments[1]) > 2:
                    break
                else:
                    seq_segments[1] += cds_child.sequence(genome_fasta,use_strand=True)

    mod_3_pre_seq = len(seq_segments[0]) % 3

    # if mod_3_pre_seq > 0:
    #     remainder_pre_seq = seq_segments[0][-mod_3_pre_seq:] 
    # else:
    #     remainder_pre_seq = ''
    
    remainder_pre_seq = seq_segments[0]

    # print('---')
    # print(len(remainder_pre_seq))
    # print(len(seq_segments[1][:2]))
    seq_adj = remainder_pre_seq + seq
    seq_adj += seq_segments[1][:2]
    seq_adj = seq_adj.upper()
    
    
    # print(seq_segments[0])
    # print(len(seq_segments[0]))
    # print(mod_3_pre_seq)
    # print(remainder_pre_seq)
    #print(seq_adj)
    #print(f'rem:{remainder_pre_seq}')

    if mRNA.strand == '+':
        
        pos_adjust = cds.start - len(remainder_pre_seq)
    else:
        # seq_adj = seq_segments[0] + seq_segments
        # seq_adj += seq_segments[1][:2]
        # seq_adj = seq_adj.upper()
        
        pos_adjust = cds.start - 1 - len(seq_segments[1][:2])
    
    return seq_adj, pos_adjust

def write_coding_variants(db, cds, gene, chromosome_variants, results_writer, found_variants, count, cds_count, cds_mismatch,syn_v,nonsyn_v):
    # Get any variants in this cds
    lower_bound_non_syn, upper_bound_non_syn = find_variants_for_feature(cds.start, cds.end, found_variants)
    cds_variants = found_variants[lower_bound_non_syn:upper_bound_non_syn]

    # remove those variants that match from the main list
    #del found_variants[lower_bound_non_syn:upper_bound_non_syn]
    
    # Make note of variants that have been found so we know by the end which variants do not appear in any coding region
    for pos in cds_variants:
        chromosome_variants['found_pos'].append(pos)

    if count < 150000000000000000000: #40000:
        # If there are variants on this cds they are written to results
        if len(cds_variants) > 0:
            cds_count += 1
            seq = cds.sequence(genome_fasta,use_strand=True)
            mRNAs = 0
            for mRNA in db.parents(cds.id,featuretype='mRNA'):
                mRNAs += 1
                # Get an adjusted sequence including the preceding codons
                seq_adj, pos_adjust = get_adjusted_cds(mRNA, seq, cds)
                
                #prot_seq1 = Seq(seq).translate()
                for pos in cds_variants:
                    relative_pos = pos - pos_adjust
                    ref = str(chromosome_variants['variants'][pos]['ref']).upper()
                    for alts in chromosome_variants['variants'][pos]['alts']:
                        for alt in alts:
                            # base at the position of the variant, to check that the reference matches
                            seq_ref = seq_adj[relative_pos]

                            if mRNA.strand == '-':
                                relative_pos = len(seq_adj) - relative_pos #+ 2
                                seq_ref = seq_adj[relative_pos]
                                seq_ref = COMP_DICT[seq_ref]
                                #print(ref)
                                #print(seq_adj[-relative_pos-5:-relative_pos+5])
                            
                            if ref == seq_ref:
                                alt = str(alt).upper()
                                
                                if mRNA.strand == '+':
                                    codons = get_alt_codon(seq_adj, relative_pos, alt)
                                else:
                                    comp_alt = COMP_DICT[alt]
                                    codons = get_alt_codon(seq_adj, relative_pos, comp_alt)
                                    relative_pos += 1
                                    #relative_pos = len(seq_adj) - relative_pos + 1

                                protein_loc = math.ceil(relative_pos / 3)
                                
                                ref_amino = str(Seq(codons[0]).translate())
                                alt_amino = str(Seq(codons[1]).translate())
                                if ref_amino == alt_amino:
                                    syn_v += 1
                                    results_writer.writerow([str(gene.chrom),str(pos),ref,alt,'Synonymous',mRNA.id,protein_loc,ref_amino,'NA'])
                                else:
                                    nonsyn_v += 1
                                    results_writer.writerow([str(gene.chrom),str(pos),ref,alt,'Non-synonymous',mRNA.id,protein_loc,ref_amino,alt_amino])
                                    
                            else:
                                # The reference does not match, skip and make note for later
                                cds_mismatch += 1
                                # print(f'pos:{pos}')print(f'pos:{pos}')
                                # print(f'start:{cds.start}')
                                # print(f'end:{cds.end}')
                                # print(f'rel_pos:{relative_pos}')
                                # print(f'c_len:{len(seq_adj)}')
                                # print(f'start:{cds.start}')
                                # print(f'end:{cds.end}')
                                # print(f'rel_pos:{relative_pos}')
                                # print(f'c_len:{len(seq_adj)}')
                                pass
            if mRNAs == 0:
                print('###')
                for parent in db.parents(cds.id):
                    print(parent.featuretype)
    return results_writer, cds_count, cds_mismatch, count, syn_v, nonsyn_v

def sort_variants_by_feature(variants_chromosome_dict, db, genome_fasta, results_filename):
    try:
        #for feature in db.all
        # if not isfile(results_filename):
        #     raise Exception('file does not exist')
        
        if not (results_filename.endswith('.csv') or results_filename.endswith('.tsv')):
            results_filename += '.tsv'
        
        genes = db.all_features(featuretype='protein_coding_gene')
        with open(results_filename,'w',newline="") as results_tsv:
            results_writer = csv.writer(results_tsv, delimiter='\t')
            results_writer.writerow(variant_results_headers)

            #for i in range(10):
            count = 0
            cds_count = 0
            syn_v = 0
            nonsyn_v = 0
            cds_mismatch = 0
            non_coding = 0
            genome_seq = Fasta(genome_fasta)
            for gene in genes:
                count += 1
                if gene.chrom in variants_chromosome_dict: 
                    chromosome_variants = variants_chromosome_dict[gene.chrom]

                    # Search for variants on corresponding chromosome 
                    lower_bound_variants, upper_bound_variants = find_variants_for_feature(gene.start, gene.end, chromosome_variants['sorted_pos_keys'])
                    found_variants = chromosome_variants['sorted_pos_keys'][lower_bound_variants:upper_bound_variants]
                    # # Make note of variants that have been found so we know by the end which variants do not appear in any coding region
                    # for pos in found_variants:
                    #     chromosome_variants['found_pos'].append(pos)

                    # Search the coding regions of this gene
                    for cds in db.children(gene.id, featuretype='CDS'):
                        results_writer, cds_count, cds_mismatch, count, syn_v, nonsyn_v = write_coding_variants(db, cds, gene, chromosome_variants, results_writer, found_variants, count, cds_count, cds_mismatch,syn_v,nonsyn_v)
                    
                    
                    # for pos in found_variants:
                    #     ref = str(chromosome_variants['variants'][pos]['ref']).upper()
                    #     seq = gene.sequence(genome_fasta,use_strand = True)
                    #     if count < 50000:
                    #         #print('bb')
                    #         #print(gene.attributes)
                    #         gene_len = gene.end - gene.start
                    #         mod_3_gene_len = gene_len % 3
                    #         num_to_append = 3 - mod_3_gene_len
                    #         relative_pos = pos - gene.start
                    #         if gene.strand == '-':
                    #             relative_pos = -relative_pos
                            
                            
                    #         seq = genome_seq[gene.chrom][gene.start:gene.end+num_to_append].seq
                    #         codon, alt_codon = get_alt_codon(seq, relative_pos, 'N')

                            # if gene.strand == '-':
                            #     codon = str(Seq(codon)reverse_complement)

                            #ref_codon = str(Seq(codon).translate())
                        

                        # for alts in chromosome_variants['variants'][pos]['alts']:
                        #     for alt in alts:
                        #         alt = str(alt).upper()
                                #results_writer.writerow([str(gene.chrom),str(pos),ref,alt,'Synonymous',gene.id,'NA',ref_codon,'NA'])
                    
            for chrom in variants_chromosome_dict:
                print(chrom)
                cnt = 0
                for found_variant in variants_chromosome_dict[chrom]['found_pos']:
                    if cnt < 2:
                        found_variant_pos = bisect.bisect_right(variants_chromosome_dict[chrom]['sorted_pos_keys'], found_variant)
                        # print('####')
                        # print(found_variant)
                        # print(variants_chromosome_dict[chrom]['sorted_pos_keys'][lower_bound_pos-1])
                        #cnt += 1
                        if variants_chromosome_dict[chrom]['sorted_pos_keys'][found_variant_pos-1] == found_variant:
                            del variants_chromosome_dict[chrom]['sorted_pos_keys'][found_variant_pos-1]
                
                highest = 0
                # if chrom == 'Pf3D7_02_v3':
                #     print(variants_chromosome_dict[chrom]['sorted_pos_keys'])

                for pos in variants_chromosome_dict[chrom]['sorted_pos_keys']:
                    if pos > highest:
                        highest = pos
                    else:
                        print(pos)
                    ref = str(variants_chromosome_dict[chrom]['variants'][pos]['ref']).upper()
                    if len(variants_chromosome_dict[chrom]['variants'][pos]['alts']) > 2:
                        print(pos)
                    for alts in variants_chromosome_dict[chrom]['variants'][pos]['alts']:
                            if len(alts) > 2:
                                print(pos)
                            for alt in alts:
                                alt = str(alt).upper()
                                non_coding += 1
                                results_writer.writerow([str(chrom),str(pos),ref,alt,'Non-Coding','NA','NA','NA','NA'])



            print(count)
            print('###')
            logger.info(f'{cds_count} variant(s) found in coding regions')
            logger.info(f'{syn_v} synonymous variant(s)')
            logger.info(f'{nonsyn_v} non-synonymous variant(s)')
            logger.info(f'{non_coding} non-coding variant(s)')
            total = syn_v + nonsyn_v + non_coding
            logger.info(f'{total} total variants')
            if cds_mismatch > 0:
                logger.error(f'{cds_mismatch} variant(s) have non-matching reference to refernence coding regions and were excluded')
           
    except Exception as err:
        logger.error(f'Error generating results: {err}; {traceback.format_exc}')
        exc_info = sys.exc_info()
        traceback.print_exception(*exc_info)
        return None

def check_results(results_filename):
    with open(results_filename,'r') as results_file:
        results_reader = csv.reader(results_file,delimiter='\t')
        pos_list = []
        for result in results_reader:
            pos = f'{result[0]}: {result[1]}'
            if pos in pos_list:
                print(pos)
            else:
                pos_list.append(pos)


if __name__ == '__main__':
    #vcf_filename = './data/Toy_Data/testData.vcf'
    vcf_filename = './data/Assessment_Data/assessmentData.vcf.gz'
    gff_filename = './data/Genome_files/PlasmoDB-54_Pfalciparum3D7.gff'
    genome_fasta = './data/Genome_files/PlasmoDB-54_Pfalciparum3D7_Genome.fasta'
    results_filename = './results.tsv'
    chromosome_filepath = './data/Genome_files/PlasmoDB-54_Pfalciparum3D7_Genome.fasta.fai'
    

    logger = logging.getLogger()
    logging.basicConfig(filename='log.txt', encoding='utf-8', filemode="w", format='%(levelname)s - %(asctime)s - %(message)s', level=logging.INFO)
    #logger.setLevel(logging.INFO)
    #handler = logging.StreamHandler()
    #handler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
    #logger.addHandler(handler)


    variants_chromosome_dict = get_variants(vcf_filename, logger)
    total_variants = 0
    for chrom in variants_chromosome_dict:
        #print(chrom)
        total_variants += len(variants_chromosome_dict[chrom]['sorted_pos_keys'])
        highest = 0
        for pos in variants_chromosome_dict[chrom]['sorted_pos_keys']:
            if pos > highest:
                highest = pos
            else:
                print(pos)
        #print(variants_chromosome_dict[chrom]['sorted_pos_keys'])

    #print(total_variants)

    #print(variants_chromosome_dict)
    db = get_feature_db(gff_filename)

    sort_variants_by_feature(variants_chromosome_dict, db, genome_fasta, results_filename)
    check_results(results_filename)
    print('a')
