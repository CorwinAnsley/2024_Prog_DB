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
from matplotlib import pyplot as plt
import seaborn as sns

# Constants
COMP_DICT = {'A': 'T',
             'T': 'A',
             'G': 'C',
             'C': 'G'}

RESULT_HEADERS = ['CHROM','POS','REF','ALT','Type','Transcript','Protein Location','Ref AA','Alt AA']

# Constants for test functions

TEST_OUTPUT_FILENAME = '__TEMP_TEST_OUTPUT__.tsv'
TEST_VARIANTS =[
    ['Chrom1','10','.','A','T','123','.','.']
    ]
TEST_POS_LIST = [1,2,3,10,12,14]
TEST_SEQ = 'AAATTTGGG'


# Loops through variants in vcf file
# Returns a dict of all variants which meet the minimum quatity threshold for each chromosome, and a sorted list of their positions which are used as keys
# Multiple variations can be stored with the same position
def get_variants(vcf_filename, logger, min_quality = 20):
    try:
        if not isfile(vcf_filename):
            raise Exception('file does not exist')
        
        if not (vcf_filename.endswith('.vcf') or vcf_filename.endswith('.vcf.gz')):
            raise Exception('Incorrect file format, should be .vcf/.vcf.gz file')
        
        vcf_reader = vcf.Reader(filename=vcf_filename)

        q_count = 0
        non_q_count = 0
        variants_chromosome_dict = {}
        record_errs = 0

        for record in vcf_reader:
            try:
                if record.QUAL > min_quality:
                    q_count += 1

                    # Make sure the relevant data is valid
                    ref = str(record.REF).upper()
                    alt =[]
                    for i in range(len(record.ALT)):
                        alt.append(str(record.ALT[i]).upper())
                    # assert ref in COMP_DICT
                    # assert alt in COMP_DICT
                    chrom = str(record.CHROM)
                    pos = int(record.POS)
                    assert pos > 0

                    # Add the variation to the dict
                    if chrom in variants_chromosome_dict:
                        if pos in variants_chromosome_dict[chrom]['variants']:
                            variants_chromosome_dict[chrom]['variants'][pos]['alts'].append(alt)
                        else:
                            variants_chromosome_dict[chrom]['variants'][pos] = {'ref':ref, 'alts':[alt]}
                            bisect.insort_right(variants_chromosome_dict[chrom]['sorted_pos_keys'], pos)
                    else:
                        variants_chromosome_dict[chrom] = {}
                        variants_chromosome_dict[chrom]['sorted_pos_keys'] = [pos]
                        variants_chromosome_dict[chrom]['found_pos'] = []
                        variants_chromosome_dict[chrom]['variants'] = {}
                        variants_chromosome_dict[chrom]['variants'][pos] = {'ref':ref, 'alts':[alt]}

                # Keep a count of variations that don't meet the quality requirement
                else:
                    non_q_count += 1
            except Exception as err:
                record_errs += 1


        logger.info(f'{non_q_count} variant(s) with quality score below {min_quality}\n')
        if record_errs > 0:
            logger.error(f'{record_errs} variant(s) were unable to be read')
        return variants_chromosome_dict, q_count, non_q_count
    
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

        if not isfile(db_filepath):
            db = gff.create_db(gff_filename, dbfn=db_filepath, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        else:
            db = gff.FeatureDB(db_filepath, keep_order=True)
        
        return db
    except Exception as err:
        logger.error(f'Error loading gene file: {err}')
        return None

# Helper function that returns upper and lower bound for a sorted list of variant positions that fall
# between a start and an end position
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
# This allows us to calculate the position the variant takes on the protein later
def get_adjusted_cds(mRNA, seq, cds, db):
    # For the negative strand we iterate backwards
    if mRNA == '+':
        mRNA_order_by = 'start'
    else:
        mRNA_order_by = 'end'
    
    # here the sections preceding and following given cds will be saved
    seq_segments = ['','']
    pre_seq = True
    # Iterate through cds children of mRNA to get sequence either side of given cds
    for cds_child in db.children(mRNA.id, featuretype='CDS', order_by=mRNA_order_by):
        if cds_child.id == cds.id:
            pre_seq = False
        else:
            if pre_seq:
                # Add this section
                seq_segments[0] += cds_child.sequence(genome_fasta,use_strand=True)
            else:
                # We only need at most 2 bases from any following cds
                if len(seq_segments[1]) > 2:
                    break
                else:
                    seq_segments[1] += cds_child.sequence(genome_fasta,use_strand=True)
    
    remainder_pre_seq = seq_segments[0]
    
    # assemble the adjusted sequence with the addition of adjacent sequence
    seq_adj = remainder_pre_seq + seq
    seq_adj += seq_segments[1][:2]
    seq_adj = seq_adj.upper()

    # Get a number that will be used to adjust the position we calculate codons from
    if mRNA.strand == '+':
        pos_adjust = cds.start - len(remainder_pre_seq)
    else:
        pos_adjust = cds.start - 1 - len(seq_segments[1][:2])
    
    return seq_adj, pos_adjust

def write_coding_variants(db, cds, gene, chromosome_variants, results_writer, found_variants, cds_count, cds_mismatch,syn_v,nonsyn_v):
    # Get any variants in this cds
    lower_bound_non_syn, upper_bound_non_syn = find_variants_for_feature(cds.start, cds.end, found_variants)
    cds_variants = found_variants[lower_bound_non_syn:upper_bound_non_syn]

    
    # Make note of variants that have been found so we know by the end which variants do not appear in any coding region
    for pos in cds_variants:
        chromosome_variants['found_pos'].append(pos)

    
    # If there are variants on this cds they are written to results
    if len(cds_variants) > 0:
        cds_count += 1
        seq = cds.sequence(genome_fasta,use_strand=True)
        for mRNA in db.parents(cds.id,featuretype='mRNA'):
            # Get an adjusted sequence including the preceding codons
            seq_adj, pos_adjust = get_adjusted_cds(mRNA, seq, cds, db)
            
            for pos in cds_variants:
                ref = chromosome_variants['variants'][pos]['ref']
                for alts in chromosome_variants['variants'][pos]['alts']:
                    for alt in alts:
                        # base at the position of the variant, to check that the reference matches
                        relative_pos = pos - pos_adjust
                        seq_ref = seq_adj[relative_pos]

                        if mRNA.strand == '-':
                            relative_pos = len(seq_adj) - relative_pos
                            seq_ref = seq_adj[relative_pos]
                            seq_ref = COMP_DICT[seq_ref]
                        
                        if ref == seq_ref:
                            if mRNA.strand == '+':
                                codons = get_alt_codon(seq_adj, relative_pos, alt)
                                protein_loc = math.ceil(relative_pos / 3)
                            else:
                                comp_alt = COMP_DICT[alt]
                                codons = get_alt_codon(seq_adj, relative_pos, comp_alt)
                                protein_loc = math.ceil((relative_pos + 1) / 3)
                            
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
                            pass

    return results_writer, cds_count, cds_mismatch, syn_v, nonsyn_v

def create_variants_report(variants_chromosome_dict, db, genome_fasta, results_filename, overwrite = True):
    try:
        # Make sure output file is correct
        if not (results_filename.endswith('.csv') or results_filename.endswith('.tsv')):
            results_filename += '.tsv'
        
        if not overwrite:
            if isfile(results_filename):
                raise Exception('Output file already exists. Use --overwrite to overwrite')
        
        # Get a generator for all the genes
        genes = db.all_features(featuretype='protein_coding_gene')
        with open(results_filename,'w',newline="") as results_tsv:
            results_writer = csv.writer(results_tsv, delimiter='\t')
            results_writer.writerow(RESULT_HEADERS)

            # counts of things we want to keep track of
            cds_count = 0
            syn_v = 0
            nonsyn_v = 0
            cds_mismatch = 0
            non_coding = 0
            
            for gene in genes:
                if gene.chrom in variants_chromosome_dict: 
                    chromosome_variants = variants_chromosome_dict[gene.chrom]

                    # Search for variants on corresponding chromosome 
                    lower_bound_variants, upper_bound_variants = find_variants_for_feature(gene.start, gene.end, chromosome_variants['sorted_pos_keys'])
                    found_variants = chromosome_variants['sorted_pos_keys'][lower_bound_variants:upper_bound_variants]
                   
                    # Search the coding regions of this gene
                    for cds in db.children(gene.id, featuretype='CDS'):
                        try:
                            results_writer, cds_count, cds_mismatch, syn_v, nonsyn_v = write_coding_variants(db, cds, gene, chromosome_variants, results_writer, found_variants, cds_count, cds_mismatch,syn_v,nonsyn_v)
                        except Exception as err:
                            try:
                                logger.error(f'Error processing cds {cds.id}: {err}')
                            except:
                                logger.error(f'Error processing cds (unable to id): {err}')
                    
            # Write the non-coding variants
            for chrom in variants_chromosome_dict:
                # Remove all the variants found on coding regions from the sorted position list 
                for found_variant in variants_chromosome_dict[chrom]['found_pos']:
                    found_variant_pos = bisect.bisect_right(variants_chromosome_dict[chrom]['sorted_pos_keys'], found_variant)
                    if variants_chromosome_dict[chrom]['sorted_pos_keys'][found_variant_pos-1] == found_variant:
                        del variants_chromosome_dict[chrom]['sorted_pos_keys'][found_variant_pos-1]

                # All the remaining variants are non-coding, write them as such
                for pos in variants_chromosome_dict[chrom]['sorted_pos_keys']:
                    ref = variants_chromosome_dict[chrom]['variants'][pos]['ref']
                    for alts in variants_chromosome_dict[chrom]['variants'][pos]['alts']:
                            for alt in alts:
                                non_coding += 1
                                results_writer.writerow([str(chrom),str(pos),ref,alt,'Non-Coding','NA','NA','NA','NA'])



            logger.info(f'{cds_count} coding regions found to contain variants')
            logger.info(f'{syn_v} synonymous variant(s)')
            logger.info(f'{nonsyn_v} non-synonymous variant(s)')
            logger.info(f'{non_coding} non-coding variant(s)')
            total = syn_v + nonsyn_v + non_coding
            logger.info(f'{total} total variants')
            if cds_mismatch > 0:
                logger.error(f'{cds_mismatch} variant(s) have non-matching reference to reference coding regions and were excluded')

            return syn_v, nonsyn_v
           
    except Exception as err:
        logger.error(f'Error generating results: {err};')
        return None

def create_variants_plot(plot_filename, non_q_count, syn_v, nonsyn_v):
    #x_axis = ['Q <= 20','Synonymous','Non-synonymous']
    data = {'Q <= 20': non_q_count, 'Synonymous': syn_v, 'Non-synonymous':nonsyn_v}
    sns.barplot(data)
    plt.savefig(plot_filename)
    plt.clf()

#### TESTING FUNCTIONS ###

def test_get_variants(vcf_filename, logger):
    try:
        variants_chromosome_dict, q_count, non_q_count = get_variants(vcf_filename, logger)
        assert 'Pf3D7_10_v3' in variants_chromosome_dict.keys()
        assert variants_chromosome_dict['Pf3D7_10_v3']['variants'][1468717]['ref'] == 'C'
        return variants_chromosome_dict
    except Exception as err:
        logger.error(f'Error occurred testing get_variants: {err}')
        return False

def test_get_feature_db(gff_filename, logger):
    try:
        db = get_feature_db(gff_filename)
        assert type(db) == gff.interface.FeatureDB
        return db
    except Exception as err:
        logger.error(f'Error occurred testing get_feature_db: {err}')
        return False

def test_find_variants_for_feature(logger):
    try:
        lower_upper_bounds = find_variants_for_feature(10,13,TEST_POS_LIST)
        assert lower_upper_bounds == (3,5)
        return True
    except Exception as err:
        logger.error(f'Error occurred testing find_variants_for_feature: {err}')
        return False

def test_get_alt_codon(logger):
    try:
        codons = get_alt_codon(TEST_SEQ,4,'C')
        assert codons == ('TTT','TCT')
        return True
    except Exception as err:
        logger.error(f'Error occurred testing get_alt_codon: {err}')
        return False

def test_create_variants_report(variants_chromosome_dict, db, genome_fasta):
    try:
        count_results = create_variants_report(variants_chromosome_dict,db,genome_fasta,TEST_OUTPUT_FILENAME,overwrite=True)
        with open( TEST_OUTPUT_FILENAME,'r') as results_file:
            results_reader = csv.reader(results_file,delimiter='\t')
            assert next(results_reader) == RESULT_HEADERS
            assert next(results_reader) == ['Pf3D7_10_v3', '1468717', 'C', 'A', 'Non-synonymous', 'PF3D7_1037000.1', '1948', 'T', 'N']
        assert count_results == (0, 2)
        return True
    except Exception as err:
        logger.error(f'Error occurred testing create_variants_report: {err}')
        return False

# Run all tests
def run_tests(test_dir,logger):
    vcf_filename = test_dir + r"/Toy Data/testData.vcf"
    genome_fasta = test_dir + r"/Genome files/PlasmoDB-54_Pfalciparum3D7_Genome.fasta"
    
    test_results = []
    test_variants_dict = test_get_variants(vcf_filename, logger)
    test_results.append(test_variants_dict)

    test_results.append(test_find_variants_for_feature(logger))
    test_results.append(test_get_alt_codon(logger))

    test_db = test_get_feature_db(gff_filename,logger)
    test_results.append(test_db)

    test_results.append(test_create_variants_report(test_variants_dict,test_db,genome_fasta))
    
    if False or None in test_results:
        raise Exception('One or more tests failed, see log for details')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Variant Parser', description='')
    parser.add_argument('--overwrite',help='overwrites output',action='store_true')
    parser.add_argument('--vcf',type=str,help='Filepath to vcf file',default=None)
    parser.add_argument('--gff',type=str,help='Filepath to gff file',default=None)
    parser.add_argument('--fasta',type=str,help='Filepath to genome fasta file',default=None)
    parser.add_argument('--output',type=str,help='Filepath to output to',default=None)
    parser.add_argument('--testdir',type=str,help='Filepath to test files, if tests are to be run',default='')
    parser.add_argument('--plotdir',type=str,help='Filepath to write plot to',default=None)

    args = parser.parse_args()
    vcf_filename = str(args.vcf)
    gff_filename = str(args.gff)
    genome_fasta = str(args.fasta) 
    results_filename = str(args.output)
    test_dir = str(args.testdir)
    plot_filename = str(args.plotdir)  

    logger = logging.getLogger()
    logging.basicConfig(filename='log.txt', encoding='utf-8', filemode="w", format='%(levelname)s - %(asctime)s - %(message)s', level=logging.INFO)

    if test_dir == '':
        variants_chromosome_dict, q_count, non_q_count = get_variants(vcf_filename, logger)
        #print(variants_chromosome_dict)
    
    
        db = get_feature_db(gff_filename)
        
        if variants_chromosome_dict != None:
            if db != None:
                counts = create_variants_report(variants_chromosome_dict, db, genome_fasta, results_filename)
                if counts != None:
                    syn_v, nonsyn_v = counts
                    if plot_filename != None:
                        create_variants_plot(plot_filename, non_q_count, syn_v, nonsyn_v)
    else:
        run_tests(test_dir,logger)
  
