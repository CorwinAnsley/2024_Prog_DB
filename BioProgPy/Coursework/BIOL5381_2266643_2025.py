import argparse
import logging
import os
import bisect
from os.path import isfile
import vcf
import gffutils as gff

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
                        variants_chromosome_dict[record.CHROM]['variants'][int(record.POS)]['alts'].append(record.ALT)
                    else:
                        variants_chromosome_dict[record.CHROM]['variants'][int(record.POS)] = {'ref':record.REF, 'alts':[record.ALT]}
                        bisect.insort_right(variants_chromosome_dict[record.CHROM]['sorted_pos_keys'], record.POS)
                else:
                    variants_chromosome_dict[record.CHROM] = {}
                    variants_chromosome_dict[record.CHROM]['sorted_pos_keys'] = [int(record.POS)]
                    #variants_chromosome_dict[record.CHROM]['found_pos'] = []
                    variants_chromosome_dict[record.CHROM]['variants'] = {}
                    variants_chromosome_dict[record.CHROM]['variants'][int(record.POS)] = {'ref':record.REF, 'alts':[record.ALT]}

                # if q_count < 10:
                #     print(record.ID)
            else:
                non_q_count += 1
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


def sort_variants_by_feature(variants_chromosome_dict, db, results_filename):
    try:
        #for feature in db.all
        features = db.all_features()
        for i in range(50):
            for child in next(features).children:
                print(child.featuretype)
    except Exception as err:
        logger.error(f'Error generatng variants: {err}')
        return None


if __name__ == '__main__':
    #vcf_filename = './data/Toy_Data/testData.vcf'
    vcf_filename = './data/Assessment_Data/assessmentData.vcf.gz'
    gff_filename = './data/Genome_files/PlasmoDB-54_Pfalciparum3D7.gff'
    

    logger = logging.getLogger()
    logging.basicConfig(filename='log.txt', encoding='utf-8', filemode="w", format='%(levelname)s - %(asctime)s - %(message)s', level=logging.INFO)
    #logger.setLevel(logging.INFO)
    #handler = logging.StreamHandler()
    #handler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
    #logger.addHandler(handler)
    
    variants_chromosome_dict = get_variants(vcf_filename, logger)
    db = get_feature_db(gff_filename)
    sort_variants_by_feature(variants_chromosome_dict, db, '')
    print('a')
