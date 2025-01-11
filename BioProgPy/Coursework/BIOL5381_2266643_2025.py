import argparse
import logging
import vcf
import bisect

def get_variants(vcf_filename, logger):
    vcf_reader = vcf.Reader(filename=vcf_filename)

    q_count = 0
    variants_chromosome_dict = {}

    for record in vcf_reader:
        if record.QUAL > 20:
            q_count += 1
            if record.CHROM in variants_chromosome_dict:
                variants_chromosome_dict[record.CHROM]['variants'][int(record.POS)] = {'ref':record.REF, 'alt':record.ALT}
                bisect.insort_right(variants_chromosome_dict[record.CHROM]['sorted_pos_keys'], record.POS)
            else:
                variants_chromosome_dict[record.CHROM] = {}
                variants_chromosome_dict[record.CHROM]['sorted_pos_keys'] = [int(record.POS)]
                variants_chromosome_dict[record.CHROM]['variants'] = {}
                variants_chromosome_dict[record.CHROM]['variants'][int(record.POS)] = {'ref':record.REF, 'alt':record.ALT}

            if q_count < 10:
                print(record.ID)
    print(q_count)

#def 

if __name__ == '__main__':
    #vcf_filename = './data/Toy_Data/testData.vcf.gz'
    vcf_filename = './data/Assessment_Data/assessmentData.vcf.gz'

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('‘%(levelname)s - %(asctime)s - %(message)s’)'))
    logger.addHandler(handler)
    
    get_variants(vcf_filename, logger)
