import csv
import argparse
import re

### Apologies for the excessive and overly verbose comments, wanted to being absolutely sure that 
### there was not too few

TEST_SAM = 'Toxo_chr8_subset.sam'
TEST_CIGAR = '1S10M7887N89M10N'

TEST_GENE_FILE = 'GenesByLocation_Summary.txt'

#  Dict of which colum in the sam file the relevant data is stored
SAM_COLUMNS = {'RNAME':2, 'POS':3, 'CIGAR':5, 'NH:i:x':-1}

# Parses sam file to extract the relevant read data and return it as a list of dictionaries
def parse_sam(filepath):
    with open(filepath) as csv_file:
        sam_reader = csv.reader(csv_file,delimiter='\t')
        read_list = []
        for row in sam_reader:
            # Skip the header rows
            if row[0].startswith('@'):
                pass
            else:
                # Obtain the relevant data and store in dictionary
                read_dict = {}
                for column in SAM_COLUMNS:
                    read_dict[column] = row[SAM_COLUMNS[column]]
                read_dict['NH:i:x'] = int(read_dict['NH:i:x'].split(':')[-1]) #Only the integer at the end of NH:i:x is needed (i.e. x)

                # Append the dict with read data to the read list
                read_list.append(read_dict)
        return read_list

# returns the junctions from cigar string and position
def get_junctions(cigar, pos):
    junctions = []
    # Make sure pos in an int
    pos = int(pos)

    # get all sections of the cigar string that are a series of digits followed by either M or N
    regex = '(\d+)([M,N])'
    sections = re.finditer(regex, cigar)
    
    # Iterate over the sections
    for section in sections:
        # for a match in the cigar string we adjust the position along
        if section.group(2) == 'M':
            pos += int(section.group(1))
        # for a skip we add the junction to the junction list and adjust the position along 
        elif section.group(2) == 'N':
            junction = (pos,pos + int(section.group(1)))
            junctions.append(junction)
            pos += int(section.group(1))
    
    return junctions



# Returns a dictionary; each key is a chromosome the values are dictionaries for 
# all the junctions in that chromosome, identified by the start position.
# Each junction is also dictionary containing the number of supporting reads and the end position
# This is a little convoluted but will make later steps more efficient as we only need to search junctions on the right chromosome
def create_chromosome_junctions_dict(read_list):
    # Creating an empty dictionary where the junctions will be saved
    chromosome_junction_dict = {}

    # Iterate over all reads
    for read in read_list:
        # Check the read only aligns once
        if read['NH:i:x'] == 1:
            # Check if the read has a split
            if 'N' in read['CIGAR']:
                # Get all the junctions from the read cigar string
                read_junctions = get_junctions(read['CIGAR'],read['POS'])
                # If the chromosome this read is from is not in the dict, add it
                if read['RNAME'] not in chromosome_junction_dict:
                    chromosome_junction_dict[read['RNAME']] = {}
                
                # Add each junction to the dictionary 
                for junction in read_junctions:
                    if junction[0] in chromosome_junction_dict[read['RNAME']].keys():
                        # If the junction is already in the chromosome dictionary simply increase the number of supporting reads for that junction
                        chromosome_junction_dict[read['RNAME']][junction[0]]['number_supporting_reads'] += 1
                    else:
                        # If the junction is not yet in the dictionary create a new dict with the key of the start position
                        chromosome_junction_dict[read['RNAME']][junction[0]] = {}
                        # Add the number of supporting reads (1 for now) and the end position
                        chromosome_junction_dict[read['RNAME']][junction[0]]['number_supporting_reads'] = 1
                        chromosome_junction_dict[read['RNAME']][junction[0]]['end_pos'] = junction[1]
    return chromosome_junction_dict

def create_gene_dict_from_file(filepath):
    with open(filepath) as csv_file:
        gene_reader = csv.reader(csv_file,delimiter='\t')
        chromosome_gene_dict = {}
        for gene in gene_reader:
            gene_location = gene[2]
            chromosome = gene_location.split(':')[0]
            if chromosome not in chromosome_gene_dict.keys():
                chromosome_gene_dict[chromosome] = {}
            
            #If the same gene is repeated multiple times we will only use the first instance
            if gene[0] not in chromosome_gene_dict[chromosome]:
                chromosome_gene_dict[chromosome][gene[0]] = {}
                position = gene_location.split(':')[1].split('..')[:-3]
                start_pos = position[0].replace(',','')
                end_pos = position[1].replace(',','')[:-3]
                chromosome_gene_dict[chromosome][gene[0]]['start_pos'] = start_pos
                chromosome_gene_dict[chromosome][gene[0]]['end_pos'] = end_pos

    return chromosome_gene_dict
# # Output a 
# def create_gene_junctions_list():


if __name__ == '__main__':
    read_list =  parse_sam(TEST_SAM)

    chromosome_junction_dict = create_chromosome_junctions_dict(read_list)
    #print(chromosome_junction_dict.keys())
    #example_chr = next(chromosome_junction_dict.keys())
    #print(chromosome_junction_dict)
    #  for i in range(5):
    #      print(read_list[i])

    #get_junctions(TEST_CIGAR, 0)

    genes = create_gene_dict_from_file(TEST_GENE_FILE)
    print(genes)
                
