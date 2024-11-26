import csv
import argparse
import re
import bisect

#  Dict of which colum in the sam file the relevant data is stored
SAM_COLUMNS = {'RNAME':2, 'POS':3, 'CIGAR':5, 'NH:i:x':-1}

# Hardcoded output filename
OUTPUT_FILE = '2266643.txt'

## Constants and filenames for testing functions ###
TEST_SAM = 'Toxo_chr8_subset.sam'
TEST_CIGAR = '1S10M7887N89M10N'
TEST_GENE_FILE = 'GenesByLocation_Summary.txt'
TEST_READ_LIST = [{'RNAME': 'TGME49_chrVIII', 'POS': '0', 'CIGAR': '1S10M20N10M20N', 'NH:i:x': 1},
                  {'RNAME': 'TGME49_chrVIII', 'POS': '0', 'CIGAR': '1S11M20N10M20N', 'NH:i:x': 2},
                  {'RNAME': 'TGME49_chrVIII', 'POS': '0', 'CIGAR': '1S10M21N10M20N', 'NH:i:x': 1},
                  ]
TEST_CHROMOSOME_JUNCTION_DICT =  {'TGME49_chrVIII': {10: {30: 1, 31: 1}, 40: {60: 1}, 41: {61: 1}}}
TEST_CHROMOSOME_GENE_DICT = {'TGME49_chrVIII':{'TEST_GENE':{'start_pos': 10, 'end_pos': 60}}}

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
# all the junction start positions in that chromosome, identified by the start position.
# Each junction start position is also dictionary containing the end positions as keys with the number of supporting reads as the corresponding value
# so multiple junctions with the same start position can be stored. 
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
                    # Check if start position is in the dict, if it has check the end position
                    if junction[0] in chromosome_junction_dict[read['RNAME']].keys():
                        # If the junction is already in the chromosome dictionary simply increase the number of supporting reads for that junction
                        if junction[1] in chromosome_junction_dict[read['RNAME']][junction[0]]:
                            chromosome_junction_dict[read['RNAME']][junction[0]][junction[1]] += 1
                        # If the junction is not yet in the dictionary create a new dict with the key of the start position
                        else:
                            chromosome_junction_dict[read['RNAME']][junction[0]][junction[1]] = 1            

                    else:
                        # If the start position is not yet in the dictionary create a new dict with the key of the start position
                        chromosome_junction_dict[read['RNAME']][junction[0]] = {}
                        # Add the number of supporting reads (1 for now) as  value with the end position as the key
                        chromosome_junction_dict[read['RNAME']][junction[0]][junction[1]] = 1
    return chromosome_junction_dict

# Returns a dict of each chromosome with each gene for the respective chromosome and its location
def create_gene_dict_from_file(filepath, headers = True):
    with open(filepath) as csv_file:
        gene_reader = csv.reader(csv_file,delimiter='\t')
        # remove headers, if the file has them
        if headers:
            headers = next(gene_reader)

        # Create the empty dict that will be populated
        chromosome_gene_dict = {}

        # Iterate over the rows in the file, i.e. the genes
        for gene in gene_reader:
            # Make sure the gene row has the right number of columns, if not we skip it
            try:
                assert len(gene) == 3
            except:
                continue

            # The location is in the 3rd column
            gene_location = gene[2]
            # The chromosome and the (start and end) position are split by a ':', here we take just the chromosome
            chromosome = gene_location.split(':')[0]

            # If the chromosome does not exist as a key in our dict we create an empty dict corresponding to it
            if chromosome not in chromosome_gene_dict.keys():
                chromosome_gene_dict[chromosome] = {}
            
            # If the same gene is repeated multiple times we will only use the first instance
            if gene[0] not in chromosome_gene_dict[chromosome]:
                # Add the gene as an empty dict into dict for the corresponding chromosome
                chromosome_gene_dict[chromosome][gene[0]] = {}
                try:
                    # Split the position from the chromosome and split that into the start and the end positions
                    position = gene_location.split(':')[1]
                    position = position.split('..')

                    # Assign the start and end positions to variables, removing extraneous characters and converting to int
                 
                    start_pos = int(position[0].replace(',',''))
                    end_pos = position[1].replace(',','')
                    end_pos = int(end_pos[:-3]) # Removing the end character which specify which strand, e.g. (+)
                except:
                    # An error occured turning the positions to ints, removing the gene and moving on to the next one
                    chromosome_gene_dict[chromosome].pop(gene[0], None)
                    continue

                # Add the start and end positions to the dict for the gene
                chromosome_gene_dict[chromosome][gene[0]]['start_pos'] = start_pos
                chromosome_gene_dict[chromosome][gene[0]]['end_pos'] = end_pos

    return chromosome_gene_dict

# Using the sorted list of junction start positions and dict of junctions, return all junctions that fall within the start and end of a gene
def get_junctions_within_gene(junction_start_list, junction_dict, gene_start, gene_end):
    
    # The junction start position list is sorted so we can use binary searches to find junction starts that fall between the start and end gene positions
    lower_bound_junctions = bisect.bisect_left(junction_start_list, gene_start)
    upper_bound_junctions = bisect.bisect_right(junction_start_list, gene_end, lo=lower_bound_junctions)

    matching_junction_start_position = junction_start_list[lower_bound_junctions:upper_bound_junctions]

    gene_junctions = {}

    # Check every junction with these start positions to make sure the end position is not greater than the end position of the gene
    for junction_start in matching_junction_start_position:
        for end_pos in junction_dict[junction_start]:
            if end_pos <= gene_end:
                # Add the junction to the gene_junctions dict with the key as a tuple of the start and end positions and the value of the number of supporting reads
                junction = (junction_start, end_pos)
                gene_junctions[junction] = junction_dict[junction_start][end_pos]            
    
    return gene_junctions

def output_gene_junctions(output_path, chromosome_gene_dict, chromosome_junction_dict):
    with open(output_path, 'w') as output_file:
        for chromosome in chromosome_gene_dict:
            junction_dict = chromosome_junction_dict[chromosome]
            junction_start_list = list(junction_dict.keys())
            junction_start_list.sort()

            for gene in chromosome_gene_dict[chromosome]:
                gene_start = chromosome_gene_dict[chromosome][gene]['start_pos']
                gene_end = chromosome_gene_dict[chromosome][gene]['end_pos']
                gene_junctions = get_junctions_within_gene(junction_start_list, junction_dict, gene_start, gene_end)

                for junction in gene_junctions:
                    junction_row = f'{gene}\t{junction[0]}\t{junction[1]}\t{gene_junctions[junction]}\n'
                    output_file.write(junction_row)

                if len(gene_junctions) > 0:
                    output_file.write('\n')


### TESTING FUNCTIONS ###

# Testing function for parse_same
def test_parse_sam():
    read_list = parse_sam(TEST_SAM)
    assert read_list[0] == {'RNAME': 'TGME49_chrVIII', 'POS': '6625858', 'CIGAR': '1S10M7887N89M', 'NH:i:x': 1}

    for read in read_list:
        for column in SAM_COLUMNS:
            assert column in read

# Testing function for get_junctions
def test_get_junctions():
    junctions = get_junctions(TEST_CIGAR, 0)
    assert junctions == [(10, 7897), (7986, 7996)]

# Testing function for create_chromosome_junctions_dict
def test_create_chromosome_junctions_dict():
    chromosome_junction_dict = create_chromosome_junctions_dict(TEST_READ_LIST)
    assert chromosome_junction_dict == {'TGME49_chrVIII': {10: {30: 1, 31: 1}, 40: {60: 1}, 41: {61: 1}}}

# Testing function for create_gene_dict_from_file
def test_create_gene_dict_from_file():
    chromosome_gene_dict = create_gene_dict_from_file(TEST_GENE_FILE)
    assert chromosome_gene_dict['TGME49_chrVIII']['TGME49_268220'] == {'start_pos': 6631349, 'end_pos': 6636865}

# Testing function for get_junctions_within_gene
def test_get_junctions_within_gene():
    junction_dict = TEST_CHROMOSOME_JUNCTION_DICT['TGME49_chrVIII']
    junction_start_list = list(junction_dict.keys())
    junction_start_list.sort()

    gene_junctions = get_junctions_within_gene(junction_start_list, junction_dict, 10, 60)

    assert gene_junctions == {(10, 30): 1, (10, 31): 1, (40, 60): 1}

# Testing function for test_output_gene_junctions
def test_output_gene_junctions():
    output_gene_junctions('test.txt', TEST_CHROMOSOME_GENE_DICT, TEST_CHROMOSOME_JUNCTION_DICT)

    with open('test.txt', 'r') as test_output_file:
        assert test_output_file.readline() == 'TEST_GENE\t10\t30\t1\n'

if __name__ == '__main__':
    read_list =  parse_sam(TEST_SAM)

    chromosome_junction_dict = create_chromosome_junctions_dict(read_list)
    #print(chromosome_junction_dict.keys())
    #example_chr = next(chromosome_junction_dict.keys())
    #print(chromosome_junction_dict)
    #  for i in range(5):
    #      print(read_list[i])

    #get_junctions(TEST_CIGAR, 0)

    chromosome_gene_dict = create_gene_dict_from_file(TEST_GENE_FILE)

    output_gene_junctions(OUTPUT_FILE, chromosome_gene_dict, chromosome_junction_dict)
    #print(genes)
    
    test_parse_sam()

    test_get_junctions()

    test_create_chromosome_junctions_dict()

    test_create_gene_dict_from_file()

    test_get_junctions_within_gene()

    test_output_gene_junctions()