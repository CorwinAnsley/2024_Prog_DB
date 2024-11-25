import csv
import argparse
import re
import bisect

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
            # Make sure the gene row has the right number of columns
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

#     # Make sure gene start is less than the gene end
#     assert gene_start < gene_end

#     indexes_to_remove = []

#     # filter out all junctions where the start position is less than the start position of the gene (the list is in ascending order so we don't need to check every entry)
#     for index, junction_start in enumerate(junction_start_list):
#         if junction_start < gene_start:
#             indexes_to_remove.append(index)
#         else:
#             break

#    # filter out junctions where the start position is greater than end position of the gene (Checked in descending order, again not every entry needs to be checked)
#     for index, junction_start in enumerate(reversed(junction_start_list)):
#         if junction_start > gene_end:
#             indexes_to_remove.append(index)
#         else:
#             break
    
#     # Check every junction remaining to make sure the end position is not greater than the end position of the gene
#     for junction_start in junction_start_list:
#         if junction_dict[junction_start]['end_pos'] > gene_end:
#             junction_start_list.pop(junction_start_list.index(junction_start))
    lower_bound_junctions = bisect.bisect_left(junction_start_list, gene_start)
    upper_bound_junctions = bisect.bisect_right(junction_start_list, gene_end, lo=lower_bound_junctions)
            
    
    return junction_start_list

def output_gene_junctions(output_path, chromosome_gene_dict, chromosome_junction_dict):
    for chromosome in chromosome_gene_dict:
        junction_dict = chromosome_junction_dict[chromosome]
        junction_start_list = list(junction_dict.keys())
        junction_start_list.sort()

        for gene in chromosome_gene_dict[chromosome]:
            gene_start = chromosome_gene_dict[chromosome][gene]['start_pos']
            gene_end = chromosome_gene_dict[chromosome][gene]['end_pos']
            gene_junctions = get_junctions_within_gene(junction_start_list, junction_dict, gene_start, gene_end)
            print(gene_junctions)
            break




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

    output_gene_junctions('', chromosome_gene_dict, chromosome_junction_dict)
    #print(genes)
                
