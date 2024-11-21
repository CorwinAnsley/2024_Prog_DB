import csv
import argparse
import re

### Apologies for the excessive and overly verbose comments, wanted to being absolutely sure that 
### there was not too few

TEST_SAM = 'Toxo_chr8_subset.sam'
TEST_CIGAR = '1S10M7887N89M10N'

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

# returns the junctions form cigar string and position
def get_junctions(cigar, pos):
    junctions = []

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



# Returns a dictionary; every key is a junction and the values are the number of reads matching the associated junction
def create_junctions_dict(read_list):
    chromosome_junction_dict = {}

    for read in read_list:
        if read['NH:i:x'] == 1:
            if 'N' in read['CIGAR']:
                read_junctions = get_junctions(read['CIGAR'],read['POS'])
                if read['RNAME'] not in chromosome_junction_dict:
                    chromosome_junction_dict[read['RNAME']] = {}
                for junction in read_junctions:

        pass


if __name__ == '__main__':
    read_list =  parse_sam(TEST_SAM)
    for i in range(5):
        print(read_list[i])

    get_junctions(TEST_CIGAR, 0)
                
