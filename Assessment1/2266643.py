import csv
import argparse

### Apologies for the excessive and overly verbose comments, wanted to being absolutely sure that 
### there was not too few

TEST_SAM = 'Toxo_chr8_subset.sam'
TEST_CIGAR = '1S10M7887N89M'

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
    junction_list = []
    
    # This placeholder string holds section of the cigar string as we iterate over it
    cigar_section = ''

    for c in cigar:
        # add digits to the placeholder
        if c.isdigit():
            cigar_section += str(c)
        # if there's an 'M' we have a match
        elif c == 'M':
            # adjust the position and erase the placeholder string
            pos += int(cigar_section)
            cigar_section = ''
            # if there's an 'N' we have a skipped region
        elif c == 'N':
            # If there's a skipped region we have a junction, using the position it is defined as string START:END 
            # and added to the list of junctions
            junction = f'{pos}:{pos + int(cigar_section)}'
            junction_list.append(junction)
            cigar_section = ''
        # If there's an 'S', 'I' or anything else it is ignored along with the preceding digits
        else:
            cigar_section = ''
        
        print(junction_list)



# Returns a dictionary; every key is a junction and the values are the number of reads matching the associated junction
def create_junctions_dict(read_list):
    junction_dict = {}

    for read in read_list:
        pass


if __name__ == '__main__':
    read_list =  parse_sam(TEST_SAM)
    #for i in range(25):
    #    print(read_list[i])

    get_junctions(TEST_CIGAR, 0)
                
