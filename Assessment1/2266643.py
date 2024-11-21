import csv
import argparse

### Apologies for the excessive and overly verbose comments, wanted to being absolutely sure that 
### there was not too few

TEST_SAM = 'Toxo_chr8_subset.sam'

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

# Returns a dictionary; every key is a junction and the values are the number of reads matching the associated junction
def get_junctions(read_list):
    junction_dict = {}

    for read in read_list:
        pass


if __name__ == '__main__':
    read_list =  parse_sam(TEST_SAM)
    for i in range(25):
        print(read_list[i])
                
