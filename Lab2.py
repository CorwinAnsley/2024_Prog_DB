import csv

WORKDIR = '/home1/bioinfo-03/repos/2024_Prog_DB/Programming-Lab2'
SEQ_FILE_PATH = f'{WORKDIR}/sequences.txt'
DATA_FILE_PATH = f'{WORKDIR}/data1.csv'
NUM = 10
EXP = 4

ADPTR_SEQ_LEN = 14

def get_trimmed_seqs(seq_filepath, adpt_seq_len = ADPTR_SEQ_LEN):
    seq_dict = {}
    with open(seq_filepath) as seq_file:

        for line in seq_file:
            seq_info = line.rstrip().split('\t')
            try:
                seq_dict[seq_info[0]] = seq_info[1][adpt_seq_len:]
            except:
                print(f'Invalid sequence: {seq_info}')
    return seq_dict

def output_trimmed_fasta(seq_filepath, output_path, output_filename="Seq.fasta", adpt_seq_len = ADPTR_SEQ_LEN):
    seq_dict = get_trimmed_seqs(seq_filepath, adpt_seq_len)
    output_filepath = f'{output_path}/{output_filename}'
    with open(output_filepath, 'w') as file:
        for seq in seq_dict:
            file.write(f'> {seq}\n{seq_dict[seq]}\n')

def create_gene_data_info(filepath, output_path, output_filename="gene_data.csv"):
    gene_data = []
    with open(filepath) as file:
        gene_reader = csv.reader(file, delimiter='\t')
        header = next(gene_reader)
        header.append('Ratio')
        header.append('Sum')
        gene_data.append(header)

        for row in gene_reader:
            print(row)
            ratio = int(row[1])/int(row[2])
            row.append(f'{ratio:.2f}')
            sum = 0
            for i in row[1:]:
                sum += float(i)
            row.append(str(sum))
            gene_data.append(row)
    
    output_filepath = f'{output_path}/{output_filename}'
    with open(output_filepath, 'w') as file:
        gene_writer = csv.writer(file, delimiter=',')
        for row in gene_data:
            gene_writer.writerow(row)

def make_powers_list(max_exp, num):
    powers_list = []
    for exp in range(max_exp+1):
        powers_list.append(num**exp)
    
    return powers_list

def cross_multiply(list1, list2):
    return_list = []
    for i in list1:
        for j in list2:
            i = i*j
        return_list.append(i)
    
    return return_list


#output_trimmed_fasta(FILE_PATH, WORKDIR)
#
# create_gene_data_info(DATA_FILE_PATH, WORKDIR)

list1 = range(1,6)
list2 = make_powers_list(EXP, NUM)
new_list = cross_multiply(list1, list2)
print(new_list)