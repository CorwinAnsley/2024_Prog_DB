FILE_PATH = '/home1/bioinfo-03/repos/2024_Prog_DB/Programming-Lab2/sequences.txt'
ADPTR_SEQ_LEN = 14

def get_trimmed_seqs(seq_filepath):
    seq_dict = {}
    with open(seq_filepath) as seq_file:

        for line in seq_file:
            seq_info = line.split('\t')
            seq_dict[seq_info[0]] = seq_info[1][ADPTR_SEQ_LEN:-1]
    return seq_dict

print(get_trimmed_seqs(FILE_PATH))