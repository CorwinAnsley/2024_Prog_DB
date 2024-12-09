import sys

def create_pos_weight_dict(fasta_filepath):
    
    with open(fasta_filepath) as fasta:
        amino_pos_dict = {}
        max_len = 0
        for line in fasta:
            if line.startswith('>'):
                id = line.rstrip().lstrip('>')
                seq = next(fasta).rstrip().upper()
                for i in range(len(seq)):
                    if i > max_len:
                        max_len = i
                    amino = seq[i]
                    if amino not in amino_pos_dict:
                        amino_pos_dict[amino] = {}
                    
                    if i not in amino_pos_dict[amino]:
                        amino_pos_dict[amino][i] = 1
                    else:
                        amino_pos_dict[amino][i] += 1
    
    return amino_pos_dict, max_len

def print_pos_weight_matrix(amino_pos_dict, max_len):
    for amino in amino_pos_dict:
        print_str = f'{amino}:'
        for i in range(max_len):
            if i in amino_pos_dict[amino]:
                print_str += f' {amino_pos_dict[amino][i]}'
            else:
                print_str += ' 0'
        print(print_str)
            
        



if __name__ == '__main__':
    # Get fasta filepath
    fasta_filepath = sys.argv[1]

    amino_pos_dict, max_len = create_pos_weight_dict(fasta_filepath)

    print_pos_weight_matrix(amino_pos_dict, max_len)