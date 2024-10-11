import traceback
import csv

#from Lab1 import find_gc_content

WORKDIR = '/home1/bioinfo-03/repos/2024_Prog_DB/Data_lab3'
SEQ_FILE_PATH = f'{WORKDIR}/sequences_lecture4.txt'
FIKK_FILEPATH = f'{WORKDIR}/FIKK_kinases.txt'
DROS_FILEPATH = f'{WORKDIR}/drosophila.csv'

CDN_LEN = 3
HEADERS = ['Species','Seq','ID','Expression']

def find_seq_content(seq, bases):
    seq = str(seq).upper()
    count = 0
    for base in seq:
        if base in bases:
            count += 1
    
    pcent = count / len(seq)

    return pcent

def find_gc_content(seq):

    return find_seq_content(seq,['G','C'])

def find_at_content(seq):

    return find_seq_content(seq,['A','T'])

def get_seqs_from_file(filepath):
    seqs = []
    with open(filepath) as seq_file:
        seq_lines = seq_file.readlines()
        for line in seq_lines:
            seqs.append(line.strip())
    return seqs

def print_gc_report(seqs, gc_hi = 0.6, gc_low = 0.4):
    for i in range(len(seqs)):
        gc_content = find_gc_content(seqs[i])
        gc_lvl = ''
        if gc_content < gc_low:
            gc_lvl = 'LOW'
        elif gc_content > gc_hi:
            gc_lvl = 'HIGH'
        else:
            gc_lvl = 'MEDIUM'
        
        print(f'Sequence {i} has {gc_lvl} GC content')

# returns list of codons
def get_codons(seq):
    seq = seq.upper()
    codons = []

    for i in range(0, len(seq), 3):
        codon = seq[i:i+CDN_LEN]
        if len(codon) == CDN_LEN:
            codons.append(codon)
    
    return codons

def test_codon_function():
    try:
        assert len(get_codons('AAA')) == 1
        assert len(get_codons('AAATTTGGG')) == 3
        assert get_codons('aaa') == ['AAA']
        assert len(get_codons('AAATTTGGGC')) == 3
    except Exception:
        print('Codons function failed test:')
        print(traceback.format_exc())

def get_seq_dict(filepath):
    seq_dict = {}
    with open(filepath) as seq_file:
        seq_lines = seq_file.readlines()
        for line in seq_lines:
            seq_info = line.rstrip().split('\t')
            try:
                seq_dict[seq_info[0]] = seq_info[1]
            except:
                print(f'Invalid sequence: {seq_info}')
    return seq_dict

def print_codons(seq_dict):
    for seq_id in seq_dict:
        print(f'{seq_id} Codons:')
        print(get_codons(seq_dict[seq_id]))
        print('###\n\n')

def get_genes_data(filepath, headers):
    genes_dict = {}
    id_pos = headers.index('ID')
    with open(filepath) as csvfile:
        genes = csv.reader(csvfile)
        for gene_info in genes:
            gene_dict = {}
            for i in range(len(headers)):
                if i != id_pos:
                    gene_dict[headers[i]] = gene_info[i]

            genes_dict[gene_info[id_pos]] = gene_dict
    
    return genes_dict

def get_genes_by_species(gene_data, species, exclude = False):
    return_genes = {}
    for gene in gene_data:
        if gene_data[gene]['Species'] == species:
            if not exclude:
                return_genes[gene] = gene_data[gene]
                

if __name__ == '__main__':
    # PROBLEM 1
    seqs = get_seqs_from_file(SEQ_FILE_PATH)
    #print_gc_report(seqs)

    # PROBLEM 2
    test_codon_function()
    fikk_kinases = get_seq_dict(FIKK_FILEPATH)
    print_codons(fikk_kinases)

    #PROBLEM 3
    gene_data = get_genes_data(DROS_FILEPATH, HEADERS)
    