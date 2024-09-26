SEQ_1 = 'ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT'
SEQ_2 = 'ATGGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAA'
TEST = 'GCA'
EXONS = [[0,62],[92,'']]
COMP_DICT = {'A': 't',
             'T': 'a',
             'G': 'c',
             'C': 'g'}

def find_gc_content(seq):
    seq = str(seq).upper()
    gc_count = 0
    for base in seq:
        if base in ['G', 'C']:
            gc_count += 1
    
    gc_pcent = gc_count / len(seq)

    return gc_pcent

def get_complement(seq):
    seq = str(seq).upper()
    for key in COMP_DICT:
        seq = seq.replace(key, COMP_DICT[key])
    
    return seq.upper()

def print_exon_info(seq):
    seq = str(seq).lower()


#print(get_complement(SEQ_1))
print(SEQ_2)