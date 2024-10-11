SEQ_1 = 'ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT'
SEQ_2 = 'ATGGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAA'
TEST = 'GCA'
EXONS = [[0,62],[92,None]]
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

def get_exon_info(seq, exon_positions):
    comp_seq = str(seq).lower()

    exons = []
    for pos in exon_positions:
        exon = seq[pos[0]:pos[1]]
        comp_seq = comp_seq[0:pos[0]] + exon + comp_seq[pos[1]:]
        exons.append(exon)

    mRNA = ''.join(exons).replace('T','U')
    coding_pc = len(mRNA)*100 / len(seq)
    return mRNA, coding_pc, comp_seq
    


#print(get_complement(SEQ_1))
#print(SEQ_2)
if __name__ == "__main__":
    mRNA, coding_pc, comp_seq = get_exon_info(SEQ_2, EXONS)
    print(f'mRNA: {mRNA}')
    print(f'%coding: {coding_pc}')
    print(f'full sequence: {comp_seq}')

