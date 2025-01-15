test = '123456'
test2 = '789'

def get_codon(seq, pos, alt, seq2=''):

    #print(seq_segments)
    return_seq = seq[:pos] + alt
    print(return_seq)
    mod_3_r_seq = len(return_seq) % 3
    
    print(mod_3_r_seq)
    #print(return_seq[pos:pos+2])
    #return_seq = return_seq[pos-1+mod_3_r_seq:]
    print('b')

    seq_len = len(return_seq)
    #print(seq_len)
    #return_seq += seq[pos+1:pos+4-seq_len]
    seq_len = len(return_seq)
    #return_seq += seq[:4-seq_len]
    return return_seq

print(get_codon(test,2,'a'))