test = '12345678AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
test2 = '90X'

# Returns the codon, and alternative codon of a sequence with SNP 
def get_alt_codon(seq, pos, alt, seq2=''):
    alt_codon = seq[:pos]
    seq_mod_3 = len(alt_codon) % 3
    if seq_mod_3 > 0:
        alt_codon = alt_codon[-seq_mod_3:]  
    else:
        alt_codon = ''
    
    codon = alt_codon
    codon += seq[pos]
    alt_codon += alt

    seq_len = len(alt_codon)
    rest_of_codon_seq = seq[pos+1:pos+4-seq_len]
    codon += rest_of_codon_seq
    alt_codon += rest_of_codon_seq
    
    # seq_len = len(alt_codon)
    # #while
    # rest_of_codon_seq2 = seq2[:3-seq_len]
    # codon += rest_of_codon_seq2
    # alt_codon += rest_of_codon_seq2

    return codon, alt_codon

# for i in range(len(test)-20):
#     print(get_alt_codon(test,i,'a',seq2=test2))
print(test2[:2])    