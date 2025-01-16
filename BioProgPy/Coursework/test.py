test = '1234567890'
test2 = '90X'

# Returns the codon, and alternative codon of a sequence with SNP
# Input sequence should already be adjusted to start at the first codon and continue until at least the last codon the snp could be located at
def get_alt_codon(seq, pos, alt):
    # Get the sequence up until the SNP
    alt_codon = seq[:pos]

    # Get rid of everything before the codon
    seq_mod_3 = len(alt_codon) % 3
    if seq_mod_3 > 0:
        alt_codon = alt_codon[-seq_mod_3:]  
    else:
        alt_codon = ''
    
    # Add the base at the site of the SNP, we need 2 codons for the reference and the alternative
    codon = alt_codon
    codon += seq[pos]
    alt_codon += alt

    # Add the rest of the codon if missing
    seq_len = len(alt_codon)
    rest_of_codon_seq = seq[pos+1:pos+4-seq_len]
    codon += rest_of_codon_seq
    alt_codon += rest_of_codon_seq

    return codon, alt_codon

# for i in range(len(test)-20):
#     print(get_alt_codon(test,i,'a',seq2=test2))
print(print(get_alt_codon(test,-3,'a')))    