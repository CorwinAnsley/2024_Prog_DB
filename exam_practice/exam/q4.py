
VALID_BASES = ['A','T','G','C']

# Takes no arguments, prompts the user for a sequence and makes sure it is valid
# If not the user is prompted again until they provide a valid sequence
def get_user_seq():
    while True:
        try:
            user_seq = input('Please input your sequence: ')
            user_seq = user_seq.rstrip().lstrip().upper()

            if len(user_seq) == 0:
                raise Exception('No sequence inputted')

            for base in user_seq:
                if base not in VALID_BASES:
                    raise Exception('Invalid base in sequence')

            break

        except Exception as err:
            print(f'{err}, try again')
    
    return user_seq

# Takes no arguments, prompts the user for a value of k and checks that it is valid
def get_k():
    while True:
        try:
            user_k = input('Please input k to find k-mers: ')
            user_k = user_k.rstrip().lstrip()
            if len(user_k) == 0:
                raise Exception('No k inputted')
            
            try:
                user_k = int(user_k)
            except:
                raise Exception('Unable to convert input to integer')
        except Exception as err:
            print(f'{err}, try again')
    
    return user_k

# With a given k and sequence as arguments get all k-mers and store in a dict
# Returns the dict with kmers
def get_kmers(seq, k):
    kmer_dict = {}
    for i in range(len(seq)):
        kmer = seq[i:i+k]
        if len(kmer) == k:
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    
    return kmer_dict

# Takes k-mer dict, sequence, and k as arguments and outputs a report for the user
def print_kmers(kmer_dict, seq, k):
    print(f'Reporting {k}-mers for sequence: {seq}\n')

if __name__ == '__main__':
    seq = get_user_seq()
    k = get_k()
    kmer_dict = get_kmers(seq, k)