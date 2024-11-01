from Lab4 import create_codon_table, translate_to_protein, get_genes_from_fasta, CODON_TABLE

WORKDIR1 = '/home1/bioinfo-03/repos/2024_Prog_DB/data_files/Lab4'
TOXO_FASTA_FILE = f'{WORKDIR1}/Toxo_Ser_Thr_kinase.fasta'

WORKDIR = '/home1/bioinfo-03/repos/2024_Prog_DB/data_files/Lab6'
OUTPUT = f'{WORKDIR}/Test.fasta'

HYDROPHOBIC_AMINOS = ['A', 'I', 'L', 'M','F', 'W', 'Y', 'V']

# create the class
class Sequence:

    # constructor
    def __init__(self, id, seq):
        self._id = id
        self._seq = seq.upper()

    # getter and setter for id
    @property
    def id(self):
        return self._id 
    
    # the setter is here as an example
    """
    @id.setter
    def id(self, id):
        self._id = id
    """

    # getter and setter for seq
    @property
    def seq(self):
        return self._seq
    
    # the setter is here as an example
    """
    @seq.setter
    def seq(self, seq):
        self._seq = seq.upper()
    """

    # gc content method
    def calc_gc_content(self, dp=2):
        c_count = self.seq.count('C')
        g_count = self.seq.count('G')
        gc_content = (c_count + g_count) / len(self.seq)
        return round(gc_content, dp)

    # Translate to protein
    def get_protein(self):
        return translate_to_protein(self.seq, CODON_TABLE)
    
    # Return a string of id and sequence
    def to_string(self):
        return f'>{self.id}\n{self.seq}\n'
    
    def __str__(self):
        return self.to_string()

class ProteinSequence(Sequence):
    
    # Returns the percentage of amino acids that are hydrophobic
    def calc_hydropho_pc(self, decpoints = 2):
        hydropho_count = 0
        for amino in HYDROPHOBIC_AMINOS:
            hydropho_count += self.seq.count(amino)
        
        hydropho_pc = hydropho_count / len(self.seq)
        
        return round(hydropho_pc, decpoints)

class FastaReader:

    def get_dna_seqs(filepath):
        seqs = get_genes_from_fasta(filepath)
        for seq in seqs:
            yield Sequence(seq, seqs[seq])

class FastaWriter:

    def write_fasta(seqs, filepath):
        with open(filepath, 'w') as file:
            for seq in seqs:
                file.write(str(seq))

if __name__ == '__main__':
    # make two instances of this class
    seq1 = Sequence('geneA', 'ATGCATG')
    seq2 = Sequence('geneB', 'ATTTTAGCGAAA')

    # get the gc content for each seq instance
    print(f'The GC content for {seq1.id} is {seq1.calc_gc_content()}')
    print(f'The GC content for {seq2.id} is {seq2.calc_gc_content(3)}')

    print(f'The protein for {seq1.id} is {seq1.get_protein()}')
    print(seq1.to_string())
    
    prot_seq1 = ProteinSequence('proteinA', seq1.get_protein())
    print(f'fasta of protein: \n{prot_seq1}')
    print(f'The hydrophobic content for {prot_seq1.id} is {prot_seq1.calc_hydropho_pc()}')

    seqs = FastaReader.get_dna_seqs(TOXO_FASTA_FILE)

    proteins = []

    for seq in seqs:
        proteins.append(ProteinSequence(seq.id, seq.get_protein()))

    FastaWriter.write_fasta(proteins, OUTPUT)
    #seq1.seq = seq1.seq.lower()
