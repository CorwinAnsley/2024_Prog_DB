import traceback
import csv

WORKDIR = '/home1/bioinfo-03/repos/2024_Prog_DB/data_files/Lab4'
TOXO_FASTA_FILE = f'{WORKDIR}/Toxo_Ser_Thr_kinase.fasta'
PEPTIDE_SIGNAL_FILE = f'{WORKDIR}/GenesWithSignalPeptide_SignalP.txt'
OUTPUT_PROTEINS = f'{WORKDIR}/Toxo_Ser_Thr_kinase_proteins.fasta'
OUTPUT_PEPTIDES = f'{WORKDIR}/Toxo_Ser_Thr_kinase_peptides.csv'

BASES = "tcag".upper()
AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'

def create_codon_table(bases, amino_acids):
    # make a codon table)
    codons = [a + b + c for a in bases for b in bases for c in bases]
    codon_table = dict(zip(codons, amino_acids))
    return codon_table

CODON_TABLE = create_codon_table(BASES, AMINO_ACIDS)

#creates a dict with a sequence for each gene
def get_genes_from_fasta(filepath):
    try:
        with open(filepath) as fasta:
            fasta_lines = fasta.readlines()
            genes_dict = {}
            current_gene = ''
            seq = ''
            for line in fasta_lines:
                if line.startswith('>'):
                    if current_gene != '':
                        genes_dict[current_gene] = seq
                    current_gene = line.lstrip('>').strip()
                else:
                    seq += line.strip().upper()
            if current_gene != '':
                genes_dict[current_gene] = seq
            return genes_dict
    except:
        print('failed to open fasta file:')
        print(traceback.format_exc())


def get_proteins(gene_dict, output_path, codon_table):
    protein_dict = {}

    for gene in gene_dict:
        protein_dict[gene] = translate_to_protein(gene_dict[gene], codon_table)
    
    with open(output_path, 'w') as fasta:
         for protein in protein_dict:
             fasta.writelines(f'>{protein}\n')
             fasta.writelines(f'{protein_dict[protein]}\n')
    
    return protein_dict
        
def translate_to_protein(seq, codon_table):
    protein = ''
    for i in range(0,len(seq),3):
        if i+2 <= len(seq):
            codon = seq[i:i+3]
            amino_a = codon_table[codon]
            protein += amino_a
    
    return protein

def get_peptide_signal_dict(filepath):
    peptide_signal_dict= {}
    try:
        with open(filepath) as csvfile:
            signal_reader = csv.reader(csvfile, delimiter='\t')
            header = next(signal_reader)
            
            for peptide_signal in signal_reader:
                data_dict = {}
                for i in range(1, len(peptide_signal)):
                    data_dict[header[i]] = peptide_signal[i]

                peptide_signal_dict[peptide_signal[0]] = data_dict
    except:
        print('failed to open peptide_signal file:')
        print(traceback.format_exc())

    return peptide_signal_dict    
        
def get_peptide_proteins(protein_dict, peptide_signal_dict, zero_indexed = False):
    peptide_protein_dict = {}    
    for protein in protein_dict:
        if protein in peptide_signal_dict:
            if zero_indexed: #everything should be zero indexed but this world is not just
                start = int(peptide_signal_dict[protein]['Start Min'])
            else:
                start = int(peptide_signal_dict[protein]['Start Min']) - 1

            end = int(peptide_signal_dict[protein]['End Max'])
            peptide_protein = protein_dict[protein][start:end]

            data_dict = {}
            data_dict['protein'] = peptide_protein
            for key in peptide_signal_dict[protein]:
                data_dict[key] = peptide_signal_dict[protein][key]

            peptide_protein_dict[protein] = data_dict
    return peptide_protein_dict

def write_peptide_proteins(peptide_protein_dict, output_path):
     with open(output_path, 'w') as file:
        headers = ['Gene ID']
        headers = headers + list(peptide_protein_dict[list(peptide_protein_dict.keys())[0]].keys())
        csv_writer = csv.writer(file, delimiter='\t') #bleh tabs
        csv_writer.writerow(headers)
        for peptide_gene in peptide_protein_dict:
            row = [peptide_gene]
            for header in headers:
                if header in peptide_protein_dict[peptide_gene]:
                    row.append(peptide_protein_dict[peptide_gene][header])
            csv_writer.writerow(row) 
            
if __name__ == '__main__':

    codon_table = create_codon_table(BASES, AMINO_ACIDS)
    gene_dict = get_genes_from_fasta(TOXO_FASTA_FILE)
    protein_dict = get_proteins(gene_dict, OUTPUT_PROTEINS, codon_table)
    peptide_signal_dict = get_peptide_signal_dict(PEPTIDE_SIGNAL_FILE)
    peptide_protein_dict = get_peptide_proteins(protein_dict, peptide_signal_dict)
    write_peptide_proteins(peptide_protein_dict, OUTPUT_PEPTIDES)
    print('a')