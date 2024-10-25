import re

from Lab4 import get_genes_from_fasta, get_peptide_signal_dict, write_peptide_proteins

WORKDIR = '/home1/bioinfo-03/repos/2024_Prog_DB/data_files/Lab5'
PROTEINS_FASTA_FILE = f'{WORKDIR}/Tgondii_Proteins.fasta'
PRODUCT_DESCRIPTIONS_FILE = f'{WORKDIR}/Tgondii_product_descriptions.txt'
OUTPUT_MOTIF_PROTIENS = f'{WORKDIR}/Motif_proteins.csv'

YXXPHI_MOTIF_REGEX = r'Y..[YFT]{2}.{,5}$'
MOTIF_LEN = 5

#
def test_regex(regex):
    print(re.search(regex, 'AAAAAAAAAAAAAAAAAAAYXXFTACTACT'))

#wrapper just so this has a better name
def load_proteins_fasta(filepath):
    return get_genes_from_fasta(filepath)

def load_descr_dict(filepath):
    return get_peptide_signal_dict(filepath)

def output_dict_csv(dict, output_path):
    write_peptide_proteins(dict, output_path)

def find_motifs(protein_dict, motif_regex):
    motif_protein_dict = {}
    for protein in protein_dict:
        motif_match = re.search(motif_regex,protein_dict[protein])
        if motif_match:
            motif_protein_dict[protein] = {'protein': protein_dict[protein], 'motif_match': motif_match}
    
    return motif_protein_dict

def create_motif_info_dict(protein_dict, motif_regex, descr_file):
    motif_protein_dict = find_motifs(protein_dict, motif_regex)
    descr_dict = load_descr_dict(descr_file)
    motif_info_dict = {}
    for protein in motif_protein_dict:
        assert protein in descr_dict
        motif_info_dict[protein] = {}
        motif_info_dict[protein]['Product Description'] = descr_dict[protein]['Product Description']
        seq_tail = motif_protein_dict[protein]['protein']#[-30:]
        match_start = motif_protein_dict[protein]['motif_match'].start()
        motif_end = seq_tail[match_start:]
        motif = motif_end[0:5].lower()
        motif_end = f'{motif}{motif_end[5:]}'
        seq_tail = f'{seq_tail[:match_start]}{motif_end}'
        seq_tail = seq_tail[-30:]
        motif_info_dict[protein]['Sequence Tail'] = seq_tail
    return motif_info_dict

if __name__ == '__main__':
    #descr_dict = load_descr_dict(PRODUCT_DESCRIPTIONS_FILE)
    protein_dict = load_proteins_fasta(PROTEINS_FASTA_FILE)
    #motif_protein_dict = find_motifs(protein_dict, YXXPHI_MOTIF_REGEX)
    motif_info_dict = create_motif_info_dict(protein_dict, YXXPHI_MOTIF_REGEX, PRODUCT_DESCRIPTIONS_FILE)
    output_dict_csv(motif_info_dict, OUTPUT_MOTIF_PROTIENS)
    test_regex(YXXPHI_MOTIF_REGEX)