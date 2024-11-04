import argparse

# Calculate GC (and optionally AT content from a DNA sequence to predict if a sequence originates from Plasmodium.

def calculateGC(seq):
    """
    Calculate GC content
    :param seq: String containing DNA sequence
    :return: Float GC content expressed as a percentage
    """
    seqLength = len(seq)
    c_count = seq.upper().count('C')
    g_count = seq.upper().count('G')
    gc_content = ((c_count + g_count) / seqLength) * 100
    return gc_content

parser = argparse.ArgumentParser(prog='GC calculater', description='Calculate GC (and optionally AT content from a DNA sequence to predict if a sequence originates from Plasmodium.')
parser.add_argument('--filename',required=True,help='Filepath to the sequnce')
parser.add_argument('--showAT',help='Show AT content as well as GC',action='store_true')
parser.add_argument('--cutoff',type=int,help='Cutoff below which to label the sequence Plasmodium',default=30)

args = parser.parse_args()
seqFile = args.filename
showAT = args.showAT
cutoff = args.cutoff

count = 0
# for each line in the file
with open(seqFile) as sequences:
    for seq in sequences:
        count += 1
        # calculate GC content and create output string
        gcContent = calculateGC(seq)
        output = 'Sequence {}: GC content {}%\n'.format(count, round(gcContent, 2))
        # if requested, calculate AT content and add to output string
        if showAT:
            atContent = 100 - gcContent
            output = output + 'Sequence {}: AT content {}%\n'.format(count, round(atContent, 2))
        # if gcContent is less than the cutoff, report that this is likely Plasmodium
        if gcContent < cutoff:
            output = output + 'Sequence {} is likely Plasmodium\n'.format(count)
        print(output)

# print the docstring for the function
# print(calculateGC.__doc__)
