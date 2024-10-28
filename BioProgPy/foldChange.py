import math
import logging
import traceback

DATADIR = '/home1/bioinfo-03/repos/2024_Prog_DB/BioProgPy/data/Lab1'

logger = logging.getLogger()
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('‘%(levelname)s - %(asctime)s - %(message)s’)'))
logger.addHandler(handler)

logger.info('Hello World\n')

def calculatePercentile(data, percentile):
    """
    calculate the chosen percentile of a list of numbers
    :param data: list of numbers (list)
    :param percentile: percentile to calculate (int)
    :return: percentile (int)
    """
    n = len(data)
    p = n * percentile / 100
    if p.is_integer():
        return sorted(data)[int(p)]
    else:
        return sorted(data)[int(math.ceil(p)) - 1]


# lists collect numeric columns to calculate percentiles
uninducedCounts = []
inducedCounts = []

# dict collects all data
countDict = {}

file = f'{DATADIR}/countData.txt'

# open the count file
try:
    with open(file) as counts:
        # skip the header row
        next(counts)
        # split line by tab and assign components
        for line in counts:
            try:
                line.rstrip()
                gene, uninduced, induced = line.split('\t')
                try:
                    if 'RNA' in gene:
                        continue
                    # append counts to lists and add data to dict
                    else:
                        uninducedCounts.append(int(uninduced))
                        inducedCounts.append(int(induced))
                        if gene in countDict:
                            raise SystemExit('Error: Duplicate for gene {} found. Please check your input\n'.format(gene))
                        else:
                            countDict[gene] = {'Uninduced': int(uninduced), 'Induced': int(induced)}
                except:
                    logger.error(f'Encountered error parsing gene: {gene}\n{traceback.format_exc()}')
                    continue
            except:
                logger.error(f'Encountered error reading line: {line}, :\n{traceback.format_exc()}')
                continue
    counts.close()
except FileNotFoundError:
    logger.error(f'File: {file} not found')
    raise SystemExit(1)

# calculate upper quartile (75th percentile) for uninduced and induced samples
uninducedUpperQuartile = calculatePercentile(uninducedCounts, 75)
inducedUpperQuartile = calculatePercentile(inducedCounts, 75)


# normalise counts using quantile values, calculate log2(foldchange) and write to file
outputFile = 'foldChange.txt'
with open(outputFile, 'wt') as output:
    output.write('GeneID\tUninducedNormalised\tInducedNormalised\tLog2FoldChange\n')
    for gene, counts in countDict.items():
        # scale the data using upper quartiles
        uninducedNormalised = (counts['Uninduced'] / uninducedUpperQuartile) * 1000
        inducedNormalised = (counts['Induced'] / inducedUpperQuartile) * 1000
        # check induced normalised reads count
        try:
            assert inducedNormalised > 0
        except:
            logger.warning(f'0 induced normalised reads for gene {gene}, setting to 0.0001 to avoid errors')
            inducedNormalised = 0.0001
        # calculate fold change
        try:
            foldChange = inducedNormalised / uninducedNormalised
        except ZeroDivisionError:
            logger.warning(f'0 uninduced normalised reads for gene {gene}, divinding by 0.0001 and moving on')
            foldChange = inducedNormalised / 0.0001
        # calculate log2 of fold change
        try:
            foldChange = math.log2(foldChange)
        except:
            logger.error(f'Error calculating foldChange for gene {gene}:\n{traceback.format_exc()}')
        # write to file
        output.write('{}\t{}\t{}\t{}\n'.format(gene, round(uninducedNormalised, 3), round(inducedNormalised, 3), round(foldChange, 3)))
output.close()
