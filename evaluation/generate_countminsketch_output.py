import sys
from countLeastSquares import countLeastSquares
import KmerSupplier as ks
import readutil


def build_countLSsketch(ksup,delta=10**-3,epsilon = 0.005):
    sketch = countLeastSquares(delta, epsilon, k=10000)  # change the entry for k
    for kmer in ksup.iterkmers():
        sketch.update(kmer,1)
    return sketch


def generate_countmin_ls_output(sketch, testFileName, cmOutFileName, lsOutFileName):
    print('countminsketch: querying from test dataset...')
    testKmers = readutil.read_kmer_counts_line_by_line(testFileName)
    testKmers = [kmertup[0] for kmertup in testKmers]
    cmEstimates = []
    for kmer in testKmers:
        cmest = sketch.get(kmer)
        cmEstimates.append(cmest)
    cmEstimates = zip(testKmers, cmEstimates)
    with open(cmOutFileName, 'w') as outfile:
        text = '\n'.join(['\t'.join([est[0], str(est[1])]) for est in cmEstimates])
        outfile.write(text)
    #
    print('lssketch: querying from test dataset...')
    uniqueTestKmers = list(set(testKmers))
    lsEstimates = sketch.lsquare(uniqueTestKmers)
    with open(lsOutFileName, 'w') as outfile:
        text = '\n'.join(['\t'.join([km, str(lsEstimates[km])])
                          for km in testKmers])
        outfile.write(text)

workdir = '/Users/princy/Documents/classes/computational_genomics/project/sketch-mer'
k = 22
datasets = ['rymv', 'tmv', 'saureus', 'ecoli', 'dmelanogaster']
#datasets = ['saureus']
batchSizes = [1,10,100,1000,10000,100000]
delta=10**-3
epsilon = 0.005
numKmers = {'rymv': 199,
            'tmv': 6374,
            'saureus': 2781271,
            'ecoli': 4563650,
            'dmelanogaster': 122701066}



for dataset in datasets:
    print(dataset + ': building sketch...')
    genomeFileName = workdir + '/data/' + dataset + '.fa'
    ksup = ks.KmerSupplier(genomeFileName, k)
    sketch = build_countLSsketch(ksup, delta, epsilon)
    for batchSize in batchSizes:
        if batchSize > numKmers[dataset]:
            continue
        print('dataset: ' + dataset + ', ' + 'batch size: ' + str(batchSize) + " ...")
        testFileName = workdir + '/data/test/test_' + dataset + '_' + str(batchSize) + '.txt'
        cmOutFileName = workdir + '/results/cmresult_' + dataset + '_' + str(batchSize) + '.txt'
        lsOutFileName = workdir + '/results/lsresult_' + dataset + '_' + str(batchSize) + '.txt'
        generate_countmin_ls_output(sketch, testFileName, cmOutFileName, lsOutFileName)

print('see output in files: results/cmresult_***_**.txt and results/lsresult_***_**.txt')