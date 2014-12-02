import random
import readutil

def generate_test_data(batchSize, numKmer, kmerFileName, outFileName):
    kmerIndexes = [random.randrange(0, numKmer) for _ in range(batchSize)]
    kmerIndexes = sorted(kmerIndexes)
    kmers = readutil.read_kmer_counts_line_by_line(kmerFileName)
    curTargetIndex = 0
    target = kmerIndexes[curTargetIndex]
    curIndex = 0
    savedKmers = []
    with open(outFileName, 'w') as of:
        while curIndex < numKmer and curTargetIndex < batchSize:
            kmer = next(kmers)
            while curIndex==target:
                of.write('\t'.join([kmer[0], str(kmer[1])])+'\n')
                savedKmers.append(kmer)
                curTargetIndex += 1
                if curTargetIndex >=  batchSize:
                    break
                target = kmerIndexes[curTargetIndex]
            curIndex += 1


workdir = '/home/ashis/work/github/sketch-mer'
k = 22
datasets = ['rymv', 'tmv', 'saureus', 'ecoli', 'dmelanogaster']
batchSizes = [1,10,100,1000,10000,100000]
datasets = datasets[2:3]
batchSizes = batchSizes[1:2]
numKmers = {'rymv': 199,
            'tmv': 6374,
            'saureus': 2781271,
            'ecoli': 4563650,
            'dmelanogaster': 122701066}

for dataset in datasets:
    for batchSize in batchSizes:
        if batchSize > numKmers[dataset]:
            continue
        ofn = workdir + '/data/test/test_'+dataset+ '_' + str(batchSize) + '.txt'
        kmerfn = workdir + '/data/22mer_exact_counts/' + dataset + '_counts.txt'
        generate_test_data(batchSize, numKmers[dataset], kmerfn, ofn)
