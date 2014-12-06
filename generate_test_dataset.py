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


#workdir = '/Users/ashis/Desktop/work/github/sketch-mer'
workdir = '.'
k = 22
#datasets = ['rymv', 'tmv', 'saureus', 'ecoli', 'dmelanogaster']
datasets = ['hbv.sim', 'rymv.sim', 'hpylori.sim', 'hiv1.sim', 'tmv.sim', 'ecoli.sim', 'saureus.sim']
batchSizes = [1,10,100, 200, 400, 600, 800, 1000]
numKmers = {'rymv': 199,
            'tmv': 6374,
            'saureus': 2781271,
            'ecoli': 4563650,
            'dmelanogaster': 122701066,
            'hbv.sim': 7600,
            'rymv.sim': 3361,
            'hpylori.sim': 1502535,
            'hiv1.sim': 25561,
            'tmv.sim': 20534,
            'ecoli.sim': 14681,
            'saureus.sim': 7721408}

for dataset in datasets:
    for batchSize in batchSizes:
        if batchSize > numKmers[dataset]:
            continue
        print(dataset + " " + str(batchSize) + ' ...')
        ofn = workdir + '/results/test_'+dataset+ '_' + str(batchSize) + '.txt'
        kmerfn = workdir + '/data/22mer_exact_counts/' + dataset + '_counts.txt'
        generate_test_data(batchSize, numKmers[dataset], kmerfn, ofn)

print('done! see output in results folder.')