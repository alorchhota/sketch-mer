import readutil

workdir = '/Users/ashis/Desktop/work/github/sketch-mer'
k = 22
#datasets = ['rymv', 'tmv', 'saureus', 'ecoli', 'dmelanogaster']
datasets = ['hbv.sim', 'rymv.sim', 'hpylori.sim', 'hiv1.sim', 'tmv.sim', 'ecoli.sim', 'saureus.sim']
#batchSizes = [1,10,100,1000,10000,100000]
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
outFileName = workdir + '/results/accuracies.txt'


def numSolidEstimates(estimated, expected):
    isSolid = [estimated[i]==expected[i] for i in range(len(estimated))]
    return sum(isSolid)

def solidError(estimated, expected):
    l = len(estimated) + 0.0
    return (l-numSolidEstimates(estimated, expected))/l

def averageError(estimated, expected):
    errSqr = [pow((estimated[i]-expected[i])/(expected[i]+0.0), 2) for i in range(len(estimated))]
    avgErr = pow( (sum(errSqr)+0.0) / len(estimated), 0.5 )
    return avgErr

def averageError2(estimated, expected):
    errAbs = [abs((estimated[i]-expected[i])/(expected[i]+0.0)) for i in range(len(estimated))]
    avgErr = (sum(errAbs)+0.0) / len(estimated)
    return avgErr

errors = []
for dataset in datasets:
    for batchSize in batchSizes:
        if batchSize > numKmers[dataset]:
            continue
        trueCountsFileName = workdir + '/data/test/test_' + dataset + '_' + str(batchSize) + '.txt'
        cmCountsFileName = workdir + '/results/cmresult_' + dataset + '_' + str(batchSize) + '.txt'
        lsCountsFileName = workdir + '/results/lsresult_' + dataset + '_' + str(batchSize) + '.txt'
        # read files
        trueKmerCounts = list(readutil.read_kmer_counts_line_by_line(trueCountsFileName))
        cmKmerCounts = list(readutil.read_kmer_counts_line_by_line(cmCountsFileName))
        lsKmerCounts = list(readutil.read_kmer_counts_line_by_line(lsCountsFileName))
        # keep only count values.
        # assumption: counts are stored in the same order.
        trueCounts = [kmcount[1] for kmcount in trueKmerCounts]
        cmCounts = [kmcount[1] for kmcount in cmKmerCounts]
        lsCounts = [kmcount[1] for kmcount in lsKmerCounts]
        # calculate solid error
        cmSolidError = solidError(cmCounts, trueCounts)
        lsSolidError = solidError(lsCounts, trueCounts)
        # calculate average error
        cmAvgError = averageError(cmCounts, trueCounts)
        lsAvgError = averageError(lsCounts, trueCounts)
        # calculate average error 2 (no square)
        cmAvgError2 = averageError2(cmCounts, trueCounts)
        lsAvgError2 = averageError2(lsCounts, trueCounts)
        # save error
        errors.append((dataset, batchSize, cmSolidError, lsSolidError, cmAvgError, lsAvgError, cmAvgError2, lsAvgError2))

# write errors in a file
with open(outFileName, 'w') as of:
    header = '\t'.join(['dataset',
                        'batchsize',
                        'countmin_solid_err',
                        'lsquare_solid_err',
                        'countmin_avg_err',
                        'lsquare_avg_err',
                        'countmin_avg_err2',
                        'lsquare_avg_err2'])
    text = '\n'.join(['\t'.join([str(item) for item in err]) for err in errors])
    of.write(header + '\n' + text)
print('done')

