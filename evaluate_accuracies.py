import readutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as p

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

def meanErr(estimated, expected):
    err = [abs(estimated[i]-expected[i]) for i in range(len(estimated))]
    stdErr = np.mean(err)
    return stdErr

def errorPlot(cmEst, lsEst, expected, plotFileName, nbins=100):
    # cmEst error
    err1 = [abs(cmEst[i]-expected[i]) for i in range(len(cmEst))]
    #y,binEdges=np.histogram(err, bins=nbins)
    #bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    #p.plot(bincenters,y,'-', label='bg', color='blue')

    # lsEst error
    err = [abs(lsEst[i]-expected[i]) for i in range(len(lsEst))]
    #y,binEdges=np.histogram(err, bins=nbins)
    #bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    #p.plot(bincenters,y,'-', label='bg', color='red')
    xmin = min(min(err1), min(err))
    xmax = max(max(err1), max(err))
    bins = np.linspace(xmin, xmax, nbins)
    plt.hist(err, bins, alpha=0.5, label='sketch-mer')
    plt.hist(err1, bins, alpha=0.5, label='countmin')
    # save
    plt.legend(loc='upper right')
    plt.savefig(plotFileName)
    plt.close()
>>>>>>> ba336893dfd46556b9269918f8ce07eb95a182fc


errors = []
for dataset in datasets:
    for batchSize in batchSizes:
        if batchSize > numKmers[dataset]:
            continue
        trueCountsFileName = workdir + '/data/test/test_' + dataset + '_' + str(batchSize) + '.txt'
        cmCountsFileName = workdir + '/results/cmresult_' + dataset + '_' + str(batchSize) + '.txt'
        lsCountsFileName = workdir + '/results/lsresult_' + dataset + '_' + str(batchSize) + '.txt'
        plotFileName = workdir + '/results/plot_' + dataset + '_' + str(batchSize) + '.png'
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
        # calculate standar deviation of error
        cmMeanError = meanErr (cmCounts, trueCounts)
        lsMeanError = meanErr(lsCounts, trueCounts)
        # save error
        errors.append((dataset, batchSize, cmSolidError, lsSolidError, cmAvgError, lsAvgError, cmAvgError2, lsAvgError2, cmMeanError, lsMeanError))
        # plot error graph
        errorPlot(cmCounts, lsCounts, trueCounts, plotFileName)

# write errors in a file
with open(outFileName, 'w') as of:
    header = '\t'.join(['dataset',
                        'batchsize',
                        'countmin_solid_err',
                        'lsquare_solid_err',
                        'countmin_avg_err',
                        'lsquare_avg_err',
                        'countmin_avg_err2',
                        'lsquare_avg_err2',
                        'countmin_mean_err',
                        'lsquare_mean_err'])
    text = '\n'.join(['\t'.join([str(item) for item in err]) for err in errors])
    of.write(header + '\n' + text)
print('done')

