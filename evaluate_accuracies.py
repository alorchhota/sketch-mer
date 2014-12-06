import readutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as p

#workdir = '/Users/ashis/Desktop/work/github/sketch-mer'
workdir = '.'
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

def errorHist(cmEst, lsEst, expected, plotFileName, nbins=100):
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
    plt.hist(err, bins, alpha=0.5, label='Sketch-mer')
    plt.hist(err1, bins, alpha=0.5, label='CountMinSketch')
    plt.xlabel('Error')
    plt.ylabel('Frequency')
    # save
    plt.legend(loc='upper right')
    plt.savefig(plotFileName)
    plt.close()

def avgErrBarChart(cmMeanErr, lsMeanErr, plotFileName):
    N = 2
    ind = np.arange(N)  # the x locations for the groups
    width = 0.05       # the width of the bars

    fig, ax = plt.subplots()
    rects =  []
    colors = ['red','green', 'blue', 'yellow', 'orange', 'black', 'gray', 'olive']
    for i in range(len(cmMeanErr)):
        errs = (cmMeanErr[i], lsMeanErr[i])
        rects.append(ax.bar(ind+i*width-.15, errs, width, color=colors[i]))

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Average Error')
    ax.set_xticks(ind+width)
    ax.set_xticklabels( ('countminsketch', 'sketch-mer'))


    ax.legend( (rects[0][0], rects[1][0], rects[2][0], rects[3][0],
                rects[4][0], rects[5][0], rects[6][0], rects[7][0]), ('1','10','100', '200', '400', '600', '800', '1000') )

    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')
    #for i in range(len(cmMeanErr)):
    #    autolabel(rects[i])

    plt.show()
    plt.savefig(plotFileName)
    plt.close()

def errComparisonBarChart(errors, batchSize, plotFileName, ylim):
    # (dataset, cmerr, lserr)
    #batchErrors = [err for err in errors if err[1]==batchSize and err[0] != 'saureus.sim' and err[0]!='hpylori.sim']
    batchErrors = [err for err in errors if err[1]==batchSize]
    batchErrors = sorted(batchErrors, key=lambda x:x[8])
    datasets = [err[0][:-4] for err in batchErrors]
    cmerr = [err[8] for err in batchErrors]
    lserr = [ err[9] for err in batchErrors]

    N = len(batchErrors)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.25       # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, cmerr, width, color='gray')
    rects2 = ax.bar(ind+width, lserr, width, color='orange')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Average Error')
    ax.set_ylim((0,ylim))
    ax.set_xticks(ind+width)
    ax.set_xticklabels( datasets)

    ax.legend( (rects1[0], rects2[0]), ('CountMinSketch', 'Sketch-mer') )

    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')
    autolabel(rects1)
    autolabel(rects2)

    plt.show()
    plt.savefig(plotFileName)
    plt.close()


def errorVsNumKmer(errors, batchSize, plotFileName, ylim):
    batchErrors = [err for err in errors if err[1]==batchSize]
    batchErrors = sorted(batchErrors, key=lambda x:numKmers[x[0]])
    datasets = [err[0] for err in batchErrors]
    cmerr = [err[8] for err in batchErrors]
    lserr = [ err[9] for err in batchErrors]
    nkm = [np.log10(numKmers[ds]) for ds in datasets]
    # plot
    p.plot(nkm,cmerr,'-', label='CountMinSketch', color='gray')
    p.plot(nkm,lserr,'-', label='Sketch-mer', color='orange')
    plt.xlabel('log(# distinct k-mers)')
    plt.ylabel('Average Error')
    plt.ylim((0,ylim))
    # save
    plt.legend(loc='upper right')
    plt.savefig(plotFileName)
    plt.close()


def hashSizeVsError(plotFileName):
    eps = [.0001, .001, .003, .005]
    w = [np.log10(int(np.ceil(np.exp(1) / e))) for e in eps]
    dataset = 'haureus'
    cmerr_hau = [734.665,7842.57,23795.075,39827.015]
    lserr_hau = [11.775,31.915,51.88,62.255]

    cmerr_hiv = [0.06,14.665,57.245,105.7]
    lserr_hiv = [0.445,2.73,3.88,4.685]

    cmerr_tmv = [0.015,12.98,57.73,104.715]
    lserr_tmv = [0.5,3.135,5.37,5.705]

    # plot
    plt.plot(w,cmerr_hiv,'-o', label='CountMinSketch (hiv)', color='gray')
    plt.plot(w,lserr_hiv,'-o', label='Sketch-mer (hiv)', color='orange')
    plt.plot(w,cmerr_tmv,'-o', label='CountMinSketch (tmv)', color='olive')
    plt.plot(w,lserr_tmv,'-o', label='Sketch-mer (tmv)', color='tomato')
    plt.xlabel('log (Hash Size)')
    plt.ylabel('Average Error')
    #plt.ylim((0,15000))
    # save
    plt.legend(loc='upper right')
    plt.savefig(plotFileName)
    plt.close()


errors = []
for dataset in datasets:
    barFileName = workdir + '/results/bar_' + dataset + '.png'
    for batchSize in batchSizes:
        if batchSize > numKmers[dataset]:
            continue
        trueCountsFileName = workdir + '/data/test/test_' + dataset + '_' + str(batchSize) + '.txt'
        cmCountsFileName = workdir + '/results/cmresult_' + dataset + '_' + str(batchSize) + '.txt'
        lsCountsFileName = workdir + '/results/lsresult_' + dataset + '_' + str(batchSize) + '.txt'
        histFileName = workdir + '/results/hist_' + dataset + '_' + str(batchSize) + '.png'
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
        errorHist(cmCounts, lsCounts, trueCounts, histFileName)
    nbat = len(batchSizes)
    dbErrors = errors[-nbat:]
    cmErr = [item[8] for item in dbErrors]
    lsErr = [item[9] for item in dbErrors]
    avgErrBarChart(cmErr, lsErr, barFileName)

for batchSize in batchSizes:
    errCompBarPlotFileName = workdir + '/results/bar_err_comp_' + str(batchSize) + '.png'
    errComparisonBarChart(errors, batchSize, errCompBarPlotFileName, 100)
    errVsNumKmerFileName = workdir + '/results/lineErrVsNumKmer_' + str(batchSize)  + '.png'
    errorVsNumKmer(errors, batchSize, errVsNumKmerFileName, 500)

errVsHashSizeFileName = workdir + '/results/lineErrVsHashSize_200.png'
hashSizeVsError(errVsHashSizeFileName)

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
print('done! see output in results folder.')

