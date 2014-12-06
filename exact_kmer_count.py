import sys
import KmerSupplier as ks
from collections import OrderedDict 
from fileinput import close
import shutil

def read_kmer_counts_line_by_line(fn):
    """ Read kmer count file, yielding (kmer, count) tuples one by one """
    with open(fn) as fh:
        for ln in fh:
            ln = ln.strip()
            if len(ln) == 0:
                continue
            items = ln.split('\t')
            yield (items[0], int(items[1]))

def copy_file(src, dst):
    shutil.copyfile(src, dst)

            
def copy_file_v0(fromfn, tofn):
    with open(fromfn, 'r') as fromfile:
        with open(tofn, 'w') as tofile:
            for ln in fromfile:
                tofile.write(ln)

class KmerCounter(object):
    """A class to count kmers. It should handle big genome.
    The idea is to save a maximum number of distinct kmer counts in memory.
    If the memory is full, it saves all info in files (filenames generated sequentially).
    Later, these files will be merged to get exact count."""
    def __init__(self, ksup, outdir, fileprefix='kmer_exact_count_', memory_limit=10):
        self._outdir = outdir
        self._fprefix = fileprefix
        self._ksup = ksup
        self._fcounter = 0
        self._memlim = memory_limit 
        self._nmem = 0
        self._kcounts = {}
    #
    def getOutputInfo(self):
        return {'outdir' : self._outdir,
                'fileprefix' : self._fprefix,
                'nfiles' : self._fcounter}
    #
    def count_and_save(self):
        """counts k-mer and save in file(s)."""
        for kmer in ksup.iterkmers():
            try:
                self._kcounts[kmer] = self._kcounts[kmer] + 1
            except:
                # store items, if memory reaches threshold.
                if(self._nmem == self._memlim):
                    self._store_counts()
                    if self._fcounter % 10 == 0:
                        print("#files: " + str(self._fcounter))
                self._kcounts[kmer] = 1
                self._nmem += 1
        # save unsaved counts
        if self._nmem > 0:
            self._store_counts()
    #
    def _store_counts(self):
        # write in file
        ofn = self._next_filename();
        with open(ofn, 'w') as of:
            text = '\n'.join(['\t'.join([km, str(self._kcounts[km])]) 
                              for km in sorted(self._kcounts)])
            of.write(text)
        # clear memory
        self._clear_memory()
    #
    def _next_filename(self):
        self._fcounter += 1
        return self._outdir + "/" + self._fprefix + str(self._fcounter);
    def _clear_memory(self):
        self._kcounts.clear()
        self._nmem = 0
    def mergeTwoCountFiles(self, fn1, fn2, ofn):
        # initialize variables for saving
        self._nmem_merge = 0
        self._kcounts_merge = []
        of = open(ofn, 'w')
        of.close()
        #
        kmers1 = read_kmer_counts_line_by_line(fn1)
        kmers2 = read_kmer_counts_line_by_line(fn2)
        #
        k1 = next(kmers1, None)
        k2 = next(kmers2, None)
        while True:
            if k1 is None:
                if k2 is not None:
                    self._store_kmer_count_merge(k2, ofn)
                    for k2 in kmers2:
                        self._store_kmer_count_merge(k2, ofn)
                break
            if k2 is None:
                if k1 is not None:
                    self._store_kmer_count_merge(k1, ofn)
                    for k1 in kmers1:
                        self._store_kmer_count_merge(k1, ofn)
                break
            if k1[0] < k2[0]:
                while k1 is not None and k1[0] < k2[0]:
                    self._store_kmer_count_merge(k1, ofn)
                    k1 = next(kmers1, None)
            elif k2[0] < k1[0]:
                while k2 is not None and k2[0] < k1[0]:
                    self._store_kmer_count_merge(k2, ofn)
                    k2 = next(kmers2, None)
            else:
                self._store_kmer_count_merge((k1[0], k1[1]+k2[1]), ofn)
                k1 = next(kmers1, None)
                k2 = next(kmers2,None)
        self._write_to_file_merge(ofn)
    #
    def _write_to_file_merge(self, ofn):
        if len(self._kcounts_merge) > 0:
            text = '\n'.join(['\t'.join([kc[0], str(kc[1])])
                    for kc in self._kcounts_merge]) + '\n'
            with open(ofn, 'a') as of:
                of.write(text)
            self._clear_memory_merge()
    def _clear_memory_merge(self):
        self._kcounts_merge = []
        self._nmem_merge = 0
    def _store_kmer_count_merge(self, kcount, ofn):
        #print(kcount)
        self._kcounts_merge.append(kcount)
        self._nmem_merge += 1
        if(self._nmem_merge == self._memlim):
            self._write_to_file_merge(ofn)
    def mergeMultipleCountFiles(self, filenames, outfilename, tmpdir='results', tmpprefix='tmp'):
        nfiles = int(len(filenames))
        nround = 1
        if nfiles == 1: copy_file(filenames[0], outfilename)
        while nfiles > 1:
            print("#files:" + str(nfiles))
            for fi in range(1, nfiles, 2):
                # input filename
                if nround == 1:
                    fn1 = filenames[fi-1]
                    fn2 = filenames[fi]
                elif nround%2==0:
                    fn1 = tmpdir+'/'+tmpprefix+'odd'+str(fi)
                    fn2 = tmpdir+'/'+tmpprefix+'odd'+str(fi+1)
                else:
                    fn1 = tmpdir+'/'+tmpprefix+'even'+str(fi)
                    fn2 = tmpdir+'/'+tmpprefix+'even'+str(fi+1)
                # output file name
                if nfiles==2:
                    tmpofn = outfilename
                else:
                    tmpofn = tmpdir+'/'+tmpprefix+('even' if nround%2==0 else 'odd') + str(int((fi+1)/2))
                # merge files
                self.mergeTwoCountFiles(fn1, fn2, tmpofn)
            # copy the last single file
            if nfiles%2==1:
                # input filename
                if nround == 1:
                    fn1 = filenames[nfiles-1]
                elif nround%2==0:
                    fn1 = tmpdir+'/'+tmpprefix+'odd'+str(int(nfiles))
                else:
                    fn1 = tmpdir+'/'+tmpprefix+'even'+str(int(nfiles))
                # output file name
                if nfiles==1:
                    tmpofn = outfilename
                else:
                    tmpofn = tmpdir+ '/' +tmpprefix + ('even' if nround%2==0 else 'odd') + str(int((nfiles+1)/2))
                # COPY FILE
                copy_file(fn1, tmpofn)
            nround += 1
            nfiles = int((nfiles+1)/2)
        

## check
# count kmers
# datasets: hbv.sim rymv.sim hpylori.sim hiv1.sim tmv.sim ecoli.sim saureus.sim dmelanogaster.sim
dataset =  'hpylori.sim' # rymv tmv saureus ecoli dmelanogaster
#workdir = '/Users/ashis/Desktop/work/github/sketch-mer'
workdir = '.'
odir = workdir + '/results/tmp'
fn = workdir + '/data/' + dataset + '.fa'
k = 22
memlim = 1000000

ksup = ks.KmerSupplier(fn, k)
kcounter = KmerCounter(ksup, odir, fileprefix= dataset + '_intermediate_', memory_limit=memlim) 

print('counting kmers and saving in intermediate files...')
kcounter.count_and_save()

# merge all intermediate files
finfo = kcounter.getOutputInfo()

print('total #intermediate files: ' + str(finfo['nfiles']))
print('merging all itermediate files...')
filenames = [finfo['outdir']+'/'+finfo['fileprefix']+str(fi) for fi in range(1,finfo['nfiles']+1)]
tmpdir = workdir + '/results/tmp'
tmpprefix = dataset + '_tmp_merge_'
ofn = workdir + '/results/' + dataset + '_counts.txt'
kcounter.mergeMultipleCountFiles(filenames, ofn, tmpdir, tmpprefix)
print('done! see output in ' + ofn)
