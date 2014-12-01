__author__ = 'princy'

import numpy as np
from countLeastSquares import countLeastSquares
import KmerSupplier as ks

def lsquare(self):
        A = np.zeros((self.d*self.w,self.n+1), dtype='int32')
        b = self.count.ravel()
        len_list = len(list(self.top_k.iterkeys()))
        b = self.count.ravel()
        flag = True
        for kmer in self.top_k.iterkeys():
            if flag:
                col = 0
                flag=False
            else:
                col += 1
            for row, hash_function in enumerate(self.hash_functions):
                column = hash_function(abs(hash(kmer)))
                A[row*self.w+column][col] = 1
        for i in xrange(A.shape[1]):
            A[i][len_list]= 1
        x = np.dot(np.linalg.pinv(A),b)
        kmer_key = list(self.top_k.iterkeys())
        result = {};
        for i in xrange(len_list):
            result[kmer_key[i]] = min(x[i], sketch.get(kmer_key[i]))
        return result

def build_countLSsketch(ksup,delta,epsilon):
    sketch = countLeastSquares(10**-7, 0.005,len(list(ksup.iterkmers()))) # change the entry for k
    for kmer in ksup.iterkmers():
        sketch.update(kmer,1)
    return sketch

ksup = ks.KmerSupplier("../DATA/g1.fa", 22)
## delta and epsilon can be taken from the user too
sketch = build_countLSsketch(ksup,delta=10**-7,epsilon=0.005)
least_est = lsquare(sketch)