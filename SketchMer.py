import sys
from countminsketch import CountMinSketch
from countLeastSquares import countLeastSquares
import KmerSupplier as ks

# configuration inputs
if len(sys.argv) < 3:
    print("parameters missing: filename, k")
    sys.exit()
fn = sys.argv[1]
k = int(sys.argv[2])

def build_countminsketch(ksup, w=1000, h=10):
    """returns a countminsketch object inserting all kmers from given KmerSupplier."""
    sketch = CountMinSketch(w, h)
    for kmer in ksup.iterkmers():
        sketch.add(kmer)
    return sketch

def build_countLSsketch(ksup,delta,epsilon):
    sketch = countLeastSquares(delta=10**-3, epsilon = 0.005, k=10)  # change the entry for k
    for kmer in ksup.iterkmers():
        sketch.update(kmer,1)
    return sketch

## for countmin
ksup = ks.KmerSupplier(fn, k)
sketch = build_countminsketch(ksup)

## for leastsquares
ksup = ks.KmerSupplier(fn, k)
## delta and epsilon can be taken from the user too
sketch_ls = build_countLSsketch(ksup,delta=10**-3,epsilon=0.005)
least_est = sketch_ls.lsquare(keys)

