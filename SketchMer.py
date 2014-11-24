import sys
from countminsketch import CountMinSketch
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

ksup = ks.KmerSupplier(fn, k)
sketch = build_countminsketch(ksup)


