def read_fasta_char_by_char(fn):
    """ Read FASTA file, yielding nucleotides one by one """
    with open(fn) as fh:
        for ln in fh:
            if ln[0] == '>':
                continue
            for c in ln.strip():
                yield c

def read_fasta_kmer_by_kmer(fn, k):
    """Creates an iterator  for all k-mers.
    Assumption: there is at least k nucleotides in the file,
    otherwise, this function will throw an exception."""
    fastaCharIterator = read_fasta_char_by_char(fn)
    curKmer = '#'
    for i in range(k-1):
        curKmer += next(fastaCharIterator)
    for c in fastaCharIterator:
        curKmer = curKmer[1:] + c
        yield curKmer

class KmerSupplier(object):
    """given filename and k, generates kmer iterators based on file type."""
    def __init__(self, fn, k, ftype='fasta'):
        if not fn or not k:
            raise ValueError("Invalid filename or k!")
        self._fn = fn
        self._k = k
        self._ftype = ftype.lower()

    def iterkmers(self):
        if self._ftype == 'fasta':
            return (read_fasta_kmer_by_kmer(self._fn, self._k))
        else:
            raise NotImplementedError('Sorry filetype ' + self._ftype + ' is not supported yet!')
        return

