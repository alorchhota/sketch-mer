def read_kmer_counts_line_by_line(fn):
    """ Read kmer count file, yielding (kmer, count) tuples one by one """
    with open(fn) as fh:
        for ln in fh:
            ln = ln.strip()
            if len(ln) == 0:
                continue
            items = ln.split('\t')
            yield (items[0], int(items[1]))
