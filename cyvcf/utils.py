

def walk_together(*readers):
    """ Walk a set of readers and return lists of records from each
        reader, with None if no record present.  Caller must check the
        inputs are sorted in the same way and use the same reference
        otherwise behaviour is undefined.
    """
    nexts = [reader.next() for reader in readers]

    while True:
        min_next = min([x for x in nexts if x is not None])

        # this line uses equality on Records, which checks the ALTs
        # not sure what to do with records that have overlapping but different
        # variation
        yield [x if x is None or x == min_next else None for x in nexts]

        # update nexts that we just yielded
        for i, n in enumerate(nexts):

            if n is not None and n == min_next:
                try:
                    nexts[i] = readers[i].next()
                except StopIteration:
                    nexts[i] = None

        if all([x is None for x in nexts]):
            break


def build_for_pickle(**kwargs):
    from .parser import Reader
    r = Reader.__new__(Reader)

    for k, v in kwargs.items():
        setattr(r, k, v)
    return r

def build_record(chrom, pos, id, ref, alt):
    from .parser import Record
    #start, end, alleles = kwargs.pop('start'), kwargs.pop('end'), kwargs.pop('alleles')
    #_sample_indexes = kwargs.pop('_sample_indexes')
    #has_genotypes = kwargs.pop('has_genotypes')
    self = Record.__new__(Record, chrom, pos, id, ref, alt, INFO={})
    # these are set in __init_ so we have to re-set them
    return self
    self.alleles = alleles
    self.start = start
    self.end = end
    self._sample_indexes = _sample_indexes
    self.has_genotypes = has_genotypes
    print self.gt_types, kwargs['gt_types']
    return self

def build_filter(*args):
    from .parser import Filter
    f = Filter(*args)
    return f

def build_format(*args):
    from .parser import Format
    f = Format(*args)
    return f

def build_info(*args):
    from .parser import Info
    i = Info(*args)
    return i
