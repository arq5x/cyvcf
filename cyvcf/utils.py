

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
