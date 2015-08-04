from .parser import Reader
from cPickle import dumps, loads
from multiprocessing import Queue, Pool


def process(args):
    rdr, line = args
    res = rdr.parse(line)
    return res

rdr = None
def set_reader(reader):
    global rdr
    rdr.reader = reader


class ParReader(object):
    def __init__(self, *args, **kwargs):
        self.ncpus = kwargs.pop('cpus', 2)
        global rdr
        rdr = self.rdr = Reader(*args, **kwargs)

    def __iter__(self):
        pool = Pool(self.ncpus)
        #pool.map(set_reader, range(self.ncpus), chunksize=1)
        for result in pool.imap(process, ((rdr, line) for line in rdr.reader)):
            yield result

    def __getattr__(self, key):
        return getattr(self.rdr, key)

