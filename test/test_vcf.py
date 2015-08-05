import unittest
import doctest
import os
import commands
from StringIO import StringIO

import cyvcf
from cyvcf import utils

suite = doctest.DocTestSuite(cyvcf.parser)


def fh(fname):
    return file(os.path.join(os.path.dirname(__file__), fname))


class TestVcfSpecs(unittest.TestCase):

    def test_vcf_4_0(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        assert reader.metadata['fileformat'] == 'VCFv4.0'

        # test we can walk the file at least
        for r in reader:

            if r.POS == 1230237:
                assert r.is_monomorphic
            else:
                assert not r.is_monomorphic

            if 'AF' in r.INFO:
                self.assertEqual(type(r.INFO['AF']),  type([]))

            for c in r:
                assert c

                # issue 19, in the example ref the GQ is length 1
                if c.called:
                    self.assertEqual(type(c.data['GQ']),  type(1))
                    if 'HQ' in c.data and c.data['HQ'] is not None:
                        self.assertEqual(type(c.data['HQ']),  type([]))



    def test_vcf_4_1(self):
        return
        reader = cyvcf.Reader(fh('example-4.1.vcf'))
        self.assertEqual(reader.metadata['fileformat'],  'VCFv4.1')

        # contigs were added in vcf4.1
        # probably need to add a reader.contigs attribute
        assert 'contig' in reader.metadata

        # test we can walk the file at least
        for r in reader:
            for c in r:
                assert c

        # asserting False while I work out what to check
        assert False

    def test_vcf_4_1_sv(self):
        return

        reader = cyvcf.Reader(fh('example-4.1-sv.vcf'))

        assert 'SVLEN' in reader.infos

        # test we can walk the file at least
        for r in reader:
            print r
            for c in r:
                print c
                assert c

        # asserting False while I work out what to check
        assert False


class TestGatkOutput(unittest.TestCase):

    filename = 'gatk.vcf'

    samples = ['BLANK', 'NA12878', 'NA12891', 'NA12892',
            'NA19238', 'NA19239', 'NA19240']
    formats = ['AD', 'DP', 'GQ', 'GT', 'PL']
    infos = ['AC', 'AF', 'AN', 'BaseQRankSum', 'DB', 'DP', 'DS',
            'Dels', 'FS', 'HRun', 'HaplotypeScore', 'InbreedingCoeff',
            'MQ', 'MQ0', 'MQRankSum', 'QD', 'ReadPosRankSum']

    n_calls = 37

    def setUp(self):
        self.reader = cyvcf.Reader(fh(self.filename))

    def testSamples(self):
        self.assertEqual(self.reader.samples, self.samples)

    def testFormats(self):
        self.assertEqual(set(self.reader.formats), set(self.formats))

    def testInfos(self):
        self.assertEqual(set(self.reader.infos), set(self.infos))


    def testCalls(self):
        n = 0

        for site in self.reader:
            n += 1
            self.assertEqual(len(site.samples), len(self.samples))


            # check sample name lookup
            for s in self.samples:
                assert site.genotype(s)

            # check ordered access
            self.assertEqual([x.sample for x in site.samples], self.samples)
        self.assertEqual(len(site.gt_phred_likelihoods), len(self.samples))
        self.assertEqual(n,  self.n_calls)


class TestFreebayesOutput(TestGatkOutput):

    filename = 'freebayes.vcf'
    formats = ['AO', 'DP', 'GL', 'GLE', 'GQ', 'GT', 'QA', 'QR', 'RO']
    infos = ['AB', 'ABP', 'AC', 'AF', 'AN', 'AO', 'BVAR', 'CIGAR',
            'DB', 'DP', 'DPRA', 'EPP', 'EPPR', 'HWE', 'LEN', 'MEANALT',
            'NUMALT', 'RPP', 'MQMR', 'ODDS', 'MQM', 'PAIREDR', 'PAIRED',
            'SAP', 'XRM', 'RO', 'REPEAT', 'XRI', 'XAS', 'XAI', 'SRP',
            'XAM', 'XRS', 'RPPR', 'NS', 'RUN', 'CpG', 'TYPE']
    n_calls = 104


    def testParse(self):
        reader = cyvcf.Reader(fh('freebayes.vcf'))
        print reader.samples
        self.assertEqual(len(reader.samples), 7)
        n = 0
        for r in reader:
            n+=1
            for x in r:
                assert x
        assert n == self.n_calls

class TestSamtoolsOutput(unittest.TestCase):

    def testParse(self):
        reader = cyvcf.Reader(fh('samtools.vcf'))

        self.assertEqual(len(reader.samples), 1)
        self.assertEqual(sum(1 for _ in reader), 11)


class Test1kg(unittest.TestCase):

    def testParse(self):
        reader = cyvcf.Reader(fh('1kg.vcf.gz'))

        self.assertEqual(len(reader.samples), 629)
        for _ in reader:
            pass


class TestWriter(unittest.TestCase):

    def testWrite(self):

        reader = cyvcf.Reader(fh('gatk.vcf'))
        out = StringIO()
        writer = cyvcf.Writer(out, reader)

        records = list(reader)

        map(writer.write_record, records)
        out.seek(0)
        reader2 = cyvcf.Reader(out)

        self.assertEquals(reader.samples, reader2.samples)
        self.assertEquals(reader.formats, reader2.formats)

        for k in reader.infos:
            self.assertEquals(reader.infos[k], reader2.infos[k], (reader.infos[k], reader2.infos[k]))

        for l, r in zip(records, reader2):
            self.assertEquals(l.samples, r.samples)

class TestRecord(unittest.TestCase):

    def test_num_calls(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            num_calls = (var.num_hom_ref + var.num_hom_alt + \
                         var.num_het + var.num_unknown)
            self.assertEqual(len(var.samples), num_calls)

    def test_call_rate(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            call_rate = var.call_rate
            if var.POS == 14370:
                self.assertEqual(3.0/3.0, call_rate)
            if var.POS == 17330:
                self.assertEqual(3.0/3.0, call_rate)
            if var.POS == 1110696:
                self.assertEqual(3.0/3.0, call_rate)
            if var.POS == 1230237:
                self.assertEqual(3.0/3.0, call_rate)
            elif var.POS == 1234567:
                self.assertEqual(2.0/3.0, call_rate)

    def test_aaf(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            aaf = var.aaf
            if var.POS == 14370:
                self.assertEqual(3.0/6.0, aaf)
            if var.POS == 17330:
                self.assertEqual(1.0/6.0, aaf)
            if var.POS == 1110696:
                self.assertEqual(None, aaf)
            if var.POS == 1230237:
                self.assertEqual(0.0/6.0, aaf)
            elif var.POS == 1234567:
                self.assertEqual(None, aaf)

    def test_pi(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            pi = var.nucl_diversity
            if var.POS == 14370:
                self.assertEqual(6.0/10.0, pi)
            if var.POS == 17330:
                self.assertEqual(1.0/3.0, pi)
            if var.POS == 1110696:
                self.assertEqual(None, pi)
            if var.POS == 1230237:
                self.assertEqual(0.0/6.0, pi)
            elif var.POS == 1234567:
                self.assertEqual(None, pi)

    def test_is_snp(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            is_snp = var.is_snp
            if var.POS == 14370:
                self.assertEqual(True, is_snp)
            if var.POS == 17330:
                self.assertEqual(True, is_snp)
            if var.POS == 1110696:
                self.assertEqual(True, is_snp)
            if var.POS == 1230237:
                self.assertEqual(False, is_snp)
            elif var.POS == 1234567:
                self.assertEqual(False, is_snp)

    def test_is_indel(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            is_indel = var.is_indel
            if var.POS == 14370:
                self.assertEqual(False, is_indel)
            if var.POS == 17330:
                self.assertEqual(False, is_indel)
            if var.POS == 1110696:
                self.assertEqual(False, is_indel)
            if var.POS == 1230237:
                self.assertEqual(True, is_indel)
            elif var.POS == 1234567:
                self.assertEqual(True, is_indel)

    def test_is_transition(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            is_trans = var.is_transition
            if var.POS == 14370:
                self.assertEqual(True, is_trans)
            if var.POS == 17330:
                self.assertEqual(False, is_trans)
            if var.POS == 1110696:
                self.assertEqual(False, is_trans)
            if var.POS == 1230237:
                self.assertEqual(False, is_trans)
            elif var.POS == 1234567:
                self.assertEqual(False, is_trans)

    def test_is_deletion(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            is_del = var.is_deletion
            if var.POS == 14370:
                self.assertEqual(False, is_del)
            if var.POS == 17330:
                self.assertEqual(False, is_del)
            if var.POS == 1110696:
                self.assertEqual(False, is_del)
            if var.POS == 1230237:
                self.assertEqual(True, is_del)
            elif var.POS == 1234567:
                self.assertEqual(False, is_del)

    def test_var_type(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            type = var.var_type
            if var.POS == 14370:
                self.assertEqual("snp", type)
            if var.POS == 17330:
                self.assertEqual("snp", type)
            if var.POS == 1110696:
                self.assertEqual("snp", type)
            if var.POS == 1230237:
                self.assertEqual("indel", type)
            elif var.POS == 1234567:
                self.assertEqual("indel", type)
        # SV tests
        reader = cyvcf.Reader(fh('example-4.1-sv.vcf'))
        for var in reader:
            type = var.var_type
            if var.POS == 2827693:
                self.assertEqual("sv", type)
            if var.POS == 321682:
                self.assertEqual("sv", type)
            if var.POS == 14477084:
                self.assertEqual("sv", type)
            if var.POS == 9425916:
                self.assertEqual("sv", type)
            elif var.POS == 12665100:
                self.assertEqual("sv", type)
            elif var.POS == 18665128:
                self.assertEqual("sv", type)


    def test_var_subtype(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            subtype = var.var_subtype
            if var.POS == 14370:
                self.assertEqual("ts", subtype)
            if var.POS == 17330:
                self.assertEqual("tv", subtype)
            if var.POS == 1110696:
                self.assertEqual("unknown", subtype)
            if var.POS == 1230237:
                self.assertEqual("del", subtype)
            elif var.POS == 1234567:
                self.assertEqual("unknown", subtype)
        # SV tests
        reader = cyvcf.Reader(fh('example-4.1-sv.vcf'))
        for var in reader:
            subtype = var.var_subtype
            if var.POS == 2827693:
                self.assertEqual("DEL", subtype)
            if var.POS == 321682:
                self.assertEqual("DEL", subtype)
            if var.POS == 14477084:
                self.assertEqual("DEL:ME:ALU", subtype)
            if var.POS == 9425916:
                self.assertEqual("INS:ME:L1", subtype)
            elif var.POS == 12665100:
                self.assertEqual("DUP", subtype)
            elif var.POS == 18665128:
                self.assertEqual("DUP:TANDEM", subtype)

    def test_is_sv(self):
        reader = cyvcf.Reader(fh('example-4.1-sv.vcf'))
        for var in reader:
            is_sv = var.is_sv
            if var.POS == 2827693:
                self.assertEqual(True, is_sv)
            if var.POS == 321682:
                self.assertEqual(True, is_sv)
            if var.POS == 14477084:
                self.assertEqual(True, is_sv)
            if var.POS == 9425916:
                self.assertEqual(True, is_sv)
            elif var.POS == 12665100:
                self.assertEqual(True, is_sv)
            elif var.POS == 18665128:
                self.assertEqual(True, is_sv)

        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            is_sv = var.is_sv
            if var.POS == 14370:
                self.assertEqual(False, is_sv)
            if var.POS == 17330:
                self.assertEqual(False, is_sv)
            if var.POS == 1110696:
                self.assertEqual(False, is_sv)
            if var.POS == 1230237:
                self.assertEqual(False, is_sv)
            elif var.POS == 1234567:
                self.assertEqual(False, is_sv)

    def test_is_sv_precise(self):
        reader = cyvcf.Reader(fh('example-4.1-sv.vcf'))
        for var in reader:
            is_precise = var.is_sv_precise
            if var.POS == 2827693:
                self.assertEqual(True, is_precise)
            if var.POS == 321682:
                self.assertEqual(False, is_precise)
            if var.POS == 14477084:
                self.assertEqual(False, is_precise)
            if var.POS == 9425916:
                self.assertEqual(False, is_precise)
            elif var.POS == 12665100:
                self.assertEqual(False, is_precise)
            elif var.POS == 18665128:
                self.assertEqual(False, is_precise)

        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            is_precise = var.is_sv_precise
            if var.POS == 14370:
                self.assertEqual(False, is_precise)
            if var.POS == 17330:
                self.assertEqual(False, is_precise)
            if var.POS == 1110696:
                self.assertEqual(False, is_precise)
            if var.POS == 1230237:
                self.assertEqual(False, is_precise)
            elif var.POS == 1234567:
                self.assertEqual(False, is_precise)

    def test_sv_end(self):
        reader = cyvcf.Reader(fh('example-4.1-sv.vcf'))
        for var in reader:
            sv_end = var.sv_end
            if var.POS == 2827693:
                self.assertEqual(2827680, sv_end)
            if var.POS == 321682:
                self.assertEqual(321887, sv_end)
            if var.POS == 14477084:
                self.assertEqual(14477381, sv_end)
            if var.POS == 9425916:
                self.assertEqual(9425916, sv_end)
            elif var.POS == 12665100:
                self.assertEqual(12686200, sv_end)
            elif var.POS == 18665128:
                self.assertEqual(18665204, sv_end)

        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            sv_end = var.sv_end
            if var.POS == 14370:
                self.assertEqual(None, sv_end)
            if var.POS == 17330:
                self.assertEqual(None, sv_end)
            if var.POS == 1110696:
                self.assertEqual(None, sv_end)
            if var.POS == 1230237:
                self.assertEqual(None, sv_end)
            elif var.POS == 1234567:
                self.assertEqual(None, sv_end)


class TestCall(unittest.TestCase):

    def test_phased(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            phases = var.gt_phases
            print var
            if var.POS == 14370:
                self.assertEqual([True, True, False], phases)
            if var.POS == 17330:
                self.assertEqual([True, True, False], phases)
            if var.POS == 1110696:
                self.assertEqual([True, True, False], phases)
            if var.POS == 1230237:
                self.assertEqual([True, True, False], phases)
            elif var.POS == 1234567:
                self.assertEqual([False, False, False], phases)

    def test_gt_bases(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            gt_bases = [s.gt_bases for s in var.samples]
            if var.POS == 14370:
                self.assertEqual(['G|G', 'A|G', 'A/A'], gt_bases)
            elif var.POS == 17330:
                self.assertEqual(['T|T', 'T|A', 'T/T'], gt_bases)
            elif var.POS == 1110696:
                self.assertEqual(['G|T', 'T|G', 'T/T'], gt_bases)
            elif var.POS == 1230237:
                self.assertEqual(['T|T', 'T|T', 'T/T'], gt_bases)
            elif var.POS == 1234567:
                self.assertEqual([None, 'GTCT/GTACT', 'G/G'], gt_bases)

    def test_gt_types(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            for s in var:
                print s.data
            gt_types = [s.gt_type for s in var.samples]
            if var.POS == 14370:
                self.assertEqual([0,1,3], gt_types)
            elif var.POS == 17330:
                self.assertEqual([0,1,0], gt_types)
            elif var.POS == 1110696:
                self.assertEqual([1,1,3], gt_types)
            elif var.POS == 1230237:
                self.assertEqual([0,0,0], gt_types)
            elif var.POS == 1234567:
                self.assertEqual([None,1,3], gt_types)

    def test_gt_depths(self):
        reader = cyvcf.Reader(fh('example-4.0.vcf'))
        for var in reader:
            for s in var:
                print s.data
            gt_depths = [s.gt_depth for s in var.samples]
            if var.POS == 14370:
                self.assertEqual([1,8,5], gt_depths)
            elif var.POS == 17330:
                self.assertEqual([3,5,3], gt_depths)
            elif var.POS == 1110696:
                self.assertEqual([6,0,4], gt_depths)
            elif var.POS == 1230237:
                self.assertEqual([7,4,2], gt_depths)
            elif var.POS == 1234567:
                self.assertEqual([4,2,3], gt_depths)

    def test_gt_ref_depths(self):

        reader = cyvcf.Reader(fh('gatk.vcf'))
        for var in reader:
            gt_ref_depths = [s.gt_ref_depth for s in var.samples]
            if var.POS == 42522392:
                self.assertEqual([6,138,169,249,248,250,250], gt_ref_depths)
            elif var.POS == 42522613:
                self.assertEqual([13,118,241,161,110,106,116], gt_ref_depths)
            elif var.POS == 42527891:
                self.assertEqual([-1,238,246,239,232,233,238], gt_ref_depths)

    def test_gt_alt_depths(self):

        reader = cyvcf.Reader(fh('gatk.vcf'))
        for var in reader:
            gt_alt_depths = [s.gt_alt_depth for s in var.samples]
            if var.POS == 42522392:
                self.assertEqual([0,107,77,0,1,0,0], gt_alt_depths)
            elif var.POS == 42522613:
                self.assertEqual([4,127,0,85,132,135,126], gt_alt_depths)
            elif var.POS == 42527891:
                self.assertEqual([-1,7,3,11,16,14,11], gt_alt_depths)

    def test_gt_quals(self):

        reader = cyvcf.Reader(fh('gatk.vcf'))
        for var in reader:
            gt_quals = [s.gt_qual for s in var.samples]
            if var.POS == 42522392:
                self.assertEqual([18.04,99,99,99,99,99,99], gt_quals)
            elif var.POS == 42522613:
                self.assertEqual([62.64,99,99,99,99,99,99], gt_quals)
            elif var.POS == 42527891:
                self.assertEqual([-1,13.70,5.97,31.42,49.09,52.10,12.71], gt_quals)


class TestTabix(unittest.TestCase):

    def setUp(self):
        self.reader = cyvcf.Reader(fh('tb.vcf.gz'))

        self.run = cyvcf.parser.pysam is not None


    def testFetchRange(self):
        if not self.run:
            return
        lines = list(self.reader.fetch('20', 14370, 14370))
        self.assertEquals(len(lines), 1)
        self.assertEqual(lines[0].POS, 14370)

        lines = list(self.reader.fetch('20', 14370, 17330))
        self.assertEquals(len(lines), 2)
        self.assertEqual(lines[0].POS, 14370)
        self.assertEqual(lines[1].POS, 17330)


        lines = list(self.reader.fetch('20', 1110695, 1234567))
        self.assertEquals(len(lines), 3)

    def testFetchSite(self):
        if not self.run:
            return
        site = self.reader.fetch('20', 14370)
        assert site.POS == 14370

        site = self.reader.fetch('20', 14369)
        assert site is None




class TestOpenMethods(unittest.TestCase):

    samples = 'NA00001 NA00002 NA00003'.split()

    def testOpenFilehandle(self):
        r = cyvcf.Reader(fh('example-4.0.vcf'))
        self.assertEqual(self.samples, r.samples)
        self.assertEqual('example-4.0.vcf', os.path.split(r.filename)[1])

    def testOpenFilename(self):
        r = cyvcf.Reader(filename='test/example-4.0.vcf')
        self.assertEqual(self.samples, r.samples)

    def testOpenFilehandleGzipped(self):
        r = cyvcf.Reader(fh('tb.vcf.gz'))
        self.assertEqual(self.samples, r.samples)

    def testOpenFilenameGzipped(self):
        r = cyvcf.Reader(filename='test/tb.vcf.gz')
        self.assertEqual(self.samples, r.samples)


class TestFilter(unittest.TestCase):


    def testApplyFilter(self):
        s, out = commands.getstatusoutput('python scripts/vcf_filter.py --site-quality 30 test/example-4.0.vcf sq')
        #print out
        assert s == 0
        buf = StringIO()
        buf.write(out)
        buf.seek(0)

        print buf.getvalue()
        reader = cyvcf.Reader(buf)


        # check filter got into output file
        assert 'sq30' in reader.filters

        print reader.filters

        # check sites were filtered
        n = 0
        for r in reader:
            if r.QUAL < 30:
                assert 'sq30' in r.FILTER
                n += 1
            else:
                assert r.FILTER is None or 'sq30' not in r.FILTER
        assert n == 2


    def testApplyMultipleFilters(self):
        s, out = commands.getstatusoutput('python scripts/vcf_filter.py --site-quality 30 '
        '--genotype-quality 50 test/example-4.0.vcf sq mgq')
        assert s == 0
        #print out
        buf = StringIO()
        buf.write(out)
        buf.seek(0)
        reader = cyvcf.Reader(buf)

        print reader.filters

        assert 'mgq50' in reader.filters
        assert 'sq30' in reader.filters


class TestRegression(unittest.TestCase):

    def test_issue_16(self):
        reader = cyvcf.Reader(fh('issue-16.vcf'))
        assert reader.next().QUAL == None

    def test_null_mono(self):
        # null qualities were written as blank, causing subsequent parse to fail
        print os.path.abspath(os.path.join(os.path.dirname(__file__),  'null_genotype_mono.vcf'))
        p = cyvcf.Reader(fh('null_genotype_mono.vcf'))
        assert p.samples
        out = StringIO()
        writer = cyvcf.Writer(out, p)
        map(writer.write_record, p)
        out.seek(0)
        print out.getvalue()
        p2 = cyvcf.Reader(out)
        rec = p2.next()
        assert rec.samples


class TestUtils(unittest.TestCase):

    def test_walk(self):
        # easy case: all same sites
        reader1 = cyvcf.Reader(fh('example-4.0.vcf'))
        reader2 = cyvcf.Reader(fh('example-4.0.vcf'))
        reader3 = cyvcf.Reader(fh('example-4.0.vcf'))

        n = 0
        for x in utils.walk_together(reader1, reader2, reader3):
            assert len(x) == 3
            assert (x[0] == x[1]) and (x[1] == x[2])
            n+= 1
        assert n == 5

        # artificial case 2 from the left, 2 from the right, 2 together, 1 from the right, 1 from the left

        expected = 'llrrttrl'
        reader1 = cyvcf.Reader(fh('walk_left.vcf'))
        reader2 = cyvcf.Reader(fh('example-4.0.vcf'))

        for ex, recs in zip(expected, utils.walk_together(reader1, reader2)):

            if ex == 'l':
                assert recs[0] is not None
                assert recs[1] is None
            if ex == 'r':
                assert recs[1] is not None
                assert recs[0] is None
            if ex == 't':
                assert recs[0] is not None
                assert recs[1] is not None


class TestAD(unittest.TestCase):
    def setUp(self):
        self.reader = cyvcf.Reader(fh('test.vcf'))

    def testRefDepth(self):
        v = self.reader.next()
        self.assertEqual(v.samples[0].gt_ref_depth, -1)

class TestGLInt(unittest.TestCase):
    def setUp(self):
        self.reader = cyvcf.Reader(fh('test-gl.vcf'))
    def testGLInt(self):
        v = next(self.reader)
        self.assertEqual(v.samples[0].gt_phred_likelihoods, None)



suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestAD))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestGatkOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFreebayesOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSamtoolsOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestWriter))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestTabix))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestOpenMethods))
#suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFilter))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(Test1kg))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestRecord))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCall))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestRegression))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestGLInt))
