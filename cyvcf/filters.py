

class Base(object):
    """ Base class for vcf_filter.py filters """

    name = 'f'
    """ name used to activate filter and in VCF headers """

    description = 'VCF filter base class'
    """ descrtiption used in vcf headers """

    @classmethod
    def customize_parser(self, parser):
        """ hook to extend argparse parser with custom arguments """
        pass

    def __init__(self, args):
        """ create the filter using argparse ``args`` """
        self.threshold = 0

    def __call__(self):
        """ filter a site, return not None if the site should be filtered """
        raise NotImplementedError('Filters must implement this method')


    def filter_name(self):
        """ return the name to put in the VCF header, default is ``name`` + ``threshold`` """
        return '%s%s' % (self.name, self.threshold)


class SiteQuality(Base):

    description = 'Filter sites by quality'
    name = 'sq'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--site-quality', type=int, default=30,
                help='Filter sites below this quality')

    def __init__(self, args):
        self.threshold = args.site_quality

    def __call__(self, record):
        if record.QUAL < self.threshold:
            return record.QUAL


class VariantGenotypeQuality(Base):

    description = 'Demand a minimum quality associated with a non reference call'
    name = 'mgq'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--genotype-quality', type=int, default=50,
                help='Filter sites with no genotypes above this quality')

    def __init__(self, args):
        self.threshold = args.genotype_quality

    def __call__(self, record):
        if not record.is_monomorphic:
            vgq = max([x['GQ'] for x in record if x.is_variant])
            if vgq < self.threshold:
                return vgq


