Filtering VCF files
===================

The filter script: vcf_filter.py
--------------------------------

Filtering a VCF file based on some properties of interest is a common enough 
operation that PyVCF offers an extensible script.  ``vcf_filter.py`` does 
the work of reading input, updating the metadata and filtering the records.


Adding a filter
---------------

You can reuse this work by providing a filter class, rather than writing your own filter.
For example, lets say I want to filter each site based on the quality of the site.
I can create a class like this::
    
    class SiteQuality(vcf.Filter):

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


This class subclasses ``vcf.Filter`` which provides the interface for VCF filters.
The ``description``` and ``name`` are metadata about the parser.
The ``customize_parser`` method allows you to add arguments to the script.
We use the ``__init__`` method to grab the argument of interest from the parser.
Finally, the ``__call__`` method processes each record and returns a value if the 
filter failed.  The base class uses the ``name`` and ``threshold`` to create
the filter ID in the VCF file.

To make vcf_filter.py aware of the filter, you need to declare a ``vcf.filters`` entry 
point in your ``setup``::

    setup(
        ...
        entry_points = {
            'vcf.filters': [
                'site_quality = module.path:SiteQuality',
            ]
        }
    )

Now when you call vcf_filter.py, you should see your filter in the list of available filters::

    >$ vcf_filter.py --help
    usage: vcf_filter.py [-h] [--no-short-circuit] [--output OUTPUT]
                         [--site-quality SITE_QUALITY]
                         [--genotype-quality GENOTYPE_QUALITY]
                         input filter [filter ...]

    Filter a VCF file

    available filters:
      sq:	Filter sites by quality

    positional arguments:
      input                 File to process (use - for STDIN)
      filter                Filters to use

    optional arguments:
      -h, --help            show this help message and exit
      --no-short-circuit    Do not stop filter processing on a site if a single
                            filter fails.
      --output OUTPUT       Filename to output (default stdout)
      --site-quality SITE_QUALITY
                            Filter sites below this quality
      --genotype-quality GENOTYPE_QUALITY
                            Filter sites with no genotypes above this quality


The filter base class: vcf.Filter
---------------------------------

.. autoclass:: vcf.Filter
   :members:

