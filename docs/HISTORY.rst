Development
===========

Please use the repository at github: https://github.com/jamescasbon/PyVCF/
Pull requests gladly accepted. 
Issues should be reported at the github issue tracker.

Changes
=======

0.4.3 Release
-------------

* Single floats in Reader._sample_parser not being converted to float #35
* Handle String INFO values when Number=1 in header #34

0.4.2 Release
-------------

* Installation problems

0.4.1 Release
-------------

* Installation problems

0.4.0 Release
-------------

* Package structure 
* add ``vcf.utils`` module with ``walk_together`` method
* samtools tests 
* support Freebayes' non standard '.' for no call
* fix vcf_melt  
* support monomorphic sites, add ``is_monomorphic`` method, handle null QUALs
* filter support for files with monomorphic calls 
* Values declared as single are no-longer returned in lists
* several performance improvements 


0.3.0 Release
-------------

* Fix setup.py for python < 2.7
* Add ``__eq__`` to ``_Record`` and ``_Call``
* Add ``is_het`` and ``is_variant`` to ``_Call``
* Drop aggressive parse mode: we're always aggressive.
* Add tabix fetch for single calls, fix one->zero based indexing
* add prepend_chr mode for ``Reader`` to add `chr` to CHROM attributes

0.2.2 Release
-------------

Documentation release

0.2.1 Release
-------------

* Add shebang to vcf_filter.py

0.2 Release 
-----------

* Replace genotype dictionary with a ``Call`` object
* Methods on ``Record`` and ``Call`` (thanks @arq5x)
* Shortcut parse_sample when genotype is None

0.1 Release 
-----------

* Added test code
* Added Writer class
* Allow negative number in ``INFO`` and ``FORMAT`` fields (thanks @martijnvermaat)
* Prefer ``vcf.Reader`` to ``vcf.VCFReader``
* Support compressed files with guessing where filename is available on fsock
* Allow opening by filename as well as filesocket
* Support fetching rows for tabixed indexed files
* Performance improvements (see ``test/prof.py``)
* Added extensible filter script (see FILTERS.md), vcf_filter.py 

Contributions
-------------

Project started by @jdoughertyii and taken over by @jamescasbon on 12th January 2011.
Contributions from @arq5x, @brentp, @martijnvermaat, @ian1roberts.


