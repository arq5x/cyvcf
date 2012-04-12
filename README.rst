CyVCF
======

A Cython port of the PyVCF library maintained by @jamescasbon.

The goal of this porject is to provide a very fast Python library for parsing and manipulating large VCF files.
Cython has been used to optimize speed.  This version is approximately 4 times faster than PyVCF,
and the parsing speed is essentially identical to that of C/C++ libraries provided by PLINKSEQ and VCFLIB.

The functionality and interface are currently the same as documented here: http://pyvcf.rtfd.org/

Installation::
============

    python setup.py build
    python setup.py install


Testing::
=======

    python setup.py test


Basic usage::
===========

    >>> import cyvcf
    >>> vcf_reader = cyvcf.Reader(open('test/example-4.0.vcf', 'rb'))
    >>> for record in vcf_reader:
    ...     print record
    20	14370	G	A	29.0	.	H2=True;NS=3;DB=True;DP=14;AF=0.5	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
    20	17330	T	A	3.0	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3:.
    20	1110696	A	G,T	67.0	.	AA=T;NS=2;DB=True;DP=10;AF=0.333,0.667	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4:.
    20	1230237	T	.	47.0	.	AA=T;NS=3;DP=13	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2:.
    20	1234567	GTCT	G,GTACT	50.0	.	AA=G;NS=3;DP=9	GT:GQ:DP	./.	0/2:17:2	1/1:40:3
