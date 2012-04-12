import ez_setup
ez_setup.use_setuptools()

import glob
import os
import sys
from setuptools import setup
from distutils.extension import Extension

# optional cython
try:
  from Cython.Distutils import build_ext
except ImportError:
  from distutils.command import build_ext as _build_ext
  class build_ext(_build_ext.build_ext):

      description = "change pyx files to corresponding .c/.cpp (fallback when cython is not installed)"

      def build_extensions(self):
          # First, sanity-check the 'extensions' list
          self.check_extensions_list(self.extensions)
          
          for extension in self.extensions:
              target_ext = '.c'

              patchedsrc = []
              for source in extension.sources:
                (root, ext) = os.path.splitext(source)
                if ext == '.pyx':
                  patchedsrc.append(root + target_ext)
                else:
                  patchedsrc.append(source)

              extension.sources = patchedsrc
              self.build_extension(extension)
  

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

version_py = os.path.join(os.path.dirname(__file__), 'cyvcf', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')

sources=["cyvcf/parser.pyx"]
exts = [ Extension("cyvcf.parser", sources=sources)]

setup(
        cmdclass= {'build_ext': build_ext},
        name="cyvcf",
        version=version,
        ext_modules=exts,
        test_suite='test.test_vcf.suite',
        packages=['cyvcf'],
        author="Aaron Quinlan, James Casbon, John Dougherty, Martin Vermaat, Brent Pedersen",
        description='A fast Python library for VCF files using Cython for speed.',
        url="none",
        package_dir = {"cyvcf": "cyvcf"},
        author_email="arq5x@virginia.edu",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']

    )