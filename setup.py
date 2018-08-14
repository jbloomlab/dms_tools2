"""Setup script for ``dms_tools2``.

Written by Jesse Bloom.
"""

import sys
import os
import subprocess
import re
import glob
try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    raise ImportError("You must install `setuptools`")

if not (sys.version_info[0] == 3 and sys.version_info[1] >= 6):
        raise RuntimeError('dms_tools2 requires Python 3.6 or higher.\n'
            'You are using Python {0}.{1}'.format(
            sys.version_info[0], sys.version_info[1]))

# get metadata, which is specified in another file
metadata = {}
with open('dms_tools2/_metadata.py') as f:
	lines = f.readlines()
for dataname in ['version', 'author', 'author_email', 'url']:
    for line in lines:
        entries = line.split('=')
        assert len(entries) == 2, "Failed parsing metadata:\n{0}".format(line)
        if entries[0].strip() == '__{0}__'.format(dataname):
            if dataname in metadata:
                raise ValueError("Duplicate metadata for {0}".format(dataname))
            else:
                metadata[dataname] = entries[1].strip()[1 : -1]
    assert dataname in metadata, "Failed to find metadata {0}".format(dataname)

with open('README.rst') as f:
    readme = f.read()

# main setup command
setup(
    name = 'dms_tools2', 
    version = metadata['version'],
    author = metadata['author'],
    author_email = metadata['author_email'],
    url = metadata['url'],
    download_url = 'https://github.com/jbloomlab/dms_tools2/tarball/{0}'.format(
		metadata['version']), # assumes tagged version is on GitHub
    description = 'Deep mutational scanning (DMS) analysis tools.',
    long_description = readme,
    license = 'GPLv3',
    install_requires = [
        'biopython>=1.68',
        'HTSeq>=0.9',
        'pysam==0.13', # got an error with later versions
        'pandas>=0.21',
        'numpy>=1.13',
        'IPython>=5.1',
        'jupyter>=1.0.0',
        'matplotlib>=2.1.1',
        'plotnine>=0.3',
        'natsort>=5.0.3',
        'pystan==2.16', # fails with 2.17
        'scipy>=1.0',
        'seaborn>=0.8',
        'phydms>=2.1.4',
        'statsmodels>=0.8',
        'regex>=2.4',
        'packaging',
        'umi_tools>=0.5.4',
        ],
    extras_require = {
        'rplot':[
                'rpy2>=2.9.1',
                'tzlocal', # required by rpy2 but not auto installed in 2.9.3
                ]
        },
    platforms = 'Linux and Mac OS X.',
    packages = ['dms_tools2'],
    package_dir = {'dms_tools2':'dms_tools2'},
    package_data = {'dms_tools2':['rplot_Rcode.R']},
    scripts = [
            'scripts/dms2_bcsubamp',
            'scripts/dms2_batch_bcsubamp',
            'scripts/dms2_prefs',
            'scripts/dms2_batch_prefs',
            'scripts/dms2_diffsel',
            'scripts/dms2_batch_diffsel',
            'scripts/dms2_fracsurvive',
            'scripts/dms2_batch_fracsurvive',
            'scripts/dms2_logoplot',
            ],
    ext_modules = [
            Extension('dms_tools2._cutils', ['dms_tools2/_cutils.c'],
                    extra_compile_args=["-Wno-error=declaration-after-statement"])
            ],
)
