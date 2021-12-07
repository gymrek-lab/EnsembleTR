import os
from setuptools import setup, find_packages

DESCRIPTION = "Toolkit for ensemble TR genotyping"
LONG_DESCRIPTION = DESCRIPTION
NAME = "EnsembleTR"
AUTHOR = "Melissa Gymrek"
AUTHOR_EMAIL = "mgymrek@ucsd.edu"
MAINTAINER = "Melissa Gymrek"
MAINTAINER_EMAIL = "mgymrek@ucsd.edu"
DOWNLOAD_URL = 'http://github.com/gymreklab/TRTools'
LICENSE = 'MIT'

# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 1
MIN = 0
REV = 0
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'ensembletr/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=DOWNLOAD_URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      python_requires='>=3.6',
      packages=find_packages(),
      include_package_data=True,
      license_file="LICENSE.txt",
      scripts=["ensembletr/testsupport/test_ensembletr.sh"],
      entry_points={
          'console_scripts': [
              'EnsembleTR=ensembletr.main:run'
          ],
      },
      install_requires=['cyvcf2',
                        'networkx',
                        'numpy',
                        'pyfaidx'],
      classifiers=['Development Status :: 4 - Beta',\
                       'Programming Language :: Python :: 3.5',\
                       'License :: OSI Approved :: MIT License',\
                       'Operating System :: OS Independent',\
                       'Intended Audience :: Science/Research',\
                       'Topic :: Scientific/Engineering :: Bio-Informatics']
     )
