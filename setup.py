from setuptools import setup
from gumpy import __version__

setup(
    name='gumpy',
    version=__version__,
    author='Philip W Fowler',
    author_email="philip.fowler@ndm.ox.ac.uk",
    packages=['gumpy'],
    package_data={'': ['../config/*']},
    install_requires=[
        "numpy >= 1.13",
        "pandas >= 0.23.1",
        "Biopython >= 1.70",
        "tqdm >= 4.19.5"
    ],
    license=None,
    long_description=open('README.md').read(),
    zip_safe=False
)

#     scripts=['bin/piezo-vcf-parse.py'],
#         "PyVCF >= 0.6.8",
