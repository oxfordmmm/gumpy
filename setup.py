from setuptools import setup
from gumpy import __version__

with open("README.md", "r") as f:
    README = f.read()

setup(
    install_requires=[
        "numpy >= 1.13",\
        "pysam >= 0.18.0",\
        "biopython >= 1.70",
        "tqdm >= 4.19.5",\
        "pytest >= 4.0.0",\
        'pytest-cov >= 3.0.0',\
        'pandas >= 1.3.0',\
        'scipy >= 1.7.0'        
        ],
    name='gumpy',
    version=__version__,
    author='Philip W Fowler',
    author_email="philip.fowler@ndm.ox.ac.uk",
    description="Genetics with Numpy",
    long_descriptions=README,
    url="https://github.com/oxfordmmm/gumpy",
    packages=['gumpy'],
    package_data={'': ['../config/*']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"],
    python_requires='>=3.8',
    license="University of Oxford, see LICENSE.md",
    scripts=['bin/gumpy-save-genome.py', 'bin/to_piezo_catalogue.py'],\
    zip_safe=False
)
