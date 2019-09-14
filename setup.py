from setuptools import setup
from gumpy import __version__

setup(
    name='gumpy',
    version=__version__,
    author='Philip W Fowler',
    author_email="philip.fowler@ndm.ox.ac.uk",
    description="Genetics with Numpy",
    url="https://github.com/philipwfowler/gumpy",
    packages=['gumpy'],
    package_data={'': ['../config/*']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"],
    python_requires='>=3.5',
    install_requires=[
        "numpy >= 1.13",
        "pytest >= 4.0.0",
        "Biopython >= 1.70",
        "tqdm >= 4.19.5",
        "pysam >= 0.15.2"
    ],
    license="University of Oxford, see LICENCE.md",
    scripts=['bin/gumpy-save-genome.py'],\
    zip_safe=False
)
