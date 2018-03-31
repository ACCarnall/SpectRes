from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='spectres', 

    version='2.0.0',

    description='Simple spectral resampling',

    long_description=long_description,

    url='https://spectres.readthedocs.io',

    author='Adam Carnall',

    author_email='adamc@roe.ac.uk',

    packages=["spectres"],

    install_requires=['numpy'],  # Optional

    project_urls={  # Optional
        "GitHub": "https://github.com/ACCarnall/spectres",
        "ArXiv": "https://arxiv.org/abs/1705.05165",
    },
)