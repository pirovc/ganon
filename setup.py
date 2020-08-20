import io
import os
import re

from setuptools import find_packages
from setuptools import setup

def read(filename):
    filename = os.path.join(os.path.dirname(__file__), filename)
    text_type = type(u"")
    with io.open(filename, mode="r", encoding='utf-8') as fd:
        return re.sub(text_type(r':[a-z]+:`~?(.*?)`'), text_type(r'``\1``'), fd.read())

setup(
    name="ganon",
    version="0.3.0",
    url="https://www.github.com/pirovc/ganon",
    license='MIT',

    author="Vitor C. Piro",
    author_email="pirovc@posteo.net",    

    description="ganon is a k-mer based read classification tool which uses Interleaved Bloom Filters in conjunction with a taxonomic clustering and a k-mer counting-filtering scheme.",
    long_description=read("README.md"),

    package_dir={'':'src'},
    packages=["ganon"],
    entry_points = {'console_scripts': ['ganon=ganon.ganon:main']},
    
    scripts=['scripts/ganon-convert-db-0.1-0.2.py',
            'scripts/ganon-convert-db-0.2-0.3.py',
            'scripts/ganon-get-len-taxid.sh'],

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
)
