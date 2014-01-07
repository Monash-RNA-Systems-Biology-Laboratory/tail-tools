#!/usr/bin/env python

from os.path import join, dirname

from setuptools import setup

directory = dirname(__file__)

# Get version
with open(join(directory,'tail_tools','__init__.py'),'rU') as f:
    exec f.readline()

setup(
    name='tail-tools',
    version=VERSION,
    description='Analyse PAT-Seq RNA expression data.',
    url = 'http://bioinformatics.net.au/software.tail-tools.shtml',
    author = 'Paul Harrison',
    author_email = 'paul.harrison@monash.edu',

    packages = [ 
        'tail_tools' 
        ],
        
    package_data = {
        'tail_tools' : [
            'web/*.html',
            'web/*.css',
            'web/*.js',
            'web/third_party/*',
            ],
        },

    entry_points = {
        'console_scripts' : [
            'tail-tools = tail_tools:main',
            ],
        },
    
    install_requires = [ 
        'nesoni' 
        ],
    
    classifiers = [
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        ],
    )

