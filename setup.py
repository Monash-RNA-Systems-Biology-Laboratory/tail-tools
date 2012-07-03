
from distutils.core import setup

import tail_tools

setup(name='tail-tools',
      version=tail_tools.VERSION,      
      packages = [ 'tail_tools' ],
      scripts = [ 'tail-tools' ],
      classifiers = [
          'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
      ],
      url = 'http://bioinformatics.net.au/software.tail-tools.shtml',
      author = 'Paul Harrison',
      author_email = 'paul.harrison@monash.edu',
)

