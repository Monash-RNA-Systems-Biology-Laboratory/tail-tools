
from distutils.core import setup

import tail_tools

setup(name='tail-tools',
      version=tail_tools.VERSION,      
      packages = [ 'tail_tools' ],
      scripts = [ 'tail-tools' ],
)

