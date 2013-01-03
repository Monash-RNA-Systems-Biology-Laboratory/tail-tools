#!/usr/bin/env pypy-bio

import sys

README = open('README','rb').read()

PAGE = """
<!--#include virtual="top.html" -->

<style>
pre { line-height: 100%%; font-size: 130%%; }
</style>

<h2>Tail tools</h2>

<h3>Download</h3>

<p>%(date)s:

<ul>
<li> <a href="%(release_tarball_name)s">%(release_tarball_name)s</a> </li>
<li><a href="https://github.com/Victorian-Bioinformatics-Consortium/tail-tools">github repository</a></li>
</ul>


<pre>
%(README)s
</pre>

<h3>Contact</h3>
<ul>
<li><a href='mailto:paul.harrison@monash.edu'>Paul Harrison</a>  
<li><a href='mailto:torsten.seemann@monash.edu'>Torsten Seemann</a>
</ul>


<!--#include virtual="bot.html" -->
"""


import os, datetime

import tail_tools

#os.environ['PATH'] = '/bio/sw/python/bin:' + os.environ['PATH']

pythonpath = os.environ['PYTHONPATH']

# RAGE
os.system('rm MANIFEST')

#assert 0 == os.system('sudo PYTHONPATH=%s pypy setup.py install --home /bio/sw/python' % pythonpath)

assert 0 == os.system('sudo /bio/sw/python/env-pypy/bin/python setup.py install')
assert 0 == os.system('sudo /bio/sw/python/env-python/bin/python setup.py install')

os.system('sudo rm -r tail_tools.egg-info')

release_tarball_name = 'tail-tools-%s.tar.gz' % tail_tools.VERSION
assert 'force' in sys.argv[1:] or not os.path.exists('dist/'+release_tarball_name), release_tarball_name + ' already exists'
date = datetime.date.today().strftime('%e %B %Y')

assert 0 == os.system('python setup.py sdist')

f = open('/home/websites/vicbioinformatics.com/software.tail-tools.shtml','wb')
f.write(PAGE % locals())
f.close()

assert 0 == os.system('cp dist/%s /home/websites/vicbioinformatics.com' % release_tarball_name)

os.system('python setup.py sdist upload')


