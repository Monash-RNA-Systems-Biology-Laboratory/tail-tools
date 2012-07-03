#!/bio/sw/python/bin/pypy

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

os.environ['PATH'] = '/bio/sw/python/bin:' + os.environ['PATH']

# RAGE
os.system('rm MANIFEST')

assert 0 == os.system('sudo -E /bio/sw/python/bin/pypy setup.py install_scripts --install-dir /bio/sw/python/bin/')
assert 0 == os.system('sudo -E /bio/sw/python/bin/pypy setup.py install_lib')
assert 0 == os.system('sudo -E /bio/sw/python/bin/python2.6 setup.py install_lib')




release_tarball_name = 'tail-tools-%s.tar.gz' % tail_tools.VERSION
#assert not os.path.exists(release_tarball_name), release_tarball_name + ' already exists'
date = datetime.date.today().strftime('%e %B %Y')

assert 0 == os.system('/bio/sw/python/bin/python2.6 setup.py sdist')

f = open('/home/torsten/public_html/vicbioinformatics.com/software.tail-tools.shtml','wb')
f.write(PAGE % locals())
f.close()

assert 0 == os.system('cp dist/%s /home/torsten/public_html/vicbioinformatics.com' % release_tarball_name)
assert 0 == os.system('cd /home/torsten/public_html/vicbioinformatics.com ; make install')
