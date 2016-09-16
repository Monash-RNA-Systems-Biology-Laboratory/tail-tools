
import os

from nesoni import config

def web_filename(filename):
    return os.path.join(os.path.dirname(__file__),'web',filename)

def style():
    with open(web_filename("style.css"),"rU") as f:
        return f.read()

def emit(f_out, filename, data={}):
    full_filename = os.path.join(os.path.dirname(__file__),'web',filename)
    with open(full_filename,'rU') as f_in:
        for line in f_in:
            if line.startswith('INCLUDE '):
                emit(f_out, line[8:].strip())
            elif line.startswith('DATA '):
                f_out.write(data[line[5:].strip()])
            else:
                f_out.write(line)

class Webapp(config.Action_with_prefix):
    def run(self):
        with open(self.prefix + '.html', 'wb') as f_out:
            emit(f_out, self.template)

@config.help(
    'Create a viewer web page for the output for "compare-peaks:".',
    )
class Geneview_webapp(Webapp):
    template = 'geneview.html'