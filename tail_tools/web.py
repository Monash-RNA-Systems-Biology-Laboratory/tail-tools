
import os

from nesoni import config

class Webapp(config.Action_with_prefix):
    def run(self):
        with open(self.prefix + '.html', 'wb') as f_out:
            self._emit(f_out, self.template)
    
    def _emit(self, f_out, filename):
        full_filename = os.path.join(os.path.dirname(__file__),'web',filename)
        with open(full_filename,'rU') as f_in:
            for line in f_in:
                if line.startswith('INCLUDE '):
                    self._emit(f_out, line[8:].strip())
                else:
                    f_out.write(line)

@config.help(
    'Create a viewer web page for the output for "compare-peaks:".',
    )
class Geneview_webapp(Webapp):
    template = 'geneview.html'