
from nesoni import config, io
from . import web

import os

@config.help(
"Create shiny report app.","""\
Create a shiny app directory that will show a report on the given pipeline output.

This can be served with Shiny Server, or viewed with runApp() from R.
""")
@config.String_flag("title", "Report title.")
@config.String_flag("species", 'Species: Currently supports Human ("Hs"), Saccharomyces cerevisiae ("Sc"), Caenorhabditis elegans ("Ce"), Mus musculus ("Mm")')
@config.Positional("pipeline", "Directory containing output of analyze-polya-batch.")
class Shiny(config.Action_with_output_dir):
    pipeline = None
    title = "Tail Tools report"
    species = ""
    
    def run(self):
        assert self.pipeline is not None, "Pipeline output directory required."
        path = os.path.abspath(self.pipeline)
        
        workspace = io.Workspace(self.output_dir, must_exist=False)
        
        with open(workspace/"app.R","wb") as f:
            print >> f, "library(tailtools)"
            print >> f, "shiny_tailtools_report(%s, species=%s, title=%s)" % (
                repr(path), 
                "NULL" if not self.species else repr(self.species), 
                repr(self.title))
            
        with open(workspace/"index.html","wb") as f:
            web.emit(f, "sorry-no-shiny.html", {})
