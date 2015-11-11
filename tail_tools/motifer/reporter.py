
from __future__ import division

from .. import env
from . import pilers, rmonkey

import os, numpy, numpy.random, random, textwrap
from os.path import join

def require_dir(dirname):
    if not os.path.exists(dirname):
       os.mkdir(dirname)

#
#def make_local_shuffle(n, sd):
#    """ Localized permutation,
#        items will be shifted a random amount, 
#        (close to) normally distributed with standard deviation sd.
#        """
#    x = numpy.arange(n) + numpy.random.normal(scale=sd,size=n)
#    return numpy.argsort(x)
#
#
#def local_shuffle(string, sd):
#    shuffler = make_local_shuffle(len(string), sd)
#    result = ''.join([ string[i] for i in shuffler ])
#    return result


def shuffle(string,n,sd):
    if len(string) <= n: 
        return string
    
    n1 = n-1
    
    string = string
    
    bins = { }
    for i in xrange(n1,len(string)):
        context = string[i-n1:i]
        if context not in bins: bins[context] = [ ]
        bins[context].append(i)
    
    for item in bins.itervalues():
        item.sort(key=lambda i: -i+random.normalvariate(0.0,sd))
    
    context = string[:n1]
    positions = range(n1)
    resets = set()
    for i in xrange(n1,len(string)):
        while context not in bins or not bins[context]:
            resets.add(i)
            context = random.choice(list(bins))
        
        j = bins[context].pop()
        positions.append(j)
        context = (context+string[j])[1:]

    #print "Resets at", sorted(resets), "of", len(string)
    offsets = positions - numpy.arange(len(positions))
    #print numpy.sqrt(numpy.mean(offsets*offsets)), "sd, aiming for", sd

    return "".join( string[i] for i in positions )


def shuffle_all(seqs,n,sd):
    """ Local-shuffle a dictionary of sequences """
    result = { }
    for name in seqs:
        result[name] = shuffle(seqs[name],n,sd)
    return result


def html_header(f, title, show_title=True):
    f.write(
        "<!doctype html>\n"+
        '<meta charset="UTF-8">\n'+
        "<title>%s</title>\n"%title)
    write(f,"""
        <style>
        body { font-family: sans-serif; }
        a { text-decoration: none; }
        </style>
        """)
    if show_title:
        f.write("<h1>%s</h1>\n"%title)


def write(f, text):
    f.write(textwrap.dedent(text))


def p_html(p_value):
    if p_value is None: return ""
    return '<div style="display: inline-block; width: %dpx; height: 10px; background: red"></div> p=%.0e' % (int(numpy.log(max(1e-100,p_value)) / numpy.log(0.1)), p_value)


class Context(object):
    def __init__(self, ref_dir, utr_filename, gene_sets):
        all_gene_set = set()
        for name, gene_set in gene_sets:
            all_gene_set.update(gene_set)
    
        self.ref = env.load_ref(ref_dir)
        
        self.genes = list(all_gene_set)
        self.n_genes = len(self.genes)
        self.gene_number = dict((b,a) for a,b in enumerate(self.genes))
        
        self.set_weights = [ ]
        for name, genes in gene_sets:
            genes = set(genes)
            weights = numpy.zeros(self.n_genes)
            for gene in genes:
                weights[ self.gene_number[gene] ] = 1.0/len(genes)
            self.set_weights.append((name, weights))
    
        self.upstrands = [ None ]*self.n_genes
        self.coding_regions = [ None ]*self.n_genes
        self.utrs = [ None ]*self.n_genes
        self.downstrands = [ None ]*self.n_genes
        
        utr_index = env.index(utr_filename, name=lambda item: item.attr["Parent"])
        
        for i, name in enumerate(self.genes):
            self.coding_regions[i] = self.ref.coding_regions[name]
            self.utrs[i] = (
                self.coding_regions[i]
                .three_prime()
                .span_with( utr_index[name].three_prime() ) 
                )
            self.upstrands[i] = (
                self.coding_regions[i]
                .five_prime()
                .shifted(-200,0)
                )
            self.downstrands[i] = (
                self.utrs[i]
                .three_prime()
                .shifted(0,200)
                )
        
        self.pilers = [
            ( "cartoon",
              pilers.Stretched_piler(
                "-200",
                [ (10, self.upstrands, False, "Coding\nstart"),
                  (50, self.coding_regions, True, "3' UTR\nstart"),
                  (20, self.utrs, True, "3' UTR\nend"),
                  (10, self.downstrands, False, "+200"), 
                  ]
                )
              ),

            ( "utr-end",
              pilers.Anchored_piler(
                 500, 500,
                 "3' UTR\nend",
                 [ item.three_prime() for item in self.utrs ],
                 stride=10
                 )
              ),
        
            ( "utr-start",
              pilers.Anchored_piler(
                500, 500,
                "3' UTR\nstart",
                [ item.five_prime() for item in self.utrs ],
                stride=10
                )
              ),
            ]
        
        self.frames = [
            ("utr",      "3' UTR", 
              self.utrs),
            ("50-peak",  "50 bases upstrand from poly(A) site", 
              [ item.three_prime().shifted(-50,0) for item in self.utrs ]),
            ("stop-50",  "50 bases downstrand from stop codon",
              [ item.five_prime().shifted(0,50) for item in self.utrs ]),
            #("500-stop", "500 bases upstrand from stop codon",
            #  [ item.five_prime().shifted(-500,0) for item in self.utrs ]),
            #("50-start", "50 bases upstrand from start codon",
            #  [ item.five_prime().shifted(-50,0) for item in self.codings ]),
            ]
    
    @env.memo_property
    def shuffle_1(self):
        return shuffle_all(self.ref.seqs, 1, 10.0)

    @env.memo_property
    def shuffle_2(self):
        return shuffle_all(self.ref.seqs, 2, 10.0)


def frame_report(
        prefix, title, features, recognizer, context
        ):
    # Construct data frame
    import pandas
    frame = pandas.DataFrame()
    frame["id"] = context.genes
    frame["gene"] = [ context.ref.genes[item].attr.get("gene","") for item in context.genes ]
    frame["in_set"] = numpy.where(context.set_weights[0][1] > 0.0, 1, 0)
    frame["length"] = [ item.get_length() for item in features ]
    frame["count"] = [
        recognizer.count( item.get_seq(context.ref.seqs).upper() )
        for item in features
        ]
    frame.to_csv(prefix+".csv")
    
    # Do stuff in R    
    env = rmonkey.R_environment()
    env["filename"] = prefix+".csv"
    env("""
        daf <- read.csv(filename)
        daf$in_set <- daf$in_set != 0
        daf$density <- daf$count / daf$length
        
        t_fit <- NULL
        try( t_fit <- t.test(density ~ in_set, data=daf) )        
        """)
    
    [all_same_length] = env("all(daf$length == daf$length[1])")
    
    if env("is.null(t_fit)")[0]:
        p_t = 1.0
    else:
        [p_t] = env("t_fit$p.value")
        [mean_neutral, mean_different] = env("t_fit$estimate")
    
    set_name = context.set_weights[0][0]
    notset_name = context.set_weights[1][0]
    
    with open(prefix+".html","wb") as f:
        html_header(f, title)
        
        if p_t < 1.0:
            write(f, """
                <h2>Density of motif</h2>
                <p>t-test p=%(p_t)g
                <p>%(set_name)s average density: %(mean_different)g
                <p>%(notset_name)s average density: %(mean_neutral)g
                """ % locals())
    
    return p_t


def recognizer_report(
        out_dir,
        recognizer,
        context,
        
        title="Motif plot"
        ):
    import pylab
    
    require_dir(out_dir)
        
    for name, piler in context.pilers:
        pylab.figure(figsize=(20,6))
        gs = pylab.GridSpec(1,1,height_ratios=[1], left=0.075, right=0.75)
        pylab.subplot(gs[0])
        piler.setup_figure()
        
        lines = [("",3.0,"-",context.ref.seqs)]
        if recognizer.length > 2:
            lines.append((", shuffle sd=10 k=2",1.0,"--",context.shuffle_2))
        if recognizer.length > 1:
            lines.append((", shuffle sd=10 k=1",1.0,":",context.shuffle_1))
    
        for i, (set_name, set_weights) in enumerate(context.set_weights):
            a = float(i)/len(context.set_weights) * numpy.pi * 2.0
            color = (
                numpy.cos(a)*0.5+0.5,
                numpy.cos(a+numpy.pi*2/3)*0.5+0.5,
                numpy.cos(a+numpy.pi*4/3)*0.5+0.5
                )
            for suffix,width,style,seqs in lines:
                pylab.plot(
                    piler.x,
                    piler.pile(seqs, recognizer, set_weights),
                    label=set_name + suffix,
                    color=color,
                    linewidth=width,
                    linestyle=style,
                    )
        
        pylab.legend(loc="upper left", bbox_to_anchor=(1,1))                
        pylab.ylim(ymin=0.0)
        pylab.ylabel("Density")
        pylab.savefig(join(out_dir,name+".png"), dpi=100)
        pylab.savefig(join(out_dir,name+".eps"))
        pylab.close()

    with open(join(out_dir,"index.html"),"wb") as f_index:
        html_header(f_index, title)
        write(f_index,"""
            <img src="cartoon.png" width="1000">
            <br/><a href="cartoon.eps">[EPS]</a>
            <table cellpadding="0" cellspacing="0"><tr>
            <td>
                <img src="utr-start.png" width="500">
                <br/><a href="utr-start.eps">[EPS]</a>
            </td>
            <td>
                <img src="utr-end.png" width="500">
                <br/><a href="utr-end.eps">[EPS]</a>
            </td>
            </table>
            """)
        
        if len(context.set_weights) != 2:
            best_p = None
        else:
            best_p = 1.0        
            for name, frame_title, features in context.frames:
                prefix = join(out_dir,name)
                p = frame_report(prefix, title+": "+frame_title, features, recognizer, context)
                write(f_index, 
                    '<br/><a href="%s.html">' % name + 
                    frame_title + "</a> " + p_html(p) + "\n"
                    )
                
                best_p = min(p, best_p)
    
    return best_p


def make_summary(out_dir, context):
    require_dir(out_dir)
    
    env = rmonkey.R_environment()
    env("sets <- list()")

    with open(join(out_dir,"index.html"),"wb") as f:
        html_header(f, "Summary")
        
        for i, (name, weights) in enumerate(context.set_weights):
            members = weights > 0.0
            n_members = members.sum()
            lengths = numpy.array([ utr.end-utr.start for utr,member in zip(context.utrs, members) if member ])
            mean_length = numpy.mean(lengths)
            write(f,"""
                <h2>%(name)s</h2>
                <p>%(n_members)d genes.
                <p>%(mean_length).1f mean 3' UTR length.
                """ % locals())
            
            env["lengths"] = lengths
            env["n"] = i+1
            env["name"] = name
            env("sets[[n]] <- data.frame(set=rep(name,length(lengths)), lengths=lengths)")
    
        env["filename"] = join(out_dir,"dist.png") 
        env("""
           dat <- do.call(rbind, sets)
           png(filename)
           print(ggplot(dat, aes(x=set,y=lengths)) + 
               geom_violin() + theme_bw() + ylab("3' UTR length")
               )
           dev.off()
           """)
        write(f, '<p><img src="dist.png">')

        

def report(
        out_dir,
        ref_dir,
        utr_filename,
        recognizers,
        gene_sets,
        
        title="Motif plots",
        ):
    """
    
    out_dir - string
        Output directory name.
    
    ref_dir - string
        Reference directory.
    
    utr_filename - string
        GFF file containing 3' UTRs.
        Only the 3' end of these features is important.
        "Parent" attribute should be gene id.
    
    recognizers - [ heading strings and (name, Recognizer)s ]
        Motifs to scan for.
    
    gene_sets - [ (name, [ locus_tag, ... ]), ... ]
        One or more gene sets.
    
    """
    context = Context(ref_dir, utr_filename, gene_sets)
    
    require_dir(out_dir)
    require_dir(join(out_dir,"motifs"))
    
    make_summary(join(out_dir,"summary"), context)
    
    n = 0
    sidebar = [ ]
    for item in recognizers:
        if isinstance(item, str):
            sidebar.append("<h2>%s</h2>\n"%item)
        else:
            name, recognizer = item
            print name
            n += 1
            item_dir = join(out_dir,"motifs",str(n))
            p = recognizer_report(
                item_dir,
                recognizer,
                context,
                title=name
                )
            sidebar.append('<a target="box" href="motifs/%d/index.html">%s</a> %s<br/>\n' % (n, name, p_html(p)))
        
        with open(join(out_dir,"index.html"),"wb") as f:
            html_header(f, title, False)
            write(f,textwrap.dedent("""
                <style>
                body { padding: 0; margin: 0; overflow: hidden; }
                </style>
                <table cellpadding="0" cellspacing="0" style="margin:0;padding0">
                    <tr>
                        <td>
                            <div style="width: 20vw; height: 100vh; overflow-y: scroll">
                                <div style="padding: 0.5em">
                                    <p><a href="summary/index.html">Summary of gene sets</a>
                                    %s
                                </div>
                            </div>
                        </td>
                        <td>
                            <iframe name="box" style="border: 0; width: 80vw; height: 100vh;"/>
                        </td>
                    </tr>
                </table>
                """) % "".join(sidebar))
            
            
            







