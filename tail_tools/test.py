
import os

from os.path import join

from nesoni import config, runr, io, selection, reporting

def term_specification(term):
    if '=' not in term: return term
    return term.split('=',1)[0]

def term_name(term):
    if '=' not in term: return term
    return term.split('=',1)[1]


TEST_R = r"""

source(SOURCE)


#toptable sanitizes names
#colnames(MODEL) <- MODEL_COLUMNS
#colnames(PAIRS_MODEL) <- PAIRS_MODEL_COLUMNS

colnames(MODEL) <- sapply(basic.seq(ncol(MODEL)), function(i) sprintf("coef%d", i))
colnames(PAIRS_MODEL) <- sapply(basic.seq(ncol(PAIRS_MODEL)), function(i) sprintf("coef%d", i))


perform_test <- function(name, elist, model, model_columns, n_alt, aveexpr.name) {
    fit <- lmFit(elist, model)
    fit <- eBayes(fit[,basic.seq(n_alt)])
    table <- topTableF(fit, number=nrow(fit))    

    for(i in basic.seq(ncol(table)))
        for(j in basic.seq(n_alt))
            if (colnames(table)[i] == sprintf("coef%d",j)) {
                colnames(table)[i] <- model_columns[j]
                break
            }

    for(i in basic.seq(ncol(table)))
        if (colnames(table)[i] == 'AveExpr') {
            colnames(table)[i] <- aveexpr.name
            break
        }
        

    write.csv(table, sprintf('%s/%s-toptable.csv',DIR,name))
    
    cat(sprintf("\n%s/%s\n",DIR,name))
    
    sink(sprintf("%s/%s.txt",DIR,name), split=TRUE)
    cat(sprintf("Raw data format: %s\n",paste(colnames(elist),collapse="; ")))
    if (!is.null(elist$info))
        cat(elist$info)
    cat(sprintf("%.1f degrees of freedom in prior\n", fit$df.prior))
    cat(sprintf("%d with fdr<=0.01\n%d with fdr<=0.05\n", 
        sum(table$adj.P.Val<0.01), 
        sum(table$adj.P.Val<0.05)))
    sink()
}

perform_tests <- function(name, counts_filename, norm_filename, select, model, model_columns, n_alt) {
    table <- read.grouped.table(counts_filename)
    counts <- as.matrix(table$Count)[,select]
    tail.counts <- as.matrix(table$Tail_count)[,select]
    tails <- as.matrix(table$Tail)[,select]
    tails[is.na(tails)] <- 0

    dgelist <- read.counts(counts_filename, norm.file=norm_filename, quiet=TRUE)
    dgelist <- dgelist[,select]

    genes <- dgelist$genes
    genes <- genes[,colnames(genes) %in% c('locus_tag','Length','gene','product')]
    genes[,'reads'] = mapply(
        function(genename) paste(counts[genename,],collapse='; '),
        rownames(genes))
    genes[,'polya.reads'] = mapply(
        function(genename) paste(tail.counts[genename,],collapse='; '),
        rownames(genes))
    genes[,'tail.lengths'] = mapply(
        function(genename) paste(sprintf('%.1f',tails[genename,]),collapse='; '),
        rownames(genes))
    dgelist$genes <- genes


    voomed <- voom(dgelist, model)    
    good <- row.apply(dgelist$counts, max) >= MIN_READS
    voomed <- voomed[good,]
    voomed$info <- sprintf(
        paste(
            '%d of %d features kept after filtering\n',
            '(required at least one sample with %d reads)\n',
            sep=''), 
        sum(good),length(good),MIN_READS)
    perform_test(sprintf("%s-voom",name), voomed, model, model_columns, n_alt, 'avg.expression')
    

    tail.elist <- elist.tails(tails, tail.counts, model, genes, MIN_READS)
    perform_test(sprintf("%s-tail",name), tail.elist, model, model_columns, n_alt, 'avg.tail')
}

perform_tests('genewise', GENEWISE_FILENAME, GENEWISE_NORM_FILENAME, SELECT, MODEL, MODEL_COLUMNS, N_ALT)
perform_tests('peakwise', PEAKWISE_FILENAME, PEAKWISE_NORM_FILENAME, SELECT, MODEL, MODEL_COLUMNS, N_ALT)
perform_tests('pairwise', PAIRWISE_FILENAME, PAIRWISE_NORM_FILENAME, PAIRS_SELECT, PAIRS_MODEL, PAIRS_MODEL_COLUMNS, PAIRS_N_ALT)

"""


@config.help(
    'Perform a set of differential expression and differential tail length tests which can be viewed using Degust.',
    'The null: and alt: sections give terms in the linear model, expressed as a selection expression (see nesoni main help text).'
    '\n\n'
    'For example if you have two conditions A and B, and have tagged your samples appropriately, you might use:'
    '\n\n'
    '  tail-tools test: result analysis-dir null: A/B alt: B'
    '\n\n'
    'The output is an HTML report containing a set of Degust pages.'
    '\n\n'
    'degust.py must be installed on the path to use this tool.'
    )
@config.String_flag('title', 'Report title.')
@config.Bool_flag('tell', 'Show R+ code instead of executing it.')
@config.Int_flag('min_reads', 
    'For expression testing, at least one sample must have this many reads. '
    'For tail length testing, sufficient samples must have this many reads to fit the linear model.'
    )
@config.Bool_flag('dedup', 'Use deduplicated counts')
@config.Positional('analysis', 'Output directory of "analyse-polya-batch:".')
@config.Section('null', 'Terms in null hypothesis (H0).')
@config.Section('alt', 'Additional terms in alternative hypothesis (H1).')
class Test(config.Action_with_output_dir):
   min_reads = 10
   tell = False
   dedup = False
   title = ''
   
   analysis = None
   null = [ ]
   alt = [ ]
   
   def get_title(self):
       title = self.title
       if not title:
           title = ', '.join(filter(term_name,self.alt)) + ' in ' + ', '.join(filter(term_name,self.null))
       if self.dedup:
          title = '[dedup] ' + title
       return title
   
   def run(self):
       title = self.get_title()
   
       n_alt = len(self.alt)
       n_null = len(self.null)
       
       suffix = '-dedup' if self.dedup else ''
   
       genewise_filename = join(self.analysis,'expression','genewise'+suffix,'counts.csv')
       genewise_norm_filename = join(self.analysis,'expression','genewise'+suffix,'norm.csv')

       peakwise_filename = join(self.analysis,'expression','peakwise'+suffix,'counts.csv')
       peakwise_norm_filename = join(self.analysis,'expression','peakwise'+suffix,'norm.csv')

       pairwise_filename = join(self.analysis,'peak-shift'+suffix,'individual-pairs.csv')
       pairwise_norm_filename = join(self.analysis,'peak-shift'+suffix,'individual-pairs-norm.csv')

   
       reader = io.Table_reader(genewise_filename, 'Count')
       reader.close()
       samples = [ item for i, item in enumerate(reader.headings) if reader.groups[i] == 'Count' ]
       tags = { }
       for item in samples:
           tags[item] = [ item ]
       for line in reader.comments:
           if line.startswith('#sampleTags='):
               parts = line[len('#sampleTags='):].split(',')
               tags[parts[0]] = parts
              
       model = [ ]
       for term in self.alt + self.null:        
           spec = term_specification(term)
           model.append([ 1 if selection.matches(spec, tags[item]) else 0 for item in samples ])
       model = zip(*model) #Transpose
       
       select = [ any(row) for row in model ]
       model = [ row for row,selected in zip(model,select) if selected ]
       model_columns = [ term_name(item) for item in self.alt + self.null ]
       
       pairs_n_alt = n_alt       
       pairs_select = select + select
       pairs_model = (
           [ (0,) * n_alt + row + (0,) for row in model ] +
           [ row[:n_alt]  + row + (1,) for row in model ] 
           )
       pairs_model_columns = (
           [ item+'-interaction' for item in model_columns[:n_alt] ] +
           model_columns +
           [ 'pair2' ]
           )
       
       workspace = self.get_workspace()
       
       runr.run_script(TEST_R, self.tell,
           SOURCE = os.path.join(os.path.dirname(__file__),'tail_tools.R'),
           DIR = workspace.working_dir,
           MIN_READS = self.min_reads,
           GENEWISE_FILENAME = genewise_filename,
           GENEWISE_NORM_FILENAME = genewise_norm_filename,
           PEAKWISE_FILENAME = peakwise_filename,
           PEAKWISE_NORM_FILENAME = peakwise_norm_filename,
           PAIRWISE_FILENAME = pairwise_filename,
           PAIRWISE_NORM_FILENAME = pairwise_norm_filename,
           
           N_ALT = n_alt,
           SELECT = select,
           MODEL = model,
           MODEL_COLUMNS = model_columns,
           PAIRS_N_ALT = pairs_n_alt,
           PAIRS_SELECT = pairs_select,
           PAIRS_MODEL = pairs_model,
           PAIRS_MODEL_COLUMNS = pairs_model_columns,
           )
       if self.tell: return
       
       reporter = reporting.Reporter(workspace.working_dir, title)
       
       if self.dedup:
           reporter.p('Read deduplication was used.')
       
       for entities, result, aveexpr, subtitle, terms in [
           ('genes', 'genewise-voom', 'avg.expression', 'Genewise expression level', model_columns[:n_alt]),
           ('genes', 'genewise-tail', 'avg.tail', 'Genewise tail length', model_columns[:n_alt]),
           ('peaks', 'peakwise-voom', 'avg.expression', 'Peakwise expression level', model_columns[:n_alt]),
           ('peaks', 'peakwise-tail', 'avg.tail', 'Peakwise tail length', model_columns[:n_alt]),
           ('peak pairs', 'pairwise-voom', 'avg.expression', 'Peak-pair expression shift', pairs_model_columns[:n_alt]),
           ('peak pairs', 'pairwise-tail', 'avg.tail', 'Peak-pair tail length shift', pairs_model_columns[:n_alt]),
           ]:
           #data = io.read_grouped_table(workspace/(result+'-toptable.csv'))['All']
           #n = 0
           #n_01 = 0
           #n_05 = 0
           #for row in data.values():
           #    fdr = float(row['adj.P.Val'])
           #    if fdr <= 0.01: n_01 += 1
           #    if fdr <= 0.05: n_05 += 1
           #    n += 1
           
           io.execute([
               'degust.py',
               '--name', title + ' : ' + subtitle,
               '--avg', aveexpr,
               '--primary', 'baseline',
               '--logFC', ','.join(terms),
               '--fdr', 'adj.P.Val',
               '--info', 'gene,locus_tag,product,reads,polya.reads,tail.lengths,'+aveexpr,
               '--notour', '1',
               '--out', workspace/(result+'.html'),
               workspace/(result+'-toptable.csv'),
               ])
            
           reporter.subheading( reporter.href(workspace/(result+'.html'), subtitle) )
           #reporter.p( '%d %s, %d with fdr&lt;=0.01, %d with fdr&lt;=0.05' % (n,entities,n_01,n_05) )
           with open(workspace/(result+'.txt'),'rU') as f:
               for line in f:
                   reporter.write(line.strip() + '<br/>\n')
        
       reporter.close()



