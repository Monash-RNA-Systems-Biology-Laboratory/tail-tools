
import os

from os.path import join

from nesoni import config, runr, io, selection, reporting


TEST_R = r"""

source(SOURCE)


#toptable sanitizes names
#colnames(MODEL) <- MODEL_COLUMNS
#colnames(PAIRS_MODEL) <- PAIRS_MODEL_COLUMNS

colnames(MODEL) <- sapply(basic.seq(ncol(MODEL)), function(i) sprintf("coef%d", i))
colnames(PAIRS_MODEL) <- sapply(basic.seq(ncol(PAIRS_MODEL)), function(i) sprintf("coef%d", i))





perform_test_fitnoise2 <- function(fitnoise_model, fitnoise_controls_model, name, elist, model, model_columns, n_alt, aveexpr.name) {
    library(fitnoise)

    colnames(model) <- model_columns
    
    noise.model <- model    
    stopifnot(ncol(model) <= nrow(model))
    comment <- ''
    if (ncol(model) == nrow(model)) {
        comment <- 'NOISE FITTED USING NULL MODEL\nFDR AND P-VALUE NOT TO BE TRUSTED, MERELY A WAY TO RANK GENES\n'
        noise.model <- model[, seq_len(ncol(model)-n_alt)+n_alt,drop=FALSE]
    }
    
    
    if (EMPIRICAL_CONTROLS) {
        cat('\nSelecting empirical controls\n')
        cfit <- fitnoise.fit(
            elist, 
            design=model,
            noise.design=noise.model, 
            model=fitnoise_controls_model,
            verbose=VERBOSE
            )
        cat(cfit$description,"\n")
        cfit <- fitnoise.test(cfit, coef=seq_len(n_alt))
        controls <- rank(cfit$p.value) > nrow(elist)*0.5
        cat(sum(controls),'control features chosen\n')
    } else {
        controls <- NULL
    }
    

    cat('\nFitting noise model\n')
    
    fit <- fitnoise.fit(
        elist, 
        design=model,
        noise.design=noise.model, 
        model=fitnoise_model, 
        controls=controls,
        control.design=model[,seq_len(ncol(model)-n_alt)+n_alt,drop=FALSE],
        verbose=VERBOSE
        )
    fit <- fitnoise.test(fit, coef=seq_len(n_alt))
    
    table <- data.frame(row.names=rownames(elist))
    
    table[,aveexpr.name] <- fit$averages
    
    for(i in seq_len(n_alt))
        table[, model_columns[i]] <- fit$coef[,i]
        
    table$P.Value <- fit$p.value
    table$adj.P.Val <- fit$q.value
    table$noise.P.Value <- fit$noise.p.values
    
    for(i in seq_len(ncol(elist$genes)))
        table[, colnames(elist$genes)[i]] <- elist$genes[,i]

    table <- table[order(table$P.Value),]

    write.csv(table, sprintf('%s/%s-toptable.csv',DIR,name))
    
    cat(sprintf("\n%s/%s\n",DIR,name))    
    sink(sprintf("%s/%s.txt",DIR,name), split=TRUE)
    cat(comment)
    cat(sprintf("Raw data format: %s\n",paste(colnames(elist),collapse="; ")))
    cat(sprintf("%s%s\n%d with fdr<=0.01\n%d with fdr<=0.05\n", 
        elist$info,
        fit$description,
        sum(table$adj.P.Val<0.01), 
        sum(table$adj.P.Val<0.05)))
    sink()
}










perform_test_fitnoise1 <- function(fitnoise_model, fitnoise_controls_model, name, elist, model, model_columns, n_alt, aveexpr.name) {
    colnames(model) <- model_columns
    
    noise.model <- model    
    stopifnot(ncol(model) <= nrow(model))
    comment <- ''
    if (ncol(model) == nrow(model)) {
        comment <- 'NOISE FITTED USING NULL MODEL\nFDR AND P-VALUE NOT TO BE TRUSTED, MERELY A WAY TO RANK GENES\n'
        noise.model <- model[, seq_len(ncol(model)-n_alt)+n_alt,drop=FALSE]
    }
    
    
    if (EMPIRICAL_CONTROLS) {
        cat('\nSelecting empirical controls\n')
        cfit <- fit.elist(
            elist, 
            design=model,
            noise.design=noise.model, 
            model=fitnoise_controls_model, 
            cores=min(detectCores(),8) #TODO: make this configurable
            )
        cat(cfit$noise.description,"\n")
        ctable <- test.fit(cfit, coefs=seq_len(n_alt), sort=F)
        controls <- rank(ctable$P.Value) > nrow(elist)*0.5
        cat(sum(controls),'control features chosen\n')
    } else {
        controls <- NULL
    }
    

    cat('\nFitting noise model\n')
    
    fit <- fit.elist(
        elist, 
        design=model,
        noise.design=noise.model, 
        model=fitnoise_model, 
        controls=controls,
        control.design=model[,seq_len(ncol(model)-n_alt)+n_alt,drop=FALSE],
        cores=min(detectCores(),8) #TODO: make this configurable
        )
    table <- test.fit(fit, coefs=seq_len(n_alt))

    for(i in basic.seq(ncol(table)))
        if (colnames(table)[i] == 'AveExpr') {
            colnames(table)[i] <- aveexpr.name
            break
        }

    write.csv(table, sprintf('%s/%s-toptable.csv',DIR,name))
    
    cat(sprintf("\n%s/%s\n",DIR,name))    
    sink(sprintf("%s/%s.txt",DIR,name), split=TRUE)
    cat(comment)
    cat(sprintf("Raw data format: %s\n",paste(colnames(elist),collapse="; ")))
    cat(sprintf("%s%s\n%d with fdr<=0.01\n%d with fdr<=0.05\n", 
        elist$info,
        fit$noise.description,
        sum(table$adj.P.Val<0.01), 
        sum(table$adj.P.Val<0.05)))
    sink()
}






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
    
    if (METHOD != "fitnoise2")
        tails[is.na(tails)] <- 0
        # fitnoise2 actually requires zero count tails to be NA

    dgelist <- read.counts(counts_filename, norm.file=norm_filename, quiet=TRUE)
    dgelist <- dgelist[,select]
    
    stopifnot(all(rownames(counts) == rownames(dgelist$genes)))

    genes <- dgelist$genes
    genes <- genes[,colnames(genes) %in% c('locus_tag','Length','gene','product')]
    genes[,'reads'] = mapply(
        function(i) paste(counts[i,],collapse='; '),
        seq_len(nrow(genes)))
    genes[,'polya.reads'] = mapply(
        function(i) paste(tail.counts[i,],collapse='; '),
        seq_len(nrow(genes)))
    genes[,'tail.lengths'] = mapply(
        function(i) paste(sprintf('%.1f',tails[i,]),collapse='; '),
        seq_len(nrow(genes)))
    dgelist$genes <- genes


    if (WEIGHT) {
        fitnoise1_voom_model <- model.t.per.sample.var
        fitnoise1_controls_voom_model <- model.t.standard
        fitnoise1_patseq_model <- model.t.patseq.per.sample.var
        fitnoise1_controls_patseq_model <- model.t.patseq
        
        fitnoise2_voom_model <- "Model_t_per_sample()"
        fitnoise2_controls_voom_model <- "Model_t()"
        fitnoise2_patseq_model <- "Model_t_patseq_per_sample()"
        fitnoise2_controls_patseq_model <- "Model_t_patseq()"
    } else {
        fitnoise1_voom_model <- model.t.standard
        fitnoise1_controls_voom_model <- model.t.standard
        fitnoise1_patseq_model <- model.t.patseq
        fitnoise1_controls_patseq_model <- model.t.patseq
        
        fitnoise2_voom_model <- "Model_t()"
        fitnoise2_controls_voom_model <- "Model_t()"
        fitnoise2_patseq_model <- "Model_t_patseq()"
        fitnoise2_controls_patseq_model <- "Model_t_patseq()"
    }
    noise.model <- model
    stopifnot(ncol(model) <= nrow(model))
    if (ncol(model) == nrow(model)) {
        noise.model <- model[, seq_len(ncol(model)-n_alt)+n_alt,drop=FALSE]
        fitnoise1_voom_model <- model.normal.standard
        fitnoise1_controls_voom_model <- model.normal.standard
        fitnoise1_patseq_model <- model.normal.patseq
        fitnoise1_controls_patseq_model <- model.normal.patseq
        
        fitnoise2_voom_model <- "Model_normal()"
        fitnoise2_controls_voom_model <- "Model_normal()"
        fitnoise2_patseq_model <- "Model_normal_patseq()"
        fitnoise2_controls_patseq_model <- "Model_normal_patseq()"
        
        if (EMPIRICAL_CONTROLS) stop("Can't find empricical controls without any degrees of freedom.")
    }

    good <- row.apply(dgelist$counts, max) >= MIN_READS
    png(sprintf("%s/%s-voom.png",DIR,name))
    if (WEIGHT && METHOD == "limma")
        voomed <- voomWithQualityWeights(dgelist[good,], noise.model, plot=T)    
    else
        voomed <- voom(dgelist[good,], noise.model, plot=T)        
    dev.off()
    voomed$info <- sprintf(
        paste(
            '%d of %d features kept after filtering\n',
            '(required at least one sample with %d reads)\n',
            sep=''), 
        sum(good),length(good),MIN_READS)
    
    if (METHOD == "fitnoise2") {    
        perform_test_fitnoise2(fitnoise2_voom_model, fitnoise2_controls_voom_model, sprintf("%s-voom",name), voomed, model, model_columns, n_alt, 'avg.expression')        
    } else if (METHOD == "fitnoise1") {
        perform_test_fitnoise1(fitnoise1_voom_model, fitnoise1_controls_voom_model, sprintf("%s-voom",name), voomed, model, model_columns, n_alt, 'avg.expression')
    } else {
        perform_test(sprintf("%s-voom",name), voomed, model, model_columns, n_alt, 'avg.expression')
    }

    if (METHOD == "fitnoise2") {    
        tail.elist <- elist.tails.for.fitnoise(tails, tail.counts, model, genes, MIN_READS)
        perform_test_fitnoise2(fitnoise2_patseq_model, fitnoise2_controls_patseq_model, sprintf("%s-tail",name), tail.elist, model, model_columns, n_alt, 'avg.tail')    
    } else if (METHOD == "fitnoise1") {
        tail.elist <- elist.tails.for.fitnoise(tails, tail.counts, model, genes, MIN_READS)
        perform_test_fitnoise1(fitnoise1_patseq_model, fitnoise1_controls_patseq_model, sprintf("%s-tail",name), tail.elist, model, model_columns, n_alt, 'avg.tail')    
    } else {
        tail.elist <- elist.tails(tails, tail.counts, model, genes, MIN_READS)
        perform_test(sprintf("%s-tail",name), tail.elist, model, model_columns, n_alt, 'avg.tail')
    }
}

perform_tests('genewise', GENEWISE_FILENAME, GENEWISE_NORM_FILENAME, SELECT, MODEL, MODEL_COLUMNS, N_ALT)
perform_tests('primarypeakwise', PRIMARYPEAKWISE_FILENAME, PRIMARYPEAKWISE_NORM_FILENAME, SELECT, MODEL, MODEL_COLUMNS, N_ALT)
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
#@config.Bool_flag('fitnoise', 'Use experimental fitnoise method.')
@config.String_flag('method', 
    'Statistical method to use:\n'
    'limma - Original method\n'
    'fitnoise1 - First (R) version of Fitnoise\n'
    'fitnoise2 - Second (Python) version of Fitnoise'
    )
@config.Bool_flag('weight', 
    'Use sample per-sample quality weights. '
    'With --fitnoise no, only expression levels are weighted using voomWithQualityWeights. With --fitnoise yes, tail lengths are also weighted.')
@config.Bool_flag('empirical_controls',
    'Fitnoise only: Select and use empirical control features, may improve sample weighting.')
@config.Int_flag('min_reads', 
    'For expression testing, at least one sample must have this many reads. '
    'For tail length testing, sufficient samples must have this many reads to fit the linear model.'
    )
@config.Bool_flag('dedup', 'Use deduplicated counts')
@config.Bool_flag('verbose', 'Extra verbosity when fitting with fitnoise2.')
@config.Positional('analysis', 'Output directory of "analyse-polya-batch:".')
@config.Section('null', 'Terms in null hypothesis (H0).')
@config.Section('alt', 'Additional terms in alternative hypothesis (H1).')
class Test(config.Action_with_output_dir):
   method = "fitnoise2"
   weight = False
   empirical_controls = False
   min_reads = 10
   tell = False
   dedup = False
   title = ''
   verbose = False
   
   analysis = None
   null = [ ]
   alt = [ ]
   
   def get_title(self):
       title = self.title
       if not title:
           title = ', '.join(filter(selection.term_name,self.alt)) + ' in ' + ', '.join(filter(selection.term_name,self.null))
       if self.dedup:
          title = '[dedup] ' + title
       return title
   
   def run(self):
       assert self.method in ("limma", "fitnoise1", "fitnoise2"), "Unknown method."
       assert self.method != "limma" or not self.empirical_controls
       
       title = self.get_title()
   
       n_alt = len(self.alt)
       n_null = len(self.null)
       
       suffix = '-dedup' if self.dedup else ''
   
       genewise_filename = join(self.analysis,'expression','genewise'+suffix,'counts.csv')
       genewise_norm_filename = join(self.analysis,'expression','genewise'+suffix,'norm.csv')

       primarypeakwise_filename = join(self.analysis,'expression','primarypeakwise'+suffix,'counts.csv')
       primarypeakwise_norm_filename = join(self.analysis,'expression','primarypeakwise'+suffix,'norm.csv')

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
           spec = selection.term_specification(term)
           model.append([ selection.weight(spec, tags[item]) for item in samples ])
       model = zip(*model) #Transpose
       
       select = [ any(row) for row in model ]
       model = [ row for row,selected in zip(model,select) if selected ]
       model_columns = [ selection.term_name(item) for item in self.alt + self.null ]
       model_rows = [ item for keep, item in zip(select, samples) if keep ]
       
       #degust complains if name starts with '-', delimits with commas
       model_columns = [ ('.' if item[:1] == '-' else '') + item.replace(',',';') for item in model_columns ]
       
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
       pairs_model_rows = [ item+'-peak1' for item in model_rows ] + [ item+'-peak2' for item in model_rows ]
       
       print
       print 'Design matrix'
       print '['+('-'*(8*n_alt-2))+'] test coefficients'
       for row, name in zip(model, model_rows):
           print ''.join('%7g ' % item for item in row), name
       print
       print 'Pair design matrix'
       print '['+('-'*(8*n_alt-2))+'] test coefficients'
       for row, name in zip(pairs_model, pairs_model_rows):
           print ''.join('%7g ' % item for item in row), name
       print
       
       
       workspace = self.get_workspace()
       
       runr.run_script(TEST_R, self.tell,
           SOURCE = os.path.join(os.path.dirname(__file__),'tail_tools.R'),
           DIR = workspace.working_dir,
           METHOD = self.method,
           WEIGHT = self.weight,
           EMPIRICAL_CONTROLS = self.empirical_controls,
           MIN_READS = self.min_reads,
           VERBOSE = self.verbose,
           
           GENEWISE_FILENAME = genewise_filename,
           GENEWISE_NORM_FILENAME = genewise_norm_filename,
           PRIMARYPEAKWISE_FILENAME = primarypeakwise_filename,
           PRIMARYPEAKWISE_NORM_FILENAME = primarypeakwise_norm_filename,
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
       
       reporter.write('<table>\n')
       for entities, result, aveexpr, subtitle, terms in [
           ('genes', 'genewise-voom', 'avg.expression', 'Genewise expression level', model_columns[:n_alt]),
           ('genes', 'genewise-tail', 'avg.tail', 'Genewise tail length', model_columns[:n_alt]),
           ('primary peaks', 'primarypeakwise-voom', 'avg.expression', 'Primary-peakwise expression level', model_columns[:n_alt]),
           ('primary peaks', 'primarypeakwise-tail', 'avg.tail', 'Primary-peakwise tail length', model_columns[:n_alt]),
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

           with open(workspace/(result+'.txt'),'rU') as f:
               lines = f.readlines()
           
           reporter.write('<tr><td valign="top" width="33%">')
           reporter.subheading( reporter.href(workspace/(result+'.html'), subtitle) )
           #reporter.p( '%d %s, %d with fdr&lt;=0.01, %d with fdr&lt;=0.05' % (n,entities,n_01,n_05) )
           line = reporter.href(workspace/(result+'-toptable.csv'), 'Spreadsheet')
           if result.endswith('voom'):
               line += ', ' + reporter.href(workspace/(result+'.png'), 'voom plot')
           reporter.p(line)
           for line in lines[-2:]:
               reporter.p(line.strip())
           reporter.write('</td><td valign="top"><br/><br/>')
           for line in lines[:-2]:
               reporter.write(line.strip() + '<br/>\n')
           reporter.write('</td></tr>')

       reporter.write('</table>\n')
        
       reporter.close()



