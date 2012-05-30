
import os, math, glob

import nesoni
from nesoni import config, working_directory, io, reporting, grace

import tail_tools
from tail_tools import clip_runs, extend_sam, proportions


class Tail_only(config.Action_filter):
    def run(self):
        fin = self.begin_input()
        fout = self.begin_output()
        
        for line in fin:
            if line.startswith('@') or '\tAA:i:' in line:
                fout.write(line)
                continue
        
        self.end_input(fin)
        self.end_output(fout)


@config.help(
'Create nesoni working directories for a single sample.',
"""\
Reads are clipped using "clip-runs:", aligned using SHRiMP, extended using "extend-sam:" and \
imported into a nesoni working directory.

Directories are created both using all reads and for reads with a poly-A tail.
""")
@config.Bool_flag('consensus', 'Look for SNPs and indels.')
@config.Positional('reference', 'Reference directory created by "nesoni make-reference:"')
@config.Positional('reads', 'Fastq file containing SOLiD reads.')
class Analyse_polya(config.Action_with_output_dir):
    reference = None
    reads = None
    
    consensus = True
    
    def run(self):
        polya_dir = os.path.normpath(self.output_dir) + '-polyA'
    
        working = working_directory.Working(self.output_dir, must_exist=False)
        working.set_reference(self.reference)
        reference = working.get_reference()
        
        polya_working = working_directory.Working(polya_dir, must_exist=False)
        polya_working.set_reference(self.reference)
        
        clipped_filename = working.object_filename('clipped_reads.csfastq')
        
        raw_filename = working/'alignments_raw.sam.gz'
        extended_filename = working/'alignments_extended.sam.gz'
        
        polya_filename = working/'alignments_filtered_polyA.sam.gz'

        if self.consensus:
            filter_tool = nesoni.Consensus
        else:
            filter_tool = nesoni.Filter
        filter_tool = filter_tool(
            monogamous=False,
            random=True,
            infidelity=0,
            userplots=False,
        )

        
        clip_runs.Clip_runs(
            self.reads,
            output=clipped_filename,
        ).make()
        
        cores = nesoni.coordinator().get_cores()
        
        nesoni.Execute(
            command = reference.shrimp_command(cs=True) + 
                      [ clipped_filename, '-N', str(cores) ],
            output=raw_filename,
            cores=cores,
        ).make()
                
        extend_sam.Extend_sam(
            input=raw_filename,
            output=extended_filename,
            reads_filename=self.reads,
            reference_filenames=[ reference.reference_fasta_filename() ],
        ).make()
        
        #Tail_only(
        #    input=extended_filename,
        #    output=polya_filename,
        #).make()
        
        

        nesoni.Import(
            input=extended_filename,
            output_dir=self.output_dir,
            reference=[ self.reference ],
        ).make()
                    
        filter_tool(
            working_dir=self.output_dir
        ).make()


        Tail_only(
            input=working/'alignments_filtered.bam',
            output=polya_filename,
        ).make()
        
        nesoni.Import(
            input=polya_filename,
            output_dir=polya_dir,
            reference=[ self.reference ],
        ).make()
        
        # This shouldn't actually filter out any alignments.
        # We do it to produce depth of coverage plots
        # and position-sorted BAM files.        
        filter_tool(
            working_dir=polya_dir
        ).make()
        
        
        
        #@nesoni.parallel_for([
        #    (extended_filename, self.output_dir),
        #    (polya_filename, polya_dir),
        #])
        #def _((sam_filename, directory)):
        #    nesoni.Import(
        #        input=sam_filename,
        #        output_dir=directory,
        #        reference=[ self.reference ],
        #    ).make()
        #                
        #    tool(
        #        working_dir=directory
        #    ).make()



@config.help(
'Analyse a set of samples and produce an HTML report, including depth of coverage plots, heatmaps, etc.',
"""\

""")
@config.String_flag('title', 'Analysis report title')
@config.String_flag('file_prefix', 'Prefix for report files')
@config.String_flag('blurb', 'Introductory HTML text for report')
@config.String_flag('genome', 'IGV .genome file, to produce IGV plots')
@config.String_flag('genome_dir', 'IGV directory of reference sequences to go with .genome file')
@config.Bool_flag('include_genome', 'Include genome in IGV plots tarball?')
@config.Int_flag('nmf', 'Perform non-negative matrix factorization up to this rank. '
                        'Defaults to the number of samples.')
@config.Positional('reference', 'Reference directory created by "nesoni make-reference:"')
@config.Main_section('reads', 'Fastq files containing SOLiD reads.')
@config.Configurable_section('analyse', 'Parameters for each "analyse-polya:"')
@config.Configurable_section('count', 'Parameters for "nesoni count:"')
@config.Section('extra_files', 'Extra files to include in report')
class Analyse_polya_batch(config.Action_with_output_dir):
    file_prefix = ''
    title = '3\'seq analysis'
    blurb = ''
    extra_files = [ ]
    
    genome = None
    genome_dir = None
    include_genome = True
    reference = None
    reads = [ ]
    
    nmf = None
    
    count = nesoni.Count(
        filter='existing',
        strand='forward',
    )
    
    analyse = Analyse_polya()

    @nesoni.stage_function
    def run(self):
        names = [
            os.path.splitext(os.path.split(item)[1])[0]
            for item in self.reads
        ]
        
        workspace = io.Workspace(self.output_dir, must_exist=False)
        plotspace = io.Workspace(workspace/'plots', must_exist=False)

        dirs = [
            workspace/item
            for item in names
        ]
        polya_dirs = [ item + '-polyA' for item in dirs ]
        
        interleaved = [ item2 for item in zip(dirs,polya_dirs) for item2 in item ]

        @nesoni.parallel_for(zip(self.reads, dirs))
        def loop((reads_filename, directory)):
            self.analyse(
                output_dir=directory,
                reference=self.reference,
                reads=reads_filename,
            ).run()

        @nesoni.parallel_for([
            ('counts', 'all', dirs),
            ('counts-polyA', 'polyA', polya_dirs),
            ('counts-both', None, dirs + polya_dirs),
        ])
        def loop((prefix, plot_name, directories)):
            prefix = workspace/prefix
            
            self.count(
                prefix = prefix,
                filenames = self.count.filenames + directories,
            ).make()
            
            nesoni.Norm_from_counts(
                prefix,
                
                tmm=False,
                # See email from Traude to Paul,
                # subject "Larger Heatmap"
                # instructions for figure 4
            ).make()  

            if plot_name is not None:
                nesoni.IGV_plots(
                    plotspace/plot_name,
                    working_dirs = directories,
                    raw = False,
                    norm = True,
                    genome = self.genome,
                    norm_file = prefix + '-norm.txt',
                    #delete_igv = False,
                ).make()

        # Make a normalization file for all
        f_in = workspace.open('counts-norm.txt', 'rb')
        f_out = workspace.open('counts-norm-all.txt', 'wb')
        f_out.write(f_in.readline())
        for line in f_in:
            parts = line.rstrip('\n').split('\t')
            print >> f_out, '\t'.join(parts)
            parts[0] += '-polyA'
            print >> f_out, '\t'.join(parts)            
        f_in.close()
        f_out.close()

        heatmaps = [ ]
        for fold, min_count in [
            (1.5, 10),
            (2, 10),
            (4, 10),
            (6, 10),
            (8, 10),
            
            (1.5, 50),
            (2, 50),
            (4, 50),
            (6, 50),
            (8, 50),
        ]:
            act = nesoni.Heatmap(
                workspace/('heatmap-%.1ffold-%dminmax'% (fold,min_count)),
                workspace/'counts.txt',
                norm_file = workspace/'counts-norm.txt',
                min_total = 0,
                min_max = min_count,
                min_span = math.log(fold)/math.log(2.0),                        
            )
            heatmaps.append(act)
            act.process_make()


        both_heatmaps = [ ]
        for deviations, min_count in [
            (1.0, 10),
            (1.0, 50),
            (1.0, 250),
            
            (2.0, 10),
            (2.0, 50),
            (2.0, 250),
        ]:
            act = nesoni.Heatmap(
                workspace/('both-heatmap-%.1fsvd-%dminmax'% (deviations,min_count)),
                workspace/'counts-both.txt',
                norm_file = workspace/'counts-both-norm.txt',
                min_total = 0,
                min_max = min_count,
                min_svd = deviations,                        
            )
            both_heatmaps.append(act)
            act.process_make()

        proportion_heatmaps = [ ]
        for min_diff in [ 0.1, 0.25, 0.5 ]:
            for min_min, min_max in [
                (10, 10),
                (10, 50),
                (50, 50),
                (50, 250),
                (250, 250),
            ]:
                act = proportions.Proportions_heatmap(
                    workspace/('proportions-heatmap-%.2fmindiff-%dminmin-%dminmax' % (min_diff, min_min, min_max)),
                    workspace/'counts.txt',
                    workspace/'counts-polyA.txt',
                    min_diff=min_diff,
                    min_min=min_min,
                    min_max=min_max,
                )
                proportion_heatmaps.append(act)
                act.process_make()


        nesoni.Stats(*self.reads, output=workspace/'stats.txt').process_make()

            
        proportions.Proportions(
            workspace/'proportions',
            workspace/'counts.txt',
            workspace/'counts-polyA.txt',
            norm_file=workspace/'counts-norm.txt'
        ).process_make()

        
        nmfs = self.nmf
        if nmfs is None:
            nmfs = len(self.reads)
        
        nmfspace = io.Workspace(workspace/'nmf', must_exist=False)
        
        order_hint = None
        nmf_actions = [ ]
        for i in xrange(2,nmfs+1):
            act = nesoni.NMF(
                nmfspace/('nmf-%d'%i),
                workspace/'counts-both.txt',
                norm_file=workspace/'counts-norm-all.txt',
                rank=i,
                order_hint = order_hint
            )
            nmf_actions.append(act)
            order_hint = act.prefix + '.rds'        

        #Need to execute in order
        @nesoni.process
        def _():
            for act in nmf_actions:
                act.make()

        nesoni.barrier() 
        #========================================================================
        
        
        #===============================================
        #                   Report        
        #===============================================

        r = reporting.Reporter(os.path.join(self.output_dir, 'report'), self.title, self.file_prefix)
        
        r.write(self.blurb)
        
        for filename in self.extra_files:
            r.p( r.get(filename) )
        
        r.heading('Alignment to reference')
        
        r.report_logs('alignment-statistics',
            [ workspace/'stats.txt' ] +
            [ os.path.join(item, 'consensus_log.txt') for item in dirs + polya_dirs ] +
            [ workspace/'counts_log.txt', workspace/'counts-polyA_log.txt' ],
            filter=lambda sample, field: (
                field not in ['fragments','fragments aligned to the reference'] and
                (not sample.endswith('-polyA') or field not in ['reads with alignments','hit multiple locations'])
            ),
        )
        
        r.heading('Normalization')
        
        r.p( 'Normalization was by total number of reads aligning to genes.' )
        
        r.p( r.get(workspace/'counts-raw.png', title='Unnormalized box and whisker plots') )
        r.p( r.get(workspace/'counts-libsize.png', title='Normalized box and whisker plots') )
        r.p( r.get(workspace/'counts-norm.txt', title='Normalization spreadsheet') )
        
        #r.heading('IGV plots')
        #
        #r.p('These files show the depth of coverage. They can be viewed with the IGV genome browser.')
        #
        #genome_files = [ ]
        #if self.include_genome:
        #    assert self.genome, '.genome file not specified.'
        #    genome_files.append(self.genome)
        #    if self.genome_dir:
        #        base = os.path.split(self.genome_dir)[1]
        #        for filename in os.listdir(self.genome_dir):
        #            genome_files.append((
        #                os.path.join(self.genome_dir, filename),
        #                os.path.join(base, filename)
        #            ))
        #
        #r.p(r.tar('igv-plots.tar.gz',
        #   genome_files +
        #   glob.glob(plotspace/'*.tdf')
        #))
        #
        #r.heading('BAM files')
        #
        #r.p('These BAM files contain the alignments of reads to the reference sequences.')
        #
        #r.p('Reads with a poly-A tail have an \'AA\' attribute.')
        #
        #bam_files = [ ]
        #for name in names:
        #    bam_files.append( (workspace/(name,'alignments_filtered_sorted.bam'),name+'.bam') )
        #    bam_files.append( (workspace/(name,'alignments_filtered_sorted.bam.bai'),name+'.bam.bai') )
        #r.p(r.tar('bam-files.tar.gz', bam_files))
        
        r.heading('Heatmaps')
        
        r.p(
            'Genes were selected on based on there being at least some number of reads '
            'in at least one of the samples (minmax), '
            'and on there being at least some fold change difference between '
            'some pair of samples (fold).'
        )
        
        for heatmap in heatmaps:
            r.report_heatmap(heatmap)
        
        r.heading('Heatmaps including poly(A) columns')
        
        r.p(
           'Below is an attempt to include the number of poly-A reads in the heatmap.'
           'There are columns for both total number of reads and for number of poly-A reads.'
        )
        r.p('Gene selection, short version: An attempt was made to select genes that '
            'represent every different pattern of expression in the data.'
        )
        r.p(
           'Gene selection, long version: Consider each gene as a point in a high-dimensional space, '
           'with coordinates being the normalized log2 number of reads from each sample, '
           'with and without poly(A) tails '
           '(2n dimensions from n samples). '
           'The distribution of points is assumed to be a multivariate normal distribution, '
           'which can be pictured as a fuzzy high-dimensional ovoid shape. '
           'Singular Value Decomposition was used to linearly transform this fuzzy ovoid '
           'into a fuzzy sphere. '
           'Points furthest from the center of the sphere were chosen.'
        )

        for heatmap in both_heatmaps:
            r.report_heatmap(heatmap)
            
        r.heading('Proportion of poly(A) reads')
        
        r.p( r.get(workspace/'proportions.csv') )
        
        r.p(
           'Genes were selected based on there being at least some number of reads '
           'in each sample (minmin), '
           'a sample with at least some number of reads (minmax), '
           'and there being a difference in the proportion of poly(A) reads '
           'between two samples of at least some amount (mindiff).'
        )
                        
        for prop in proportion_heatmaps:
            r.report_heatmap(prop)

        if nmf_actions:
            r.heading('Non-negative matrix factorization')
            
            r.p(
                'This is an experimental method of clustering genes and samples.'
            )
            
            r.p(
                'The best number of clusters to divide the genes and samples into is '
                'a matter of judgement. Choose the clustering you think is best.'
            )
            
            for act in nmf_actions:
                for filename in glob.glob(act.prefix + '*'):
                    g = r.get(filename, prefix='')
                r.p( r.get(act.prefix+'.html', title='%d clusters' % act.rank) )

        r.write('<p/><hr>\n')
        r.p('tail_tools version '+tail_tools.VERSION)
        r.p('nesoni version '+nesoni.VERSION)
        r.p('SHRiMP version '+grace.get_shrimp_2_version())
        
        r.close()




