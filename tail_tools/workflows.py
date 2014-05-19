
import os, math, glob
from os.path import join

import nesoni
from nesoni import config, workspace, working_directory, reference_directory, io, reporting, grace, annotation

import tail_tools
from . import clip_runs, extend_sam, proportions, tail_lengths, web, alternative_tails


@config.help(
    'Call peaks in depth of coverage indicating transcription end sites.',
    'Note that you must specify --select (the annotation type to use from the annotation file), '
    '--shift-start, and --shift-end.'
    '\n\n'
    'Operates as follows:'
    '\n\n'
    '- Call end points using "nesoni modes: --what 3prime".\n'
    '- Extend called modes back by --peak-length.\n'
    '- Collapse any overlapping annotations down to a single annotation.\n'
    '- Relate resultant peaks to the given (collapsed) annotations.\n'
    '\n'
    'Note: As this is a workflow, you may need to specify "--make-do all" to force everything to recompute if an input file is changed. '
    'Use "--make-do -modes" to recompute everything but the peak calling.'
    )
@config.Int_flag('lap', '--lap value for "nesoni modes:". How fuzzy the pileup of 3\' ends can be when calling a peak.')
@config.Int_flag('radius', '--radius value for "nesoni modes:". How close peaks can be to one another.')
@config.Int_flag('min_depth', '--min-depth value for "nesoni modes:".')
@config.Int_flag('peak_length', 'Number of bases to extend peak back from 3\' end point.')
@config.String_flag('annotations', 'Annotation file.')
@config.String_flag('types', 'What feature types from annotation file to relate peaks to (comma separated list).')
@config.Int_flag('shift_start', 'How far downstrand of the start of features do peaks need to end to be related to a feature.')
@config.Int_flag('shift_end', 'How far downstrand of the end of the feature can a peak start in order to be related to a feature.')
@config.Main_section('polyas', 'List of ...-polyA directories as produced by "analyse-polya:" or "analyse-polya-batch:".')
class Call_peaks(config.Action_with_output_dir):
    lap = 10
    radius = 50
    min_depth = 50
    peak_length = 100
    
    annotations = None
    types = None # Must specify
    shift_start = None # Must specify
    shift_end = None
    
    
    polyas = [ ]
    
    def run(self):
        assert self.shift_start is not None, '--shift-start must be specified'
        assert self.shift_end is not None, '--shift-end must be specified'
        assert self.types is not None, '--select must be specified'
        assert self.annotations is not None, '--annotations must be specified'
    
        outspace = self.get_workspace()
        working = workspace.Workspace(outspace / 'working', must_exist=False)
        
        nesoni.Modes(
            working/'modes',
            filenames = self.polyas,
            what = '3prime',
            lap = self.lap,
            radius = self.radius,
            min_depth = self.min_depth,
            ).make()
        
        nesoni.Modify_features(
            working/'peaks',
            working/'modes.gff',
            shift_start = str(-self.peak_length),
            ).make()
        
        nesoni.Collapse_features(
            working/'collapsed',
            self.annotations,
            select = '/'.join(self.types.split(',')),
            ).make()
        
        nesoni.Relate_features(
            outspace/'relation',
            parent = working/'collapsed.gff',
            child = working/'peaks.gff',
            upstrand = -self.shift_start,
            downstrand = self.shift_end,
            use = 'in/upstrand/downstrand',
            to_child = 'Name/Product',
            ).make()
        



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
@config.Section('tags', 'Tags for this sample. (See "nesoni tag:".)')
@config.Positional('reference', 'Reference directory created by "nesoni make-reference:"')
@config.Section('reads', 'Fastq files containing SOLiD reads.')
class Analyse_polya(config.Action_with_output_dir):
    reference = None
    tags = [ ]
    reads = [ ]    
    consensus = False
    
    _workspace_class = working_directory.Working
    
    def get_polya_dir(self):
        return os.path.normpath(self.output_dir) + '-polyA'
    
    def get_filter_tool(self):
        if self.consensus:
            filter_tool = nesoni.Consensus
        else:
            filter_tool = nesoni.Filter
        return filter_tool(
            monogamous=False,
            random=True,
            infidelity=0,
            userplots=False,
        )
    
    def get_filter_action(self):
        return self.get_filter_tool()(working_dir = self.output_dir)

    def get_polya_filter_action(self):
        return self.get_filter_tool()(working_dir = self.get_polya_dir())
    
    def run(self):
        assert self.reads, 'No read files given.'
        colorspace = [ io.is_colorspace(item) for item in self.reads ]
        assert len(set(colorspace)) == 1, 'Mixture of colorspace and basespace reads is not currently supported.'
        colorspace = colorspace[0]
        
        polya_dir = self.get_polya_dir()
    
        working = working_directory.Working(self.output_dir, must_exist=False)
        working.set_reference(self.reference)
        reference = working.get_reference()
        
        polya_working = working_directory.Working(polya_dir, must_exist=False)
        polya_working.set_reference(self.reference)
        
        clipped_prefix = working/'clipped_reads'
        clipped_filename = clipped_prefix+('.csfastq.gz' if colorspace else '.fastq.gz')
        
        raw_filename = working/'alignments_raw.sam.gz'
        extended_filename = working/'alignments_extended.sam.gz'
        
        polya_filename = working/'alignments_filtered_polyA.sam.gz'

        if colorspace:
            clip_runs.Clip_runs_colorspace(
                filenames=self.reads,
                prefix=clipped_prefix,
                sample=working.name,
            ).make()
        else:
            clip_runs.Clip_runs_basespace(
                filenames=self.reads,
                prefix=clipped_prefix,
                sample=working.name,
            ).make()        

        cores = min(nesoni.coordinator().get_cores(), 8)
        
        if colorspace:
            nesoni.Execute(
                command = reference.shrimp_command(cs=colorspace, parameters=[ clipped_filename ]) + [ '--qv-offset', '33' ],
                execution_options = [ '-N', str(cores) ],
                output=raw_filename,
                cores=cores,
                ).make()
        
        else:
            nesoni.Execute(
                command = [ 
                    'bowtie2', 
                    '--rg-id', '1',
                    '--rg', 'SM:'+working.name,
                    '-k', '10', #Up to 10 alignments per read
                    '-x', reference.get_bowtie_index_prefix(),
                    '-U', clipped_filename,
                    ],
                execution_options = [ '--threads', str(cores) ],
                output=raw_filename,
                cores=cores,
                ).make()
                
        if colorspace:
            extend_sam.Extend_sam_colorspace(
                input=raw_filename,
                output=extended_filename,
                reads=self.reads,
                reference_filenames=[ reference.reference_fasta_filename() ],
            ).make()
        else:    
            extend_sam.Extend_sam_basespace(
                input=raw_filename,
                output=extended_filename,
                clips=[ clipped_prefix+'.clips.gz' ],
                reference_filenames=[ reference.reference_fasta_filename() ],
            ).make()
        
        nesoni.Import(
            input=extended_filename,
            output_dir=self.output_dir,
            reference=[ self.reference ],
        ).make()
        
        self.get_filter_action().make()

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
        self.get_polya_filter_action().make()
        
        nesoni.Tag(self.output_dir, tags=self.tags).make()
        nesoni.Tag(polya_dir, tags=self.tags).make()



@config.help(
'Analyse a set of samples and produce an HTML report, including depth of coverage plots, heatmaps, etc.',
"""\

""")
@config.String_flag('title', 'Analysis report title.')
@config.String_flag('file_prefix', 'Prefix for report filenames.')
@config.Int_flag('peak_min_depth', 
    'Number of poly(A) reads ending at nearly the same position required in order to call a peak.'
    )
@config.Int_flag('extension', 'How far downstrand of the given annotations a read or peak belonging to a gene might be.')
#@config.String_flag('blurb', 'Introductory HTML text for report')
#@config.String_flag('genome', 'IGV .genome file, to produce IGV plots')
#@config.String_flag('genome_dir', 'IGV directory of reference sequences to go with .genome file')
@config.Bool_flag('include_plots', 'Include plots in report?')
@config.Bool_flag('include_genome', 'Include genome in IGV plots tarball?')
@config.Bool_flag('include_bams', 'Include BAM files in report?')
@config.Positional('reference', 'Reference directory created by "nesoni make-reference:"')
@config.Configurable_section('template', 
    'Common options for each "sample:". '
    'Put this section before the actual "sample:" sections.',
    presets = [ ('default', lambda obj: Analyse_polya(), 'default "analyse-polya:" options') ],
    )
@config.Configurable_section_list('samples',
    'Samples for analysis. Give one "sample:" section for each sample.',
    templates = [ ],
    sections = [ ('sample', lambda obj: obj.template, 'A sample, parameters as per "analyse-polya:".') ],
    )
#@config.Configurable_section('analyse_tail_lengths', 'Parameters for "analyse-tail-lengths:"')
#@config.Section('extra_files', 'Extra files to include in report')
@config.Section('groups', 'Sample groups, given as a list of <selection>=<name>.')
@config.Configurable_section_list('tests',
    'Differential tests to perform.',
    templates = [ ],
    sections = [ ('test', lambda obj: tail_tools.Test(), 'A differential test, parameters as per "test:".') ],
    )
class Analyse_polya_batch(config.Action_with_output_dir):
    file_prefix = ''
    title = 'PAT-Seq analysis'
    #blurb = ''
    #extra_files = [ ]
    
    extension = 1000
    
    peak_min_depth = 50
    
    #genome = None
    #genome_dir = None
    include_plots = True
    include_genome = False
    include_bams = True
    reference = None
    samples = [ ]
    
    groups = [ ]
    
    tests = [ ]
    
    #analyse_tail_lengths = tail_lengths.Analyse_tail_lengths()
    
    def run(self):
        names = [ sample.output_dir for sample in self.samples ]
            #os.path.splitext(os.path.split(item)[1])[0]
            #for item in self.reads
            #]
        
        reference = reference_directory.Reference(self.reference, must_exist=True)
        
        workspace = io.Workspace(self.output_dir, must_exist=False)
        samplespace = io.Workspace(workspace/'samples', must_exist=False)
        plotspace = io.Workspace(workspace/'plots', must_exist=False)
        expressionspace = io.Workspace(workspace/'expression', must_exist=False)
        testspace = io.Workspace(workspace/'test', must_exist=False)
                
        file_prefix = self.file_prefix
        if file_prefix and not file_prefix.endswith('-'):
            file_prefix += '-'


        #dirs = [
        #    workspace/item
        #    for item in names
        #]

        samples = [ ]
        for sample in self.samples:
            samples.append(sample(
                samplespace / sample.output_dir,
                reference = self.reference,
                ))
        
        dirs = [ item.output_dir for item in samples ]
        polya_dirs = [ item + '-polyA' for item in dirs ]        
        interleaved = [ item2 for item in zip(dirs,polya_dirs) for item2 in item ]
        
        clipper_logs = [ join(item.output_dir, 'clipped_reads_log.txt') for item in samples ]
        filter_logs = [ join(item.output_dir, 'filter_log.txt') for item in samples ]
        filter_polya_logs = [ join(item.output_dir + '-polyA', 'filter_log.txt') for item in samples ]                
        #filter_logs = [ item.get_filter_action().log_filename() for item in samples ]
        #filter_polya_logs = [ item.get_polya_filter_action().log_filename() for item in samples ]

        analyse_template = tail_lengths.Analyse_tail_counts(
            working_dirs = dirs,
            saturation = 0,
            extension = self.extension,
            annotations = reference/'reference.gff',
            types = 'gene',
            )
        

        with nesoni.Stage() as stage:        
            for item in samples:
                item.process_make(stage)

        
        
        nesoni.Norm_from_samples(
            workspace/'norm',
            working_dirs = dirs
            ).make()

        def writer():
            for row in io.read_table(workspace/'norm.csv'):
                row['Name'] = row['Name']+'-polyA'
                yield row
        io.write_csv(workspace/'norm-polyA.csv', writer(), comments=['Normalization'])


        with nesoni.Stage() as stage:
            if self.include_plots:        
                for plot_name, directories, norm_filename in [
                      ('all',   dirs,       workspace/'norm.csv'),
                      ('polyA', polya_dirs, workspace/'norm-polyA.csv'),
                      ]:
                    nesoni.IGV_plots(
                        plotspace/plot_name,
                        working_dirs = directories,
                        label_prefix = plot_name+' ',
                        raw = True,
                        norm = True,
                        genome = reference.get_genome_filename(),
                        norm_file = norm_filename,
                        #delete_igv = False,
                        ).process_make(stage)

            analyse_gene_counts_0 = analyse_template(
                output_dir = expressionspace/'genewise',
                saturation = 0,
                extension = self.extension,
                title = 'Genewise expression - ' + self.title,
                file_prefix = file_prefix+'genewise-',
                )
            analyse_gene_counts_0.process_make(stage)
            
            analyse_gene_counts_1 = analyse_template(
                output_dir = expressionspace/'genewise-dedup',
                saturation = 1,
                title = 'Genewise expression with read deduplication - ' + self.title,
                file_prefix = file_prefix+'genewise-dedup-',
                )
            analyse_gene_counts_1.process_make(stage)
            
            stage.process(self._run_peaks, 
                workspace=workspace, expressionspace=expressionspace, reference=reference, 
                polya_dirs=polya_dirs, analyse_template=analyse_template, file_prefix=file_prefix,
                )
            
        with nesoni.Stage() as stage:
            for test in self.tests:
                test(
                    output_dir = testspace/test.output_dir,
                    analysis = self.output_dir
                    ).process_make(stage)
        
        #===============================================
        #                   Report        
        #===============================================

        r = reporting.Reporter(os.path.join(self.output_dir, 'report'), self.title, self.file_prefix)
                    
        r.heading('Alignment to reference')
        
        r.report_logs('alignment-statistics',
            #[ workspace/'stats.txt' ] +
            clipper_logs + filter_logs + #filter_polya_logs +
            [ expressionspace/('genewise-dedup','aggregate-tail-counts_log.txt') ],
            filter=lambda sample, field: (
                field not in [
                    
                    'fragments','fragments aligned to the reference','reads kept',
                    'average depth of coverage, ambiguous',
                    'average depth of coverage, unambiguous',
                    ]
            ),
        )


        if self.include_plots:        
            r.heading('IGV plots')
            
            r.p('These files show the depth of coverage. They can be viewed with the IGV genome browser.')
            
            genome_files = [ ]
            if self.include_genome:
                genome_files.append(reference.get_genome_filename())
                genome_dir = reference.get_genome_dir()
                base = os.path.split(self.genome_dir)[1]
                for filename in os.listdir(genome_dir):
                    genome_files.append((
                        os.path.join(genome_dir, filename),
                        os.path.join(base, filename)
                        ))
            
            r.p(r.tar('igv-plots',
                genome_files +
                glob.glob(plotspace/'*.tdf')
                ))
        

        if self.include_bams:
            r.heading('BAM files')
            
            r.p('These BAM files contain the alignments of reads to the reference sequences.')
            
            r.p('Reads with a poly(A) tail have an \'AN\' attribute giving the length of non-templated poly(A) sequence. '
                'Tail-tools only treats a read as having a tail if this length is at least 4.')
            
            bam_files = [ ]
            for name in names:
                bam_files.append( (samplespace/(name,'alignments_filtered_sorted.bam'),name+'.bam') )
                bam_files.append( (samplespace/(name,'alignments_filtered_sorted.bam.bai'),name+'.bam.bai') )
            r.p(r.tar('bam-files', bam_files))


        r.heading('Genewise expression')
        
        io.symbolic_link(source=expressionspace/('genewise','report'),link_name=r.workspace/'genewise')
        r.subheading('<a href="genewise/index.html">&gt; Genewise expression</a>')

        io.symbolic_link(source=expressionspace/('genewise-dedup','report'),link_name=r.workspace/'genewise-dedup')
        r.subheading('<a href="genewise-dedup/index.html">&gt; Genewise expression with read deduplication</a>')


        r.heading('Peakwise expression')

        web.Geneview_webapp(r.workspace/'view').run()        
        
        peak_filename = expressionspace/('peakwise','features-with-data.gff')
        n_peaks = len(list(annotation.read_annotations(peak_filename)))
        r.p('%d peaks called (%d poly(A) reads were required to call a peak).' % (n_peaks, self.peak_min_depth))
        
        r.p(r.get(peak_filename, name='peaks.gff') + ' - peaks called')        

        if self.groups:
            r.subheading('Peak shift between groups')
            r.p(r.get(workspace/('peak-shift','grouped.csv')) + ' - genes with a potential peak shift')        
            r.get(workspace/('peak-shift','grouped.json'))
            r.p('<a href="view.html?json=%sgrouped.json">Gene viewer</a>' % r.file_prefix)

        r.subheading('Peak shift between samples')
        r.p(r.get(workspace/('peak-shift','individual.csv')) + ' - genes with a potential peak shift')        
        r.get(workspace/('peak-shift','individual.json'))
        r.p('<a href="view.html?json=%sindividual.json">Gene viewer</a>' % r.file_prefix)
       
        
        io.symbolic_link(source=expressionspace/('peakwise','report'),link_name=r.workspace/'peakwise')
        r.subheading('<a href="peakwise/index.html">&gt; Peakwise expression</a>')

        io.symbolic_link(source=expressionspace/('peakwise-dedup','report'),link_name=r.workspace/'peakwise-dedup')
        r.subheading('<a href="peakwise-dedup/index.html">&gt; Peakwise expression with read deduplication</a>')
        
        
        if self.tests:
            r.heading('Differential tests')
            for test in self.tests:
                io.symbolic_link(source=testspace/test.output_dir,link_name=r.workspace/('test-'+test.output_dir))
                r.p('<a href="test-%s">%s</a>' % (test.output_dir, test.get_title()))


        r.write('<p/><hr>\n')
        
        r.p('This set of genes was used in the analysis:')
        
        r.p(r.get(reference/'reference.gff') + ' - Reference annotations in GFF3 format')
        r.p(r.get(reference/'utr.gff') + ' - 3\' UTR regions')

        r.p('tail-tools version '+tail_tools.VERSION)
        r.p('nesoni version '+nesoni.VERSION)
        #r.p('SHRiMP version '+grace.get_shrimp_2_version())
        
        r.close()


    def _run_peaks(self, workspace, expressionspace, reference, polya_dirs, analyse_template, file_prefix):
        shiftspace = io.Workspace(workspace/'peak-shift')

        Call_peaks(
            workspace/'peaks',
            annotations = reference/'reference.gff',
            shift_start = 0,
            shift_end = self.extension,
            min_depth = self.peak_min_depth,
            types = 'gene',
            polyas = polya_dirs,
            ).make()

        peak_template = analyse_template(
            annotations=workspace/('peaks','relation-child.gff'), 
            extension=0,
            types='peak',
            )

        with nesoni.Stage() as stage:
            peak_template(
                expressionspace/'peakwise',
                saturation=0,
                title='Peakwise expression - ' + self.title,
                file_prefix=file_prefix+'peakwise-',
                ).process_make(stage)
                
            peak_template(
                expressionspace/'peakwise-dedup',
                saturation=1,
                title='Peakwise expression with read deduplication - ' + self.title,
                file_prefix=file_prefix+'peakwise-dedup-',
                ).process_make(stage)

        alternative_tails.Compare_peaks(
            shiftspace/'individual',
            norm_file=expressionspace/('peakwise','norm.csv'),
            utrs=reference/'utr.gff',
            utr_only=True,
            top=2,
            reference=reference/'reference.fa',
            parents=workspace/('peaks','relation-parent.gff'),
            children=workspace/('peaks','relation-child.gff'),
            counts=expressionspace/('peakwise','counts.csv'),
            ).make()
        
        if self.groups:
            tail_lengths.Collapse_counts(
                shiftspace/'grouped-counts',
                counts=expressionspace/('peakwise','counts.csv'),
                groups=self.groups
                ).make()
            
            alternative_tails.Compare_peaks(
                shiftspace/'grouped',
                utrs=reference/'utr.gff',
                utr_only=True,
                top=2,
                reference=reference/'reference.fa',
                parents=workspace/('peaks','relation-parent.gff'),
                children=workspace/('peaks','relation-child.gff'),
                counts=shiftspace/'grouped-counts.csv',
                ).make()
            
            
            
                




