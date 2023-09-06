
import os, math, glob, json, collections
from os.path import join

import nesoni
from nesoni import config, workspace, working_directory, reference_directory, io, reporting, grace, annotation, selection, span_index

import tail_tools
from . import clip_runs, extend_sam, proportions, tail_lengths, web, alternative_tails, bigwig, web, peaks

def _do_nothing():
    pass

def _make_each(actions):
    for item in actions:
        item.make()

def _make(action):
    action.make()

class _call(object):
    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs
    
    def __call__(self):
        self.func(*self.args, **self.kwargs)

def _serial(*items):
    for item in items:
        item()
        
def _parallel(*items):
    with nesoni.Stage() as stage:
        for item in items:
            stage.process(item)


        



@config.help(
'Create nesoni working directories for a single sample.',
"""\
Reads are clipped using "clip-runs:", aligned using SHRiMP, extended using "extend-sam:" and \
imported into a nesoni working directory.

Directories are created both using all reads and for reads with a poly-A tail.
""")
@config.Bool_flag('consensus', 'Look for SNPs and indels.')
@config.Bool_flag('discard_multimappers', 'If "yes", multimapping reads are discarded. If "no", multimapping reads are assigned at random.')
@config.Bool_flag('delete_files', 'Delete files not needed by downstream processing.')
@config.Section('tags', 'Tags for this sample. (See "nesoni tag:".)')
@config.Positional('reference', 'Reference directory created by "nesoni make-reference:"')
@config.Section('reads', 'Fastq files containing SOLiD or Illumina reads.')
@config.Configurable_section('clip_runs_colorspace', 'Configuration options for clip-runs-colorspace:. Will be used with colorspace reads.')
@config.Configurable_section('clip_runs_basespace', 'Configuration options for clip-runs-basespace:. Will be used with basespace reads.')
@config.String_flag("aligner", "Aligner to use, basespace only. Options are 'bowtie2' or 'STAR'.")
@config.Int_flag("min_match", "STAR only: minimum number of matches required for alignment.")
@config.Float_flag('extension_prop_a', 'Extending alignments over genomic "A"s, what is the lowest proportion of "A"s allowed? (Basespace only.)')
@config.Int_flag('min_tail',
    'Minimum number of As to be called a poly(A) tail. Note if using Analyse_polya_batch (as usual), the setting in Analyse_poly_batch will override this value.'
    )
class Analyse_polya(config.Action_with_output_dir):
    reference = None
    tags = [ ]
    reads = [ ]    
    consensus = False
    discard_multimappers = False
    delete_files = True
    
    aligner = "star"
    min_match = 0
    
    min_tail = 4
    
    clip_runs_colorspace = clip_runs.Clip_runs_colorspace()
    clip_runs_basespace = clip_runs.Clip_runs_basespace()
    extension_prop_a = 0.6
    
    _workspace_class = working_directory.Working
    
    #def get_polya_dir(self):
    #    return os.path.normpath(self.output_dir) + '-polyA'
    
    def get_filter_tool(self):
        if self.consensus:
            filter_tool = nesoni.Consensus
        else:
            filter_tool = nesoni.Filter
        return filter_tool(
            monogamous=self.discard_multimappers,
            random=True,
            infidelity=0,
            userplots=False,
        )
    
    def get_filter_action(self):
        return self.get_filter_tool()(working_dir = self.output_dir)

    #def get_polya_filter_action(self):
    #    return self.get_filter_tool()(working_dir = self.get_polya_dir())
    
    def run(self):
        assert self.reads, 'No read files given.'
        colorspace = [ io.is_colorspace(item) for item in self.reads ]
        assert len(set(colorspace)) == 1, 'Mixture of colorspace and basespace reads is not currently supported.'
        colorspace = colorspace[0]
        
        #polya_dir = self.get_polya_dir()
        
        working = working_directory.Working(self.output_dir, must_exist=False)
        working.set_reference(self.reference)
        reference = working.get_reference()
        
        #polya_working = working_directory.Working(polya_dir, must_exist=False)
        #polya_working.set_reference(self.reference)
        
        clipped_prefix = working/'clipped_reads'
        clipped_filename = clipped_prefix+('.csfastq.gz' if colorspace else '.fastq.gz')
        
        raw_filename = working/'alignments_raw.sam.gz'
        extended_filename = working/'alignments_extended.sam.gz'
        
        aligner = self.aligner.lower()
        assert aligner in ("bowtie2", "star")
        
        #polya_filename = working/'alignments_filtered_polyA.sam.gz'
        
        if colorspace:
            self.clip_runs_colorspace(
                filenames=self.reads,
                prefix=clipped_prefix,
                sample=working.name,
            ).make()
        else:
            self.clip_runs_basespace(
                filenames=self.reads,
                prefix=clipped_prefix,
                sample=working.name,
            ).make()        

        cores = nesoni.coordinator().get_cores()
        
        if colorspace:
            nesoni.Execute(
                command = reference.shrimp_command(cs=colorspace, parameters=[ clipped_filename ]) + [ '--qv-offset', '33' ],
                execution_options = [ '-N', str(cores) ],
                output=raw_filename,
                cores=cores,
                prefix=working/'run_alignment'
                ).make()
        
        elif aligner == "bowtie2":
            nesoni.Execute(
                command = [ 
                    'bowtie2', 
                    '--rg-id', '1',
                    '--rg', 'SM:'+working.name,
                    '--sensitive-local',
                    '-k', '10', #Up to 10 alignments per read
                    '-x', reference.get_bowtie_index_prefix(),
                    '-U', clipped_filename,
                    ],
                execution_options = [ '--threads', str(cores) ],
                output=raw_filename,
                cores=cores,
                prefix=working/'run_alignment'
                ).make()
        else:
            nesoni.Execute(
                command = [
                    'STAR',
                    '--genomeDir', reference/'star',
                    '--outFileNamePrefix', working/'star',
                    '--outSAMtype', 'SAM', #'Unsorted',
                    '--outStd', 'SAM',
                    '--readFilesIn', clipped_filename,
                    '--readFilesCommand', 'zcat',
                    # If we wanted no minimum proportion alignment could use this:
                    # (reads are clipped, so safe to leave default of 2/3 alignment)
                    #'--outFilterScoreMinOverLread', '0',
                    #'--outFilterMatchNminOverLread', '0',
                    '--outMultimapperOrder', 'Random',
                    # Only output 1 alignment for multmappers (NH still set)
                    #'--outSAMmultNmax', '1',
                    # Require alignment from start (not needed for clipped reads)
                    #'--alignEndsType', 'Extend5pOfRead1',
                    # No de novo introns, annotated introns will still be used
                    '--alignIntronMax', '20',
                    ] +
                    ([ '--outFilterMatchNmin', str(self.min_match) ] if self.min_match else [ ]),
                execution_options = [ '--runThreadN', str(cores) ],
                output=raw_filename,
                cores=cores,
                prefix=working/'run_alignment'
                ).make()
                
        if colorspace:
            extend_sam.Extend_sam_colorspace(
                input=raw_filename,
                output=extended_filename,
                reads=self.reads,
                reference_filenames=[ reference.reference_fasta_filename() ],
                tail=self.min_tail,
            ).make()
        else:    
            extend_sam.Extend_sam_basespace(
                input=raw_filename,
                output=extended_filename,
                clips=[ clipped_prefix+'.clips.gz' ],
                reference_filenames=[ reference.reference_fasta_filename() ],
                prop_a = self.extension_prop_a,
                tail=self.min_tail,
            ).make()
        
        nesoni.Import(
            input=extended_filename,
            output_dir=self.output_dir,
            reference=[ self.reference ],
        ).make()
        
        self.get_filter_action().make()

        #Tail_only(
        #    input=working/'alignments_filtered.bam',
        #    output=polya_filename,
        #).make()
        
        #nesoni.Import(
        #    input=polya_filename,
        #    output_dir=polya_dir,
        #    reference=[ self.reference ],
        #).make()
        
        # This shouldn't actually filter out any alignments.
        # We do it to produce depth of coverage plots
        # and position-sorted BAM files.
        #self.get_polya_filter_action().make()
        
        nesoni.Tag(self.output_dir, tags=self.tags).make()
        #nesoni.Tag(polya_dir, tags=self.tags).make()
        
        if self.delete_files:
            # Delete unneeded files
            os.unlink(clipped_prefix+'.state')
            os.unlink(clipped_filename)
            os.unlink(working/'alignments.bam')
            os.unlink(working/'alignments_filtered.bam')
            os.unlink(working/'run_alignment.state')
            os.unlink(raw_filename)
            os.unlink(extended_filename)
            #os.unlink(polya_filename)


@config.help(
    "Use a sample from a previous pipeline run.",
    "")
@config.String_flag("abductee", "Path to sample.")
@config.Positional('reference', 'Reference directory created by "nesoni make-reference:". Filled in by Analyse_polya_batch.')
class Abduct_polya(config.Action_with_output_dir):
    abductee = None
    reference = None

    # No reads (don't participate in duplicate filename checking)
    reads = [ ]

    def __init__(self, *args, **argd):
        super(Abduct_polya, self).__init__(*args, **argd)
        if self.abductee:
            abductee = working_directory.Working(self.abductee, must_exist=True)
            self.tags = abductee.param["tags"]

    def run(self):
        working = working_directory.Working(self.output_dir, must_exist=False)
        abductee = working_directory.Working(self.abductee, must_exist=True)

        for filename in [
                "parameters",
                "alignments_filtered_sorted.bam",
                "alignments_filtered_sorted.bam.bai",
                "depths.pickle.gz",
                "clipped_reads_log.txt",
                "filter_log.txt",
                ]:
            io.symbolic_link(abductee/filename, working/filename)


@config.help(
'Analyse a set of samples and produce an HTML report, including depth of coverage plots, heatmaps, etc.',
"""\

""")
@config.String_flag('title', 'Analysis report title.')
@config.String_flag('file_prefix', 'Prefix for report filenames.')
@config.Bool_flag('do_fragile', 'Do steps that might fail? Some steps fail if too few peaks are called.')
@config.Bool_flag('do_bigwigs', 'Do bigwig file creation?')
@config.Int_flag('adaptor',
    'Minimum number of adaptor bases required, 0 for no filtering.'
    )
@config.Int_flag('min_tail',
    'Minimum number of non-templated As for a read to be called as having a poly(A) tail. Affects both peak calling and estimation of average tail lengths.'
    )
@config.Int_flag('clip_tail',
    'Tails longer than this will be reduced to this length. 0 for no clipping.'
    )
@config.Int_flag('peak_min_depth', 
    'Number of poly(A) reads ending at nearly the same position required in order to call a peak.'
    )
@config.Bool_flag('peak_polya', 
    'Call peaks using poly(A) reads only.'
    )
@config.Int_flag('peak_length', 'Length of peak features, ie how far back a read may be from a peak and still be counted. Should be around read length or a little shorter.')
@config.Float_flag('peak_min_tail', 'Minimum average tail length required for a peak.')
@config.Bool_flag('peak_pair_utr', 'For the older peak-pair tests, consider only peaks in the 3\'UTR?')
@config.Int_flag('extension', 'How far downstrand of the given annotations a read or peak belonging to a gene might be. However tail-tools (with a recently created reference directory) will not extend over coding sequence.')
@config.String_flag('types', 'Comma separated list of feature types to use as genes. Default is "gene".')
@config.String_flag('parts', 'Comma separated list of feature types that make up features. Default is "exon". Alternatively you might use "three_prime_utr" for a stricter definition of where we expect reads or "gene" for a broad definition including introns.')
@config.Bool_flag('umis',
    'Should counting be based on UMIs (contained in read names)?')
@config.String_flag("species", 'Species for GO analysis. Currently supports Human ("Hs"), Saccharomyces cerevisiae ("Sc"), Caenorhabditis elegans ("Ce"), Mus musculus ("Mm")')
@config.String_flag("shiny_report_url", 'The main report contains a link to a shiny-app section. What URL will the shiny-app be served from?')
#@config.String_flag('blurb', 'Introductory HTML text for report')
#@config.String_flag('genome', 'IGV .genome file, to produce IGV plots')
#@config.String_flag('genome_dir', 'IGV directory of reference sequences to go with .genome file')
#@config.Bool_flag('include_plots', 'Include plots in report?')
#@config.Bool_flag('include_genome', 'Include genome in IGV plots tarball?')
#@config.Bool_flag('include_bams', 'Include BAM files in report?')
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
    
    extension = 1000
    
    do_fragile = True
    do_bigwigs = True

    adaptor = 0
    min_tail = 4
    clip_tail = 0
    
    peak_min_depth = 50
    peak_polya = True
    peak_length = 300
    peak_min_tail = 15.0
    peak_pair_utr = True
    
    types = "gene"
    parts = "exon"
    
    umis = False
    
    species = ""
    shiny_report_url = "shiny/"
    
    reference = None
    samples = [ ]
    
    groups = [ ]
    
    tests = [ ]
    
    def run(self):
        #===============================================
        #                Sanity checks
        #===============================================
        
        assert len(set([ item.output_dir for item in self.samples ])) == len(self.samples), "Duplicate sample name."
        
        all_inputs = [ ]
        for sample in self.samples:
            all_inputs.extend(sample.reads)
        assert len(set(all_inputs)) == len(all_inputs), "Duplicate read filename."
        
        assert len(set([ item.output_dir for item in self.tests ])) == len(self.tests), "Duplicate test name."
        
        for test in self.tests:
            assert not test.analysis, "analysis parameter for tests should not be set, will be filled in automatically"
        
        #===============================================
        #                Run pipeline
        #===============================================
        
        names = [ sample.output_dir for sample in self.samples ]
        
        reference = reference_directory.Reference(self.reference, must_exist=True)
        
        workspace = io.Workspace(self.output_dir, must_exist=False)
        samplespace = io.Workspace(workspace/'samples', must_exist=False)
        expressionspace = io.Workspace(workspace/'expression', must_exist=False)
        testspace = io.Workspace(workspace/'test', must_exist=False)
        
        self._create_json()
                
        file_prefix = self.file_prefix
        if file_prefix and not file_prefix.endswith('-'):
            file_prefix += '-'


        samples = [ ]
        for sample in self.samples:
            samples.append(sample(
                samplespace / sample.output_dir,
                reference = self.reference,
                min_tail = self.min_tail,
                ))
        
        dirs = [ item.output_dir for item in samples ]
        
        clipper_logs = [ join(item.output_dir, 'clipped_reads_log.txt') for item in samples ]
        filter_logs = [ join(item.output_dir, 'filter_log.txt') for item in samples ]
        
        analyse_template = tail_lengths.Analyse_tail_counts(
            working_dirs = dirs,
            extension = self.extension,
            annotations = reference/'reference.gff',
            types = self.types,
            parts = self.parts,
            adaptor = self.adaptor,
            tail = self.min_tail,
            clip_tail = self.clip_tail,
            umis = self.umis,
            )
        
        with nesoni.Stage() as stage:        
            for item in samples:
                item.process_make(stage)

        job_gene_counts = analyse_template(
            output_dir = expressionspace/'genewise',
            extension = self.extension,
            title = 'Genewise expression - ' + self.title,
            file_prefix = file_prefix+'genewise-',
            ).make
        
        job_peaks = _call(self._run_peaks, 
            workspace=workspace, 
            expressionspace=expressionspace, 
            reference=reference, 
            dirs = dirs,
            analyse_template = analyse_template,
            file_prefix=file_prefix,
            )
        
        job_norm = nesoni.Norm_from_samples(
            workspace/'norm',
            working_dirs = dirs
            ).make
            
        job_bigwig = bigwig.Polya_bigwigs(
            workspace/'bigwigs', 
            working_dirs = dirs, 
            norm_file = workspace/"norm.csv",
            peaks_file = workspace/("peaks", "relation-child.gff"),
            title = "IGV tracks - "+self.title
            ).make if self.do_bigwigs else _do_nothing
        
        job_norm_bigwig = _call(_serial, job_norm, job_bigwig)

        job_utrs = tail_tools.Call_utrs(
            workspace/('peaks','primary-peak'),
            self.reference,
            self.output_dir,
            extension=self.extension
            ).make if self.do_fragile else _do_nothing
            
        job_primpeak_counts = analyse_template(
            expressionspace/'primarypeakwise',
            annotations=workspace/('peaks','primary-peak-peaks.gff'), 
            extension=0,
            types='peak',
            parts='peak',
            title='Primary-peakwise expression - ' + self.title,
            file_prefix=file_prefix+'primarypeakwise-',
            ).make if self.do_fragile else _do_nothing
        
        job_primpeak = _call(_serial, job_utrs, job_primpeak_counts)
        
        job_peak_primpeak_bigwig = _call(_serial, 
            job_peaks, 
            _call(_parallel, job_norm_bigwig, job_primpeak))
        
        job_count = _call(_parallel, job_gene_counts, job_peak_primpeak_bigwig)
            
        test_jobs = [ ]
        for test in self.tests:
            test_jobs.append(test(
                output_dir = testspace/test.output_dir,
                analysis = self.output_dir,
                ).make)

        job_test = _call(_parallel, *test_jobs)

        job_raw = self._extract_raw

        job_all = _call(_serial, job_count, _call(_parallel, job_raw, job_test))        
        
        job_all()



        #===============================================
        #                   Report        
        #===============================================

        r = reporting.Reporter(workspace/'report', self.title, self.file_prefix, style=web.style())
        
        if self.do_bigwigs:
            io.symbolic_link(source=workspace/'bigwigs', link_name=r.workspace/'bigwigs')
            r.write('<div style="font-size: 150%; margin-top: 1em; margin-bottom: 1em;"><a href="bigwigs/index.html">&rarr; Load tracks into IGV</a></div>')

        tail_tools.Shiny(workspace/('report','shiny'), self.output_dir, title=self.title, species=self.species).run()
        r.write('<div style="font-size: 150%; margin-top: 1em; margin-bottom: 1em;"><a href="'+self.shiny_report_url+'" target="_blank">&rarr; Interactive report (shiny)</a></div>')
        
        r.heading('Alignment to reference')
        
        r.report_logs('alignment-statistics',
            #[ workspace/'stats.txt' ] +
            clipper_logs + filter_logs + #filter_polya_logs +
            [ expressionspace/('genewise','aggregate-tail-counts_log.txt') ],
            filter=lambda sample, field: (
                field not in [
                    'fragments','fragments aligned to the reference','reads kept',
                    'average depth of coverage, ambiguous',
                    'average depth of coverage, unambiguous',
                    ]
            ),
        )
        

        r.heading('Genewise expression')
        
        r.p("This is based on all reads within each gene (possibly from multiple peaks, or decay products).")
        
        io.symbolic_link(source=expressionspace/('genewise','report'),link_name=r.workspace/'genewise')
        r.p('<a href="genewise/index.html">&rarr; Genewise expression</a>')


        r.heading('Peakwise expression')
        
        r.p("This shows results from all called peaks.")
        
        peak_filename = expressionspace/('peakwise','features-with-data.gff')
        r.p(r.get(peak_filename, name='peaks.gff') + ' - peaks called')        

        self._describe_peaks(r)
        
        io.symbolic_link(source=expressionspace/('peakwise','report'),link_name=r.workspace/'peakwise')
        r.p('<a href="peakwise/index.html">&rarr; Peakwise expression</a>')

        if self.do_fragile:
            r.subheading('Primary-peakwise expression')
            
            r.p("This is based on the most prominent peak in the 3'UTR for each gene. (Peak can be up to %d bases downstrand of the annotated 3'UTR end, but not inside another gene on the same strand.)" % self.extension)
            
            io.symbolic_link(source=expressionspace/('primarypeakwise','report'),link_name=r.workspace/'primarypeakwise')
            r.p('<a href="primarypeakwise/index.html">&rarr; Primary-peakwise expression</a>')

            r.p(r.get(workspace/('peaks','primary-peak-peaks.gff')) + ' - primary peaks for each gene.')
            r.p(r.get(workspace/('peaks','primary-peak-utrs.gff')) + ' - 3\' UTR regions, based on primary peak call.')
            r.p(r.get(workspace/('peaks','primary-peak-genes.gff')) + ' - full extent of gene, based on primary peak call.')


        if self.tests:
            r.heading('Differential tests')
            for test in self.tests:
                io.symbolic_link(source=testspace/test.output_dir,link_name=r.workspace/('test-'+test.output_dir))
                r.p('<a href="test-%s">&rarr; %s</a> '
                    % (test.output_dir, test.get_title()))


        web.Geneview_webapp(r.workspace/'view').run()        
                
        r.heading('Gene viewers')
        r.p('Having identified interesting genes from heatmaps and differential tests above, '
            'these viewers allow specific genes to be examined in detail.')
        
        if self.groups and os.path.exists(workspace/('peak-shift','grouped.json')):
            r.get(workspace/('peak-shift','grouped.json'))
            r.p('<a href="view.html?json=%sgrouped.json">&rarr; Gene viewer, grouped samples</a>' % r.file_prefix)
        
        if os.path.exists(workspace/('peak-shift','individual.json')):
            r.get(workspace/('peak-shift','individual.json'))
            r.p('<a href="view.html?json=%sindividual.json">&rarr; Gene viewer, individual samples</a>' % r.file_prefix)
        
        r.heading('Raw data')
        
        r.p(r.tar('csv-files',glob.glob(workspace/('raw','*.csv'))))
        
        r.write('<ul>\n')
        r.write('<li> -info.csv = gene name and product, etc\n')
        r.write('<li> -count.csv = read count\n')
        r.write('<li> -mlog2-RPM.csv = moderated log2 Reads Per Million\n')
        r.write('<li> -tail.csv = average poly(A) tail length\n')
        r.write('<li> -tail-count.csv = poly(A) read count\n')
        r.write('<li> -proportion.csv = proportion of reads with poly(A)\n')
        r.write('<li> -norm.csv = read count normalization used for log2 transformation, heatmaps, differential tests, etc etc\n')
        r.write('</ul>\n')

        r.p('This set of genes was used in the analysis:')
        
        r.p(r.get(reference/'reference.gff') + ' - Reference annotations in GFF3 format')
        r.p(r.get(reference/'utr.gff') + ' - 3\' UTR regions')

        r.p('<b>%d further bases 3\' extension was allowed</b> beyond the GFF files above (but not extending into the next gene on the same strand).' % self.extension)

        r.write('<p/><hr>\n')
        r.subheading('About normalization and log transformation')

        r.p('Counts are converted to '
            'log2 Reads Per Million using Anscombe\'s variance stabilizing transformation '
            'for the negative binomial distribution, implemented in '
            'R package "varistran".')
                
        r.write('<p/><hr>\n')

        r.p('Reference directory '+self.reference)
        r.p('Tail Tools version '+tail_tools.VERSION)
        r.p('Nesoni version '+nesoni.VERSION)
        
        r.close()


    def _run_peaks(self, workspace, expressionspace, reference, dirs, analyse_template, file_prefix):
        shiftspace = io.Workspace(workspace/'peak-shift')

        peaks.Call_peaks(
            workspace/'peaks',
            annotations = reference/'reference.gff',
            extension = self.extension,
            min_depth = self.peak_min_depth,
            polya = self.peak_polya,
            min_tail = self.peak_min_tail,
            peak_length = self.peak_length,
            samples = dirs,
            ).make()

        analyse_template(
            expressionspace/'peakwise',
            annotations=workspace/('peaks','relation-child.gff'), 
            extension=0,
            types='peak',
            parts='peak',
            title='Peakwise expression - ' + self.title,
            file_prefix=file_prefix+'peakwise-',
            ).make()
        
        if self.do_fragile:    
            alternative_tails.Compare_peaks(
                shiftspace/'individual',
                norm_file=expressionspace/('peakwise','norm.csv'),
                utrs=reference/'utr.gff',
                utr_only=self.peak_pair_utr,
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


    def _extract_raw(self):            
        work = io.Workspace(self.output_dir, must_exist=False)
        raw = io.Workspace(work/'raw', must_exist=False)
        
        for name, counts, norms in [
                ('genewise',
                    work/('expression','genewise','counts.csv'),
                    work/('expression','genewise','norm.csv'),
                    ),
                ('primarypeakwise',
                    work/('expression','primarypeakwise','counts.csv'),
                    work/('expression','primarypeakwise','norm.csv'),
                    ),
                ('peakwise',
                    work/('expression','peakwise','counts.csv'),
                    work/('expression','peakwise','norm.csv'),
                    ),
                ('pairwise',
                    work/('peak-shift','individual-pairs.csv'),
                    work/('peak-shift','individual-pairs-norm.csv'),
                    ),
                ]:
            
            if not self.do_fragile and name in ('primarypeakwise','pairwise'):
                continue
            
            nesoni.Vst(
                raw/(name+'-mlog2-RPM'),
                counts,
                norm_file = norms
                ).make()
            
            counts_table = io.read_grouped_table(counts)
            io.write_csv_2(raw/(name+'-info.csv'), counts_table['Annotation'])
            io.write_csv_2(raw/(name+'-count.csv'), counts_table['Count'])
            io.write_csv_2(raw/(name+'-tail.csv'), counts_table['Tail'])
            io.write_csv_2(raw/(name+'-tail-count.csv'), counts_table['Tail_count'])
            io.write_csv_2(raw/(name+'-proportion.csv'), counts_table['Proportion'])
            
            norm_table = io.read_grouped_table(norms)
            io.write_csv_2(raw/(name+'-norm.csv'), norm_table['All'])


    def _create_json(self):                    
        workspace = io.Workspace(self.output_dir, must_exist=False)
        
        samples = [ ]
        groups = [ ]
        for sample in self.samples:
            this_groups = [ ]
            for item in self.groups:
                if selection.matches(
                        selection.term_specification(item),
                        sample.tags + [ sample.output_dir ]
                        ):
                    this_groups.append(selection.term_name(item))
            group = ','.join(this_groups) if this_groups else 'ungrouped'
            if group not in groups: groups.append(group)
            
            item = {
                'name' : sample.output_dir,
                'bam' : os.path.abspath( 
                    workspace/('samples',sample.output_dir,'alignments_filtered_sorted.bam')
                    ),
                'group' : group,
                'tags' : sample.tags,
                }
            samples.append(item)
            
        obj = collections.OrderedDict()
        obj['reference'] = os.path.abspath( self.reference )
        obj['extension'] = self.extension
        obj['genes'] = os.path.abspath( workspace/('peaks','relation-parent.gff') )
        obj['peaks'] = os.path.abspath( workspace/('peaks','relation-child.gff') )
        obj['groups'] = groups
        obj['samples'] = samples
        
        with open(workspace/"plotter-config.json","wb") as f:
            json.dump(obj, f, indent=4)
    
    
    def _describe_peaks(self, r):        
        workspace = io.Workspace(self.output_dir, must_exist=False)
        
        counts = io.read_grouped_table(workspace/("expression","peakwise","counts.csv"))["Count"]
        
        peak_counts = collections.defaultdict(int)
        read_counts = collections.defaultdict(int)
        total = 0
        for item in annotation.read_annotations(workspace/("peaks","relation-child.gff")):
            peak_counts[item.attr.get("Relation","None")] += 1
            read_counts[item.attr.get("Relation","None")] += sum(int(c) for c in counts[item.get_id()].values())            
            total += 1
        
        total_reads = sum(read_counts.values())
        
        r.write("<p>\n")
        r.write("%d peaks\n" % total)
        
        for name, desc in [
                ("3'UTR", "in a 3' UTR"),
                ("Exon", "otherwise in an exon"),
                ("Downstrand", "otherwise downstrand of a non-coding RNA"),
                ("Intron", "otherwise in an intron"),
                ("Antisense", "otherwise antisense to a gene"),
                ("None", "couldn't be related to annotated genes"),
                ]:
            r.write("<br/>%d peaks and %.1f%% of reads %s\n" % (peak_counts[name], read_counts[name]*100.0/total_reads, desc))
        r.write("</p>\n")
        






