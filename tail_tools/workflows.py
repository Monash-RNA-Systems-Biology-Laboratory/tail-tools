
import os, math, glob

import nesoni
from nesoni import config, workspace, working_directory, reference_directory, io, reporting, grace

import tail_tools
from tail_tools import clip_runs, extend_sam, proportions, tail_lengths


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
    min_depth = 10
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
            to_child = 'gene/product',
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
    consensus = True
    
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
        
        nesoni.Execute(
            command = reference.shrimp_command(cs=colorspace, parameters=[ clipped_filename ]),
            execution_options = [ '-N', str(cores) ] + [ '--qv-offset', '33' ] if not colorspace else [ ],
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
@config.String_flag('title', 'Analysis report title')
@config.String_flag('file_prefix', 'Prefix for report files')
@config.String_flag('blurb', 'Introductory HTML text for report')
@config.String_flag('genome', 'IGV .genome file, to produce IGV plots')
@config.String_flag('genome_dir', 'IGV directory of reference sequences to go with .genome file')
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
@config.Configurable_section('analyse_tail_lengths', 'Parameters for "analyse-tail-lengths:"')
@config.Section('extra_files', 'Extra files to include in report')
class Analyse_polya_batch(config.Action_with_output_dir):
    file_prefix = ''
    title = '3\'seq analysis'
    blurb = ''
    extra_files = [ ]
    
    genome = None
    genome_dir = None
    include_plots = True
    include_genome = True
    include_bams = True
    reference = None
    samples = [ ]
    
    analyse_tail_lengths = tail_lengths.Analyse_tail_lengths()
    
    def run(self):
        stage = nesoni.Stage()

        names = [ sample.output_dir for sample in self.samples ]
            #os.path.splitext(os.path.split(item)[1])[0]
            #for item in self.reads
            #]
        
        reference = reference_directory.Reference(self.reference, must_exist=True)
        
        workspace = io.Workspace(self.output_dir, must_exist=False)
        samplespace = io.Workspace(workspace/'samples', must_exist=False)
        plotspace = io.Workspace(workspace/'plots', must_exist=False)
        expressionspace = io.Workspace(workspace/'expression', must_exist=False)

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
        
        filter_logs = [ item.get_filter_action().log_filename() for item in samples ]
        filter_polya_logs = [ item.get_polya_filter_action().log_filename() for item in samples ]
        
        for item in samples:
            item.process_make(stage)

        #for reads_filename, directory in zip(self.reads, dirs):
        #    action = self.analyse(
        #        output_dir=directory,
        #        reference=self.reference,
        #        reads=reads_filename,
        #    )
        #    stage.process(action.run)
        #    filter_logs.append(action.get_filter_action().log_filename())
        #    filter_polya_logs.append(action.get_filter_action().log_filename())

        stage.barrier() #====================================================
        
        
        nesoni.Norm_from_samples(
            workspace/'norm',
            working_dirs = dirs
            ).make()

        def writer():
            for row in io.read_table(workspace/'norm.csv'):
                row['Name'] = row['Name']+'-polyA'
                yield row
        io.write_csv(workspace/'norm-polyA.csv', writer(), comments=['Normalization'])


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
                    genome = self.genome,
                    norm_file = norm_filename,
                    #delete_igv = False,
                ).process_make(stage)

        
        #nesoni.Stats(
        #    *self.reads, 
        #    output=workspace/'stats.txt'
        #).process_make(stage)
        
        stage.barrier() #====================================================

        analyse_tail_lengths_0 = self.analyse_tail_lengths(
            prefix = expressionspace/'nodedup',
            working_dirs = dirs,
            saturation = 0,
        )
        analyse_tail_lengths_0.process_make(stage)

        analyse_tail_lengths_1 = self.analyse_tail_lengths(
            prefix = expressionspace/'dedup',
            working_dirs = dirs,
            saturation = 1,
        )
        analyse_tail_lengths_1.process_make(stage)
        
        stage.barrier() #====================================================
        
        #===============================================
        #                   Report        
        #===============================================

        r = reporting.Reporter(os.path.join(self.output_dir, 'report'), self.title, self.file_prefix)
        
        r.write(self.blurb)
        
        for filename in self.extra_files:
            r.p( r.get(filename) )
        
        r.heading('Alignment to reference')
        
        r.report_logs('alignment-statistics',
            #[ workspace/'stats.txt' ] +
            filter_logs + filter_polya_logs +
            [ expressionspace/'dedup_aggregate_log.txt' ],
            filter=lambda sample, field: (
                field not in ['fragments','fragments aligned to the reference'] and
                (not sample.endswith('-polyA') or field not in ['reads with alignments','hit multiple locations'])
            ),
        )


        if self.include_plots:        
            r.heading('IGV plots')
            
            r.p('These files show the depth of coverage. They can be viewed with the IGV genome browser.')
            
            genome_files = [ ]
            if self.include_genome:
                assert self.genome, '.genome file not specified.'
                genome_files.append(self.genome)
                if self.genome_dir:
                    base = os.path.split(self.genome_dir)[1]
                    for filename in os.listdir(self.genome_dir):
                        genome_files.append((
                            os.path.join(self.genome_dir, filename),
                            os.path.join(base, filename)
                        ))
            
            r.p(r.tar('igv-plots',
               genome_files +
               glob.glob(plotspace/'*.tdf')
            ))
        

        if self.include_bams:
            r.heading('BAM files')
            
            r.p('These BAM files contain the alignments of reads to the reference sequences.')
            
            r.p('Reads with a poly(A) tail have an \'AA\' attribute.')
            
            bam_files = [ ]
            for name in names:
                bam_files.append( (samplespace/(name,'alignments_filtered_sorted.bam'),name+'.bam') )
                bam_files.append( (samplespace/(name,'alignments_filtered_sorted.bam.bai'),name+'.bam.bai') )
            r.p(r.tar('bam-files.tar.gz', bam_files))
        
        r.write('<div style="background: #ddddff; padding: 1em;">\n')
        r.heading('Expression levels and tail lengths <b>without</b> read deduplication')        
        analyse_tail_lengths_0.report_tail_lengths(r)
        r.write('</div>')

        r.write('<div style="background: #ddffdd; padding: 1em;">\n')
        r.heading('Expression levels and tail lengths <b>with</b> read deduplication')        
        analyse_tail_lengths_1.report_tail_lengths(r)
        r.write('</div>\n')

        r.write('<p/><hr>\n')
        r.p('tail_tools version '+tail_tools.VERSION)
        r.p('nesoni version '+nesoni.VERSION)
        r.p('SHRiMP version '+grace.get_shrimp_2_version())
        
        r.close()




