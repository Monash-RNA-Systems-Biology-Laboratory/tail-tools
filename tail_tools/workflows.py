
import os, math, glob

import nesoni
from nesoni import config, working_directory, reference_directory, io, reporting, grace

import tail_tools
from tail_tools import clip_runs, extend_sam, proportions, tail_lengths


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
        polya_dir = self.get_polya_dir()
    
        working = working_directory.Working(self.output_dir, must_exist=False)
        working.set_reference(self.reference)
        reference = working.get_reference()
        
        polya_working = working_directory.Working(polya_dir, must_exist=False)
        polya_working.set_reference(self.reference)
        
        clipped_prefix = working/'clipped_reads'
        clipped_filename = clipped_prefix+'.csfastq.gz'
        
        raw_filename = working/'alignments_raw.sam.gz'
        extended_filename = working/'alignments_extended.sam.gz'
        
        polya_filename = working/'alignments_filtered_polyA.sam.gz'

        
        clip_runs.Clip_runs(
            filenames=self.reads,
            prefix=clipped_prefix,
            sample=working.name,
        ).make()
        
        cores = min(nesoni.coordinator().get_cores(), 8)
        
        nesoni.Execute(
            command = reference.shrimp_command(cs=True, parameters=[ clipped_filename ]),
            execution_options = [ '-N', str(cores) ],
            output=raw_filename,
            cores=cores,
        ).make()
                
        extend_sam.Extend_sam(
            input=raw_filename,
            output=extended_filename,
            reads=self.reads,
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

        #tail_tools.Tail_lengths(
        #    working_dir=self.output_dir, 
        #    types=self.types,
        #    ).make()    
        
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

        if self.include_plots:
            for plot_name, directories in [
                ('all',   dirs),
                ('polyA', polya_dirs),
            ]:
                nesoni.IGV_plots(
                    plotspace/plot_name,
                    working_dirs = directories,
                    label_prefix = plot_name+' ',
                    raw = False,
                    norm = True,
                    genome = self.genome,
                    #norm_file = workspace/norm_filename,
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




