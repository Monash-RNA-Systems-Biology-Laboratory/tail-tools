VERSION = '1.2dev'
#^ Note: this first line is read by the setup.py script to get the version

import nesoni

from .fasta_qual_merge import Fasta_qual_merge
from .clip_runs import Clip_runs_colorspace, Clip_runs_basespace
from .extend_sam import Extend_sam_colorspace, Extend_sam_basespace
from .proportions import Proportions, Proportions_heatmap 
from .tail_lengths import Tail_count, Aggregate_tail_counts, Plot_pooled, Plot_comparison, Analyse_tail_counts
from .alternative_tails import Compare_peaks
from .call_utrs import Call_utrs
from .test import Test
from .web import Geneview_webapp
from .reference_directory import Make_tt_reference, Make_ucsc_reference
from .reference_directory_ensembl import Make_ensembl_reference
from .primer_gff import Primer_gff
from .rnaseq import Make_rnaseq_reference
from .bigwig import Bam_to_bigwig, Bam_ambiguity, Polya_bigwigs
from .shiny import Shiny
from .peaks import Call_peaks
from .workflows import Analyse_polya, Analyse_polya_batch


def main():
    nesoni.run_toolbox([
            'tail-tools ' + VERSION,
            'Tools:',
            'These tools are used by the workflows below, you generally won\'t need to use these directly.',
            Fasta_qual_merge,
            Clip_runs_colorspace,
            Clip_runs_basespace,
            Extend_sam_colorspace,
            Extend_sam_basespace,
            Proportions,
            Proportions_heatmap,
            Tail_count,
            Aggregate_tail_counts,
            Plot_pooled,
            #Plot_comparison,
            Compare_peaks,
            Call_utrs,
            Test,
            Geneview_webapp,
            Shiny,
            
            'Component workflows:',
            Call_peaks,
            Polya_bigwigs,   
            Analyse_polya,
            Analyse_tail_counts,
            
            'Primary workflow:',
            Analyse_polya_batch,

            'Reference directory creation:',
            Make_tt_reference,
            Make_ucsc_reference,
            Make_ensembl_reference,

            'RNA-Seq analysis:',
            Make_rnaseq_reference,
            Bam_to_bigwig,
            Bam_ambiguity,

            'Utilities:',
            Primer_gff,            
        ], 
        'tail-tools',
        )
