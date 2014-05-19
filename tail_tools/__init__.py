VERSION = '0.30'
#^ Note: this first line is read by the setup.py script to get the version

import nesoni

from .fasta_qual_merge import Fasta_qual_merge
from .clip_runs import Clip_runs_colorspace, Clip_runs_basespace
from .extend_sam import Extend_sam_colorspace, Extend_sam_basespace
from .proportions import Proportions, Proportions_heatmap 
from .tail_lengths import Tail_count, Aggregate_tail_counts, Plot_pooled, Plot_comparison, Collapse_counts, Analyse_tail_counts
from .alternative_tails import Compare_peaks
from .test import Test
from .web import Geneview_webapp
from .workflows import Call_peaks, Analyse_polya, Analyse_polya_batch
from .reference_directory import Make_tt_reference, Make_ucsc_reference

def main():
    nesoni.run_toolbox([
            'tail-tools ' + VERSION,
            'Tools:',
            'These tools are used by the workflows below, you generally won\'t need to use these directly.',
            Make_tt_reference,
            Make_ucsc_reference,
            Fasta_qual_merge,
            Clip_runs_colorspace,
            Clip_runs_basespace,
            Extend_sam_colorspace,
            Extend_sam_basespace,
            Proportions,
            Proportions_heatmap,
            Tail_count,
            Aggregate_tail_counts,
            #Tail_stats,
            Plot_pooled,
            Plot_comparison,
            Collapse_counts,
            Compare_peaks,
            Test,
            Geneview_webapp,
            
            'Workflows:',
            Analyse_polya,
            Analyse_tail_counts,
            Analyse_polya_batch,
            Call_peaks,
        ], 
        'tail-tools',
        )
