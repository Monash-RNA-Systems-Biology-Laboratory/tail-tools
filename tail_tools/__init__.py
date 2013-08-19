VERSION = '0.19'
#^ Note: this first line is read by the setup.py script to get the version

import nesoni

from fasta_qual_merge import Fasta_qual_merge
from clip_runs import Clip_runs_colorspace, Clip_runs_basespace
from extend_sam import Extend_sam_colorspace, Extend_sam_basespace
from proportions import Proportions, Proportions_heatmap 
from tail_lengths import Tail_lengths, Aggregate_tail_lengths, Plot_pooled, Plot_comparison, Analyse_tail_lengths
from alternative_tails import Compare_peaks
from web import Geneview_webapp
from workflows import Call_peaks, Analyse_polya, Analyse_polya_batch

def main():
    nesoni.run_toolbox([
            'tail-tools ' + VERSION,
            'Tools:',
            Fasta_qual_merge,
            Clip_runs_colorspace,
            Clip_runs_basespace,
            Extend_sam_colorspace,
            Extend_sam_basespace,
            Proportions,
            Proportions_heatmap,
            Tail_lengths,
            Aggregate_tail_lengths,
            #Tail_stats,
            Plot_pooled,
            Plot_comparison,

            Compare_peaks,
            Geneview_webapp,
            
            'Workflows:',
            Analyse_polya,
            Analyse_tail_lengths,
            Analyse_polya_batch,
            Call_peaks,
        ], 
        'tail-tools',
        )
