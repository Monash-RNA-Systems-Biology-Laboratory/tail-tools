VERSION = '0.17'
#^ Note: this first line is read by the setup.py script to get the version

import nesoni

from fasta_qual_merge import Fasta_qual_merge
from clip_runs import Clip_runs
from extend_sam import Extend_sam
from proportions import Proportions, Proportions_heatmap 
from tail_lengths import Tail_lengths, Aggregate_tail_lengths, Plot_pooled, Plot_comparison, Analyse_tail_lengths
from alternative_tails import Compare_peaks
from workflows import Analyse_polya, Analyse_polya_batch

def main():
    nesoni.run_toolbox([
            'tail-tools ' + VERSION,
            'Tools:',
            Fasta_qual_merge,
            Clip_runs,
            Extend_sam,
            Proportions,
            Proportions_heatmap,
            Tail_lengths,
            Aggregate_tail_lengths,
            #Tail_stats,
            Plot_pooled,
            Plot_comparison,

            Compare_peaks,
            
            'Workflows:',
            Analyse_polya,
            Analyse_tail_lengths,
            Analyse_polya_batch,
        ], 
        'tail-tools',
        )
