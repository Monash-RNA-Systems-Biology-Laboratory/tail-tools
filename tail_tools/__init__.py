VERSION = '0.10'

import nesoni

from fasta_qual_merge import Fasta_qual_merge
from clip_runs import Clip_runs
from extend_sam import Extend_sam
from proportions import Proportions, Proportions_heatmap 
from tail_lengths import Tail_lengths, Aggregate_tail_lengths, Tail_stats, Plot_pooled, Plot_comparison, Analyse_tail_lengths
from workflows import Analyse_polya, Analyse_polya_batch

__all__ = [
    'Fasta_qual_merge',
    'Fasta_qual_merge',
    'Clip_runs',
    'Extend_sam',
    'Proportions',
    'Proportions_heatmap',
    'Analyse_polya',
    'Analyse_polya_batch',
    ]

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
            Tail_stats,
            Plot_pooled,
            Plot_comparison,
            
            'Workflows:',
            Analyse_polya,
            Analyse_tail_lengths,
            Analyse_polya_batch,
        ], 
        'tail-tools',
        )
