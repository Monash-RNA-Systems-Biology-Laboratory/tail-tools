
import nesoni

import tail_tools

nesoni.run_toolbox(
    [
      'tail-tools ' + tail_tools.VERSION,
      'Tools:',
      tail_tools.Fasta_qual_merge,
      tail_tools.Clip_runs,
      tail_tools.Extend_sam,
      tail_tools.Proportions,
      tail_tools.Proportions_heatmap,
      tail_tools.Tail_lengths,
      tail_tools.Aggregate_tail_lengths,
      tail_tools.Tail_stats,
      tail_tools.Plot_pooled,
      tail_tools.Plot_comparison,
      
      'Workflows:',
      tail_tools.Analyse_polya,
      tail_tools.Analyse_tail_lengths,
      tail_tools.Analyse_polya_batch,
    ], 
    'tail-tools',
)

