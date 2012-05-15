
import nesoni

import tail_tools

nesoni.run_toolbox([
    tail_tools.Fasta_qual_merge,
    tail_tools.Clip_runs,
    tail_tools.Extend_sam,
    tail_tools.Ratios,
    tail_tools.Analyse_polya,
    tail_tools.Analyse_polya_batch,
])

