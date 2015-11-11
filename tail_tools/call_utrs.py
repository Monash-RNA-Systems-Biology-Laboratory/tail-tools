
import nesoni
from nesoni import config, annotation
from . import env
import collections

def _extend(feature, extension):
    return feature.shifted(0, min(extension,int(feature.attr.get("max_extension","10000"))))


@config.help(
    'Call 3\' UTR regions based on the most prominent peak.'
    )
@config.Int_flag('extension',
    'Peak may be this far downstrand of the annotated 3\'UTR '
    '(but not into another gene on the same strand).')
@config.Positional('reference_dir')
@config.Positional('analysis_dir')
@config.Section('samples', 'Samples to use, defaults to all.')
class Call_utrs(config.Action_with_prefix):
    reference_dir = None
    analysis_dir = None
    extension = 1000
    samples = [ ]
    
    def run(self):
        ref = env.load_ref(self.reference_dir)
        analysis = env.load_analysis(self.analysis_dir)
        
        samples = self.samples
        if not samples:
            samples = analysis.peak_counts['Count'].value_type().keys()
        
        def total_count(peak_id):
            return sum(analysis.peak_counts['Count'][peak_id][sample] for sample in samples)
        
        #for each gene
        #determine a cutoff region
        #get the most highly expressed peak in that region
        #output something useful
        
        peaks_of = collections.defaultdict(list)
        for peak in analysis.peaks.values():
            if "Parent" in peak.attr:
                peaks_of[peak.attr["Parent"]].append(peak)
        
        called_peaks = [ ]
        called_utrs = [ ]
        called_genes = [ ]
        
        n_good = 0
        n_bad = 0
        for utr in ref.utrs.values():
            query = _extend(utr, self.extension)
            
            candidates = [ ]
            for peak in peaks_of[utr.attr["Parent"]]:
            	if peak.attr.get("Relation") != "3'UTR": continue
            
                relpeak = peak.relative_to(query)
                if relpeak.end <= 0: continue
                
                candidates.append((-total_count(peak.get_id()), peak.get_id()))
            candidates.sort()
            
            if not candidates:
                #print 'no peak for ', utr.attr['Parent']
                n_bad += 1
                continue
            
            peak = analysis.peaks_asis[candidates[0][1]].copy()
            peak.attr['Parent'] = utr.attr['Parent']
            peak.attr['Name'] = utr.attr.get('Name','')
            peak.attr['Product'] = utr.attr.get('Product','')        
            called_peaks.append(peak)
            
            called_utr = utr.five_prime().span_with(peak.three_prime())
            called_utr.attr['Peak'] = peak.get_id()
            called_utrs.append(called_utr)

            called_gene = ref.genes[utr.attr['Parent']].five_prime().span_with(peak.three_prime())
            called_gene.attr['Peak'] = peak.get_id()
            called_genes.append(called_gene)

            n_good += 1
        
        print n_good, 'UTR called'
        print n_bad, 'no UTR called'
        
        annotation.write_gff3(self.prefix + '-peaks.gff', called_peaks)
        annotation.write_gff3(self.prefix + '-utrs.gff', called_utrs)
        annotation.write_gff3(self.prefix + '-genes.gff', called_genes)
            
            
            
            
            
            
if __name__ == '__main__':
    nesoni.run_tool(Call_utrs)
