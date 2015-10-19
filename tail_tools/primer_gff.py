
import csv, re

from nesoni import config, io, bio, annotation
from . import env

@config.help("Create a GFF file for regions amplified in a REPAT experiment.","""\
<csv-file> should be a CSV file with headings "ID" and "Primer", for example:

ID,Primer
gene1,ACGTACGTACGT
gene2,CGTACGTAGCAT

""")
@config.Int_flag("skip",
    "Discard first n bases of primer sequences.")
@config.Int_flag("length",
    "Length of features to produce, should be at least intended read length.")
@config.Positional("reference",
    "Reference directory."
    )
@config.Positional("csv_file",
    "CSV file of primers."
    )
class Primer_gff(config.Action_with_prefix):
    reference = None
    csv_file = None

    skip = 0
    length = 400

    def run(self):
        seqs = env.load_ref(self.reference).seqs
        
        result = [ ]
        errors = [ ]
        with open(self.csv_file, "rU") as f:
            reader = csv.reader(f)
            headings = reader.next()
            headings = [ item.lower() for item in headings ]
            assert "id" in headings
            assert "primer" in headings
            id_col = headings.index("id")
            primer_col = headings.index("primer")
            
            for row in reader:
                if len(row) == 0 or (not row[id_col].strip() and not row[primer_col].strip()):
                    continue
                
                id = row[id_col].strip()
                assert " " not in id, "ID contains space: "+id
                primer = row[primer_col].strip().upper()
                assert len(primer) > self.skip, "Primer too short: "+id
                assert [ char in "ACGT" for char in primer ], "Primer not ACGT: "+id
                primer = primer[self.skip:]
                rprimer = bio.reverse_complement(primer)
                
                hits = [ ]
                for seq_name in seqs:
                    for match in re.finditer( 
                            primer, seqs[seq_name], re.IGNORECASE):
                        hits.append( (seq_name, 1, match.start(), match.start()+self.length) )
                    for match in re.finditer(
                            rprimer, seqs[seq_name], re.IGNORECASE):
                        hits.append( (seq_name, -1, match.end()-self.length, match.end()) )
                    if len(hits) > 100:
                        raise config.Error("Many many hits for "+id+".")
                
                if not hits:
                    errors.append("No hits for "+id+".")
                    continue

                if len(hits) > 1:
                    self.log.log("Warning: %d hits for %s.\n" % (len(hits),id))
                    
                for i, hit in enumerate(hits):
                    hit_name = id
                    if len(hits) > 1:
                        hit_name += "-%dof%d" % (i+1,len(hits))
                    result.append(annotation.Annotation(
                        seqid = hit[0],
                        source = "tail-tools",
                        type = "region",
                        start = hit[2],
                        end = hit[3],
                        strand = hit[1],
                        attr = dict(
                            ID=hit_name,
                            Primer=primer
                            )
                        ))
        
        if errors:
            raise config.Error("\n".join(errors))
        
        annotation.write_gff3(self.prefix+".gff", result)
        
        
        
        
        