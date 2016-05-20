
#
# Tools associated with RNA-Seq alternate 3' UTR usages
#
# These support the end_shift_rnaseq test, which is implemented in R
#



from nesoni import annotation, span_index, io, config

def write_gff3(filename, items, header):
    with io.open_possibly_compressed_writer(filename) as f:
        f.write(header)
        for item in items:
            print >> f, item.as_gff()

@config.help("Process annotations suitable for RNA-Seq end shift test.","""\
features file should contain features of type:

exon            - exon locations
CDS             - CDS locations
three_prime_UTR - 3' UTR locations (configurable with --what)

Note: Use "--gene-level yes --what exon" to apply the end-shift test to whole genes!

""")
@config.Bool_flag("gene_level", "Hack to use genes rather than transcripts.")
@config.Int_flag("support", "Maximum transcript_support_level to use (this is an attribute in Ensembl gff files).")
@config.String_flag("what", "Feature type for UTR.")
@config.Positional("features", "GFF file containing genome annotation.")
class Make_rnaseq_reference(config.Action_with_output_dir):
    what = "three_prime_UTR"
    support = None
    features = None
    gene_level = False

    def run(self):
        workspace = self.get_workspace()
                
        header = [ "##gff-version 3\n" ]
        lengths = { }
        with io.open_possibly_compressed_file(self.features) as f:
            f.next()
            for line in f:
                if not line.startswith("#"): break
                if line.startswith("##gff-version"): continue
                header.append(line)
                parts = line.strip().split()
                if parts[0] == "##sequence-region":
                    lengths[parts[1]] = int(parts[3])
                    
        header = "".join(header)
                
        items = list(annotation.read_gff(self.features, "/"))
        annotation.link_up_annotations(items)
        for item in items:
            assert len(item.parents) < 2
            if "ID" in item.attr:
                item.attr["ID"] = item.attr["ID"].split(":")[1]
            if "Parent" in item.attr:
                item.attr["Parent"] = item.attr["Parent"].split(":")[1]
            if item.parents:
                item.parent = item.parents[0]
            
        
        def well_supported(item):
            if self.support is None: return True
            return int(item.attr.get("transcript_support_level","NA").split()[0]) <= self.support
        
        exons = [ item for item in items if item.type == "exon" and well_supported(item.parent) ]
        exon_index = span_index.index_annotations(exons)
        
        utrs = [ ]
        extended_utrs = [ ]
        utr_parts = [ ]
        exons_kept = [ ]
        cds_kept = [ ]
        transcripts_kept = [ ]
        for item in items:
            this_exons = [ item2 for item2 in item.children if item2.type == "exon" ]
            if this_exons and well_supported(item):    
                 transcripts_kept.append(item)
                 exons_kept.extend(this_exons)
                 cds_kept.extend([ item2 for item2 in item.children if item2.type == "CDS" ])
        
            if self.gene_level:
                utr_bits = [ item3 for item2 in item.children  if well_supported(item2)
                                   for item3 in item2.children if item3.type == self.what ] 
            else:
                if not well_supported(item): continue
                utr_bits = [ item2 for item2 in item.children if item2.type == self.what ] 
            
            if not utr_bits:
                continue
            
            utr = utr_bits[0].copy()
            for item2 in utr_bits[1:]:
                utr = utr.span_with(item2)
            
            gene = item if self.gene_level else item.parent
            
            utr.attr = dict(
                ID=item.get_id(),
                Name=item.attr["Name"],
                gene_id=gene.get_id(),
                gene=gene.attr["Name"],
                description=gene.attr.get("description",""),
                biotype=item.attr["biotype"]
                )
        
            max_extension = 10000
            if item.strand < 0:
                max_extension = min(max_extension, utr.start)
            else:
                max_extension = min(max_extension, lengths[utr.seqid] - utr.end)
            assert max_extension >= 0, utr
            
            end = utr.three_prime()
            for hit in exon_index.get(end.shifted(0,max_extension), same_strand=True):
                #if hit.parent.get_id() == item.get_id():
                #    continue
                rel = hit.relative_to(end).start
                if rel >= 0:
                    max_extension = min(max_extension, rel)
        
            extended_utr = utr.shifted(0,max_extension)
            extended_utr.start = max(extended_utr.start, 0)
            utr.attr["max_extension"] = str(max_extension)
        
            utrs.append(utr)
            extended_utrs.append(extended_utr)
            
            for item2 in utr_bits:
                part = item2.copy()
                part.attr = dict(Parent=item.get_id())
                part.type = "part"
                utr_parts.append(part)
                            
        write_gff3(workspace/"utr.gff",utrs,header)
        write_gff3(workspace/"utr_extended.gff",extended_utrs,header)
        write_gff3(workspace/"utr_part.gff",utr_parts,header)
        write_gff3(workspace/"transcript.gff",transcripts_kept,header)
        write_gff3(workspace/"exon.gff",exons_kept,header)
        write_gff3(workspace/"cds.gff",cds_kept,header)

     




        