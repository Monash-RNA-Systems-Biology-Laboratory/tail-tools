

from nesoni import config, workspace, io, annotation, span_index
from . import reference_directory

import collections, re

attr_renaming = {
    "biotype" : "Biotype",
    "description" : "Product",    
}

source_levels = { "ensembl_havana":0, "havana":1, "ensembl":2 }

def translate_attr(attr):
    result = { }
    for name in attr:
        result[ attr_renaming.get(name,name) ] = attr[name]
    return result


def matches(pattern, value):
    if not pattern: return True
    return re.match("^"+pattern+"$", value) is not None


def natural_sorted(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def show_tree(item, indent=0):
    print "  "*indent, item.type
    for item2 in item.children:
        show_tree(item2, indent+1)

def all_within(item):
    result = [item]
    for item2 in item.children:
        result.extend(all_within(item2))
    return result


def translate_gene(item, transcript_type, transcript_biotype):
    if "ID" not in item.attr: return [ ]
    
    result = item.copy()
    assert result.attr["ID"].count(":") == 1
    result.attr["ID"] = result.attr["ID"].split(":")[1]
    result.attr = translate_attr(result.attr)
    result.type = "gene"
    result.children = [ ]
    for child in item.children:
        if not matches(transcript_type, child.type): continue
        if not matches(transcript_biotype, child.attr.get("biotype","")): continue
        result.children.extend(translate_transcript(child))
    if result.children:
        return [result]
    else:
        return [ ]


def translate_transcript(item):
    result = item.copy()
    #assert result.attr["ID"].count(":") == 1
    #result.attr["ID"] = result.attr["ID"].split(":")[1]
    result.attr = translate_attr(result.attr)
    result.type = "mRNA" #"transcript" would be better but this is what tail tools expects
    result.children = [ ]
    for child in item.children:
        if child.type not in ["exon","CDS"]: continue
        child = child.copy()
        child.attr = translate_attr(child.attr)
        child.attr["Parent"] = result.attr["ID"]
        child.children = [ ]
        result.children.append(child)
    if result.children:
        return [result]
    else:
        return [ ]


def extract_and_translate(items, type, biotype, transcript_type, transcript_biotype):
    result = [ ]
    for item in items:
        if not matches(type, item.type): continue
        if not matches(biotype, item.attr.get("biotype","")): continue
        result.extend(translate_gene(item, transcript_type, transcript_biotype))
    return result



def ablate(items, ablators, radius):
    result = [ ]
    ablated_any = False
    for ablator in ablators:
        ablator = ablator.shifted(-radius, radius)
        todo = items[:]
        items = [ ]
        while todo:
            item = todo.pop()
            if not item.overlaps(ablator):
                items.append(item)
                continue

            ablated_any = True

            item1 = item.copy()
            item1.children = [ ]
            item1.end = ablator.start
            if item1.start < item1.end:
                items.append(item1)

            item2 = item.copy()
            item2.children = [ ]
            item2.start = ablator.end
            if item2.start < item2.end:
                items.append(item2)

    return items, ablated_any


def dominate(genes, log, radius=20):
    """ Reduce overlaps between exons in different genes. 
        Priority is given to ensemb+havana then havana then ensembl, 
        and then to furthest upstrand-ending transcript (differences less than radius are ignored),
        then to protein_coding genes """
    transcripts = [ ]
    for gene in genes:
        for transcript in gene.children:
            transcript.parent = gene
        transcripts.extend(gene.children)

    transcript_index = span_index.index_annotations(transcripts)

    result = [ ]
    for gene in genes:
        new_gene = gene.copy()
        new_gene.children = [ ]
        for transcript in gene.children:
            new_exons = [ item for item in transcript.children if item.type == "exon" ]
            for hit in transcript_index.get(transcript, same_strand=True):
                if hit.parent == transcript.parent: continue
                if hit.strand > 0:
                    offset = hit.end - transcript.end
                else:
                    offset = transcript.start - hit.start
                
                hit_score = (
                    source_levels.get(hit.source, 0),
                    min(offset+radius,0) if offset < 0 else max(0,offset-radius),
                    hit.attr["Biotype"] != "protein_coding"
                )
                transcript_score = (
                    source_levels.get(transcript.source, 0),
                    0,
                    transcript.attr["Biotype"] != "protein_coding"
                )

                if hit_score < transcript_score:
                    new_exons, ablated_any = ablate(new_exons, 
                        [ item for item in hit.children if item.type == "exon" ],
                        radius)
                    #if ablated_any:
                    #    print hit.attr["Biotype"], transcript.attr["Biotype"], transcript.attr["Name"]

            if new_exons:
                new_transcript = transcript.copy()
                new_transcript.children = [ item for item in transcript.children if item.type != "exon" ] + new_exons
                new_gene.children.append(new_transcript)

        if new_gene.children:
            result.append(new_gene)

    return result


class Union_find(object):
    def __init__(self, items):
        self.items = items
        self.parent = { }
    
    def root(self, item):
        if item not in self.parent: return item
        original = item
        while item in self.parent:
            item = self.parent[item]
        self.parent[original] = item
        return item 

    def join(self, a,b):
        a = self.root(a)
        b = self.root(b)
        if a != b:
            self.parent[a] = b

    def get_sets(self):
        sets = collections.defaultdict(set)
        for item in self.items:
            sets[self.root(item)].add(item)
        return sets.values()


def join_up(items, brief=False):
    items = natural_sorted(items)
    if brief and len(items) > 3:
        items = items[:2] + ["etc"]
    return  "/".join(items)


def merge(genes, log):
    good_genes = [ ]
    exons = [ ]
    gene_exons = { }
    exon_gene = { }
    for gene in genes:
        this_all = all_within(gene)
        seqid = gene.seqid
        strand = gene.strand
        sane = True
        for item in this_all:
            sane = sane and item.seqid == seqid
            sane = sane and item.strand == strand
        
        if not sane:
            log.log("Skipping gene "+str(gene)+" due to inconsistent chromosome or strand.\n")
            continue
        
        this_exons = [ item for item in this_all if item.type == "exon" ]

        good_genes.append(gene)
        gene_exons[gene] = this_exons 
        for exon in this_exons:
            exon_gene[exon] = gene
        exons.extend(this_exons)
        
        
    index = span_index.index_annotations(exons)
    union = Union_find(good_genes)
    for exon in exons:
        gene1 = exon_gene[exon]
        for hit in index.get(exon, same_strand=True):
            gene2 = exon_gene[hit]
            if gene1 != gene2: union.join(gene1,gene2)
    
    
    sets = union.get_sets()    
    counts = collections.defaultdict(int)
    for item in sets:
        counts[len(item)] += 1
    for n in sorted(counts):
        log.log("Merging produced %d sets of %d genes\n" % (counts[n], n))
    
    result = [ ]
    for gene_set in sets:
        gene_list = list(gene_set)
        seqid = gene_list[0].seqid
        strand = gene_list[0].strand        
        source = gene_list[0].source
        type = gene_list[0].type
        
        for gene in gene_set:
            assert gene.seqid == seqid
            assert gene.strand == strand
        
        start = min(gene.start for gene in gene_list)
        end = max(gene.end for gene in gene_list)        
        attr = collections.defaultdict(set)
        for gene in gene_set:
            for name in gene.attr:
                attr[name].add(gene.attr[name])
        for name in attr:
            attr[name] = join_up(attr[name], brief = name == "Name")
        attr = dict(attr)

        if len(gene_set) >= 2:
            print "Merged gene:", attr.get("Name",""), attr.get("Biotype","")
        
        merged_gene = annotation.Annotation(
            seqid = seqid,
            source = source,
            type = type,
            start = start,
            end = end,
            strand = strand,
            attr = attr
            )
        merged_gene.children = [ ]
        for gene in gene_list:
            for transcript in gene.children:
                new_transcript = transcript.copy()
                new_transcript.children = transcript.children
                new_transcript.attr["Parent"] = merged_gene.get_id()
                merged_gene.children.append(new_transcript)
        result.append(merged_gene) 
    
    return result
    



def show_tree(items, log):
    def stack(item):
        if not item.parents: result = ()
        else:
            [parent] = item.parents
            result = stack(parent)
        return result + (item.type+":"+item.attr.get("biotype",""),)
    
    stacks = collections.defaultdict(int)
    for item in items:
        stacks[" ".join(stack(item))] += 1
    
    for item in sorted(stacks):
        log.log("% 8d %s\n"%(stacks[item],item))


def get_genes(items, extractions, log):
    annotation.link_up_annotations(items)

    log.log("Content of downloaded GFF file (type:biotype ... paths):\n")
    show_tree(items, log)
    log.log("\n")        
    
    genes = list()
    for gene_type,gene_biotype,transcript_type,transcript_biotype in extractions:
        this_genes = extract_and_translate(items, gene_type, gene_biotype, transcript_type, transcript_biotype)
        log.log("%d genes from %s/%s/%s/%s\n" % (len(this_genes), gene_type,gene_biotype,transcript_type,transcript_biotype))
        genes.extend(this_genes)
    
    log.log("\n")
    genes = dominate(genes,log)
    log.log("\n")
    genes = merge(genes,log)
    features = [ ]
    for gene in genes:
        features.extend(all_within(gene))
    features.sort(key=lambda i: (i.seqid,i.start))
    return features





# ================ Actual tool =================

@config.help("""\
Create a tail-tools reference directory using Ensembl files.
""", """\
Genome and annotation should be downloaded from the Ensembl servers.

Find the appropriate fasta file in ftp.ensembl.org. The "primary_assembly" is preferred (otherwise "toplevel"), tail-tools is not able to cope with the variant contigs in Ensembl.

For example for release 82 of homo sapiens, you would download the files:

ftp://ftp.ensembl.org/pub/release-82/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
ftp://ftp.ensembl.org/pub/release-82/gff3/homo_sapiens/Homo_sapiens.GRCh38.82.gff3.gz

""")
@config.Bool_flag('index', 'Generate bowtie2 and shrimp indexes. Only disable if re-building reference directory.')
@config.Bool_flag('shrimp', 'Generate SHRiMP colorspace index.')
@config.Bool_flag('bowtie', 'Generate bowtie2 index.')
@config.Bool_flag('star', 'Generate STAR index.')
@config.String_flag('genes',
    'Comma separated list of gene_type/gene_biotype/transcript_type/transcript_biotype. Any of these can be left blank. Regular expressions can be used.'
    )
@config.String_flag('rename',
    'Comma separated list of old=new chromosome renamings.'
    )
@config.Positional('genome', 'Genome FASTA file.')
@config.Positional('annotation', 'Genome annotation in GFF3 format.')
@config.Section('extra', 'Extra genomes and annotations, as for nesoni make-reference.')
class Make_ensembl_reference(config.Action_with_output_dir):
    _workspace_class = reference_directory.Tailtools_reference
    
    index = True
    shrimp = False
    bowtie = True
    star = True
    rename = ''    
    genes = '///'
    
    genome = None
    annotation = None
    extra = [ ]
    
    def run(self):
        extractions = [ ]
        for item in self.genes.split(','):
            extraction = item.split('/')
            assert len(extraction) == 4
            extractions.append(extraction)
            
        rename = { }
        if self.rename:
            for item in self.rename.split(','):
                old,new = item.split('=')
                rename[old] = new

        work = self.get_workspace()        
        
        with workspace.tempspace() as temp:
            items = list(annotation.read_annotations(self.annotation))
            for item in items:
                item.seqid = rename.get(item.seqid, item.seqid)
            annotation.write_gff3(temp/'temp.gff', get_genes(items, extractions, self.log))
            del items
            
            with open(temp/'temp.fa','wb') as f:
                for name,seq in io.read_sequences(self.genome):
                    name = name.split()[0]
                    name = rename.get(name,name)
                    io.write_fasta(f, name, seq)
            
            reference_directory.Make_tt_reference(
                self.output_dir,
                filenames = [ temp/'temp.fa', temp/'temp.gff' ] + self.extra,
                index = self.index, shrimp = self.shrimp, 
                bowtie = self.bowtie, star = self.star
                ).run()
        
        
        
        
        