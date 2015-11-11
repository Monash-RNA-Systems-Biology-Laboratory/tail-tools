

from nesoni import config, workspace, io, annotation, span_index
from . import reference_directory

import collections, re

attr_renaming = {
    "biotype" : "Biotype",
    "description" : "Product",    
}

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


def merge(genes, log):
    exons = [ ]
    gene_exons = { }
    exon_gene = { }
    for gene in genes:
        gene_exons[gene] = [ item for item in all_within(gene) if item.type == "exon" ]
        for exon in gene_exons[gene]:
            exon_gene[exon] = gene
        exons.extend(gene_exons[gene])
        
    index = span_index.index_annotations(exons)
    union = Union_find(genes)
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
            attr[name] = "/".join(natural_sorted(attr[name]))
        attr = dict(attr)
        
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
    genes = merge(genes,log)
    features = [ ]
    for gene in genes:
        features.extend(all_within(gene))
    features.sort(key=lambda i: (i.seqid,i.start))
    return features





# ================ Actual tool =================

@config.help("""\
Create a tail-tools reference directory using Ensembl.

Genome and annotation is downloaded from the Ensembl servers.

Find the appropriate fasta file in ftp.ensembl.org to work out correct flags. The "primary_assembly" is preferred (otherwise "toplevel"), tail-tools is not able to cope with the variant contigs in Ensembl.

For example if the fasta file url was:

ftp://ftp.ensembl.org/pub/release-82/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

use flags

--release 82 --species Homo_sapiens --assembly GRCh38 --dna dna.primary_assembly

""")
@config.Bool_flag('download', 'Download latest Ensembl files. Only disable if re-building reference directory.')
@config.Bool_flag('index', 'Generate bowtie2 and shrimp indexes. Only disable if re-building reference directory.')
@config.String_flag('release',
    'Ensemble release number.'
    )
@config.String_flag('species',
    'Species. First letter is uppercase.'
    )
@config.String_flag('assembly',
    'Name of the genome assembly.'
    )
@config.String_flag('dna', 
    'If available specify "dna.primary_assembly", otherwise "dna.toplevel".'
    )
@config.String_flag('genes',
    'Comma separated list of gene_type/gene_biotype/transcript_type/transcript_biotype. biotypes can be left blank.'
    )
@config.String_flag('rename',
    'Comma separated list of old=new chromosome renamings.'
    )
class Make_ensembl_reference(config.Action_with_output_dir):
    _workspace_class = reference_directory.Tailtools_reference
    
    download = True
    index = True
    
    release = None
    species = None
    assembly = None
    dna = None
    
    rename = ''
    
    genes = "gene/protein_coding/transcript/protein_coding"
    
    def run(self):
        assert self.release
        assert self.species
        assert self.assembly
        assert self.dna
        
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
        ensembl = workspace.Workspace(work/'ensembl')
        
        genome_filename = self.species+"."+self.assembly+"."+self.dna+".fa.gz"
        genome_url = "rsync://ftp.ensembl.org/ensembl/pub/release-"+self.release+"/fasta/"+self.species.lower()+"/dna/"+genome_filename
        
        gff_filename = self.species+"."+self.assembly+"."+self.release+".gff3.gz"
        gff_url = "rsync://ftp.ensembl.org/ensembl/pub/release-"+self.release+"/gff3/"+self.species.lower()+"/"+gff_filename
        
        
        if self.download:
            self.log.log("Fetching "+genome_url+"\n")
            io.execute(['rsync','-aP',genome_url, ensembl/genome_filename])
            self.log.log("Fetching "+gff_url+"\n")
            io.execute(['rsync','-aP',gff_url, ensembl/gff_filename])
        
        with workspace.tempspace() as temp:
            items = list(annotation.read_annotations(ensembl/gff_filename))
            for item in items:
                item.seqid = rename.get(item.seqid, item.seqid)
            annotation.write_gff3(temp/'temp.gff', get_genes(items, extractions, self.log))
            del items
            
            with open(temp/'temp.fa','wb') as f:
                for name,seq in io.read_sequences(ensembl/genome_filename):
                    name = name.split()[0]
                    name = rename.get(name,name)
                    io.write_fasta(f, name, seq)
            
            reference_directory.Make_tt_reference(
                self.output_dir,
                filenames = [ temp/'temp.fa', temp/'temp.gff' ],
                index = self.index,
                ).run()
        
        
        
        
        