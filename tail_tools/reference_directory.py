
import os, collections, gzip, re, sys

import nesoni
from nesoni import reference_directory, workspace, io, config, annotation, annotation_tools, span_index, grace

def natural_sorted(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)



class Tailtools_reference(reference_directory.Reference):
    VERSION = 0

    def __init__(self, working_dir, must_exist):
        super(Tailtools_reference,self).__init__(working_dir, must_exist)
        
        if must_exist:
            param = self.param
            assert 'tail_tools_reference_version' in param, 'Not a tail-tools reference directory.'
            assert param['tail_tools_reference_version'] >= self.VERSION, 'Reference directory from previous version of tail-tools, please rebuild.'
        
        #make-tt should do this
        #else:
        #    self.update_param(tail_tools_reference_version = VERSION)



def _max_extension(exon, exon_index, mrna_end_index):
    # Get maximum extension until hit a CDS, arbitrarily max out at 10000
    # Also do not extend through the end of another mRNA
    max_extension = 10000
    exon_end = exon.three_prime()
    for hit in exon_index.get(exon_end.shifted(0,10000), same_strand=True):
        max_extension = min(max_extension, hit.relative_to(exon_end).start)
    for hit in mrna_end_index.get(exon_end.shifted(2,10000), same_strand=True):
        max_extension = min(max_extension, hit.relative_to(exon_end).start)
    return max(0,max_extension)


@config.help(
    'Create a tail tools reference directory. '
    'See the requirements for this in the README, they are quite strict.',
    'You need a FASTA file containing sequences and a GFF file containing annotations.'
    )
@config.Bool_flag('index', 'Generate indexes. Only disable if re-building reference directory.')
@config.Bool_flag('shrimp', 'Generate SHRiMP colorspace index.')
@config.Bool_flag('bowtie', 'Generate bowtie2 index.')
@config.Bool_flag('star', 'Generate STAR index.')
@config.Main_section('filenames', 'Sequence and annotation files.')
class Make_tt_reference(config.Action_with_output_dir):
    _workspace_class = Tailtools_reference

    index = True
    shrimp = False
    bowtie = False
    star = True
    filenames = [ ]
    
    def run(self):    
        work = self.get_workspace()
        work.update_param(remove=['tail_tools_reference_version'])
        
        nesoni.Make_reference(
            self.output_dir,
            filenames = self.filenames,
            snpeff = False,
            cs = 'ifavailable' if self.index and self.shrimp else False,
            ls = False,
            bowtie = 'ifavailable' if self.index and self.bowtie else False,
            ).run()
            
        annotations = list(annotation.read_annotations(work/'reference.gff'))
        annotation.link_up_annotations(annotations)
        
        exon_index = span_index.index_annotations([
            item for item in annotations if item.type == "exon"
            ])
        mrna_end_index = span_index.index_annotations([
            item.three_prime() for item in annotations if item.type == "mRNA"
            ])
        
        mrna_utrs = [ ]
        gene_utrs = [ ]
        
        for gene in annotations:
            if gene.type != 'gene': continue

            mrnas = [ item for item in gene.children if item.type == 'mRNA' ]
            assert mrnas, "Gene without any mRNAs: "+gene.get_id()

            gene.attr['color'] = '#880088'
            gene.start = min(item.start for item in mrnas)
            gene.end = max(item.end for item in mrnas)
            gene.attr["max_extension"] = str(_max_extension(gene, exon_index, mrna_end_index))
        
            gene_utr_5primes = [ ]
            
            for mrna in mrnas:
                assert mrna.strand == gene.strand, mrna
                assert mrna.seqid == gene.seqid, mrna
                
                mrna.attr["max_extension"] = str(_max_extension(mrna, exon_index, mrna_end_index))
            
                cdss = [ item for item in mrna.children if item.type == 'CDS' ]
                exons = [ item for item in mrna.children if item.type == 'exon' ]
                
                if not exons: continue
                
                #link up annotations sorts children, so final is really final
                for item in exons[:-1]:
                    item.attr["max_extension"] = "0"
                exons[-1].attr["max_extension"] = mrna.attr["max_extension"]
                
                if not cdss: continue
                
                mrna_utr_5primes = [ ]
                if gene.strand >= 0:
                   cds_3prime = max(item.end for item in cdss)
                   for item in exons:
                       if item.end >= cds_3prime:
                           mrna_utr_5primes.append(max(item.start,cds_3prime))
                else:
                   cds_3prime = min(item.start for item in cdss)
                   for item in exons:
                       if item.start <= cds_3prime:
                           mrna_utr_5primes.append(min(item.end,cds_3prime))
                
                if mrna.strand >= 0:
                    utr_start = min(mrna_utr_5primes) if mrna_utr_5primes else mrna.end
                    utr_end = max(utr_start+1,mrna.end)
                    gene_utr_5primes.append(utr_start)
                else:
                    utr_end = max(mrna_utr_5primes) if mrna_utr_5primes else mrna.start
                    utr_start = min(mrna.start,utr_end-1)
                    gene_utr_5primes.append(utr_end)
                
                attr = mrna.attr.copy()
                attr['Parent'] = attr['ID']
                attr['ID'] = attr['ID']+'-3UTR'
                attr['color'] = '#008888'
                utr = annotation.Annotation(
                    source = 'tt',
                    type = 'three_prime_utr',
                    seqid = mrna.seqid,
                    strand = mrna.strand,
                    start = utr_start,
                    end = utr_end,
                    attr = attr,
                    )
                max_ext = _max_extension(utr, exon_index, mrna_end_index)
                utr.attr["max_extension"] = str(max_ext)
                #Only include if there is an annotated 3' UTR or end is not in the middle of some other isoform's exon, and it has position > 0
                if utr_end-utr_start+max_ext > 1 and utr_start >= 0:
                    mrna_utrs.append(utr)
            
            if gene.strand >= 0:
                utr_start = max(gene_utr_5primes) if gene_utr_5primes else gene.end
                utr_end = max(utr_start+1,gene.end)
            else:
                utr_end = min(gene_utr_5primes) if gene_utr_5primes else gene.start
                utr_start = min(gene.start,utr_end-1)
            
            attr = gene.attr.copy()
            attr['Parent'] = attr['ID']
            attr['ID'] = attr['ID']+'-3UTR'
            attr['color'] = '#008888'
            utr = annotation.Annotation(
                source = 'tt',
                type = 'three_prime_utr',
                seqid = gene.seqid,
                strand = gene.strand,
                start = utr_start,
                end = utr_end,
                attr = attr,
                )
            utr.attr["max_extension"] = str(_max_extension(utr, exon_index, mrna_end_index))
            gene_utrs.append(utr)
        
        annotation.write_gff3(work/'reference.gff', annotations + mrna_utrs)
        annotation.write_gff3(work/'utr.gff', gene_utrs)
        
        if self.index and self.star and grace.can_execute("STAR"):
            star_work = workspace.Workspace(work/'star')
            io.execute([
                'STAR','--runMode','genomeGenerate',
                '--outFileNamePrefix',star_work.working_dir+'/',
                '--genomeDir',star_work.working_dir,
                '--genomeFastaFiles',work/'reference.fa',
                '--sjdbGTFfile',work/'reference.gff',
                '--sjdbGTFtagExonParentTranscript','Parent',
                '--sjdbOverhang','199',
                ])
            
        work.update_param(tail_tools_reference_version=work.VERSION)
        
        
@config.help(
    'Create a tail-tools reference directory using a UCSC-browser genome.',
    'Any gene prediction track may be used, consult the UCSC table browser for options:\n'
    'http://genome.ucsc.edu/cgi-bin/hgTables\n'
    '\n'
    'The default is to use the refGene table containing genes from the NCBI RefSeq database, '
    'which is available for most genomes, '
    'but notably not for yeast genomes (see below).\n'
    '\n'    
    'mRNAs are named after the "name" field in the gene prediction table, '
    'forcibly uniquified where this field is non-unique.\n'
    '\n'
    'For yeast we suggest using:\n'
    '\n'
    '  tail-tools make-ucsc-reference: sacCer3 sacCer3 \\\n'
    '      --table sgdGene \\\n'
    '      --name name,sgdToName,name,value \\\n'
    '      --product name,sgdDescription,name,description\n'
    )
@config.Bool_flag('download', 'Download latest UCSC files. Only disable if re-building reference directory.')
@config.Bool_flag('index', 'Generate bowtie2 and shrimp indexes. Only disable if re-building reference directory.')
@config.String_flag('table', 'UCSC table containing mRNAs.')
@config.String_flag('name', 'How to obtain a friendly gene name. Blank for none, otherwise "field" or "join_field,table_name,joined_field,result_field".')
@config.String_flag('product', 'How to obtain a gene description. Blank for none, otherwise "field" or "join_field,table_name,joined_field,result_field".')
@config.Positional('ucsc_name', 'UCSC genome, eg "hg19".')
class Make_ucsc_reference(config.Action_with_output_dir):
    _workspace_class = Tailtools_reference
    
    download = True
    index = True
    table = 'refGene'
    name = 'name,refLink,mrnaAcc,name'
    product = 'name,refLink,mrnaAcc,product'

    ucsc_name = None

    def run(self):
        assert self.ucsc_name, 'Need a UCSC genome name'
        
        scratch = _ucsc_scratch(self)
        
        # Load annotations
        
        source = 'tt-ucsc-%s-%s' % (self.ucsc_name, self.table)
        
        table = scratch.get_table(self.table)
        get_name = scratch.getter(self.name)
        get_product = scratch.getter(self.product)

        mrnas = [ ]
        
        for item in table:
            ann = annotation.Annotation(
                seqid = item.chrom,
                source = source,
                type = 'mRNA',
                strand = {'+':1, '-':-1}[item.strand],
                start = int(item.txStart),
                end = int(item.txEnd),
                attr = {
                    'ID' : item.name,
                    'Name' : get_name(item),
                    'Product' : get_product(item),
                    #'UCSC_name2' : item.name2,
                    }
                )
            
            ann.record = item
            mrnas.append(ann)

        _uniquify_ids(mrnas)
        
        annotations = [ ]
        
        for group in _grouped_features(mrnas):
            ID = '/'.join(item.attr['ID'] for item in group)
            for item in group:
                item.attr['Parent'] = ID
                item.attr['ID'] = item.attr['ID'] + '-mRNA'
            
            annotations.append(annotation.Annotation(
                source = source,
                type = 'gene',
                seqid = group[0].seqid,
                strand = group[0].strand,
                start = min(item.start for item in group),
                end = max(item.end for item in group),
                attr = {
                    'ID' : ID,
                    'Name' : annotation_tools.join_descriptions([ item.attr['Name'] for item in group ], '/'),
                    'Product' : annotation_tools.join_descriptions([ item.attr['Product'] for item in group ], '/'),
                    #'UCSC_name2' : annotation_tools.join_descriptions([ item.attr['UCSC_name2'] for item in group ], '/'),
                    }
                ))
            for item in group:
                annotations.append(item)
                
                exonStarts = _parse_ints(item.record.exonStarts)
                exonEnds = _parse_ints(item.record.exonEnds)
                cdsStart = int(item.record.cdsStart)
                cdsEnd = int(item.record.cdsEnd)
                for start,end in zip(exonStarts,exonEnds):
                    annotations.append(annotation.Annotation(
                        source = source,
                        type = 'exon',
                        seqid = item.seqid,
                        strand = item.strand,
                        start = start,
                        end = end,
                        attr = {
                            'Parent' : item.attr['ID'],
                            }
                        ))
                    if max(cdsStart,start) < min(cdsEnd,end):
                        annotations.append(annotation.Annotation(
                            source = source,
                            type = 'CDS',
                            seqid = item.seqid,
                            strand = item.strand,
                            start = max(cdsStart,start),
                            end = min(cdsEnd,end),
                            #TODO: phase
                            attr = {
                                'Parent' : item.attr['ID'],
                                }
                            ))

        # Load sequence
        
        if self.download:
            io.execute(['rsync','-P','rsync://hgdownload.cse.ucsc.edu/goldenPath/'+self.ucsc_name+'/bigZips/chromFa.tar.gz',scratch.ucsc/'chromFa.tar.gz'])
        
        with workspace.tempspace() as temp:
            io.execute(['tar','-C',temp.working_dir,'-zxf',scratch.ucsc/'chromFa.tar.gz'])
            sequences = [ temp/item for item in natural_sorted(os.listdir(temp.working_dir)) ]
            
            with open(temp/'reference.gff','wb') as f:
                annotation.write_gff3_header(f)
                for item in annotations:
                    print >> f, item.as_gff()
            
            Make_tt_reference(
                self.output_dir,
                filenames = sequences + [ temp/'reference.gff' ],
                index = self.index,
                ).run()
            


class _ucsc_scratch(object):
    def __init__(self, action):
        self.action = action
        self.work = action.get_workspace()        
        self.ucsc = workspace.Workspace(self.work/'ucsc')
        self.tables = { }

    def get_table(self, table_name):
        if table_name not in self.tables:
            if self.action.download: 
                for filename in [ table_name+'.txt.gz', table_name+'.sql' ]:
                    io.execute(['rsync','-P','rsync://hgdownload.cse.ucsc.edu/goldenPath/'+self.action.ucsc_name+'/database/'+filename, self.ucsc/filename])
            
            fields = [ ]
            with open(self.ucsc/table_name+'.sql','rU') as f:
                for line in f:
                    if line.startswith('  `'):
                        parts = line.strip().split()
                        assert parts[0][0] == '`' and parts[0][-1] == '`'
                        fields.append(parts[0][1:-1])

            tup_class = collections.namedtuple(table_name, fields)
            data = [ ]
            with gzip.open(self.ucsc/table_name+'.txt.gz','rb') as f:
                for line in f:
                    data.append(tup_class(* line.rstrip('\n').split('\t') ))
            self.tables[table_name] = data
         
        return self.tables[table_name]
    
    def getter(self, spec):
        if not spec:
            return lambda record: ''
        
        if ',' not in spec:
            return lambda record: getattr(record,spec)
    
        join_field, table_name, table_join_field, result_field = spec.split(',')
        index = collections.defaultdict(list)
        for item in self.get_table(table_name):
            index[getattr(item,table_join_field)].append(getattr(item,result_field))
        def get(record):
            return annotation_tools.join_descriptions(index[getattr(record,join_field)],'/')
        return get


def _uniquify_ids(features):
    count = collections.defaultdict(int)
    for item in features: 
        count[item.attr['ID']] += 1
    result = [ ]
    count2 = collections.defaultdict(int)
    for item in features:
        ID = item.attr['ID']
        count2[ID] += 1
        if count[ID] != 1:
            item.attr['ID'] = ID + '--%d-of-%d'%(count2[ID],count[ID])


def _grouped_features(features):    
    def get_key(item):
        return (item.seqid, item.strand)
    features = sorted(features, key=lambda item: (get_key(item), item.start))
            
    group = [ ]
    groups = [ ]
    def emit():
        if not group: return
        groups.append(group[:])
        del group[:]        
    key = None
    end = 0
    for item in features:
        if get_key(item) != key or item.start >= end:
            emit()
            key = get_key(item)
            end = item.end
        group.append(item)
        end = max(item.end, end)
    emit()
    
    return groups
  
    
def _parse_ints(text):
    return [ int(item) for item in  text.split(',')[:-1] ]


