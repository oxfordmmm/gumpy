import gzip, os, pickle, time, copy, re
from collections import defaultdict

import numpy

import pysam
from Bio import SeqIO
from tqdm import tqdm

from gumpy import Gene
from gumpy import Genotype


class Genome(object):

    def __init__(self,genbank_file=None,fasta_file=None,show_progress_bar=False,default_promoter_length=100,name=None):

        '''
        Instantiates a genome object by loading a VCF file and storing the whole genome as a numpy array

        Args:
            genbank_file (str): path to the GenBank file to build the reference genome
            fasta_file (str): path to the FASTA file to build the reference genome
        '''

        assert ((genbank_file is not None) or (fasta_file is not None)), "one of a GenBank file or a FASTA file  must be specified!"
        assert default_promoter_length>=0, "the promoter length must be a positive integer!"
        assert isinstance(default_promoter_length,int), "the promoter length must be a positive integer!"
        assert name is not None, "you must specify a name for this genome!"
        self.name=name

        self.id=""
        self.organism=""
        self.sample_metadata={}

        # load the specified GenBank file
        if genbank_file is not None:

            # create the genbank file and store in a BioPython object
            reference_genome=SeqIO.read(genbank_file,'genbank')

            # convert to a numpy array at the first opportunity since slicing BioPython is between 10 and 50,000 times slower!
            self.genome_coding_strand=numpy.array([i.lower() for i in str(reference_genome.seq)])

            # create the complementary strand upfront
            self.genome_noncoding_strand=self._complement(self.genome_coding_strand)

            # store the length of the genome
            self.genome_length=len(self.genome_coding_strand)

            # create an array of the genome indices
            self.genome_index=numpy.arange(1,self.genome_length+1)

            self.version=reference_genome.annotations['accessions'][0]+"."+str(reference_genome.annotations['sequence_version'])

            # store some of the metadata, if it is present
            self.id=reference_genome.id

            if 'organism' in reference_genome.annotations.keys():
                self.organism=reference_genome.annotations['organism']
            if 'sequence_version' in reference_genome.annotations.keys():
                self.sample_metadata['SEQUENCE_VERSION']=reference_genome.annotations['sequence_version']
            if 'source' in reference_genome.annotations.keys():
                self.sample_metadata['SOURCE']=reference_genome.annotations['source']
            if 'taxonomy' in reference_genome.annotations.keys():
                self.sample_metadata['TAXONOMY']=reference_genome.annotations['taxonomy']

            self.genome_feature_name=numpy.zeros(self.genome_length,dtype="<U10")
            self.genome_feature_type=numpy.zeros(self.genome_length,dtype="<U5")
            self.genome_is_cds=numpy.zeros(self.genome_length,dtype=bool)
            self.genome_is_promoter=numpy.zeros(self.genome_length,dtype=bool)
            self.genome_on_noncoding_strand=numpy.zeros(self.genome_length,dtype=bool)

            self._gene_type=defaultdict(str)
            self._gene_codes_protein=defaultdict(str)
            self.genes=defaultdict(str)

            self.genome_sequence=copy.deepcopy(self.genome_coding_strand)
            self.genome_numbering=numpy.zeros(self.genome_length,int)
            self.genome_positions=numpy.zeros(self.genome_length,int)

            self.is_indel=numpy.zeros(self.genome_length,dtype=bool)
            self.indel_length=numpy.zeros(self.genome_length,int)

            previous_gene_reversed=False

            # go through the GenBank file record-by-record
            for record in tqdm(reference_genome.features,disable=not(show_progress_bar)):

                if record.type in ['CDS','rRNA']:

                    gene_name=None
                    type=None

                    if 'gene' in record.qualifiers.keys():

                        gene_name=record.qualifiers['gene'][0]

                        if record.type=='rRNA':
                            type="RNA"
                        else:
                            type="GENE"

                    elif 'locus_tag' in record.qualifiers.keys():

                        gene_name=record.qualifiers['locus_tag'][0]

                        if record.type=='rRNA':
                            type="RNA"
                        else:
                            type="LOCUS"

                    else:
                        continue

                    if gene_name is not None:

                        gene_start=int(record.location.start)
                        gene_end=int(record.location.end)
                        if type in ["GENE","LOCUS"]:
                            codes_protein=True
                        else:
                            codes_protein=False

                        if record.strand==1:

                            # This is a bit hacky; we are assuming that the genes are being read from the Genbank file
                            # in sequential order, hence for +ve strand genes the start of the next may overwrite the end of the previous,
                            # which is what we want. But the next gene after a -ve strand gene may overwrite the _start_ of that gene, which we don't want
                            # hence the logic to remember if the previous gene is on the reverse strand or not!

                            if previous_gene_reversed:
                                gene_mask=(self.genome_index>gene_start) & (self.genome_index<=gene_end) & (~self.genome_is_cds)
                            else:
                                gene_mask=(self.genome_index>gene_start) & (self.genome_index<=gene_end)

                            previous_gene_reversed=False

                            promoter_mask=(self.genome_index>gene_start-default_promoter_length) & (self.genome_index<=gene_start) & (~self.genome_is_cds) & (~self.genome_is_promoter)

                            mask=gene_mask+promoter_mask

                            # be paranoid and set even though the default value is False
                            self.genome_on_noncoding_strand[mask]=False

                            self.genome_sequence[mask]=self.genome_coding_strand[mask]

                            promoter_coding_numbering=self.genome_index[promoter_mask]-gene_start-1

                            if codes_protein:
                                gene_coding_numbering=numpy.floor_divide(self.genome_index[gene_mask]-gene_start+2,3)
                            else:
                                gene_coding_numbering=self.genome_index[gene_mask]-gene_start

                            gene_coding_position=self.genome_index[gene_mask]-gene_start

                        elif record.strand==-1:

                            if previous_gene_reversed:
                                gene_mask=(self.genome_index>gene_start) & (self.genome_index<=gene_end) & (~self.genome_is_cds)
                            else:
                                gene_mask=(self.genome_index>gene_start) & (self.genome_index<=gene_end)

                            previous_gene_reversed=True

                            promoter_mask=(self.genome_index>gene_end) & (self.genome_index<=gene_end+default_promoter_length) & (~self.genome_is_cds) & (~self.genome_is_promoter)

                            mask=gene_mask+promoter_mask

                            # the default value is False, so only need to set those which are reversed
                            self.genome_on_noncoding_strand[mask]=True

                            # replace the coding sequence with the complement
                            self.genome_sequence[mask]=self.genome_noncoding_strand[mask]

                            promoter_coding_numbering=-1*(self.genome_index[promoter_mask]-gene_end)

                            if codes_protein:
                                gene_coding_numbering=-1*(numpy.floor_divide(self.genome_index[gene_mask]-gene_end-1,3))
                            else:
                                gene_coding_numbering=-1*(self.genome_index[gene_mask]-gene_end)

                            gene_coding_position=-1*(self.genome_index[gene_mask]-gene_end)
                        else:
                            raise TypeError("gene in GenBank file has strand that is not 1 or -1")


                        self.genome_feature_type[mask]=type
                        self.genome_feature_name[mask]=gene_name

                        self.genome_numbering[gene_mask]=gene_coding_numbering
                        self.genome_numbering[promoter_mask]=promoter_coding_numbering

                        promoter_coding_position=promoter_coding_numbering

                        self.genome_positions[gene_mask]=gene_coding_position
                        self.genome_positions[promoter_mask]=promoter_coding_position

                        self.genome_is_cds[gene_mask]=True
                        self.genome_is_promoter[gene_mask]=False
                        self.genome_is_promoter[promoter_mask]=True

                        self._gene_type[gene_name]=type
                        self._gene_codes_protein[gene_name]=codes_protein

            # store a list of all the gene names
            self.gene_names=numpy.unique(self.genome_feature_name[self.genome_feature_name!=""])

            # # pass ALL the gene names to create all the Gene objects for the first time
            self._recreate_genes(self.gene_names,show_progress_bar=show_progress_bar)

        # otherwise there must be a FASTA file so load that instead
        elif fasta_file is not None:

            header,nucleotide_sequence=self._load_fastafile(fasta_file)

            nucleotide_sequence=nucleotide_sequence.lower()

            cols=header[1:].split("|")
            if len(cols)>1:
                self.id=cols[0]
                self.organism=cols[1]
                self.name=cols[2]
            # if len(cols)>3:
            #     self.additional_metadata=cols[3]

            self.genome_coding_strand=numpy.array(list(nucleotide_sequence))

            self.genome_noncoding_strand=self._complement(self.genome_coding_strand)

            # store the length of the genome
            self.genome_length=len(self.genome_coding_strand)

            # create an array of the genome indices
            self.genome_index=numpy.arange(1,self.genome_length+1)

        # insist that bases are lower case
        self.genome_coding_strand=numpy.char.lower(self.genome_coding_strand)

        # store the sequence as integers 0,1,2,3 with which bases they refer to in bases_integer_lookup
        self.bases_to_integer, self.genome_coding_integers = numpy.unique(self.genome_coding_strand, return_inverse=True)

    def _recreate_genes(self,list_of_genes,show_progress_bar=False):
        """
        Private method to re-instantiate the passed list Genes.

        This translates the nucleotide sequence into amino acids (if the gene codes protein) and is
        hence necessary after applying a vcf file, albeit only for those genes whose sequence has been altered.

        Args:
            list_of_genes (list)            list of genes to recreate
            show_progress_bar (True/False)  whether to show the (tqdm) progress bar

        Returns:
            None
        """
        # pass ALL the gene names to create all the Gene objects for the first time
        for gene in tqdm(list_of_genes,disable=not(show_progress_bar)):

            mask=self.genome_feature_name==gene

            assert numpy.count_nonzero(mask)>0, "gene not found in genome!"

            self.genes[gene]=Gene(  gene_name=gene,\
                                    sequence=self.genome_sequence[mask],\
                                    index=self.genome_index[mask],\
                                    numbering=self.genome_numbering[mask],\
                                    positions=self.genome_positions[mask],
                                    is_indel=self.is_indel[mask],
                                    indel_length=self.indel_length[mask],
                                    codes_protein=self._gene_codes_protein[gene],\
                                    feature_type=self._gene_type[gene]  )


    def __repr__(self):

        '''
        Overload the print function to write a summary of the genome.
        '''

        output=""
        if hasattr(self,'id'):
            output+=self.id+"\n"
        if hasattr(self,'organism'):
            output+=self.organism+"\n"
        if hasattr(self,'name'):
            output+=self.name+"\n"
        output+=str(self.genome_length)+" bases\n"
        output+=''.join(i for i in self.genome_coding_strand[0:3])
        output+="..."
        output+=''.join(i for i in self.genome_coding_strand[-3:])

        return(output)

    def list_variants_wrt(self,other):

        assert self.genome_length==other.genome_length, "genomes must have the same length!"

        mask=self.genome_sequence!=other.genome_sequence

        ref=other.genome_sequence[mask]
        idx=self.genome_index[mask]
        alt=self.genome_sequence[mask]

        variants=[]
        for (r,i,a) in zip(ref,idx,alt):
            variants.append(str(i)+r+">"+a)

        mask=self.is_indel
        idx=self.genome_index[mask]
        length=self.indel_length[mask]

        for (i,l) in zip(idx,length):
            variants.append(str(i)+"_indel_"+str(l))

        return(variants)


    def __sub__(self,other):

        """
        Overload the subtraction operator so it returns a tuple of the differences between the two genomes
        """

        assert self.genome_length==other.genome_length, "genomes must have the same length!"

        mask=self.genome_coding_strand!=other.genome_coding_strand

        mask+=self.indel_length!=0

        return(self.genome_index[mask])

    @staticmethod
    def _complement(nucleotides_array):
        """
        Simple private method for returning the complement of an array of bases.

        Note that takes account of HET and NULL calls via z and x, respectively
        """

        complementary_bases = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'x':'x', 'z':'z'}

        complement=[complementary_bases[i] for i in nucleotides_array]

        return(numpy.array(complement))

    def contains_gene(self,gene_name):
        '''
        Simply checks to see if the specified gene exists in the Genome object.

        Args:
            gene_name (str) e.g. katG

        Returns:
            True/False
        '''

        # convert from numpy array to Python list to avoid FutureWarning about comparing numpy arrays and Python objects
        # https://stackoverflow.com/questions/40659212/futurewarning-elementwise-comparison-failed-returning-scalar-but-in-the-futur
        list_of_genes=list(self.gene_names)

        if gene_name in list_of_genes:
            return(True)
        else:
            return(False)


    def at_index(self,index):

        assert index > 0, "index must be a positive integer!"
        assert isinstance(index,int), "index must be a positive integer!"
        assert index <= self.genome_length, "index must be less than the length of the genome!"

        mask=self.genome_index==index

        putative_gene=self.genome_feature_name[mask]

        if putative_gene!="":
            return putative_gene[0]
        else:
            return None

    def snp_distance(self,other):
        return (numpy.count_nonzero(self.genome_coding_strand!=other.genome_coding_strand))

    def apply_vcf_file(self,vcf_file=None,ignore_filter=False, ignore_status=False,show_progress_bar=False,metadata_fields=None):
        """
        Load a VCF file and apply the variants to the whole genome sequence.

        Args:
            vcf_file (str): path to the VCF file to be loaded and applied
            ignore_filter (bool): whether to ignore the FILTER column in the VCF file (Clockwork hasn't always written it correctly)
            ignore_status (bool): ditto
            show_progress_bar (bool): whether to draw a nice tqdm progress bar (False by default)
        """

        # if we are showing a TQDM progress bar, count the rows in the VCF file (0.3s overhead)
        if show_progress_bar:
            try:
                lines_in_vcf=len(open(vcf_file).readlines())
            except:
                lines_in_vcf=None
        else:
            lines_in_vcf=None

        # since we are now applying a VCF file, it makes sense to create these numpy arrays
        self.coverage=numpy.zeros(self.genome_length,int)
        self.indel_ref=numpy.zeros(self.genome_length,dtype='<U50')
        self.indel_alt=numpy.zeros(self.genome_length,dtype='<U50')

        # create a set of mutually exclusive Boolean arrays that tell you what the 'single sequence' result is
        self.is_ref=numpy.zeros(self.genome_length,dtype=bool)
        self.is_null=numpy.zeros(self.genome_length,dtype=bool)
        self.is_het=numpy.zeros(self.genome_length,dtype=bool)
        self.is_snp=numpy.zeros(self.genome_length,dtype=bool)


        # set up a dictionary for the metadata since the field names will vary between calling calling codes
        # for Clockwork these will be GT_CONF and GT_CONF_PERCENTILE
        self.metadata_fields=metadata_fields
        if self.metadata_fields is not None:
            self.genome_sequence_metadata={}
            for field in self.metadata_fields:
                self.genome_sequence_metadata[field]=numpy.zeros(self.genome_length,float)

        # to deal with HET calls we need to setup some diploid arrays
        self.het_variations=numpy.zeros((self.genome_length,2),str)
        self.het_coverage=numpy.zeros((self.genome_length,2),int)
        self.het_indel_length=numpy.zeros((self.genome_length,2),int)
        self.het_ref=numpy.zeros(self.genome_length,dtype='<U50')
        self.het_alt=numpy.zeros((self.genome_length,2),dtype='<U50')

        # split and remember the path, filename and stem of the VCF file
        (self.vcf_folder,self.vcf_file_name)=os.path.split(vcf_file)
        self.vcf_file_stem, file_extension = os.path.splitext(self.vcf_file_name)

        # it may have been compressed, in which case there will be TWO fileextensions to remove
        if self.vcf_file_stem[-4:]==".vcf":
            self.vcf_file_stem, file_extension = os.path.splitext(self.vcf_file_stem)

        # assume that the sample name is the filestem and remember
        self.name=self.vcf_file_stem

        self.genes_mutated=[]

        # open the supplied VCF file
        # note that this will read in bgzip compressed vcf files (from htslib) but not gzip compressed files
        # even though the file extension is the same
        vcf_reader = pysam.VariantFile(vcf_file.rstrip())

        # now iterate through the records found in the VCF file
        for record in tqdm(vcf_reader,disable=not(show_progress_bar),total=lines_in_vcf):

            # check to see the filter is ok (or we are ignoring it)
            if self._is_record_invalid(ignore_filter,record):
                continue

            # cope with multiple entries in a row
            for sample_idx, (sample_name, sample_info) in enumerate(
                record.samples.items()
            ):

                # check to see if the status is ok (or we are ignoring it)
                if not ignore_status and sample_info["STATUS"] == "FAIL":
                    continue

                # ugly; deals with a problem with Minos/Clockwork getting the GT in the wrong place
                try:
                    genotype = Genotype(*sample_info["GT"])
                except TypeError as err:
                    genotype = self._minos_gt_in_wrong_position_fix(record, sample_idx)
                    if genotype is None:
                        raise err

                # return the call
                ref_bases,index,alt_bases = self._get_variant_for_genotype_in_vcf_record(genotype, record)

                # bypass (for speed) if this is a REF call
                if alt_bases=="":
                    continue

                # deal with everything except HET calls
                if not isinstance(alt_bases,tuple):

                    # one or more SNPs (this will naturally catch NULLs as well)
                    if len(ref_bases)==len(alt_bases):

                        for before,after in zip(ref_bases,alt_bases):

                            # only make a change if the ALT is different to the REF
                            if before!=after:

                                # find out the coverage
                                coverage=sample_info['COV'][genotype.call1]

                                # record any additional metadata
                                self._set_sequence_metadata(index,sample_info)

                                # make the mutation
                                self._permute_sequence(index,coverage,after=after)

                            # increment the position in the genome
                            index+=1

                    # an INDEL
                    else:

                        # calculate the length of the indel
                        indel_length=len(alt_bases)-len(ref_bases)

                        assert indel_length!=0, "REF: "+ref_bases+" and ALT: "+alt_bases+" same length?"

                        # find out the coverage
                        coverage=sample_info['COV'][genotype.call1]

                        # record any additional metadata
                        self._set_sequence_metadata(index,sample_info)

                        # make the mutation
                        self._permute_sequence(index,coverage,indel_length=indel_length,indel_bases=(ref_bases,alt_bases))

                # HET calls
                else:

                    # alt_bases is now a 2-tuple, so iterate
                    for strand,alt in enumerate(alt_bases):

                        # one or more SNPs
                        if len(ref_bases)==len(alt):

                            # have to create a copy of index so it is unaltered for the next strand...
                            idx=index

                            # walk down the bases
                            for before,after in zip(ref_bases,alt):

                                # calculate a Boolean mask identifying where we are in the genome
                                mask=self.genome_index==idx

                                # record any additional metadata
                                self._set_sequence_metadata(idx,sample_info)

                                # remember the coverage in the diploid representation for this het
                                self.het_coverage[(mask,strand)]=sample_info['COV'][genotype.call()[strand]]

                                # only record a SNP if there is a change
                                if (before!=after):
                                    self.het_variations[(mask,strand)]=after
                                    self.het_ref[mask]=ref_bases
                                    self.het_alt[(mask,0)]=alt_bases[0]
                                    self.het_alt[(mask,1)]=alt_bases[1]
                                idx+=1
                        else:

                            # calculate the length of the indel
                            indel_length=len(alt)-len(ref_bases)

                            # calculate a Boolean mask identifying where we are in the genome
                            mask=self.genome_index==index

                            # record any additional metadata
                            self._set_sequence_metadata(index,sample_info)

                            # remember the het indel
                            self.het_coverage[(mask,strand)]=sample_info['COV'][genotype.call()[strand]]
                            self.het_indel_length[(mask,strand)]=indel_length
                            self.het_variations[(mask,strand)]="i"
                            self.het_ref[mask]=ref_bases
                            self.het_alt[(mask,0)]=alt_bases[0]
                            self.het_alt[(mask,1)]=alt_bases[1]

        # now that we've parsed the VCF file, and hence all the HETs, we need to update the main sequence
        # to show that there are HETs

        # pick out all genome locations where one of the diploid sequences has been altered
        all_hets_mask=(self.het_variations[:,0]!="") | (self.het_variations[:,1]!="")

        # iterate through the genome positions
        for idx in self.genome_index[all_hets_mask]:

            # where are we?
            mask=self.genome_index==idx

            # define the total coverage as the sum of the two HET calls
            coverage=numpy.sum(self.het_coverage[mask])

            # make the mutation, identifying this as a HET call
            self._permute_sequence(idx,coverage,after='z')

        # first recompute the complementary strand
        self.genome_noncoding_strand=self._complement(self.genome_coding_strand)

        # reintialise the coding sequence
        self.genome_sequence=copy.deepcopy(self.genome_coding_strand)

        # ..and then replace with any relevant sections that code from the complementary strand
        self.genome_sequence[self.genome_on_noncoding_strand]=self.genome_noncoding_strand[self.genome_on_noncoding_strand]

        self._recreate_genes(self.genes_mutated,show_progress_bar=show_progress_bar)

    def _set_sequence_metadata(self,idx,sample_info):

        mask=self.genome_index==idx
        assert numpy.count_nonzero(mask)>0, idx

        altered_gene=self.genome_feature_name[mask][0]

        if altered_gene!="" and altered_gene not in self.genes_mutated:
            self.genes_mutated.append(altered_gene)

        if self.metadata_fields is not None:
            for field in self.metadata_fields:
                if field in sample_info.keys():
                    self.genome_sequence_metadata[field][mask]=sample_info[field]

    def _permute_sequence(self,idx,coverage,after=None,indel_length=0,indel_bases=(None,None)):

        # calculate a Boolean mask identifying where we are in the genome
        mask=self.genome_index==idx

        # use to assign the coverage
        self.coverage[mask]=coverage

        # permuate the sequence
        if after is not None:
            self.genome_coding_strand[mask]=after
            if after=='x':
                self.is_null[mask]=True
            elif after=='z':
                self.is_het[mask]=True
            elif after in ['a','c','t','g']:
                self.is_snp[mask]=True
            else:
                raise TypeError("passed base "+after+" not one of a,t,c,g,z,x")

        # .. and remember the indel length
        if indel_length!=0:
            self.indel_length[mask]=indel_length
            self.is_indel[mask]=True
            self.indel_ref[mask]=indel_bases[0]
            self.indel_alt[mask]=indel_bases[1]

    def save_pickle(self,filename=None,compression=True,compresslevel=1):

        assert compression in [True,False]
        assert filename is not None
        assert compresslevel in range(1,10), "compresslevel must be in range 1-9!"

        if compression:
            OUTPUT=gzip.open(filename+".gz",'wb',compresslevel=compresslevel)
        else:
            OUTPUT=open(filename,'wb')

        pickle.dump(self,OUTPUT)
        OUTPUT.close()

    def save_sequence(self,filename=None):

        '''
        Save the genome as a compressed NPZ file (compressed internally using gzip).

        This is purely done because loading an NPZ file back into memory is FAST (~200µs) so this could allow future analyses

        Args:
            filename (str): path of the output file without the file extension
        '''

        numpy.savez_compressed(filename,sequence=self.genome_sequence)

    def save_fasta(self,filename=None,compression=False,compresslevel=2,chars_per_line=70,nucleotides_uppercase=True):

        '''
        Save the genome as a FASTA file.

        Args:
            filename (str): path of the output file
            compression (bool): If True, save compressed using gzip. (bzip2 is too slow)
            compresslevel (0-9): the higher the number, the harder the algorithm tries to compress but it takes longer. Default is 2.
            # additional_metadata (str): will be added to the header of the FASTA file
            chars_per_line (int): the number of characters per line. Default=70. Must be either a positive integer or None (i.e. no CRs)
        '''

        # check the arguments are well formed
        assert compression in [True,False]
        assert nucleotides_uppercase in [True,False]
        assert compresslevel in range(1,10), "compresslevel must be in range 1-9!"
        if chars_per_line is not None:
            assert chars_per_line > 0, "number of characters per line in the FASTA file must be an integer!"

        # check the specified fileextension to see if the FASTA file needs compressing
        if compression:
            OUTPUT=gzip.open(filename+".gz",'wb',compresslevel=compresslevel)
        else:
            OUTPUT=open(filename,'w')

        # create the header line for the FASTA file using "|" as delimiters
        header=">"
        if hasattr(self,'id'):
            header+=self.id+"|"
        if hasattr(self,'organism'):
            header+=self.organism+"|"
        header+=self.name
        # if additional_metadata is not None:
        #     header+="|" + additional_metadata
        header+="\n"

        # create a string of the genome
        genome_string=''.join(self.genome_coding_strand)

        # insert carriage returns so it looks pretty in the file...
        output_string=self._insert_newlines(genome_string,every=chars_per_line)
        output_string+="\n"

        # set the case accordingly
        if nucleotides_uppercase:
            output_string=output_string.upper()
        else:
            output_string=output_string.lower()

        # write out the FASTA files
        if compression:
            OUTPUT.write(str.encode(header))
            OUTPUT.write(str.encode(output_string))
        else:
            OUTPUT.write(header)
            OUTPUT.write(output_string)

        OUTPUT.close()

    @staticmethod
    def _get_variant_for_genotype_in_vcf_record(
        genotype: Genotype, record: pysam.VariantRecord
    ) -> str:
        """Retrieves the variant a genotype maps to for a given record.

        Args:
            genotype: The genotype call for the sample
            record: A VCF record object.
        Returns:
            str: A hyphen if the call is null (ie ./.) or the alt variant if
            the call is alt. Returns an empty string if the call is ref or heterozygous.
        """

        #  find out what the reference bases are
        ref_bases=record.ref.lower()

        if genotype.is_reference():
            variant = ""
        elif genotype.is_heterozygous():
            variant = record.alleles[genotype.call1].lower(),record.alleles[genotype.call2].lower()
        elif genotype.is_alt():
            variant = record.alleles[genotype.call1].lower()
        elif genotype.is_null():
            variant = "x"*len(ref_bases)
        else:
            raise UnexpectedGenotypeError(
                """Got a genotype for which a Ref/Alt/Null call could not be
                    determined: {}.\nPlease raise this with the developers.""".format(
                    genotype.call()
                )
            )


        return ref_bases,int(record.pos),variant

    def _is_record_invalid(self, ignore_filter: bool, record: pysam.VariantRecord) -> bool:
        """
        Simple private method for parsing VCF record.

        Args:
            ignore_filter (bool):  ignore the FILTER column in the VCF file?
            record (pysam.VariantRecord): the record object from the VCF file to consider

        Returns:
            True/False
        """
        return ( not ignore_filter and "PASS" not in record.filter.keys() )

    def _minos_gt_in_wrong_position_fix(self,record, sample_idx):
        """
        Hacky private method to fix a minos mistake

        (A version of minos had GT in the second column instead of the first)
        """

        info = str(record).strip().split("\t")[9 + sample_idx]
        for field in info.split(":"):
            if "/" in field:
                return Genotype.from_string(field)


    @staticmethod
    def _load_fastafile(fasta_file):
        """
        Loads the fasta file whether uncompressed or compressed by gzip

        Args:
            fasta_file(path): the path to the fasta file

        Returns:
            header (str): the first line of the FASTA file which will hopefully contain some metadata
            nucleotide_sequence (str): the nucleotide sequence as a string
        """

        # check if it is compressed and load it accordingly
        if fasta_file.endswith(".gz"):
            INPUT = gzip.open(fasta_file,'rb')
            header=INPUT.readline().decode()
            nucleotide_sequence=INPUT.read().decode()
        else:
            INPUT = open(fasta_file,'r')
            header=INPUT.readline()
            nucleotide_sequence=INPUT.read()

        nucleotide_sequence=nucleotide_sequence.replace('\n','')

        return(header,nucleotide_sequence)

    @staticmethod
    def _insert_newlines(string: str, every=70):
        '''
        Simple private method for inserting a carriage return every N characters into a long string.

        Args:
            string (str): the string to insert carriage returns
            every (int): how many characters between each carriage return
        '''

        assert every>0, "every must be an integer greater than zero"

        assert len(string)>1, "string is too short!"

        return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

    def valid_genome_variant(self,genome_variant: str):
        '''
        Check whether the supplied genome_variant is valid and correctly formed.

        Args:
            genome_variant (str): in the format a100t (at genome index 100 mutate the adenosine to thymine )

        Returns:
            either fails an assertion or is True
        '''

        assert genome_variant is not None, "a genome_variant must be specified! e.g. FIXME"

        # since the smallest variant is '1a>c'
        assert len(genome_variant)>=4, "a genome variant must have at least 4 characters e.g. 3a>c"
        assert genome_variant[-2]=='>', "a genome variant must be of the form 3a>c where 3 is the 1-based position in the genome"

        before=genome_variant[-3]
        after=genome_variant[-1]

        assert before!=after, "before and after are identical hence this is not a variant!"

        try:
            position=int(genome_variant[:-3])
        except:
            raise TypeError("position "+genome_variant[:-3]+" is not an integer!")

        # check that the specified base is actually a base!
        assert before in ['c','t','g','a'], before+" is not a valid nucleotide!"

        # having checked that position is an integer, make a mask
        mask=self.genome_index==position

        # if the mask contains anything other than a single True, the position is not in the genome
        assert numpy.count_nonzero(mask)==1, "position specified not in the genome!"

        # confirm that the given base matches what is in the genome
        assert before==self.genome_sequence[mask][0], "base in genome is "+self.genome_sequence[mask][0]+" but specified base is "+before

        # check that the base to be mutated to is valid (z=het, ?=any other base according to the grammar)
        assert after in ['c','t','g','a','?','z'], after+" is not a valid nucleotide!"

        return True

    def convert_variant_to_mutation(self,gene_variant: str):
        """
        Function for converting a nucleotide gene variant to an amino acid mutation

        This is necessary since the TB genetic catalogues "speak" gene_mutations not gene_variants

        E.g. rpoB_c1349t -> rpoB_S450L

        Args:
            gene_variant (str): must be in the form e.g. rpoB_c1349t

        Returns:
            gene_mutation (str): e.g. rpoB_S450L
        """

        # break apart the gene variant
        cols=gene_variant.split("_")

        # check it is in the format we expect
        assert len(cols)==2, "a gene_variant can only contain two elements"

        # retrieve the gene name
        gene_name=cols[0]

        # check the gene is found in this genome
        assert self.contains_gene(gene_name), "gene not in the genome!"

        # make a copy of the gene so we can mutate it to get at the answer
        tmp_gene=copy.deepcopy(self.genes[gene_name])

        # split apart the variant
        ref_base=cols[1][0]
        assert ref_base in ['a','c','t','g'], 'reference base is not a, t, c or g!'

        alt_base=cols[1][-1]
        assert alt_base in ['a','c','t','g','x','z'], 'alt base is not a, t, c, g, z or x!'

        try:
            position=int(cols[1][1:-1])
        except:
            raise TypeError("the position "+cols[1][1:-1]+" is not an integer!")

        # now we can make the mask based on the nucleotide position
        mask=tmp_gene.positions==position

        # more paranoia and check the provided reference base matches what is in the gene
        assert tmp_gene.sequence[mask]==ref_base, 'provided reference base '+ref_base+' does not match '+tmp_gene.sequence[mask][0]+" which is what is in the genome"

        # find out the amino acid numbering
        amino_acid_position=tmp_gene.numbering[mask][0]

        # find out the reference amino acid
        ref=tmp_gene.amino_acid_sequence[tmp_gene.amino_acid_numbering==amino_acid_position][0]

        # now mutate the base
        tmp_gene.sequence[mask]=alt_base

        # retranslate the amino acid sequence
        tmp_gene._translate_sequence()

        # now find out the new amino acid
        alt=tmp_gene.amino_acid_sequence[tmp_gene.amino_acid_numbering==amino_acid_position][0]

        # and return the gene_mutation
        return(gene_name+"_"+ref+str(amino_acid_position)+alt)


    def valid_gene_mutation(self,mutation):
        '''
        Parse the mutation and return a collection of variables and Booleans.

        Args:
            mutation (str) e.g. katG_S315T, katG_c-15t, katG_200_ins_3

        Returns:
            gene_name (str): the name of the gene
            before (str): the amino acid or base in the reference genomes
            position (int or str): the numerical position of either the amino acid or base
            after (str): the mutated amino acid or base
            wildcard (bool): whether the mutation applies to all positions or not
            promoter (bool): whether the mutation applies to the promoter of a gene
            mutation_type (str): either PROMOTER, SNP or INDEL
        '''

        # first, parse the mutation
        cols=mutation.split("_")

        before=None
        after=None
        wildcard=False

        # check the mutation at least comprises the expected number of sections
        assert len(cols) in [2,3,4], "mutation "+mutation+" not in correct format!"

        # the gene/locus name should always be the first component
        gene_name=cols[0]

        assert self.contains_gene(gene_name), "gene not found in Genome! "+gene_name

        # find out if it is a GENE, LOCUS or RNA
        gene_type=self._gene_type[gene_name]

        if gene_type=="RNA":
            nucleotide_mutation=True
        else:
            nucleotide_mutation=False

        # determine if this is a CDS or PROMOTER SNP mutation
        if len(cols)==2:

            if '*' in cols[1]:

                # there can be no 'before' amino acid if there is a wildcard at position
                wildcard=True

                if cols[1][0]=='-':
                    nucleotide_mutation=True
                    position=str(cols[1][1:-1])
                else:
                    nucleotide_mutation=False
                    position=str(cols[1][0:-1])

                # all the positions should be a wildcard otherwise something is wrong..
                assert position=='*', mutation+' has a * but not formatted like a wildcard'


            else:

                wildcard=False

                before=cols[1][0]
                position=int(cols[1][1:-1])

                if position<0 or before.islower():
                    nucleotide_mutation=True

            # they all have the after amino acid in the same place
            after=cols[1][-1]

            if after.islower():
                nucleotide_mutation=True

            # if it is a promoter mutation
            if nucleotide_mutation:

                assert after in ['c','t','g','a','?','z'], after+" is not a nucleotide!"

                if wildcard:
                    return(True)

                else:
                    assert before in ['c','t','g','a'], before+" is not a nucleotide!"

                    return(self.genes[gene_name].valid_variant(cols[1]))

            # ..otherwise it is an amino acid SNP
            else:

                assert after in ['=','?',"!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid!"

                if not wildcard:

                    assert before in ["!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], before+" is not an amino acid!"

                    try:
                        position=int(cols[1][1:-1])
                    except:
                        raise TypeError("position "+cols[1]+" is not an integer")

                    mask=self.genes[gene_name].amino_acid_numbering==position

                    assert numpy.count_nonzero(mask)==1, "position "+str(position)+" is not in the genome"

                    assert before==self.genes[gene_name].amino_acid_sequence[mask][0], "specified amino acid is "+before+" but is "+self.genes[gene_name].amino_acid_sequence[mask][0]+" in the reference genome"

                return(True)

        # otherwise it must be an INDEL, which is always nucleotide based
        else:

            # deal with wildcards in the position
            if cols[1] in ["*","-*"]:

                assert cols[2] in ["ins","del","indel","fs"], "INDEL must be on the format katG_*_fs i.e. the third element must be ins or del, not "+cols[2]

                return(True)

            else:

                assert "*" not in cols[1], "mutation "+mutation+" contains a wildcard (*) but is badly formed"

                try:
                    position=int(cols[1])
                except:
                    raise TypeError("the position "+cols[1]+" is not an integer!")

                mask=self.genes[gene_name].positions==position

                assert numpy.count_nonzero(mask)==1, "specified position "+cols[1]+" not in the gene"

                # be defensive here also!
                assert cols[2] in ["ins","del","indel","fs"], "INDEL must be on the format rpoB_1300_ins_1 i.e. the third element must be ins or del, not "+cols[2]

                if len(cols)==4:
                    try:
                        number_nucleotides=int(cols[3])
                        assert number_nucleotides!=0, "an INDEL must be a non-zero number of bases"
                    except:
                        assert bool(re.match('^[catg]+$', cols[3])), cols[3]+" INDEL contains bases other than a,t,c,g"

                    # if cols[3].isnumeric():
                    #     assert int(cols[3])>0, "number of nucleotides inserted or deleted must be >0"
                    # else:
                    #     assert cols[2] in ["ins","indel"], cols[2]+" can only specify precise bases for an insertion!"


                return True
