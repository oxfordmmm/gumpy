import gzip, os, pickle, time, copy
from collections import defaultdict

import numpy, h5py

import pysam
from Bio import SeqIO
from tqdm import tqdm

import gumpy
from gumpy import Gene
from gumpy import Genotype


class Genome(object):

    def __init__(self,genbank_file=None,fasta_file=None,show_progress_bar=False,default_promoter_length=100):

        '''
        Instantiates a genome object by loading a VCF file and storing the whole genome as a numpy array

        Args:
            genbank_file (str): path to the GenBank file to build the reference genome
            fasta_file (str): path to the FASTA file to build the reference genome
        '''

        assert ((genbank_file is not None) or (fasta_file is not None)), "one of a GenBank file or a FASTA file must be specified!"

        self.id=""
        self.organism=""
        self.sample_name=""
        self.sample_metadata={}

        # load the specified GenBank file
        if genbank_file is not None:

            # create the genbank file and store in a BioPython object
            reference_genome=SeqIO.read(genbank_file,'genbank')

            self.sample_name="Reference Genome"

            # convert to a numpy array at the first opportunity since slicing BioPython is between 10 and 50,000 times slower!
            self.sequence=numpy.array([i.lower() for i in str(reference_genome.seq)])

            # create the complementary strand upfront
            self.complementary_bases = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'x':'x', 'z':'z'}
            self.complementary_sequence=numpy.array([self.complementary_bases[i] for i in self.sequence])

            # store the length of the genome
            self.length=len(self.sequence)

            # create an array of the genome indices
            self.index=numpy.arange(1,self.length+1)

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

            self.gene=numpy.zeros(self.length,dtype="<U10")
            self.is_gene=numpy.zeros(self.length,dtype=bool)
            self.promoter=numpy.zeros(self.length,dtype="<U10")
            self.is_promoter=numpy.zeros(self.length,dtype=bool)
            self.reverse=numpy.zeros(self.length,dtype=bool)
            self.type=numpy.zeros(self.length,dtype="<U5")
            self.gene_or_promoter=numpy.zeros(self.length,dtype="<U10")

            self.gene_type=defaultdict(str)
            self.gene_is_reverse=defaultdict(str)
            self.gene_index_start=defaultdict(str)
            self.gene_index_end=defaultdict(str)
            self.gene_codes_protein=defaultdict(str)
            self.overlaps_existing_gene=defaultdict(str)
            self.moo=defaultdict(str)

            self.coding_sequence=copy.deepcopy(self.sequence)
            self.coding_position=numpy.zeros(self.length,int)

            previous_gene_reversed=False

            # go through the GenBank file record-by-record
            for record in tqdm(reference_genome.features,disable=not(show_progress_bar)):

                if record.type in ['CDS','rRNA']:

                    gene_name=None
                    gene_type=None

                    if 'gene' in record.qualifiers.keys():

                        gene_name=record.qualifiers['gene'][0]

                        if record.type=='rRNA':
                            gene_type="RNA"
                        else:
                            gene_type="GENE"

                    elif 'locus_tag' in record.qualifiers.keys():

                        gene_name=record.qualifiers['locus_tag'][0]

                        if record.type=='rRNA':
                            gene_type="RNA"
                        else:
                            gene_type="LOCUS"

                    else:
                        continue

                    if gene_name is not None:

                        gene_start=int(record.location.start)
                        gene_end=int(record.location.end)
                        if gene_type in ["GENE","LOCUS"]:
                            gene_codes_protein=True
                        else:
                            gene_codes_protein=False

                        if record.strand==1:

                            self.gene_is_reverse[gene_name]=False

                            # This is a bit hacky; we are assuming that the genes are being read from the Genbank file
                            # in sequential order, hence for +ve strand genes the start of the next may overwrite the end of the previous,
                            # which is what we want. But the next gene after a -ve strand gene may overwrite the _start_ of that gene, which we don't want
                            # hence the logic to remember if the previous gene is on the reverse strand or not!

                            if previous_gene_reversed:
                                gene_mask=(self.index>gene_start) & (self.index<=gene_end) & (~self.is_gene)
                            else:
                                gene_mask=(self.index>gene_start) & (self.index<=gene_end)

                            previous_gene_reversed=False

                            promoter_mask=(self.index>gene_start-default_promoter_length) & (self.index<=gene_start) & (~self.is_gene) & (~self.is_promoter)

                            mask=gene_mask+promoter_mask

                            # be paranoid and set even though the default value is True
                            self.reverse[mask]=True

                            promoter_coding_position=self.index[promoter_mask]-gene_start-1

                            if gene_codes_protein:
                                gene_coding_position=numpy.floor_divide(self.index[gene_mask]-gene_start+2,3)
                            else:
                                gene_coding_position=self.index[gene_mask]-gene_start

                        elif record.strand==-1:

                            self.gene_is_reverse[gene_name]=True

                            if previous_gene_reversed:
                                gene_mask=(self.index>gene_start) & (self.index<=gene_end) & (~self.is_gene)
                            else:
                                gene_mask=(self.index>gene_start) & (self.index<=gene_end)

                            previous_gene_reversed=True

                            promoter_mask=(self.index>gene_end) & (self.index<=gene_end+default_promoter_length) & (~self.is_gene) & (~self.is_promoter)

                            mask=gene_mask+promoter_mask

                            # the default value is False, so only need to set those which are reversed
                            self.reverse[mask]=True

                            # replace the coding sequence with the complement
                            self.coding_sequence[mask]=self.complementary_sequence[mask]

                            promoter_coding_position=-1*(self.index[promoter_mask]-gene_end)

                            if gene_codes_protein:
                                gene_coding_position=-1*(numpy.floor_divide(self.index[gene_mask]-gene_end-1,3))
                            else:
                                gene_coding_position=-1*(self.index[gene_mask]-gene_end)

                        else:
                            raise TypeError("gene in GenBank file has strand that is not 1 or -1")


                        self.type[mask]=gene_type
                        self.gene[gene_mask]=gene_name
                        self.promoter[gene_mask]=''
                        self.promoter[promoter_mask]=gene_name
                        self.gene_or_promoter[mask]=gene_name

                        self.coding_position[gene_mask]=gene_coding_position
                        self.coding_position[promoter_mask]=promoter_coding_position

                        self.is_gene[gene_mask]=True
                        self.is_promoter[gene_mask]=False
                        self.is_promoter[promoter_mask]=True

                        self.gene_type[gene_name]=gene_type
                        self.gene_index_start[gene_name]=gene_start
                        self.gene_index_end[gene_name]=gene_end
                        self.gene_codes_protein=gene_codes_protein

            # store a list of all the gene names
            self.gene_names=numpy.unique(self.gene[self.gene!=""])
            self.is_gene_or_promoter=self.is_gene+self.is_promoter
            self.gene_or_promoter=numpy.core.defchararray.add(self.gene,self.promoter)

            for i in tqdm(self.gene_names):
                if self.gene_type[i] in ["GENE",'LOCUS']:
                    cds=True
                else:
                    cds=False
                mask=self.gene_or_promoter==i
                self.moo[i]=Gene(gene_name=i,sequence=self.coding_sequence[mask],index=self.index[mask],position=self.coding_position[mask],codes_protein=cds)


        # otherwise there must be a FASTA file so load that instead
        elif fasta_file is not None:

            header,nucleotide_sequence=self._load_fastafile(fasta_file)

            cols=header[1:].split("|")
            if len(cols)>1:
                self.id=cols[0]
                self.organism=cols[1]
                self.sample_name=cols[2]
            # if len(cols)>3:
            #     self.additional_metadata=cols[3]

            self.sequence=numpy.array(list(nucleotide_sequence))

            # store the length of the genome
            self.length=len(self.sequence)

            # create an array of the genome indices
            self.index=numpy.arange(1,self.length+1)

        # insist that bases are lower case
        self.sequence=numpy.char.lower(self.sequence)

        # store the sequence as integers 0,1,2,3 with which bases they refer to in bases_integer_lookup
        self.bases_integer_lookup, self.integers = numpy.unique(self.sequence, return_inverse=True)

    def __repr__(self):

        '''
        Overload the print function to write a summary of the genome.
        '''

        output=""
        if hasattr(self,'id'):
            output+=self.id+"\n"
        if hasattr(self,'organism'):
            output+=self.organism+"\n"
        if hasattr(self,'sample_name'):
            output+=self.sample_name+"\n"
        output+=str(self.length)+" bases\n"
        output+=''.join(i for i in self.sequence[0:3])
        output+="..."
        output+=''.join(i for i in self.sequence[-3:])

        return(output)

    def __sub__(self,other):

        """
        Overload the subtraction operator so it returns a tuple of the differences between the two genomes
        """

        assert self.length==other.length, "genomes must have the same length!"

        mask=self.sequence!=other.sequence

        return(self.index[mask])

    @staticmethod
    def _complement(nucleotides_array):

        complementary_bases = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'x':'x', 'z':'z'}

        complement=[complementary_bases[i] for i in nucleotides_array]

        return(complement)

    def contains_gene(self,gene_name):

        if gene_name in self.gene_names:
            return True
        else:
            return False

    def at_index(self,index):

        assert index > 0, "index must be a positive integer!"
        assert index <= self.length, "index must be less than the length of the genome!"

        mask=self.index==index

        putative_gene=self.gene[mask]

        if putative_gene!="":
            return (putative_gene[0],self.gene_type[putative_gene[0]])
        else:
            putative_promoter=self.promoter[mask]
            if putative_promoter!="":
                return (putative_promoter[0],"PROM")
            else:
                return None

    def snp_distance(self,other):
        return (numpy.count_nonzero(self.sequence!=other.sequence))

    def apply_vcf_file(self,vcf_file=None,ignore_filter=False, ignore_status=False,show_progress_bar=False,metadata_fields=None):
        """
        Load a VCF file and apply the variants to the whole genome sequence.

        Args:
            vcf_file (str): path to the VCF file to be loaded and applied
            ignore_filter (bool): whether to ignore the FILTER column in the VCF file (Clockwork hasn't always written it correctly)
            ignore_status (bool): ditto
            show_progress_bar (bool): whether to draw a nice tqdm progress bar (False by default)
        """

        # since we are now applying a VCF file, it makes sense to create these numpy arrays
        self.coverage=numpy.zeros(self.length,int)
        self.indel_length=numpy.zeros(self.length,int)
        self.indel_ref=numpy.zeros(self.length,dtype='<U50')
        self.indel_alt=numpy.zeros(self.length,dtype='<U50')

        # create a set of mutually exclusive Boolean arrays that tell you what the 'single sequence' result is
        self.is_ref=numpy.zeros(self.length,dtype=bool)
        self.is_null=numpy.zeros(self.length,dtype=bool)
        self.is_het=numpy.zeros(self.length,dtype=bool)
        self.is_snp=numpy.zeros(self.length,dtype=bool)
        self.is_indel=numpy.zeros(self.length,dtype=bool)

        # set up a dictionary for the metadata since the field names will vary between calling calling codes
        # for Clockwork these will be GT_CONF and GT_CONF_PERCENTILE
        self.metadata_fields=metadata_fields
        if self.metadata_fields is not None:
            self.sequence_metadata={}
            for field in self.metadata_fields:
                self.sequence_metadata[field]=numpy.zeros(self.length,float)

        # to deal with HET calls we need to setup some diploid arrays
        self.het_variations=numpy.zeros((self.length,2),str)
        self.het_coverage=numpy.zeros((self.length,2),int)
        self.het_indel_length=numpy.zeros((self.length,2),int)
        self.het_ref=numpy.zeros(self.length,dtype='<U50')
        self.het_alt=numpy.zeros((self.length,2),dtype='<U50')

        # split and remember the path, filename and stem of the VCF file
        (self.vcf_folder,self.vcf_file_name)=os.path.split(vcf_file)
        self.vcf_file_stem, file_extension = os.path.splitext(self.vcf_file_name)

        # it may have been compressed, in which case there will be TWO fileextensions to remove
        if self.vcf_file_stem[-4:]==".vcf":
            self.vcf_file_stem, file_extension = os.path.splitext(self.vcf_file_stem)

        # assume that the sample name is the filestem and remember
        self.sample_name=self.vcf_file_stem

        # open the supplied VCF file
        # note that this will read in bgzip compressed vcf files (from htslib) but not gzip compressed files
        # even though the file extension is the same
        vcf_reader = pysam.VariantFile(vcf_file.rstrip())

        # now iterate through the records found in the VCF file
        for record in tqdm(vcf_reader,disable=not(show_progress_bar)):

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



                # bypass for speed if this is a REF call
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
                                mask=self.index==idx

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
                            mask=self.index==index

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
        for idx in self.index[all_hets_mask]:

            # where are we?
            mask=self.index==idx

            # define the total coverage as the sum of the two HET calls
            coverage=numpy.sum(self.het_coverage[mask])

            # make the mutation, identifying this as a HET call
            self._permute_sequence(idx,coverage,after='z')

        # first recompute the complementary strand
        self.complementary_sequence=numpy.array([self.complementary_bases[i] for i in self.sequence])

        # reintialise the coding sequence
        self.coding_sequence=copy.deepcopy(self.sequence)

        # ..and then replace with any relevant sections that code from the complementary strand
        self.coding_sequence[self.reverse]=self.complementary_sequence[self.reverse]

    def _set_sequence_metadata(self,idx,sample_info):

        if self.metadata_fields is not None:
            mask=self.index==idx
            for field in self.metadata_fields:
                if field in sample_info.keys():
                    self.sequence_metadata[field][mask]=sample_info[field]

    def _permute_sequence(self,idx,coverage,after=None,indel_length=0,indel_bases=(None,None)):

        # calculate a Boolean mask identifying where we are in the genome
        mask=self.index==idx

        # use to assign the coverage
        self.coverage[mask]=coverage

        # permuate the sequence
        if after is not None:
            self.sequence[mask]=after
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

    def save_sequence(self,filename=None):

        '''
        Save the genome as a compressed NPZ file (compressed internally using gzip).

        This is purely done because loading an NPZ file back into memory is FAST (~200µs) so this could allow future analyses

        Args:
            filename (str): path of the output file without the file extension
        '''

        numpy.savez_compressed(filename,sequence=self.sequence)

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
        header+=self.sample_name
        # if additional_metadata is not None:
        #     header+="|" + additional_metadata
        header+="\n"

        # create a string of the genome
        genome_string=''.join(self.sequence)

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
        Simple private method for parsing VCF record
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
    def _insert_newlines(string, every=70):
        '''
        Simple private method for inserting a carriage return every N characters into a long string.

        Args:
            string (str): the string to insert carriage returns
            every (int): how many characters between each carriage return
        '''

        assert every>0, "every must be an integer greater than zero"

        assert len(string)>1, "string is too short!"

        return '\n'.join(string[i:i+every] for i in range(0, len(string), every))
