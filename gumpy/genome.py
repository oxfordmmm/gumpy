'''
Genome object
'''
import base64
import copy
import gzip
import json
import pathlib
import time
from collections import defaultdict

import numpy
from Bio import SeqIO
from tqdm import tqdm

from gumpy import Gene, GenomeDifference, VCFFile
from gumpy.variantfile import VCFRecord


class Genome(object):

    """Genome object"""

    def __init__(self, *args, **kwargs):
        '''Constructor for the Genome object.

        Args:
            genbank_file (str) : The path to the genbank file.
            show_progress_bar (bool, optional) : Boolean as whether to show a progress bar when building Gene objects. Defaults to False.
            gene_subset (list, optional) : List of gene names used to extract just a subset of genes. Defaults to None
            max_promoter_length (int, optional) : Size of the default maximum number of upstream bases to consider the promoter of a gene. Defaults to 100
            max_gene_name_length (int, optional) : Length of the longest gene name. Defaults to 20
            verbose (bool, optional) : Give verbose statements? Defaults to False
            is_reference (bool, optional) : Is this a reference genome? i.e. mutations can be derived with respect to it? Defaults to False
        '''
        if len(args) != 1:
            if "reloading" not in kwargs.keys():
                #A genbank file has not been given so warn user
                print("No genbank file was given, setting up minimal detail for this Genome")
            #Setup the kwargs
            #UPDATE THIS AS REQUIRED
            allowed_kwargs = ['show_progress_bar','verbose', 'gene_subset', 'max_gene_name_length', 'is_reference']
            seen = set()
            for key in kwargs.keys():
                '''
                Use of a whitelist of kwargs allowed to stop overriding of functions from loading malicious file
                If a malicious user overwrites any of these variable names with a function, the function **should** never be called as none
                  of these are called within this code. I'm sure there will be a way to break this though, such as custom objects within a dict which should
                  have objects - and then overwriting functions to allow arbitrary execution.
                  A possible way to circumvent this kind of behaviour would be to cryptographically sign and verify each object loaded
                  when saving/loading respectively. However, this may make the saves less portable than would be desired.
                However, all of this depends on the use case. If it is not possible for a user to upload/download the saved files, this should not be a risk.
                    This would mean that this kind of caching would only be used internally
                '''
                if key in allowed_kwargs:
                    setattr(self, key, kwargs[key])
                    seen.add(key)
                for key in set(allowed_kwargs).difference(seen):
                    #Default values to None if the kwarg has not been passed
                    setattr(self, key, None)

            return
        else:
            #Set values for kwargs to defaults if not provided
            self.show_progress_bar = kwargs.get("show_progress_bar", False)
            self.gene_subset = kwargs.get("gene_subset")
            self.max_promoter_length = kwargs.get("max_promoter_length", 100)
            self.max_gene_name_length = kwargs.get("max_gene_name_length", 20)
            self.verbose = kwargs.get("verbose", False)
            self.is_reference = kwargs.get("is_reference", False)
            #Set the args value to the genbank file
            genbank_file = pathlib.Path(args[0])

        assert genbank_file.is_file(), 'GenBank file does not exist!'
        assert isinstance(self.max_promoter_length,int) and self.max_promoter_length>=0, "the promoter length must be zero or a positive integer!"

        assert isinstance(self.verbose,bool)
        assert isinstance(self.is_reference,bool)
        assert isinstance(self.show_progress_bar,bool)
        assert isinstance(self.max_gene_name_length,int) and self.max_gene_name_length>0
        if self.gene_subset is not None:
            # first check it is a list
            assert isinstance(self.gene_subset,list)
            # then check all elements in the list are strings
            assert all(isinstance(i,str) for i in self.gene_subset)

        if self.verbose:
            timings=defaultdict(list)
            start_time=time.time()

        self.__parse_genbank_file(pathlib.Path(genbank_file))

        if self.verbose:
            timings['parse genbank'].append(time.time()-start_time)
            start_time=time.time()

        self.__setup_arrays()

        if self.verbose:
            timings['define arrays'].append(time.time()-start_time)
            start_time=time.time()

        if self.max_promoter_length>0:
            self.__assign_promoter_regions()

        if self.verbose:
            timings['promoter'].append(time.time()-start_time)
            start_time=time.time()

        self.__convert_references()

    def __convert_references(self):
        '''Convert BIOPython Reference objects to normal dictionaries. They do not
            appear to have any greater application than storing structured data, so
            removing the object wrappers appears to be a clean way to combat the object's
            issues with serialization.
        '''
        for (i, reference) in enumerate(self.annotations["references"]):
            new_ref = {}
            for key in vars(reference):
                #This key contains unhelpfully structured data
                if key == "location":
                    new_loc = []
                    for item in vars(reference)[key]:
                        loc = {}
                        for item_key in vars(item):
                            if item_key == "_start" or item_key == "_end":
                                loc[item_key] = getattr(item, item_key).position
                            elif getattr(item, item_key) is not None:
                                loc[item_key] = getattr(item, item_key)
                        new_loc.append(loc)
                    new_ref[key] = new_loc
                else:
                    new_ref[key] = vars(reference)[key]
            self.annotations["references"][i] = new_ref

    def __repr__(self):
        '''Overload the print function to write a summary of the genome.

        Returns:
            str : String including attributes for the genome
        '''

        output=""
        if hasattr(self,'name'):
            output+=self.name+"\n"
        if hasattr(self,'id'):
            output+=self.id+"\n"
        if hasattr(self,'description'):
            output+=self.description+"\n"
        output+=str(self.length)+" bases\n"
        output+=''.join(i for i in self.nucleotide_sequence[0:6])
        output+="..."
        output+=''.join(i for i in self.nucleotide_sequence[-6:])+'\n'
        if self.gene_subset is None:
            output+='metadata for all genes/loci have been included\n'
        elif len(self.gene_subset)<10:
            output+='the following '+str(len(self.gene_subset))+' genes have been included: '
            for i in self.gene_subset:
                output+=str(i)+', '
        else:
            output+=str(len(self.gene_subset))+' gene/loci have been included.'
        return(output)

    def __sub__(self, other):
        '''Generate a GenomeDifference object for a in-depth difference of the two Genomes

        Args:
            other (gumpy.Genome) : The other genome used in the subtraction

        Returns:
            GenomeDifference: object containing numpy arrays of the differences (variants)
        '''

        assert isinstance(other,Genome), 'RHS must be a gumpy.Genome object'

        return(GenomeDifference(self,other))

    def __eq__(self, other):
        '''Overloading the equality operator so two Genome objects can be compared directly
        Checks for the equality based on fields, but does not check for filename equality
        Args:
            other (gumpy.Genome) : The other Genome object to compare to
        Returns:
            bool : Boolean showing equality of the objects
        '''
        assert isinstance(other,Genome)

        check = True
        check = check and self.genes == other.genes
        check = check and self.name == other.name
        check = check and self.id == other.id
        check = check and self.description == other.description
        check = check and numpy.all(self.nucleotide_sequence == other.nucleotide_sequence)
        check = check and numpy.all(self.nucleotide_index == other.nucleotide_index)
        check = check and self.length == other.length
        check = check and numpy.all(self.stacked_gene_name.tolist() == other.stacked_gene_name.tolist())

        return check

    def __len__(self):
        '''Adding len functionality - len(genome) returns the length of the genome

        Returns:
            int : Length of the genome
        '''
        return self.length

    def contains_gene(self,gene_name):
        '''
        Simply checks to see if the specified gene exists in the Genome object.

        Args:
            gene_name (str) : Name of the gene e.g. katG

        Returns:
            bool : Boolean showing if the genome contains a gene with that name
        '''
        assert type(gene_name) == str, "Gene name must be string. Gene name provided was of type: "+str(type(gene_name))
        #Use of dict.get(obj) returns an object or None if obj does not exist in dict
        #bool(None) = False, bool(obj) = True
        return bool(self.genes.get(gene_name))


    def at_index(self,index):
        '''
        Returns the name of any genome features (genes, loci) at a specified genome index (1-based).

        Args:
            index (int): Genome index to check for genes at.

        Returns:
            list : list of gene_names or locus_tags at that index in the genome

        '''
        assert isinstance(index,int), "index must be an integer!"
        assert index > 0, "index must be a positive integer!"
        assert index <= self.length, "index must be less than the length of the genome!"

        mask=self.stacked_nucleotide_index==index

        foo=self.stacked_gene_name[mask]

        putative_genes=list(foo[foo!=''])

        if not putative_genes:
            return(None)
        else:
            return(putative_genes)

    def save(self, filename, compression_level=None):
        '''Experimental way to save the entire object (and Gene objects) based on
            json and base64 encoding rather than relying on pickles (due to security implications of pickle)
            Based on numpy serialisation detailed here by daniel451:
                https://stackoverflow.com/questions/30698004/how-can-i-serialize-a-numpy-array-while-preserving-matrix-dimensions
            This is definitely space inefficient (~1.7GB for a TB genome) but faster than re-instanciation. Using compression, this can be reduced to <100MB
            If this causes issues (possibly due to transferring saved files between machines), it may be possible to
                change numpy serialization to convert to/from lists, although it is likely that this will be significantly more
                computationally expensive

        Args:
            filename (str): Path to the file to save in
            compression_level (int, optional): Level of compression to use (1-9), when None is given, no compression is used. Defaults to None.
        '''
        output = self.__save(self, None)
        #Write the output to a json file
        if compression_level is not None and compression_level in range(1,10):
            json.dump(output, gzip.open(filename, "wt", compresslevel=compression_level))
        else:
            json.dump(output, open(filename, "w"), indent=2)

    def __save(self, obj, output, name=None):
        '''Helper function to recursively convert and save each object

        Args:
            obj (object): Any object
            output (list/dict): Aggregator for outputs
            name (str, optional): Name of the attribute to be stored as in the dict. Defaults to None.

        Returns:
            dict/list/tuple: Either aggregated output or a single output depending on if aggregation was required
        '''
        if type(obj) in [bool, int, str, float, complex, bytes, bytearray] or obj is None:
            #Fundamental data types which need no conversions
            to_return = obj
        elif type(obj) == type(numpy.array([])):
            #Convert numpy arrays to 3 item lists
            if obj.flags["C_CONTIGUOUS"] == False:
                #Some arrays are not contiguous, so make them contiguous as base64 requires it
                obj = numpy.ascontiguousarray(obj, obj.dtype)
            to_return = [str(obj.dtype), base64.b64encode(obj).decode("utf-8"), obj.shape]
        elif type(obj) == list:
            #Convert items in a list
            to_return = [self.__save(x, []) for x in obj]
        elif type(obj) == tuple:
            #Convert items in a tuple
            to_return = tuple([self.__save(x, []) for x in obj])
        elif type(obj) == dict:
            #Keys should be hashable so ignore them but convert the values
            to_return = {key: self.__save(obj[key], {}) for key in obj.keys()}
        elif type(obj) in [Gene, Genome, VCFFile, VCFRecord]:
            #Convert items to dicts
            attributes = [(attr, getattr(obj, attr)) for attr in vars(obj)]
            to_return = {}
            for (a_name, attr) in attributes:
                to_return[a_name] = self.__save(attr, {})
        elif str(type(obj)) == "<class 'numpy.str_'>":
            #Not sure where these come from, but some numpy.str_ items exist so treat them as strings
            to_return = str(obj)
        else:
            #Other types are not allowed.
            assert False, "Object of weird type: "+str(obj)+str(type(obj))

        if name is not None:
            #This attribute has a name so add it to the output and return it
            if type(output) == dict:
                output[name] = (str(type(obj)), to_return)
            elif type(output) == list:
                output.append(str(type(obj)), to_return)
            return output
        else:
            #Unnamed, so just return it
            return (str(type(obj)), to_return)

    @staticmethod
    def load(filename):
        '''Load the object using base64 encoding and JSON

        Args:
            filename (str): Path to the saved file
        '''
        def _load(type_, obj, output, name=None):
            '''Helper function to recursively load an object

            Args:
                type_ (str): String of the type which the object should be
                obj (object): Any input object
                output (list/dict): Either a list or a dictionary used as an accumulator
            '''
            if type_ in [str(type(t)) for t in [bool(), int(), str(), float(), complex(), bytes(), bytearray(), None]]:
                #Fundamental data types which need no conversions
                to_return = obj
            elif type_ == str(type(numpy.array([]))):
                #Convert numpy arrays to 3 item lists
                d_type = numpy.dtype(obj[0])
                d_array = numpy.frombuffer(base64.decodebytes(bytes(obj[1], 'utf-8')), d_type)
                to_return = d_array.reshape(obj[2]).copy()
            elif type_ == str(type(list())):
                #Convert items in a list
                to_return = [_load(t, o, []) for (t, o) in obj]
            elif type_ == str(type(tuple())):
                #Convert items in a tuple
                to_return = tuple([_load(t, o, []) for (t, o) in obj])
            elif type_ == str(type(dict())):
                #Keys should be hashable so ignore them but convert the values
                to_return = {key: _load(*obj[key], {}) for key in obj.keys()}
            elif type_ == str(Gene):
                to_return = Gene(**{key: _load(*obj[key], {}) for key in obj.keys()}, reloading=True)
            elif type_ == str(Genome):
                to_return = Genome(**{key: _load(*obj[key], {}) for key in obj.keys()}, reloading=True)
            elif type_ == str(VCFFile):
                to_return = VCFFile(**{key: _load(*obj[key], {}) for key in obj.keys()}, reloading=True)
            elif type_ == str(VCFRecord):
                to_return = VCFRecord(**{key: _load(*obj[key], {}) for key in obj.keys()}, reloading=True)
            elif type_ == "<class 'numpy.str_'>":
                to_return = numpy.str_(obj)
            else:
                #Other types are not allowed.
                assert False, "Object of weird type: "+str(obj)+str(type(obj))

            if name is not None:
                #This attribute has a name so add it to the output and return it
                if type(output) == dict:
                    output[name] = to_return
                elif type(output) == list:
                    output.append(to_return)
                return output
            else:
                #Unnamed, so just return it
                return to_return
        try:
            inp = json.load(gzip.open(filename, "rt"))
        except gzip.BadGzipFile:
            #Not a gzipped file so try with normal open
            inp = json.load(open(filename))
        return _load(*inp, {})



    def save_sequence(self,filename=None):

        '''
        Save the genome as a compressed NPZ file (compressed internally using gzip).

        This is purely done because loading an NPZ file back into memory is FAST (~200µs) so this could allow future analyses

        Args:
            filename (str): path of the output file without the file extension
        '''
        numpy.savez_compressed(filename,sequence=self.nucleotide_sequence)

    def __build_genome_variable_length_string(self,indices):
        genome_string=''
        # work backwards as easier to deal with insertions/deletions when you've already gone past them
        for i in indices[::-1]:
            mask=self.nucleotide_index==i
            base=self.nucleotide_sequence[mask][0]
            genome_string=base+genome_string
            if self.is_indel[mask]:
                indel_length=self.indel_length[mask][0]
                if indel_length>0:
                    genome_string=self.indel_nucleotides[mask][0]+genome_string
                elif indel_length<0:
                    genome_string=genome_string[abs(indel_length):]
        return(genome_string)

    def build_genome_string(self,fixed_length=False,nucleotide_index_range=None):
        '''
        Generate a string of the nucleotides in the genome (positive strand if DNA).

        Args:
            fixed_length (bool): if True, then do not add insertions and deletions. Default False.
            nucleotide_index_range (tuple, ints): the 1-based positions of the sequence to return with start<=index<end.

        Returns:
            (str): the genome as a string.
        '''
        # create a string of the genome
        if fixed_length:
            if nucleotide_index_range is not None:
                assert isinstance(nucleotide_index_range, tuple)
                start,end=nucleotide_index_range
                genome_string=''.join(self.nucleotide_sequence[start-1:end-1])
            else:
                genome_string=''.join(self.nucleotide_sequence)
        else:
            if nucleotide_index_range is not None:
                start,end=nucleotide_index_range
                genome_string=self.__build_genome_variable_length_string(self.nucleotide_index[start-1:end-1])
            else:
                genome_string=self.__build_genome_variable_length_string(self.nucleotide_index)

        return(genome_string)


    def save_fasta(self,filename,fixed_length=False,nucleotide_index_range=None,compression=False,compresslevel=2,chars_per_line=70,nucleotides_uppercase=True,description=None, overwrite_existing=True):

        '''
        Save the genome as a FASTA file.

        Args:
            filename (str): path of the output file
            fixed_length (bool): If True, ignore indels and only output a genome the same length as the reference but with SNPs. This is useful for phylogeny analyses and relatedness. If false, a genome including indels is produced. Default is false.
            nucleotide_index_range (tuple, optional): A tuple of (start,end) genome indices
            compression (bool): If True, save compressed using gzip. (bzip2 is too slow)
            compresslevel (0-9): the higher the number, the harder the algorithm tries to compress but it takes longer. Default is 2.
            chars_per_line (int): the number of characters per line. Default=70. Must be either a positive integer or None (i.e. no CRs)
            nucleotide_uppercase (bool): If True, provide the nucleotides in UPPER CASE. Default is True.
            description (str, optional): what to write on the header line of the FASTA file. If not provided, then a description will be automatically generated from the GenBank file metadata.
            overwrite_existing (bool): If False, then the code will refuse to overwrite a FASTA file already on disc. Default is True. 
        '''

        # check the arguments are well formed
        if not overwrite_existing:
            assert not pathlib.Path(filename).is_file(), 'filename already exists! '+filename
        assert isinstance(compression,bool)
        assert isinstance(fixed_length,bool)
        assert isinstance(nucleotides_uppercase,bool)
        assert isinstance(chars_per_line,int)
        if nucleotide_index_range is not None:
            assert isinstance(nucleotide_index_range,tuple)
            assert isinstance(nucleotide_index_range[0],int)
            assert isinstance(nucleotide_index_range[1],int)
            assert nucleotide_index_range[0]>=1, 'genomes are 1-based so the first base must be >=1'
            assert nucleotide_index_range[1]<self.length, 'longer than the genome!'
        assert compresslevel in range(1,10), "compresslevel must be in range 1-9!"
        assert chars_per_line > 0, "number of characters per line in the FASTA file must be a positive integer!"
        if description is not None:
            assert isinstance(description,str)

        # check the specified fileextension to see if the FASTA file needs compressing
        if compression:
            OUTPUT=gzip.open(filename+".gz",'wb',compresslevel=compresslevel)
        else:
            OUTPUT=open(filename,'w')

        # create the header line for the FASTA file using "|" as delimiters
        header=">"
        if description is None:
            if hasattr(self,'name'):
                header+=self.name+"|"
            if hasattr(self,'id') and isinstance(self.id,str) and len(self.id)>0:
                header+=self.id+"|"
            if hasattr(self,'description') and isinstance(self.description,str) and len(self.description)>0:
                header+=self.description+"|"
            header=header[:-1]
        else:
            header+=description
        header+="\n"

        genome_string=self.build_genome_string(fixed_length,nucleotide_index_range)

        # insert carriage returns so it looks pretty in the file...
        output_string=self.__insert_newlines(genome_string,every=chars_per_line)
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


    def __add_empty_row(self,array):
        '''
        Private function to add an empty row of the correct type to a numpy array
        Args:
            array (numpy.array) : Array to add an empty row to
        Returns:
            (numpy.array): The same array with an empty row of the same length and dtype appended
        '''

        empty_row=numpy.zeros((1,array.shape[1]),dtype=array.dtype)

        return(numpy.vstack((array,empty_row)))

    def __parse_genbank_file(self,genbank_file):
        '''
        Private function to parse a genbank file
        Args:
            genbank_file (Path) : pathlib.Path object of the genbank file
        '''

        if genbank_file.suffix == '.gz':
            file_handle  = gzip.open(genbank_file,'rt')
        elif genbank_file.suffix == '.gbk':
            file_handle  = open(genbank_file,'rt')

        reference_genome=SeqIO.read(file_handle,'genbank')

        # convert to a numpy array at the first opportunity since slicing BioPython is between 10 and 50,000 times slower!
        self.nucleotide_sequence=numpy.array([i.lower() for i in str(reference_genome.seq)])

        self.name=reference_genome.name
        self.id=reference_genome.id
        self.description=reference_genome.description

        # store the length of the genome
        self.length=len(self.nucleotide_sequence)

        assert self.length>0, 'genome length zero!'

        # create an array of the genome indices
        self.nucleotide_index=numpy.arange(1,self.length+1,dtype="int")

        self.stacked_gene_name=numpy.zeros((1,self.length),dtype='<U'+str(int(self.max_gene_name_length)))
        self.stacked_is_cds=numpy.zeros((1,self.length),dtype=bool)
        self.stacked_is_promoter=numpy.zeros((1,self.length),dtype=bool)
        self.stacked_nucleotide_number=numpy.zeros((1,self.length),dtype='int')
        self.stacked_is_reverse_complement=numpy.zeros((1,self.length),dtype=bool)

        self.is_indel=numpy.zeros(self.length,dtype=bool)
        self.indel_length=numpy.zeros(self.length,int)
        self.indel_nucleotides=numpy.empty(self.length,dtype=object)

        assert len(reference_genome.annotations['accessions'])==1, 'only GenBank files with a single accessions currently allowed'

        self.annotations={}
        for i in reference_genome.annotations.keys():
            self.annotations[i]=reference_genome.annotations[i]

        self.genes={}

        # loop through the features listed in the GenBank File
        if self.verbose:
            print("Iterating through features in GenBank file...")

        for record in tqdm(reference_genome.features, disable=(not self.show_progress_bar)):

            # only parse coding sequences and rRNA features
            if record.type not in ['CDS','rRNA']:
                continue

            gene_name=None
            type_=None
            codes_protein=True

            # try and use the gene name if available, otherwise use the locus
            if 'gene' in record.qualifiers.keys():
                gene_name=record.qualifiers['gene'][0]
                type_='GENE'

            elif 'locus_tag' in record.qualifiers.keys():
                gene_name=record.qualifiers['locus_tag'][0]
                type_="LOCUS"

            if gene_name is None or (self.gene_subset is not None and gene_name not in self.gene_subset):
                continue

            # if this is ribosomal RNA, then record as such
            if record.type=='rRNA':
                type_="RNA"
                codes_protein=False

            # determine if this is a reverse complement gene (only relevant to dsDNA genomes)
            rev_comp=True if record.strand==-1 else False

            # sigh, you can't assume that a gene_name is unique in a GenBank file
            # this only allows for duplicates though.
            # duplicates of duplicates will be foo_2_2
            gene_name+="_2" if gene_name in self.genes.keys() else ''

            # check the gene_name will fit in the max gene name length
            assert len(gene_name)<=self.max_gene_name_length, "Gene "+gene_name+" is too long at "+str(len(gene_name))+" chars; need to specify max_gene_name_length"

            # note that BioPython "helpfully" turns these from 1-based into 0-based coordinates, hence the +1
            # gene_end has also been incremented by 1 so that slicing naturally works
            gene_start=int(record.location.start)+1
            gene_end=int(record.location.end)+1

            #Check for ribosomal shift
            #This happens when a start position < end position
            positions = [(loc.start.position, loc.end.position) for loc in record.location.parts]
            shifts = []
            #Check for -1 PFS
            if len(positions) > 1 and positions[0][1] > positions[1][0]:
                shifts.append(positions[1][0] - gene_start + 1)

            # record feature metadata in a dict
            self.genes[gene_name]={  'reverse_complement':rev_comp,\
                                        'type':type_,\
                                        'codes_protein':codes_protein,\
                                        'start':gene_start,\
                                        'end':gene_end,
                                        'ribosomal_shifts': shifts }

        # now we can check that all the genes in the gene_subset exist in the GenBank file!
        if self.gene_subset is not None:
            for i in self.gene_subset:
                assert self.contains_gene(i), 'Gene '+i+' not found in the Genbank file!'

    def __handle_rev_comp(self, rev_comp, start, end, i):
        '''
        Private function to handle the rev-comp changes required
        Args:
            rev_comp (bool) : Boolean to show if rev-comp is required
            start (int) : Start index of the gene
            end (int) : End index of the gene
            i (int) : The index of the stacked row
        '''
        #Check if the arrays have the correct number of rows and add as required
        while len(self.stacked_nucleotide_number) <= i:
            self.stacked_is_cds=self.__add_empty_row(self.stacked_is_cds)
            self.stacked_is_reverse_complement=self.__add_empty_row(self.stacked_is_reverse_complement)
            self.stacked_is_promoter=self.__add_empty_row(self.stacked_is_promoter)
            self.stacked_nucleotide_number=self.__add_empty_row(self.stacked_nucleotide_number)
        #Update the required items for rev_comp
        #Use of array slicing here introduces speedups as not all of the array should be considered
        #(Previous method applied a mask which requires full array consideration rather than direct access)
        if rev_comp:
            self.stacked_nucleotide_number[i][start-1:end-1] = numpy.mod(-1*(self.nucleotide_index[start-1:end-1]-end),self.length)
            self.stacked_is_reverse_complement[i][start-1:end-1] = True
        else:
            self.stacked_nucleotide_number[i][start-1:end-1] = numpy.mod(1+self.nucleotide_index[start-1:end-1]-start,self.length)

    def __fit_gene(self, mask, genes, genes_mask, start, end, gene_name, rev_comp):
        '''
        Private function to fit a gene into the genes based on the dot product of the masks
            numpy.dot([bool], [bool])-> bool showing if there are collisions of True values within args
            This takes 10^-5 seconds which is a significant improvement on use of numpy.all() iteration of 10^-2 seconds for TB length genome
        Args:
            mask (numpy.array) : Boolean array showing positions where the gene lies
            genes (numpy.array) : 2D numpy array of the format used for all stacked values
            genes_mask (numpy.array) : The corresponding boolean mask arrays for the `genes` arg
            start (int) : Start index of the gene
            end (int) : End index of the gene
            gene_name (str) : Name of the gene
            rev_comp (bool) : Boolean to show whether the gene required a reverse complement
        Returns:
            (numpy_array) : Updated genes_mask array
            (numpy_array) : Updated genes array
        '''

        for (i, row) in enumerate(genes_mask):

            # use of the dot product of masks allows determining if a gene will fit in an array
            if numpy.dot(row, mask) == False:

                # there is no collision with this row so add the row
                # start/end have to be adjusted to account for 0 indexing of arrays and 1 indexing of genetics
                row[start-1:end-1] = True
                genes[i][start-1:end-1] = gene_name
                self.__handle_rev_comp(rev_comp, start, end, i)
                return genes_mask, genes, i

        # if this point is reached, there has been no rows without collisions, so add one
        genes_mask = numpy.vstack((genes_mask, mask))
        new_row = numpy.array([gene_name if m else '' for m in mask])
        genes = numpy.vstack((genes, new_row))
        i += 1
        self.__handle_rev_comp(rev_comp, start, end, i)
        return genes_mask, genes, i

    def __find_overlaps(self):
        '''
        Private function to find the sections of the genome in which there are overlapping genes
        This should be more efficient than the older version as it avoids consistent genome iteration
        Use of the dot product on boolean arrays returns a single boolean showing collisions in almost constant time (10^-5 secs for TB size)
            This can be used to determine which row the gene should be in
        '''

        # Boolean mask to show gene presence at indicies (default to all False values)
        genes_mask = numpy.array([numpy.array([False for x in range(self.length)])])

        # gene names
        genes = numpy.array([numpy.array(['' for x in range(self.length)])], dtype="U"+str(self.max_gene_name_length))

        # dict to pull out row indicies for each gene in the stacked arrays
        self.gene_rows = dict()

        if self.verbose:
            print("Finding overlaps...")

        for gene_name in tqdm(self.genes, disable=(not self.show_progress_bar)):

            # get the start/end/rev_comp values
            start = self.genes[gene_name]["start"]
            end = self.genes[gene_name]["end"]
            rev_comp=self.genes[gene_name]['reverse_complement']

            # determine the boolean mask for this gene
            if end<start:
                mask = numpy.logical_or((self.nucleotide_index>=start), (self.nucleotide_index<end))
                end += self.length
            else:
                mask=(self.nucleotide_index>=start) & (self.nucleotide_index<end)

            # fit the gene into the stacked arrays
            genes_mask, genes, row = self.__fit_gene(mask, genes, genes_mask, start, end, gene_name, rev_comp)
            self.gene_rows[gene_name] = row

        #Singular array to determine if there are genes in places within the genome
        self.genes_mask = numpy.any(genes_mask, axis=0)

        #Set instance variable for the gene names
        self.stacked_gene_name = genes

    def __setup_arrays(self):
        '''
        Private function to initalise all of the required arrays, fitting the gene names into the
            correct places within stacked arrays
        '''
        self.__find_overlaps()

        # do as many assignments outside the loop, i.e. in one go to improve performance
        self.stacked_is_cds=self.stacked_gene_name!=''

        self.n_rows=self.stacked_gene_name.shape[0]

        self.stacked_nucleotide_index=numpy.tile(self.nucleotide_index,(self.n_rows,1))

        self.stacked_nucleotide_sequence=numpy.tile(self.nucleotide_sequence,(self.n_rows,1))

    def __assign_promoter_regions(self):
        '''
        Private function to assign promoter regions to genes
        '''
        assert isinstance(self.max_promoter_length,int), 'max_promoter_length must be an integer!'

        assert self.max_promoter_length>0, 'max_promoter_length must be greater than zero'

        # labelling promoters is a difficult problem since
        #  (i)  it is arbitrary and
        #  (ii) we need to ensure that only unassigned bases can be labelled as promoters and each should only 'belong' to a single feature
        # the latter is especially difficult when you have two genes next to one another, one reverse complement, since their promoters can
        # 'fight' for space. It is this problem that means we have to grow each promoter out one base at a time

        # populate a dictionary to store the starts/ends of genes as they grow with promoters
        start_end = {gene_name : {
                                "start": self.genes[gene_name]["start"],
                                "end": self.genes[gene_name]["end"]}
                    for gene_name in self.genes}
        if self.verbose:
            print("Assigning promoters...")

        for promoter in tqdm(range(1,self.max_promoter_length+1), disable=(not self.show_progress_bar)):

            #Replacement `start_end` because dictionaries can't be changed during iteration
            new_start_end = dict()
            for gene_name in start_end:
                #Get the associated start/end
                start = start_end[gene_name]["start"]
                end = start_end[gene_name]["end"]
                rev_comp = self.genes[gene_name]["reverse_complement"]
                #Check if the region which the gene would grow into is empty
                if rev_comp:
                    #Indexing is weird so stacked_array[i][end-2] is the end of the gene
                    #   making stacked_array[i][end-1] the next item on the right
                    if end == len(self.nucleotide_sequence):
                        #If the end would be out of range, loop back around to position 0
                        end = 0
                    pos = end -1
                else:
                    #Similar indexing issue except indexing starts on start-1
                    #   so start-2 is the next item on the left
                    pos = start -2
                if self.genes_mask[pos] == True:
                    #There is a gene here already so skip it
                    continue
                else:
                    #This position is free so set the appropriate values
                    new_start_end[gene_name] = start_end[gene_name] #Retain gene for future expansions
                    #Get the row index for stacked arrays
                    row = self.gene_rows[gene_name]

                    #Set appropriate values
                    self.stacked_gene_name[row][pos] = gene_name
                    self.stacked_nucleotide_number[row][pos] = -1*promoter
                    self.stacked_is_reverse_complement[row][pos] = rev_comp
                    self.stacked_is_promoter[row][pos] = True
                    self.genes_mask[pos] = True

                    # move the start/end values appropriately
                    if rev_comp:
                        new_start_end[gene_name]["end"] = end + 1
                    else:
                        new_start_end[gene_name]["start"] = start - 1
            start_end = new_start_end

    def __insert_newlines(self, string: str, every=70):
        '''
        Simple private method for inserting a carriage return every N characters into a long string.

        Args:
            string (str): the string to insert carriage returns
            every (int): how many characters between each carriage return
        '''

        assert every>0, "every must be an integer greater than zero"

        assert len(string)>1, "string is too short!"

        return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

    def build_gene(self, gene):
        '''
        Public function to build the gumpy.Gene object

        Args:
            gene (str) : The name of the gene

        Returns:
            gumpy.Gene : The instanciated gene object.
        '''

        # The mask for all stacked arrays (N-dim)
        stacked_mask=self.stacked_gene_name==gene

        # The mask for singular arrays (1-dim) by collapsing stacked mask to 1-dim
        mask=numpy.any(stacked_mask,axis=0)

        # check that the genome does contain the gene name
        assert numpy.count_nonzero(mask)>0, "gene ("+gene+") not found in genome!"

        # instantiate a Gene object
        g = Gene(name=gene,
                    nucleotide_sequence=self.nucleotide_sequence[mask],
                    nucleotide_index=self.nucleotide_index[mask],
                    nucleotide_number=self.stacked_nucleotide_number[stacked_mask],
                    is_cds=self.stacked_is_cds[stacked_mask],
                    is_promoter=self.stacked_is_promoter[stacked_mask],
                    is_indel=self.is_indel[mask],
                    indel_length=self.indel_length[mask],
                    indel_nucleotides=self.indel_nucleotides[mask],
                    codes_protein=self.genes[gene]['codes_protein'],
                    reverse_complement=self.genes[gene]['reverse_complement'],
                    feature_type=self.genes[gene]['type'],
                    ribosomal_shifts=self.genes[gene]['ribosomal_shifts'] )

        return g

    def __add__(self, vcf):
        '''Function to apply a VCF file to the genome  - producing a replica genome with the specified changes

        Args:
            vcf (gumpy.VCFFile): The VCFFile object for the VCF

        Returns:
            gumpy.Genome: The resulting Genome object
        '''

        assert isinstance(vcf,VCFFile), 'RHS must be a gumpy.VCFFile object!'

        assert isinstance(vcf.calls,dict), 'something wrong with the gumpy.VCFFile object!'

        indices=[i[0] for i in vcf.calls.keys()]

        assert max(indices) <= self.length, "The VCF file details changes outside of this genome!"

        if self.verbose:
            print("Copying the genome...")

        # Replicate this Genome object
        genome = copy.deepcopy(self)

        '''
        Using numpy's fancy array indexing may provide neat code, and provides some speed
            in some cases, the constant time access of a standard dictionary results in
            faster code when the mask only contains a few True values.
        For TB length arrays with a ~0.1% True mask, applying a mask takes ~10^-3s
            Applying a dictionary to the same array takes ~10^-4s
            ~0.1% True mask is a reasonable amount for this task as a VCF file is ~4000 entries
        '''

        if self.verbose:
            print("Updating the genome...")

        # use the calls dict to change the nucleotide indicies in the copy of the genome
        for (idx,type_) in tqdm(vcf.calls.keys(), disable=(not self.show_progress_bar)):

            # deal with changes at a single nucleotide site
            if type_ in ['snp','null','het']:

                # only set values if the idx is to a single nucleotide
                genome.nucleotide_sequence[idx-1] = vcf.calls[(idx,type_)]['call']

            # deal with insertions and deletions
            elif type_ in ['indel']:

                genome.is_indel[idx-1] = True
                genome.indel_nucleotides[idx-1] = vcf.calls[(idx,type_)]['call'][1]

                if vcf.calls[(idx,type_)]['call'][0]=='ins':
                    genome.indel_length[idx-1] = len(vcf.calls[(idx,type_)]['call'][1])
                else:
                    genome.indel_length[idx-1] = -1*len(vcf.calls[(idx,type_)]['call'][1])

            else:
                raise Exception('variant type not recognised!', vcf.calls[(idx,type_)])

        # the genome has been altered so not a reference genome
        genome.is_reference = False

        return genome
