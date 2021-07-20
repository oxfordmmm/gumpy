import base64
import copy
import gzip
import json
import math
import multiprocessing
import pickle
import time
from collections import defaultdict

import numpy
from Bio import SeqIO
from tqdm import tqdm

from gumpy import Gene


class Genome(object):

    def __init__(self, *args, **kwargs):
        if len(args) != 1:
            if "reloading" not in kwargs.keys():
                #A genbank file has not been given so warn user
                print("No genbank file was given, setting up minimal detail for this Genome")
            #Setup the kwargs
            #UPDATE THIS AS REQUIRED
            allowed_kwargs = ['verbose', 'gene_subset', 'max_gene_name_length', 'nucleotide_sequence', 'name', 'id', 
                            'description', 'length', 'nucleotide_index', 'stacked_gene_name', 'stacked_is_cds', 'stacked_is_promoter', 
                            'stacked_nucleotide_number', 'stacked_is_reverse_complement', 'is_indel', 'indel_length', 'annotations', 
                            'genes_lookup', 'gene_rows', 'genes_mask', 'n_rows', 'stacked_nucleotide_index', 'stacked_nucleotide_sequence', 'genes']
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
            return
        else:
            #Set values for kwargs to defaults if not provided
            show_progress_bar = kwargs.get("show_progress_bar", False)
            gene_subset = kwargs.get("gene_subset")
            default_promoter_length = kwargs.get("default_promoter_length", 100)
            max_gene_name_length = kwargs.get("max_gene_name_length", 20)
            verbose = kwargs.get("verbose", False)
            #Set the args value to the genbank file
            genbank_file = args[0]
            

        assert isinstance(default_promoter_length,int) and default_promoter_length>=0, "the promoter length must be a positive integer!"

        self.verbose=verbose
        self.gene_subset=gene_subset
        self.max_gene_name_length=max_gene_name_length

        if self.verbose:
            timings=defaultdict(list)
            start_time=time.time()

        self.__parse_genbank_file(genbank_file,gene_subset)

        if verbose:
            timings['parse genbank'].append(time.time()-start_time)
            start_time=time.time()

        self.__setup_arrays()

        if verbose:
            timings['define arrays'].append(time.time()-start_time)
            start_time=time.time()

        if default_promoter_length>0:
            self.__assign_promoter_regions(default_promoter_length)

        if verbose:
            timings['promoter'].append(time.time()-start_time)
            start_time=time.time()

        self.__recreate_genes(show_progress_bar=True)

        if verbose:
            timings['create genes'].append(time.time()-start_time)
            for i in timings:
                print("%20s %6.3f s" % (i, numpy.sum(timings[i])))
        
        self.__convert_references()

        #Set default attributes for items populated when a VCF is applied
        self.indels = None
        self.changes = None
        self.original = None
        self.calls = None
    
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

        '''
        Overload the print function to write a summary of the genome.
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
            output+='all genes/loci have been included\n'
        elif len(self.gene_subset)<10:
            output+='the following '+str(len(self.gene_subset))+' genes have been included: '
            for i in self.gene_subset:
                output+=str(i)+', '
        else:
            output+=str(len(self.gene_subset))+' gene/loci have been included.'
        return(output)

    def __sub__(self,other):

        """
        Overload the subtraction operator so it returns a tuple of the indices where there differences between the two genomes
        Args:
            other (gumpy.Genome) : The other genome used in the subtraction
        Returns:
            numpy.array : Numpy array of the indicies where the two genomes differ in nucleotides
        """

        assert self.length==other.length, "genomes must have the same length!"

        mask=self.nucleotide_sequence!=other.nucleotide_sequence

        mask+=self.indel_length!=0

        return(self.nucleotide_index[mask])
    
    def __eq__(self, other):
        '''
        Overloading the equality operator so two Genome objects can be compared directly
        Checks for the equality based on fields, but does not check for filename equality
        Args:
            other (gumpy.Genome) : The other Genome object to compare to
        Returns:
            bool : Boolean showing equality of the objects
        '''
        check = True
        check = check and self.genes == other.genes
        # print("genes", self.genes == other.genes)
        check = check and self.name == other.name
        # print("Name", self.name == other.name)
        check = check and self.id == other.id
        # print("ID", self.id == other.id)
        check = check and self.description == other.description
        # print("Desc", self.description == other.description)
        check = check and numpy.all(self.nucleotide_sequence == other.nucleotide_sequence)
        # print("NS", numpy.all(self.nucleotide_sequence == other.nucleotide_sequence))
        check = check and numpy.all(self.nucleotide_index == other.nucleotide_index)
        # print("NI", numpy.all(self.nucleotide_index == other.nucleotide_index))
        check = check and self.genes_lookup == other.genes_lookup
        # print("gene lookup", self.genes_lookup == other.genes_lookup)
        check = check and self.length == other.length
        # print("len", self.length == other.length)
        check = check and numpy.all(self.stacked_gene_name.tolist() == other.stacked_gene_name.tolist())
        # print("SGN", numpy.all(self.stacked_gene_name.tolist() == other.stacked_gene_name.tolist()))
        # print(self.stacked_gene_name, self.stacked_gene_name.shape)
        # print(other.stacked_gene_name, other.stacked_gene_name.shape)
        # print()

        return check
    
    def __len__(self):
        '''
        Adding len functionality - len(genome) will now return length of the genome
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
        #Use of dict.get(obj) returns an object or None if obj does not exist in dict
        #bool(None) = False, bool(obj) = True
        return bool(self.genes_lookup.get(gene_name))


    def at_index(self,index):
        '''
        Returns the name of any genome features (genes, loci) at a specified genome index (1-based).

        Args:
            index : (int)

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


    def save_pickle(self,filename=None,compression=True,compresslevel=1):
        '''
        Save the genome as a pickle.

        The main use for this is to save time by having to recreate a reference Genome object each time, which
        for bacterial genomes can take 5-10 mins.

        Args:
            filename (str): path of the output file without the file extension
            compression (bool): whether to compress the pkl or not using gzip (default: True)
            compresslevel (int): the gzip compression level, lower=faster (default:1)
        '''

        assert isinstance(compression,bool)
        assert filename is not None
        assert compresslevel in range(1,10), "compresslevel must be in range 1-9!"

        if compression:
            OUTPUT=gzip.open(filename+".gz",'wb',compresslevel=compresslevel)
        else:
            OUTPUT=open(filename,'wb')

        pickle.dump(self,OUTPUT)
        OUTPUT.close()
    
    
    @staticmethod
    def load(filename):
        '''Experimental loading utilising JSON and base64 encoding of numpy arrays

        Args:
            filename (str): Path to the file to load

        Returns:
            gumpy.Genome: The resulting Genome object
        '''        

        #Define a helper function to load a numpy array
        def load_numpy(dtype, data, shape):
            d_type = numpy.dtype(dtype)
            d_array = numpy.frombuffer(base64.decodestring(bytes(data, 'utf-8')), d_type)
            d_array = d_array.reshape(shape)
            return d_array

        #Define a helper function to load a Gene object
        def load_gene(name, value):
            loaded_gene = {"reloading": True}
            for g_attr in value.keys():
                g_a_type = value[g_attr][0]
                g_a_val = value[g_attr][1]
                if (g_a_type == "<class 'bool'>" or g_a_type == "<class 'NoneType'>" or
                    g_a_type == "<class 'str'>" or g_a_type == "<class 'int'>" or g_a_type == "<class 'float'>"
                    or g_a_type == "dict" or g_a_type == "<class 'dict'>"):
                    loaded_gene[g_attr] = g_a_val
                elif g_a_type == "numpy.array":
                    val = load_numpy(*g_a_val)
                    loaded_gene[g_attr] = val
                else:
                    assert False, name+" has a weird value of type: "+g_a_type
            return Gene(**loaded_gene)

        inp = json.load(open(filename))

        to_load = {"reloading": True}
        for field in tqdm(inp.keys()):
            # print(field)
            f_type = inp[field][0]
            f_value = inp[field][1]
            if (f_type == "<class 'bool'>" or f_type == "<class 'NoneType'>" or
                f_type == "<class 'str'>" or f_type == "<class 'int'>" or f_type == "<class 'float'>"
                or f_type == "dict"):
                to_load[field] = f_value
            else:
                #These fields require some wrangling to re-load
                if f_type == "genes":
                    loaded_genes = {}
                    for gene_name in f_value.keys():
                        loaded_genes[gene_name] = load_gene(gene_name, f_value[gene_name])
                    to_load["genes"] = loaded_genes
                if f_type == "numpy.array":
                    to_load[field] = load_numpy(*f_value)
            # print()
        return Genome(**to_load)
    
    def save(self, filename):
        '''Experimental way to save the entire object (and Gene objects) based on
            json and base64 encoding rather than relying on pickles
            Based on numpy serialisation detailed here by daniel451:
                https://stackoverflow.com/questions/30698004/how-can-i-serialize-a-numpy-array-while-preserving-matrix-dimensions
            This is definitely space inefficient (~1.7GB for a TB genome) but faster than re-instanciation
            If this causes issues (possibly due to transferring saved files between machines), it may be possible to 
                change numpy serialization to convert to/from lists, although it is likely that this will be significantly more
                computationally expensive

        Args:
            filename (str): Filename for the output file
        '''
        def save_numpy(numpy_obj):
            '''Helper function to save a numpy object

            Args:
                numpy_obj (numpy.array): Numpy array to save
            Returns:
                list : List of (dtype, base64 encoded object, shape)
            '''            
            return [str(numpy_obj.dtype), base64.b64encode(numpy_obj).decode("utf-8"), numpy_obj.shape]
        def save_gene(gene):
            '''Helper function to save a Gene object

            Args:
                gene (gumpy.gene): Gene object to save
            Returns:
                dict : Dictionary mapping attribute names to serialized numpy arrays
            '''
            #Output dict
            output = {}
            #Get all attributes (variables)
            attributes = [(attr, getattr(gene, attr)) for attr in vars(gene)]
            for (name, attr) in attributes:
                #Find the attributes which require attention
                if type(attr) not in [bool, str, int, float, list, type(None)]:
                    #Check specifically for numpy arrays
                    if type(attr) == type(numpy.array([])):
                        val = ("numpy.array", save_numpy(attr))
                    #And dicts as the genes dict stores Gene objects
                    elif type(attr) == dict:
                        if Gene in [type(value) for value in attr.values()]:
                            val = {}
                            for key in attr.values():
                                #Convert each Gene object
                                val[key] = save_gene(attr[key])
                            val = ("genes", val)
                        else:
                            val = (str(type(attr)), attr)
                    else:
                        val = (str(type(attr)), attr)
                else:
                    val = (str(type(attr)), attr)
                output[name] = val
            return output
        #Output dict
        output = {}
        #Get all attributes (variables)
        attributes = [(attr, getattr(self, attr)) for attr in vars(self)]
        for (name, attr) in tqdm(attributes):
            #Find the attributes which require attention
            if type(attr) not in [bool, str, int, float, list, type(None)]:
                #Check specifically for numpy arrays
                if type(attr) == type(numpy.array([])):
                    val = ("numpy.array", save_numpy(attr))
                #And dicts as the genes dict stores Gene objects
                elif type(attr) == dict:
                    if Gene in [type(value) for value in attr.values()]:
                        val = {}
                        for key in attr.keys():
                            #Convert each Gene object
                            val[key] = save_gene(attr[key])
                        val = ("genes", val)
                    else:
                        val = ("dict", attr)
                else:
                    val = (str(type(attr)), attr)
            else:
                val = (str(type(attr)), attr)
            output[name] = val
        #Write the output to a json file
        json.dump(output, open(filename, "w"))


    def save_sequence(self,filename=None):

        '''
        Save the genome as a compressed NPZ file (compressed internally using gzip).

        This is purely done because loading an NPZ file back into memory is FAST (~200µs) so this could allow future analyses

        Args:
            filename (str): path of the output file without the file extension
        '''
        numpy.savez_compressed(filename,sequence=self.nucleotide_sequence)

    def save_fasta(self,filename,compression=False,compresslevel=2,chars_per_line=70,nucleotides_uppercase=True):

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
        assert isinstance(compression,bool)
        assert isinstance(nucleotides_uppercase,bool)
        assert compresslevel in range(1,10), "compresslevel must be in range 1-9!"
        assert chars_per_line > 0, "number of characters per line in the FASTA file must be a positive integer!"

        # check the specified fileextension to see if the FASTA file needs compressing
        if compression:
            OUTPUT=gzip.open(filename+".gz",'wb',compresslevel=compresslevel)
        else:
            OUTPUT=open(filename,'w')

        # create the header line for the FASTA file using "|" as delimiters
        header=">"
        if hasattr(self,'name'):
            header+=self.name+"|"
        if hasattr(self,'id') and isinstance(self.id,str) and len(self.id)>0:
            header+=self.id+"|"
        if hasattr(self,'description') and isinstance(self.description,str) and len(self.description)>0:
            header+=self.description+"|"
        header=header[:-1]
        header+="\n"

        # create a string of the genome
        genome_string=''.join(self.nucleotide_sequence)

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

    def __parse_genbank_file(self,genbank_file,gene_subset):
        '''
        Private function to parse a genbank file
        Args:
            genbank_file (str) : Filename of the genbank file
            gene_subset (list) : List of gene names used to extract just a subset of genes
        '''

        reference_genome=SeqIO.read(genbank_file,'genbank')

        # convert to a numpy array at the first opportunity since slicing BioPython is between 10 and 50,000 times slower!
        self.nucleotide_sequence=numpy.array([i.lower() for i in str(reference_genome.seq)])

        self.name=reference_genome.name
        self.id=reference_genome.id
        self.description=reference_genome.description

        # store the length of the genome
        self.length=len(self.nucleotide_sequence)

        # create an array of the genome indices
        self.nucleotide_index=numpy.arange(1,self.length+1,dtype="int")

        self.stacked_gene_name=numpy.zeros((1,self.length),dtype='<U'+str(int(self.max_gene_name_length)))
        self.stacked_is_cds=numpy.zeros((1,self.length),dtype=bool)
        self.stacked_is_promoter=numpy.zeros((1,self.length),dtype=bool)
        self.stacked_nucleotide_number=numpy.zeros((1,self.length),dtype='int')
        self.stacked_is_reverse_complement=numpy.zeros((1,self.length),dtype=bool)

        self.is_indel=numpy.zeros(self.length,dtype=bool)
        self.indel_length=numpy.zeros(self.length,int)

        assert len(reference_genome.annotations['accessions'])==1, 'only GenBank files with a single accessions currently allowed'

        self.annotations={}
        for i in reference_genome.annotations.keys():
            self.annotations[i]=reference_genome.annotations[i]

        self.genes_lookup={}

        # loop through the features listed in the GenBank File
        for record in tqdm(reference_genome.features):

            # only parse coding sequences and rRNA features
            if record.type not in ['CDS','rRNA']:
                continue

            gene_name=None
            type=None
            codes_protein=True

            # try and use the gene name if available, otherwise use the locus
            if 'gene' in record.qualifiers.keys():
                gene_name=record.qualifiers['gene'][0]
                type='GENE'

            elif 'locus_tag' in record.qualifiers.keys():
                gene_name=record.qualifiers['locus_tag'][0]
                type="LOCUS"

            # if this is ribosomal RNA, then record as such
            if record.type=='rRNA':
                type="RNA"
                codes_protein=False

            # determine if this is a reverse complement gene (only relevant to dsDNA genomes)
            rev_comp=True if record.strand==-1 else False

            if gene_name is None or (gene_subset is not None and gene_name not in gene_subset):
                continue

            # sigh, you can't assume that a gene_name is unique in a GenBank file
            gene_name+="_2" if gene_name in self.genes_lookup.keys() else ''

            # since we've defined the feature_name array to be 20 chars, check the gene_name will fit
            assert len(gene_name)<=self.max_gene_name_length, "Gene "+gene_name+" is too long at "+str(len(gene_name))+" chars; need to change numpy.zeros definiton U20 to match"

            # note that BioPython "helpfully" turns these from 1-based into 0-based coordinates, hence the +1
            # gene_end has also been incremented by 1 so that slicing naturally works
            gene_start=int(record.location.start)+1
            gene_end=int(record.location.end)+1

            # record feature metadata in a dict
            self.genes_lookup[gene_name]={  'reverse_complement':rev_comp,\
                                        'type':type,\
                                        'codes_protein':codes_protein,\
                                        'start':gene_start,\
                                        'end':gene_end }
    
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
            #Use of the dot product of masks allows determining if a gene will fit in an array
            if numpy.dot(row, mask) == False:
                #There is not a collision with this row
                #Add the row
                #Start/End have to be adjusted to account for 0 indexing of arrays and 1 indexing of genetics
                row[start-1:end-1] = True
                genes[i][start-1:end-1] = gene_name
                self.__handle_rev_comp(rev_comp, start, end, i)
                return genes_mask, genes, i
        #If this point is reached, there has been no rows without collisions, so add one
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
        Use of the dot product on boolean arrays returns a single boolean showing collisions in almost linear time (10^-5 secs for TB size)
            This can be used to determine which row the gene should be in
        '''
        #Default to all False values
        genes_mask = numpy.array([numpy.array([False for x in range(self.length)])]) #Boolean mask to show gene presence at indicies
        genes = numpy.array([numpy.array(['' for x in range(self.length)])], dtype="U"+str(self.max_gene_name_length)) #Gene names
        self.gene_rows = dict()#Dict to pull out row indicies for each gene in the stacked arrays
        for gene_name in tqdm(self.genes_lookup):
            #Get the start/end/rev_comp values
            start = self.genes_lookup[gene_name]["start"]
            end = self.genes_lookup[gene_name]["end"]
            rev_comp=self.genes_lookup[gene_name]['reverse_complement']

            #Determine the boolean mask for this gene
            if end<start:
                mask = numpy.logical_or((self.nucleotide_index>=start), (self.nucleotide_index<end))
                end += self.length
            else:
                mask=(self.nucleotide_index>=start) & (self.nucleotide_index<end)
            
            #Fit the gene into the stacked arrays
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


    def __assign_promoter_regions(self,default_promoter_length):
        '''
        Private function to assign promoter regions to genes
        Args:
            default_promoter_length (int) : The default length a promoter for a gene should be
        '''
        assert isinstance(default_promoter_length,int), 'default_promoter_length must be an integer!'

        assert default_promoter_length>0, 'default_promoter_length must be greater than zero'

        # labelling promoters is a difficult problem since
        #  (i)  it is arbitrary and
        #  (ii) we need to ensure that only unassigned bases can be labelled as promoters and each should only 'belong' to a single feature
        # the latter is especially difficult when you have two genes next to one another, one reverse complement, since their promoters can
        # 'fight' for space. It is this problem that means we have to grow each promoter out one base at a time
        
        #Populate a dictionary to store the starts/ends of genes as they grow with promoters
        start_end = {gene_name : {
                                "start": self.genes_lookup[gene_name]["start"],
                                "end": self.genes_lookup[gene_name]["end"]}
                    for gene_name in self.genes_lookup}
        for promoter in tqdm(range(1,default_promoter_length+1)):

            #Replacement `start_end` because dictionaries can't be changed during iteration
            new_start_end = dict()
            for gene_name in start_end:
                #Get the associated start/end
                start = start_end[gene_name]["start"]
                end = start_end[gene_name]["end"]
                rev_comp = self.genes_lookup[gene_name]["reverse_complement"]
                #Check if the region which the gene would grow into is empty
                if rev_comp:
                    #Indexing is weird so stacked_array[i][end-1] is the end of the gene
                    #   making stacked_array[i][end] the next item on the right
                    if end == len(self.nucleotide_sequence):
                        #If the end would be out of range, loop back around to position 0
                        end = 0
                    pos = end
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
                    #Move the start/end values appropriately
                    if rev_comp:
                        new_start_end[gene_name]["end"] = end + 1
                    else:
                        new_start_end[gene_name]["start"] = start - 1
            start_end = new_start_end

    @staticmethod
    def __insert_newlines(string: str, every=70):
        '''
        Simple private method for inserting a carriage return every N characters into a long string.

        Args:
            string (str): the string to insert carriage returns
            every (int): how many characters between each carriage return
        '''

        assert every>0, "every must be an integer greater than zero"

        assert len(string)>1, "string is too short!"

        return '\n'.join(string[i:i+every] for i in range(0, len(string), every))
    
    def _build_gene(self, gene, conn=None):
        '''
        Private function to build the gumpy.Gene object
        (Should be private with leading `__` but a cross-platform bug means this causes crashes)
        Args:
            gene (str) : The name of the gene
            conn (multiprocessing.connection=None) : Connection object from a multiprocessing.Pipe(), default to None
                                                        for cases where it is faster to single thread
        Returns:
            gumpy.Gene : The instanciated gene object. Returns None in cases where a Connection object is passed.
        '''
        # FIXME: this gene has ribosomal slippage in it...
        if gene=='ORF1ab':
            if conn:
                conn.send(None)
            return None

        #The mask for all stacked arrays (N-dim)
        stacked_mask=self.stacked_gene_name==gene
        #The mask for singular arrays (1-dim) by collapsing stacked mask to 1-dim
        mask=numpy.any(stacked_mask,axis=0)

        assert numpy.count_nonzero(mask)>0, "gene ("+gene+") not found in genome!"

        g = Gene(  name=gene,\
                                nucleotide_sequence=self.nucleotide_sequence[mask],\
                                index=self.nucleotide_index[mask],\
                                nucleotide_number=self.stacked_nucleotide_number[stacked_mask],\
                                is_cds=self.stacked_is_cds[stacked_mask],\
                                is_promoter=self.stacked_is_promoter[stacked_mask],\
                                is_indel=self.is_indel[mask],
                                indel_length=self.indel_length[mask],
                                codes_protein=self.genes_lookup[gene]['codes_protein'],\
                                reverse_complement=self.genes_lookup[gene]['reverse_complement'],\
                                feature_type=self.genes_lookup[gene]['type'])
        if conn:
            conn.send(g)
        else:
            return g

    def __recreate_genes(self,show_progress_bar=False):
        """
        Private method to re-instantiate the passed list Genes.

        This translates the nucleotide sequence into amino acids (if the gene codes protein) and is
        hence necessary after applying a vcf file, albeit only for those genes whose sequence has been altered.

        Args:
            show_progress_bar (True/False)  whether to show the (tqdm) progress bar

        Returns:
            None
        """
        list_of_genes = list(self.genes_lookup.keys())

        self.genes={}

        #Limit of how many processes can be open at once
        limit = 50

        #As there is some overhead for multithreading, there are cases where this is actually slower
        #So add a check to use a single thread if the number of required iterations is less than the limit
        if len(list_of_genes) <= limit:
            #Single threaded
            for gene in tqdm(list_of_genes, disable=not(show_progress_bar)):
                self.genes[gene] = self._build_gene(gene, conn=None)
            return

        #Multithreading
        for thread_index in tqdm(range(0, math.ceil(len(list_of_genes)/limit)),disable=not(show_progress_bar)):
            #Get the communication pipes required
            pipes = [multiprocessing.Pipe() for i in range(limit)]
            #Get some threads
            threads = [(multiprocessing.Process(target=self._build_gene, args=(gene, child)), gene, parent) 
                        for (gene, (parent, child)) in zip(list_of_genes[thread_index*limit:thread_index*limit+limit], pipes)]
            #Start some threads
            for i in range(limit):
                if i < len(threads):
                    threads[i][0].start()
            #Wait for the threads to finish
            for i in range(limit):
                if i < len(threads):
                    #Get data from the pipe
                    recieved = threads[i][2].recv()
                    if recieved is None:
                        continue
                    #Rejoin the main thread
                    threads[i][0].join()
                    #Set the appropriate value in the genes dict
                    self.genes[threads[i][1]] = recieved
                    #Close the threads and connections
                    threads[i][0].close()
    
    def apply_variant_file(self, vcf):
        '''Function to apply a variant file to the genome  - producing a replica genome with the specified changes

        Args:
            vcf (gumpy.VariantFile): The VariantFile object for the VCF

        Returns:
            gumpy.Genome: The resulting Genome object
        '''
        assert max(vcf.changes.keys()) <= self.length, "The VCF file details changes outside of this genome!"
        #Replicate this Genome object
        print("Copying the genome...")
        genome = copy.deepcopy(self)

        '''
        Using numpy's fancy array indexing may provide neat code, and provides some speed
            in some cases, the constant time access of a standard dictionary results in 
            faster code when the mask only contains a few True values. 
        For TB length arrays with a ~0.1% True mask, applying a mask takes ~10^-3s
            Applying a dictionary to the same array takes ~10^-4s
            ~0.1% True mask is a reasonable amount for this task as a VCF file is ~4000 entries
        '''
        genome.indels = dict()
        genome.changes = dict()
        genome.original = dict()
        #Change the nucleotide indicies
        print("Updating the genome...")
        for change in tqdm(vcf.changes.keys()):
            genome.original[change] = genome.nucleotide_sequence[change]
            if type(vcf.changes[change][0]) == str:
                #Only set values if the change is to a single nucleotide
                genome.nucleotide_sequence[change] = vcf.changes[change][0]
                #Record the changes in the format (old_base, new_base)
                genome.changes[change] = (genome.original[change], vcf.changes[change][0])
            else:
                #It was an indel, so add the indel call to the indels dict
                genome.indels[change] = vcf.changes[change][0]

        #Rebuild the genes with this new information
        print("Rebuilding the Gene objects with the updated genome...")
        genome.__recreate_genes(show_progress_bar=True)

        #Save all of the calls in the format {arr_index: (n_reads, call)}
        genome.calls = {index: vcf.changes[index][1] for index in vcf.changes.keys()}

        return genome
