import gzip
import math
import multiprocessing
import pickle
import time
from collections import defaultdict

import numpy
from Bio import SeqIO
from scipy import ndimage
from tqdm import tqdm

from gumpy import Gene


class Genome(object):

    def __init__(self,\
                genbank_file,\
                show_progress_bar=False,\
                gene_subset=None,\
                default_promoter_length=100,\
                max_gene_name_length=20,
                verbose=False):

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
        check = check and self.name == other.name
        check = check and self.id == other.id
        check = check and self.description == other.description
        check = check and numpy.all(self.nucleotide_sequence == other.nucleotide_sequence)
        check = check and numpy.all(self.nucleotide_index == other.nucleotide_index)
        check = check and self.genes_lookup == other.genes_lookup
        check = check and self.length == other.length
        check = check and numpy.all(self.stacked_gene_name.tolist() == other.stacked_gene_name.tolist())
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
                row[start-1:end-1] = True
                genes[i][start-1:end-1] = gene_name
                self.__handle_rev_comp(rev_comp, start, end, i)
                return genes_mask, genes
        #If this point is reached, there has been no rows without collisions, so add one
        genes_mask = numpy.vstack((genes_mask, mask))
        new_row = numpy.array([gene_name if m else '' for m in mask])
        genes = numpy.vstack((genes, new_row))
        i += 1
        self.__handle_rev_comp(rev_comp, start, end, i)
        return genes_mask, genes
    
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
            genes_mask, genes = self.__fit_gene(mask, genes, genes_mask, start, end, gene_name, rev_comp)

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

        for promoter in tqdm(range(1,default_promoter_length+1)):

            # awkward logic to cope with genetic numbering without zeros i.e. -2, -1, 1, 2
            if promoter==1:
                mask=self.stacked_nucleotide_number==1
            else:
                mask=self.stacked_nucleotide_number==-1*(promoter-1)

            # pick out the forward features first
            before=(mask) & (~self.stacked_is_reverse_complement)
            gaps=numpy.tile(~numpy.any(self.stacked_gene_name!='',axis=0),(self.n_rows,1))

            # use a scikit-image function to grow out
            after=ndimage.grey_dilation(before,footprint=[[1,0,0]],mode='wrap') & gaps
            before=ndimage.grey_dilation(after,footprint=[[0,0,1]],mode='wrap')

            if numpy.sum(after)>0:
                self.stacked_gene_name[after]=self.stacked_gene_name[before]
                self.stacked_nucleotide_number[after]=-1*promoter

            before=(mask) & (self.stacked_is_reverse_complement)
            gaps=numpy.tile(~numpy.any(self.stacked_gene_name!='',axis=0),(self.n_rows,1))
            after=ndimage.grey_dilation(before,footprint=[[0,0,1]],mode='wrap') & gaps
            before=ndimage.grey_dilation(after,footprint=[[1,0,0]],mode='wrap')

            if numpy.sum(after)>0:
                self.stacked_gene_name[after]=self.stacked_gene_name[before]
                self.stacked_nucleotide_number[after]=-1*promoter
                self.stacked_is_reverse_complement[after]=True

        mask=self.stacked_nucleotide_number<0
        self.stacked_is_promoter[mask]=True

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
    
    def __build_gene(self, gene, conn=None):
        '''
        Private function to build the gumpy.Gene object
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
        limit = 25

        #As there is some overhead for multithreading, there are cases where this is actually slower
        #So add a check to use a single thread if the number of required iterations is less than the limit
        if len(list_of_genes) <= limit:
            #Single threaded
            for gene in tqdm(list_of_genes, disable=not(show_progress_bar)):
                self.genes[gene] = self.__build_gene(gene, conn=None)
            return

        #Multithreading
        for thread_index in tqdm(range(0, math.ceil(len(list_of_genes)/limit)),disable=not(show_progress_bar)):
            #Get the communication pipes required
            pipes = [multiprocessing.Pipe() for i in range(limit)]
            #Get some threads
            threads = [(multiprocessing.Process(target=self.__build_gene, args=(gene, child)), gene, parent) 
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
        




'''
Older version used for a comparison - checking that the multithreaded method produces the same results
'''
class Genome2(object):

    def __init__(self,\
                genbank_file,\
                show_progress_bar=False,\
                gene_subset=None,\
                default_promoter_length=100,\
                max_gene_name_length=20,
                verbose=False):

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

    def __repr__(self):

        '''
        Overload the print function to write a summary of the genome.
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
        """

        assert self.length==other.length, "genomes must have the same length!"

        mask=self.nucleotide_sequence!=other.nucleotide_sequence

        mask+=self.indel_length!=0

        return(self.nucleotide_index[mask])

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
        # list_of_genes=list(self.genes_lookup.keys())
        list_of_genes = self.genes_lookup.keys()

        if gene_name in list_of_genes:
            return(True)
        else:
            return(False)


    def at_index(self,index):
        '''
        Returns the name of any genome features (genes, loci) at a specified genome index (1-based).

        Args:
            index (int)

        Returns:
            list of gene_names or locus_tags at that index in the genome

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
            compresslevel (1-10): the gzip compression level, lower=faster (default:1)
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

        empty_row=numpy.zeros((1,array.shape[1]),dtype=array.dtype)

        return(numpy.vstack((array,empty_row)))

    def __parse_genbank_file(self,genbank_file,gene_subset):

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
    def make_arrays(self):
        self.stacked_gene_name=self.__add_empty_row(self.stacked_gene_name)
        self.stacked_is_cds=self.__add_empty_row(self.stacked_is_cds)
        self.stacked_is_reverse_complement=self.__add_empty_row(self.stacked_is_reverse_complement)
        self.stacked_is_promoter=self.__add_empty_row(self.stacked_is_promoter)
        self.stacked_nucleotide_number=self.__add_empty_row(self.stacked_nucleotide_number)
    def check(self, row, mask):
        return ~(numpy.all(self.stacked_gene_name[(row,mask)]==''))
    def mask(self, gene_start, gene_end):
        return numpy.logical_or((self.nucleotide_index>=gene_start), (self.nucleotide_index<gene_end))
    def _rev_comp(self, row, mask, gene_start, gene_end, rev_comp):
        if rev_comp:
            self.stacked_nucleotide_number[(row,mask)]=numpy.mod(-1*(self.nucleotide_index[mask]-gene_end),self.length)
            self.stacked_is_reverse_complement[(row,mask)]=True
        else:
            self.stacked_nucleotide_number[(row,mask)]=numpy.mod(1+self.nucleotide_index[mask]-gene_start,self.length)
    def tile_check(self):
        self.stacked_nucleotide_index=numpy.tile(self.nucleotide_index,(self.n_rows,1))

        self.stacked_nucleotide_sequence=numpy.tile(self.nucleotide_sequence,(self.n_rows,1))
    def __setup_arrays(self):

        for gene_name in tqdm(self.genes_lookup):

            gene_start=self.genes_lookup[gene_name]['start']
            gene_end=self.genes_lookup[gene_name]['end']
            rev_comp=self.genes_lookup[gene_name]['reverse_complement']

            # deal with features that overlap the 'start' of the genome
            if gene_end<gene_start:

                mask=self.mask()
                gene_end+=self.length

            else:

                mask=(self.nucleotide_index>=gene_start) & (self.nucleotide_index<gene_end)

            # start in the top row of the arrays to see if the feature will "fit"
            row=0

            # below is True if it does fit, otherwise it will try the next row
            while self.check(row, mask):

                row+=1

                n_rows=self.stacked_gene_name.shape[0]-1

                # if we need a new row to accommodate overlapping genes, add one to all the essential numpy arrays
                if row>n_rows:
                    self.make_arrays()

            self.stacked_gene_name[(row,mask)]=gene_name

            self._rev_comp(row, mask, gene_start, gene_end, rev_comp)

        # do as many assignments outside the loop, i.e. in one go to improve performance
        self.stacked_is_cds=self.stacked_gene_name!=''

        self.n_rows=self.stacked_gene_name.shape[0]

        self.tile_check()


    def __assign_promoter_regions(self,default_promoter_length):

        assert isinstance(default_promoter_length,int), 'default_promoter_length must be an integer!'

        assert default_promoter_length>0, 'default_promoter_length must be greater than zero'

        # labelling promoters is a difficult problem since
        #  (i)  it is arbitrary and
        #  (ii) we need to ensure that only unassigned bases can be labelled as promoters and each should only 'belong' to a single feature
        # the latter is especially difficult when you have two genes next to one another, one reverse complement, since their promoters can
        # 'fight' for space. It is this problem that means we have to grow each promoter out one base at a time

        for promoter in tqdm(range(1,default_promoter_length+1)):

            # awkward logic to cope with genetic numbering without zeros i.e. -2, -1, 1, 2
            if promoter==1:
                mask=self.stacked_nucleotide_number==1
            else:
                mask=self.stacked_nucleotide_number==-1*(promoter-1)

            # pick out the forward features first
            before=(mask) & (~self.stacked_is_reverse_complement)
            gaps=numpy.tile(~numpy.any(self.stacked_gene_name!='',axis=0),(self.n_rows,1))

            # use a scikit-image function to grow out
            after=ndimage.grey_dilation(before,footprint=[[1,0,0]],mode='wrap') & gaps
            before=ndimage.grey_dilation(after,footprint=[[0,0,1]],mode='wrap')

            if numpy.sum(after)>0:
                self.stacked_gene_name[after]=self.stacked_gene_name[before]
                self.stacked_nucleotide_number[after]=-1*promoter

            before=(mask) & (self.stacked_is_reverse_complement)
            gaps=numpy.tile(~numpy.any(self.stacked_gene_name!='',axis=0),(self.n_rows,1))
            after=ndimage.grey_dilation(before,footprint=[[0,0,1]],mode='wrap') & gaps
            before=ndimage.grey_dilation(after,footprint=[[1,0,0]],mode='wrap')

            if numpy.sum(after)>0:
                self.stacked_gene_name[after]=self.stacked_gene_name[before]
                self.stacked_nucleotide_number[after]=-1*promoter
                self.stacked_is_reverse_complement[after]=True

        mask=self.stacked_nucleotide_number<0
        self.stacked_is_promoter[mask]=True

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
    
    def build_gene(self, gene):
        # FIXME: this gene has ribosomal slippage in it...
        if gene=='ORF1ab':
            return None

        stacked_mask=self.stacked_gene_name==gene

        mask=numpy.any(stacked_mask,axis=0)
        # print("Got mask: ", gene)

        assert numpy.count_nonzero(mask)>0, "gene not found in genome!"
        # print("Passed check: ", gene)

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
        # print("Built gene: ", gene)
        self.genes[gene] = g
        # print("Sent gene: ", gene)

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

        list_of_genes=list(self.genes_lookup.keys())
        # print(list_of_genes)

        self.genes={}
        # pool = multiprocessing.Pool(5)
        # out = pool.map(self.build_gene, list_of_genes)
        # pool.close()
        # threads = [multiprocessing.Process(target=self.build_gene, args=(gene,)) for gene in list_of_genes]
        # limit = 20

        # # pass ALL the gene names to create all the Gene objects for the first time
        # for thread_index in tqdm(range(0, len(list_of_genes), limit),disable=not(show_progress_bar)):
        #     threads = [(multiprocessing.Process(target=self.build_gene, args=(gene, child)), gene, parent) for (gene, (parent, child)) in zip(list_of_genes[thread_index:thread_index+limit], [multiprocessing.Pipe() for i in range(limit)])]
        #     # print(thread_index+limit)
        #     #Start some threads
        #     for i in range(limit):
        #         if thread_index+i < limit:
        #             threads[thread_index+i][0].start()
        #     #Wait for the threads to finish
        #     for i in range(limit):
        #         if thread_index + 1 < limit:
        #             threads[thread_index+i][0].join()
        #             recieved = threads[thread_index+i][2].recv()
        #             print(recieved)
        #             self.genes[threads[thread_index+i][1]] = recieved
        #             threads[thread_index+i][0].close()
            # thread.join()
        for gene in tqdm(list_of_genes, disable=not(show_progress_bar)):
            self.build_gene(gene)
