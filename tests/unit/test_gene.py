import pytest, copy, numpy
from gumpy.genome import Genome

test_nucleotide_string="aaatttcccggg"
test_promoter_string="atcg"

reference=Genome(genbank_file="config/NC_004148.2.gbk",name="HMPV")

test_gene=reference.genes["F"]

# def test_Gene_instantiation():
#
#     assert test_gene.total_number_nucleotides==1720
#     assert test_gene.gene_name=="F"
#     assert test_gene.on_noncoding_strand==False
#     assert test_gene.numbering[0]==-100
#     assert test_gene.numbering[-1]==540
#     assert test_gene.numbering[-2]==540
#     assert test_gene.index[test_gene.positions==1]==3067
#     assert test_gene.index[test_gene.positions==2]==3068
#     assert test_gene.index[test_gene.positions==-1]==3066
#     assert test_gene.sequence[0]=='a'
#     assert test_gene.sequence[-1]=='g'
#     assert test_gene.sequence[-2]=='a'
#
#     F_gene="MSWKVVIIFSLLITPQHGLKESYLEESCSTITEGYLSVLRTGWYTNVFTLEVGDVENLTCSDGPSLIKTELDLTKSALRELKTVSADQLAREEQIENPRQSRFVLGAIALGVATAAAVTAGVAIAKTIRLESEVTAIKNALKTTNEAVSTLGNGVRVLATAVRELKDFVSKNLTRAINKNKCDIDDLKMAVSFSQFNRRFLNVVRQFSDNAGITPAISLDLMTDAELARAVSNMPTSAGQIKLMLENRAMVRRKGFGILIGVYGSSVIYMVQLPIFGVIDTPCWIVKAAPSCSGKKGNYACLLREDQGWYCQNAGSTVYYPNEKDCETRGDHVFCDTAAGINVAEQSKECNINISTTNYPCKVSTGRHPISMVALSPLGALVACYKGVSCSIGSNRVGIIKQLNKGCSYITNQDADTVTIDNTVYQLSKVEGEQHVIKGRPVSSSFDPIKFPEDQFNVALDQVFENIENSQALVDQSNRILSSAEKGNTGFIIVIILIAVLGSSMILVSIFIIIKKTKKPTGAPPELSGVTNNGFIPHS!"
#
#     assert F_gene==''.join(test_gene.amino_acid_sequence)
#
#
# # def test_Gene___repr__():
# #
# #     expected_output="F gene\n1720 nucleotides, codes for protein\nactac...taaaa\n-100 -99 -98 -97 -96 ...-5 -4 -3 -2 -1\nMSWKV...IPHS!\n1 2 3 4 5 ...536 537 538 539 540 \n"
# #
# #     assert test_gene.__repr__()==expected_output
#
# new_gene=copy.deepcopy(test_gene)
# new_gene.sequence[new_gene.positions==1]='t'
# new_gene._translate_sequence()
#
# def test_Gene_list_mutations_wrt_1():
#
#     assert new_gene.list_mutations_wrt(test_gene)==['M1L']
#
# def test_Gene_table_mutations_wrt_1():
#
#     MUTATIONS=new_gene.table_mutations_wrt(test_gene)
#
#     assert list(MUTATIONS['MUTATION'])==['M1L']
#     assert list(MUTATIONS['REF'])==['atg']
#     assert list(MUTATIONS['ALT'])==['ttg']
#     assert list(MUTATIONS['POSITION'])==[1]
#     assert list(MUTATIONS['AMINO_ACID_NUMBER'])==[1]
#     assert numpy.isnan(MUTATIONS['NUCLEOTIDE_NUMBER'])
#     assert list(MUTATIONS['IS_SNP'])==[True]
#     assert list(MUTATIONS['IS_INDEL'])==[False]
#     assert list(MUTATIONS['IN_CDS'])==[True]
#     assert list(MUTATIONS['IN_PROMOTER'])==[False]
#     assert numpy.isnan(MUTATIONS['INDEL_LENGTH'])
#     assert list(MUTATIONS['MUTATION_TYPE'])==['SNP']
#
# def test_Gene___sub__1():
#
#     assert new_gene-test_gene==[1]
#
# new_gene2=copy.deepcopy(test_gene)
# new_gene2.sequence[new_gene2.positions==1]='c'
# new_gene2.sequence[new_gene2.positions==4]='c'
# new_gene2._translate_sequence()
#
# def test_Gene_list_mutations_wrt_2():
#
#     assert new_gene2.list_mutations_wrt(test_gene)==['M1L','S2P']
#
# def test_Gene___sub__2():
#
#     assert list(new_gene2-test_gene)==[1,4]
#
# def test_Gene_table_mutations_wrt_2():
#
#     MUTATIONS=new_gene2.table_mutations_wrt(test_gene)
#
#     assert list(MUTATIONS['GENE'])==['F','F']
#     assert list(MUTATIONS['MUTATION'])==['M1L','S2P']
#     assert list(MUTATIONS['REF'])==['atg','tct']
#     assert list(MUTATIONS['ALT'])==['ctg','cct']
#     assert list(MUTATIONS['POSITION'])==[1,2]
#     assert list(MUTATIONS['AMINO_ACID_NUMBER'])==[1,2]
#     assert list(MUTATIONS['NUCLEOTIDE_NUMBER'])==[None,None]
#     assert list(MUTATIONS['IS_SNP'])==[True,True]
#     assert list(MUTATIONS['IS_INDEL'])==[False,False]
#     assert list(MUTATIONS['IN_CDS'])==[True,True]
#     assert list(MUTATIONS['IN_PROMOTER'])==[False,False]
#     assert list(MUTATIONS['INDEL_LENGTH'])==[None,None]
#     assert list(MUTATIONS['MUTATION_TYPE'])==['SNP','SNP']
#
# # mutate it away, and then back again
# new_gene3=copy.deepcopy(test_gene)
# new_gene3.sequence[new_gene3.positions==1]='t'
# new_gene3._translate_sequence()
# new_gene3.sequence[new_gene3.positions==1]='a'
# new_gene3._translate_sequence()
#
# def test_Gene_list_mutations_wrt_3():
#
#     assert new_gene3.list_mutations_wrt(test_gene) is None
#
# new_gene4=copy.deepcopy(test_gene)
# new_gene4.is_indel[new_gene4.positions==32]=True
# new_gene4.indel_length[new_gene4.positions==32]=-12
# new_gene4._translate_sequence()
#
# def test_Gene_table_mutations_wrt_4():
#
#     MUTATIONS=new_gene4.table_mutations_wrt(test_gene)
#
#     assert list(MUTATIONS['GENE'])==['F']
#     assert list(MUTATIONS['MUTATION'])==['32_indel']
#     assert list(MUTATIONS['REF'])==[None]
#     assert list(MUTATIONS['ALT'])==[None]
#     assert list(MUTATIONS['POSITION'])==[32]
#     assert list(MUTATIONS['AMINO_ACID_NUMBER'])==[11]
#     assert list(MUTATIONS['NUCLEOTIDE_NUMBER'])==[32]
#     assert list(MUTATIONS['IS_SNP'])==[False]
#     assert list(MUTATIONS['IS_INDEL'])==[True]
#     assert list(MUTATIONS['IN_CDS'])==[True]
#     assert list(MUTATIONS['IN_PROMOTER'])==[False]
#     assert list(MUTATIONS['INDEL_LENGTH'])==[-12]
#     assert list(MUTATIONS['INDEL_1'])==['32_del']
#     assert list(MUTATIONS['INDEL_2'])==['32_del_12']
#     assert list(MUTATIONS['MUTATION_TYPE'])==['INDEL']
#
# new_gene5=copy.deepcopy(test_gene)
# new_gene5.is_indel[new_gene5.positions==-10]=True
# new_gene5.indel_length[new_gene5.positions==-10]=4
# new_gene5._translate_sequence()
#
# def test_Gene_table_mutations_wrt_5():
#
#     MUTATIONS=new_gene5.table_mutations_wrt(test_gene)
#
#     assert list(MUTATIONS['GENE'])==['F']
#     assert list(MUTATIONS['MUTATION'])==['-10_indel']
#     assert list(MUTATIONS['REF'])==[None]
#     assert list(MUTATIONS['ALT'])==[None]
#     assert list(MUTATIONS['POSITION'])==[-10]
#     assert list(MUTATIONS['AMINO_ACID_NUMBER'])==[None]
#     assert list(MUTATIONS['NUCLEOTIDE_NUMBER'])==[-10]
#     assert list(MUTATIONS['IS_SNP'])==[False]
#     assert list(MUTATIONS['IS_INDEL'])==[True]
#     assert list(MUTATIONS['IN_CDS'])==[False]
#     assert list(MUTATIONS['IN_PROMOTER'])==[True]
#     assert list(MUTATIONS['INDEL_LENGTH'])==[4]
#     assert list(MUTATIONS['INDEL_1'])==['-10_ins']
#     assert list(MUTATIONS['INDEL_2'])==['-10_ins_4']
#     assert list(MUTATIONS['MUTATION_TYPE'])==['INDEL']
#
# new_gene6=copy.deepcopy(test_gene)
# new_gene6.sequence[new_gene6.positions==-1]='t'
# new_gene6._translate_sequence()
#
# def test_Gene_table_mutations_wrt_6():
#
#     MUTATIONS=new_gene6.table_mutations_wrt(test_gene)
#
#     assert list(MUTATIONS['GENE'])==['F']
#     assert list(MUTATIONS['MUTATION'])==['a-1t']
#     assert list(MUTATIONS['REF'])==['a']
#     assert list(MUTATIONS['ALT'])==['t']
#     assert list(MUTATIONS['POSITION'])==[-1]
#     assert list(MUTATIONS['AMINO_ACID_NUMBER'])==[None]
#     assert list(MUTATIONS['NUCLEOTIDE_NUMBER'])==[-1]
#     assert list(MUTATIONS['IS_SNP'])==[True]
#     assert list(MUTATIONS['IS_INDEL'])==[False]
#     assert list(MUTATIONS['IN_CDS'])==[False]
#     assert list(MUTATIONS['IN_PROMOTER'])==[True]
#     assert list(MUTATIONS['INDEL_LENGTH'])==[None]
#     assert list(MUTATIONS['INDEL_1'])==[None]
#     assert list(MUTATIONS['INDEL_2'])==[None]
#     assert list(MUTATIONS['MUTATION_TYPE'])==['SNP']
#
#
# def test_Gene_valid_variant():
#
#     # badly formed variant
#     with pytest.raises(Exception):
#         assert test_gene.valid_variant("1a>g")
#
#     assert test_gene.valid_variant("a1g")
#     assert test_gene.valid_variant("t2g")
#     assert test_gene.valid_variant("g3a")
#     assert test_gene.valid_variant("g3z")
#     assert test_gene.valid_variant("g3?")
#     assert test_gene.valid_variant("a-100g")
#     assert test_gene.valid_variant("g1620a")
#
#     # out of range position
#     with pytest.raises(Exception):
#         assert test_gene.valid_variant("a1621g")
#     with pytest.raises(Exception):
#         assert test_gene.valid_variant("a-101g")
#
#     # badly formed reference base
#     with pytest.raises(Exception):
#         assert test_gene.valid_variant("A1g")
#     with pytest.raises(Exception):
#         assert test_gene.valid_variant("p2g")
#     with pytest.raises(Exception):
#         assert test_gene.valid_variant("z2g")
#
#     # badly formed positions
#     with pytest.raises(Exception):
#         assert test_gene.valid_variant("zag")
#     with pytest.raises(Exception):
#         assert test_gene.valid_variant("z1o1g")
#
#
# # use the subset argument to speed up the genome creation
# h37rv=Genome(genbank_file="config/H37rV_v3.gbk",name="H37rV_v3",gene_subset=['katG','rpoB','rrs'])
#
#
#
# def test_katG_instantiation():
#
#     test_gene=h37rv.genes['katG']
#
#     assert test_gene.total_number_nucleotides==2323
#     assert test_gene.gene_name=="katG"
#     assert test_gene.on_noncoding_strand==True
#     assert test_gene.numbering[0]==-100
#     assert test_gene.numbering[-1]==741
#     assert test_gene.numbering[-2]==741
#     assert test_gene.index[test_gene.positions==1]==2156111
#     assert test_gene.index[test_gene.positions==2]==2156110
#     assert test_gene.index[test_gene.positions==-1]==2156112
#     assert test_gene.sequence[0]=='g'
#     assert test_gene.sequence[-1]=='a'
#     assert test_gene.sequence[-2]=='g'
#
#     amino_acids="VPEQHPPITETTTGAASNGCPVVGHMKYPVEGGGNQDWWPNRLNLKVLHQNPAVADPMGAAFDYAAEVATIDVDALTRDIEEVMTTSQPWWPADYGHYGPLFIRMAWHAAGTYRIHDGRGGAGGGMQRFAPLNSWPDNASLDKARRLLWPVKKKYGKKLSWADLIVFAGNCALESMGFKTFGFGFGRVDQWEPDEVYWGKEATWLGDERYSGKRDLENPLAAVQMGLIYVNPEGPNGNPDPMAAAVDIRETFRRMAMNDVETAALIVGGHTFGKTHGAGPADLVGPEPEAAPLEQMGLGWKSSYGTGTGKDAITSGIEVVWTNTPTKWDNSFLEILYGYEWELTKSPAGAWQYTAKDGAGAGTIPDPFGGPGRSPTMLATDLSLRVDPIYERITRRWLEHPEELADEFAKAWYKLIHRDMGPVARYLGPLVPKQTLLWQDPVPAVSHDLVGEAEIASLKSQIRASGLTVSQLVSTAWAAASSFRGSDKRGGANGGRIRLQPQVGWEVNDPDGDLRKVIRTLEEIQESFNSAAPGNIKVSFADLVVLGGCAAIEKAAKAAGHNITVPFTPGRTDASQEQTDVESFAVLEPKADGFRNYLGKGNPLPAEYMLLDKANLLTLSAPEMTVLVGGLRVLGANYKRLPLGVFTEASESLTNDFFVNLLDMGITWEPSPADDGTYQGKDGSGKVKWTGSRVDLVFGSNSELRALVEVYGADDAQPKFVQDFVAAWDKVMNLDRFDVR!"
#
#     assert amino_acids==''.join(test_gene.amino_acid_sequence)
#
#
# def test_rpoB_instantiation():
#
#     test_gene=h37rv.genes['rpoB']
#
#     assert test_gene.total_number_nucleotides==3619
#     assert test_gene.gene_name=="rpoB"
#     assert test_gene.on_noncoding_strand==False
#     assert test_gene.numbering[0]==-100
#     assert test_gene.numbering[-1]==1173
#     assert test_gene.index[test_gene.positions==1]==759807
#     assert test_gene.index[test_gene.positions==2]==759808
#     assert test_gene.index[test_gene.positions==-1]==759806
#     assert test_gene.sequence[0]=='c'
#     assert test_gene.sequence[-1]=='a'
#     assert test_gene.sequence[-2]=='a'
#
#     amino_acids="LADSRQSKTAASPSPSRPQSSSNNSVPGAPNRVSFAKLREPLEVPGLLDVQTDSFEWLIGSPRWRESAAERGDVNPVGGLEEVLYELSPIEDFSGSMSLSFSDPRFDDVKAPVDECKDKDMTYAAPLFVTAEFINNNTGEIKSQTVFMGDFPMMTEKGTFIINGTERVVVSQLVRSPGVYFDETIDKSTDKTLHSVKVIPSRGAWLEFDVDKRDTVGVRIDRKRRQPVTVLLKALGWTSEQIVERFGFSEIMRSTLEKDNTVGTDEALLDIYRKLRPGEPPTKESAQTLLENLFFKEKRYDLARVGRYKVNKKLGLHVGEPITSSTLTEEDVVATIEYLVRLHEGQTTMTVPGGVEVPVETDDIDHFGNRRLRTVGELIQNQIRVGMSRMERVVRERMTTQDVEAITPQTLINIRPVVAAIKEFFGTSQLSQFMDQNNPLSGLTHKRRLSALGPGGLSRERAGLEVRDVHPSHYGRMCPIETPEGPNIGLIGSLSVYARVNPFGFIETPYRKVVDGVVSDEIVYLTADEEDRHVVAQANSPIDADGRFVEPRVLVRRKAGEVEYVPSSEVDYMDVSPRQMVSVATAMIPFLEHDDANRALMGANMQRQAVPLVRSEAPLVGTGMELRAAIDAGDVVVAEESGVIEEVSADYITVMHDNGTRRTYRMRKFARSNHGTCANQCPIVDAGDRVEAGQVIADGPCTDDGEMALGKNLLVAIMPWEGHNYEDAIILSNRLVEEDVLTSIHIEEHEIDARDTKLGAEEITRDIPNISDEVLADLDERGIVRIGAEVRDGDILVGKVTPKGETELTPEERLLRAIFGEKAREVRDTSLKVPHGESGKVIGIRVFSREDEDELPAGVNELVRVYVAQKRKISDGDKLAGRHGNKGVIGKILPVEDMPFLADGTPVDIILNTHGVPRRMNIGQILETHLGWCAHSGWKVDAAKGVPDWAARLPDELLEAQPNAIVSTPVFDGAQEAELQGLLSCTLPNRDGDVLVDADGKAMLFDGRSGEPFPYPVTVGYMYIMKLHHLVDDKIHARSTGPYSMITQQPLGGKAQFGGQRFGEMECWAMQAYGAAYTLQELLTIKSDDTVGRVKVYEAIVKGENIPEPGIPESFKVLLKELQSLCLNVEVLSSDGAAIELREGEDEDLERAAANLGINLSRNESASVEDLA!"
#
#     assert amino_acids==''.join(test_gene.amino_acid_sequence)
#
#
# def test_rrs_instantiation():
#
#     test_gene=h37rv.genes['rrs']
#     assert test_gene.total_number_nucleotides==1637
#     assert test_gene.gene_name=="rrs"
#     assert test_gene.on_noncoding_strand==False
#     assert test_gene.numbering[0]==-100
#     assert test_gene.numbering[-1]==1537
#     assert test_gene.index[test_gene.positions==1]==1471846
#     assert test_gene.index[test_gene.positions==2]==1471847
#     assert test_gene.index[test_gene.positions==-1]==1471845
#
#     sequence="ggccatgctcttgatgccccgttgtcgggggcgtggccgtttgttttgtcaggatatttctaaatacctttggctcccttttccaaagggagtgtttgggttttgtttggagagtttgatcctggctcaggacgaacgctggcggcgtgcttaacacatgcaagtcgaacggaaaggtctcttcggagatactcgagtggcgaacgggtgagtaacacgtgggtgatctgccctgcacttcgggataagcctgggaaactgggtctaataccggataggaccacgggatgcatgtcttgtggtggaaagcgctttagcggtgtgggatgagcccgcggcctatcagcttgttggtggggtgacggcctaccaaggcgacgacgggtagccggcctgagagggtgtccggccacactgggactgagatacggcccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagcgacgccgcgtgggggatgacggccttcgggttgtaaacctctttcaccatcgacgaaggtccgggttctctcggattgacggtaggtggagaagaagcaccggccaactacgtgccagcagccgcggtaatacgtagggtgcgagcgttgtccggaattactgggcgtaaagagctcgtaggtggtttgtcgcgttgttcgtgaaatctcacggcttaactgtgagcgtgcgggcgatacgggcagactagagtactgcaggggagactggaattcctggtgtagcggtggaatgcgcagatatcaggaggaacaccggtggcgaaggcgggtctctgggcagtaactgacgctgaggagcgaaagcgtggggagcgaacaggattagataccctggtagtccacgccgtaaacggtgggtactaggtgtgggtttccttccttgggatccgtgccgtagctaacgcattaagtaccccgcctggggagtacggccgcaaggctaaaactcaaaggaattgacgggggcccgcacaagcggcggagcatgtggattaattcgatgcaacgcgaagaaccttacctgggtttgacatgcacaggacgcgtctagagataggcgttcccttgtggcctgtgtgcaggtggtgcatggctgtcgtcagctcgtgtcgtgagatgttgggttaagtcccgcaacgagcgcaacccttgtctcatgttgccagcacgtaatggtggggactcgtgagagactgccggggtcaactcggaggaaggtggggatgacgtcaagtcatcatgccccttatgtccagggcttcacacatgctacaatggccggtacaaagggctgcgatgccgcgaggttaagcgaatccttaaaagccggtctcagttcggatcggggtctgcaactcgaccccgtgaagtcggagtcgctagtaatcgcagatcagcaacgctgcggtgaatacgttcccgggccttgtacacaccgcccgtcacgtcatgaaagtcggtaacacccgaagccagtggcctaaccctcgggagggagctgtcgaaggtgggatcggcgattgggacgaagtcgtaacaaggtagccgtaccggaaggtgcggctggatcacctcctttct"
#
#     assert sequence==''.join(test_gene.sequence)
#     assert test_gene.amino_acid_sequence==None
#
#     new_gene=copy.deepcopy(test_gene)
#     new_gene.sequence[new_gene.positions==1]='c'
#     new_gene._translate_sequence()
#
#     MUTATIONS=new_gene.table_mutations_wrt(test_gene)
#
#     assert list(MUTATIONS['GENE'])==['rrs']
#     assert list(MUTATIONS['MUTATION'])==['t1c']
#     assert list(MUTATIONS['REF'])==['t']
#     assert list(MUTATIONS['ALT'])==['c']
#     assert list(MUTATIONS['POSITION'])==[1]
#     assert numpy.isnan(MUTATIONS['AMINO_ACID_NUMBER'])
#     assert list(MUTATIONS['NUCLEOTIDE_NUMBER'])==[1]
#     assert list(MUTATIONS['IS_SNP'])==[True]
#     assert list(MUTATIONS['IS_INDEL'])==[False]
#     assert list(MUTATIONS['IN_CDS'])==[True]
#     assert list(MUTATIONS['IN_PROMOTER'])==[False]
#     assert list(MUTATIONS['INDEL_LENGTH'])==[None]
#     assert list(MUTATIONS['INDEL_1'])==[None]
#     assert list(MUTATIONS['INDEL_2'])==[None]
#     assert list(MUTATIONS['MUTATION_TYPE'])==['SNP']
