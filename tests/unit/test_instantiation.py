import numpy, gumpy, pytest, math


#As BioPython thinks that the locus line of the TEST-RNA.gbk is malformed, it gives a warning
#So ignore it to stop failing tests...
pytestmark = pytest.mark.filterwarnings("ignore")

@pytest.mark.tb
def test_instanciate_genome_tb():

    tb_reference = gumpy.Genome('config/NC_000962.3.gbk',gene_subset=['katG','rpoB','pncA','Rv2042c','rrs'])

    assert len(tb_reference)==4411532
    assert tb_reference.name=='NC_000962'
    assert tb_reference.id=='NC_000962.3'

    # check to see if the stacking can cope with two genes overlapping
    mask=tb_reference.stacked_nucleotide_index==2288681
    assert set(tb_reference.stacked_gene_name[mask]) == {'pncA','Rv2042c'}, 'stacking failing at 2288681 in H37rV v3'

    # ..and the at_index function
    assert set(tb_reference.at_index(2288681)) == {'pncA','Rv2042c'}, 'not correctly detecting that pncA and Rv2042c both include 2288681 in NC_000962.3'

    assert tb_reference.contains_gene('katG')
    assert tb_reference.contains_gene('rpoB')
    assert tb_reference.contains_gene('pncA')
    assert tb_reference.contains_gene('Rv2042c')
    #The ~ operator uses bitwise NOT rather than logical NOT, ~True == -2, ~False == -1 which always pass assertations
    #The ~ operator only works this way for numpy arrays: ~[True, False] == [False, True] for numpy arrays
    assert not tb_reference.contains_gene('rpoC')

    # check the first and last dozen bases of the whole sequence
    assert ''.join(i for i in tb_reference.nucleotide_sequence[:12]) == 'ttgaccgatgac'
    assert ''.join(i for i in tb_reference.nucleotide_sequence[-12:]) == 'ggagatacgtcg'

@pytest.mark.tb
def test_instanciate_genes_tb():

    tb_reference = gumpy.Genome('config/NC_000962.3.gbk',gene_subset=['katG','rpoB','pncA','Rv2042c','rrs'])

    truth_gene_sequence={}
    truth_gene_sequence['katG']='VPEQHPPITETTTGAASNGCPVVGHMKYPVEGGGNQDWWPNRLNLKVLHQNPAVADPMGAAFDYAAEVATIDVDALTRDIEEVMTTSQPWWPADYGHYGPLFIRMAWHAAGTYRIHDGRGGAGGGMQRFAPLNSWPDNASLDKARRLLWPVKKKYGKKLSWADLIVFAGNCALESMGFKTFGFGFGRVDQWEPDEVYWGKEATWLGDERYSGKRDLENPLAAVQMGLIYVNPEGPNGNPDPMAAAVDIRETFRRMAMNDVETAALIVGGHTFGKTHGAGPADLVGPEPEAAPLEQMGLGWKSSYGTGTGKDAITSGIEVVWTNTPTKWDNSFLEILYGYEWELTKSPAGAWQYTAKDGAGAGTIPDPFGGPGRSPTMLATDLSLRVDPIYERITRRWLEHPEELADEFAKAWYKLIHRDMGPVARYLGPLVPKQTLLWQDPVPAVSHDLVGEAEIASLKSQIRASGLTVSQLVSTAWAAASSFRGSDKRGGANGGRIRLQPQVGWEVNDPDGDLRKVIRTLEEIQESFNSAAPGNIKVSFADLVVLGGCAAIEKAAKAAGHNITVPFTPGRTDASQEQTDVESFAVLEPKADGFRNYLGKGNPLPAEYMLLDKANLLTLSAPEMTVLVGGLRVLGANYKRLPLGVFTEASESLTNDFFVNLLDMGITWEPSPADDGTYQGKDGSGKVKWTGSRVDLVFGSNSELRALVEVYGADDAQPKFVQDFVAAWDKVMNLDRFDVR!'
    truth_gene_sequence['rpoB']='LADSRQSKTAASPSPSRPQSSSNNSVPGAPNRVSFAKLREPLEVPGLLDVQTDSFEWLIGSPRWRESAAERGDVNPVGGLEEVLYELSPIEDFSGSMSLSFSDPRFDDVKAPVDECKDKDMTYAAPLFVTAEFINNNTGEIKSQTVFMGDFPMMTEKGTFIINGTERVVVSQLVRSPGVYFDETIDKSTDKTLHSVKVIPSRGAWLEFDVDKRDTVGVRIDRKRRQPVTVLLKALGWTSEQIVERFGFSEIMRSTLEKDNTVGTDEALLDIYRKLRPGEPPTKESAQTLLENLFFKEKRYDLARVGRYKVNKKLGLHVGEPITSSTLTEEDVVATIEYLVRLHEGQTTMTVPGGVEVPVETDDIDHFGNRRLRTVGELIQNQIRVGMSRMERVVRERMTTQDVEAITPQTLINIRPVVAAIKEFFGTSQLSQFMDQNNPLSGLTHKRRLSALGPGGLSRERAGLEVRDVHPSHYGRMCPIETPEGPNIGLIGSLSVYARVNPFGFIETPYRKVVDGVVSDEIVYLTADEEDRHVVAQANSPIDADGRFVEPRVLVRRKAGEVEYVPSSEVDYMDVSPRQMVSVATAMIPFLEHDDANRALMGANMQRQAVPLVRSEAPLVGTGMELRAAIDAGDVVVAEESGVIEEVSADYITVMHDNGTRRTYRMRKFARSNHGTCANQCPIVDAGDRVEAGQVIADGPCTDDGEMALGKNLLVAIMPWEGHNYEDAIILSNRLVEEDVLTSIHIEEHEIDARDTKLGAEEITRDIPNISDEVLADLDERGIVRIGAEVRDGDILVGKVTPKGETELTPEERLLRAIFGEKAREVRDTSLKVPHGESGKVIGIRVFSREDEDELPAGVNELVRVYVAQKRKISDGDKLAGRHGNKGVIGKILPVEDMPFLADGTPVDIILNTHGVPRRMNIGQILETHLGWCAHSGWKVDAAKGVPDWAARLPDELLEAQPNAIVSTPVFDGAQEAELQGLLSCTLPNRDGDVLVDADGKAMLFDGRSGEPFPYPVTVGYMYIMKLHHLVDDKIHARSTGPYSMITQQPLGGKAQFGGQRFGEMECWAMQAYGAAYTLQELLTIKSDDTVGRVKVYEAIVKGENIPEPGIPESFKVLLKELQSLCLNVEVLSSDGAAIELREGEDEDLERAAANLGINLSRNESASVEDLA!'
    truth_gene_sequence['pncA']='MRALIIVDVQNDFCEGGSLAVTGGAALARAISDYLAEAADYHHVVATKDFHIDPGDHFSGTPDYSSSWPPHCVSGTPGADFHPSLDTSAIEAVFYKGAYTGAYSGFEGVDENGTPLLNWLRQRGVDEVDVVGIATDHCVRQTAEDAVRNGLATRVLVDLTAGVSADTTVAALEEMRTASVELVCSS!'
    truth_gene_sequence['Rv2042c']='MAPPNRDELLAAVERSPQAAAAHDRAGWVGLFTGDARVEDPVGSQPQVGHEAIGRFYDTFIGPRDITFHRDLDIVSGTVVLRDLELEVAMDSAVTVFIPAFLRYDLRPVTGEWQIAALRAYWELPAMMLQFLRTGSGATRPALQLSRALLGNQGLGGTAGFLTGFRRAGRRHKKLVETFLNAASRADKSAAYHALSRTATMTLGEDELLDIVELFEQLRGASWTKVTGAGSTVAVSLASDHRRGIMFADVPWRGNRINRIRYFPA!'
    truth_gene_sequence['rrs']='ttttgtttggagagtttgatcctggctcaggacgaacgctggcggcgtgcttaacacatgcaagtcgaacggaaaggtctcttcggagatactcgagtggcgaacgggtgagtaacacgtgggtgatctgccctgcacttcgggataagcctgggaaactgggtctaataccggataggaccacgggatgcatgtcttgtggtggaaagcgctttagcggtgtgggatgagcccgcggcctatcagcttgttggtggggtgacggcctaccaaggcgacgacgggtagccggcctgagagggtgtccggccacactgggactgagatacggcccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagcgacgccgcgtgggggatgacggccttcgggttgtaaacctctttcaccatcgacgaaggtccgggttctctcggattgacggtaggtggagaagaagcaccggccaactacgtgccagcagccgcggtaatacgtagggtgcgagcgttgtccggaattactgggcgtaaagagctcgtaggtggtttgtcgcgttgttcgtgaaatctcacggcttaactgtgagcgtgcgggcgatacgggcagactagagtactgcaggggagactggaattcctggtgtagcggtggaatgcgcagatatcaggaggaacaccggtggcgaaggcgggtctctgggcagtaactgacgctgaggagcgaaagcgtggggagcgaacaggattagataccctggtagtccacgccgtaaacggtgggtactaggtgtgggtttccttccttgggatccgtgccgtagctaacgcattaagtaccccgcctggggagtacggccgcaaggctaaaactcaaaggaattgacgggggcccgcacaagcggcggagcatgtggattaattcgatgcaacgcgaagaaccttacctgggtttgacatgcacaggacgcgtctagagataggcgttcccttgtggcctgtgtgcaggtggtgcatggctgtcgtcagctcgtgtcgtgagatgttgggttaagtcccgcaacgagcgcaacccttgtctcatgttgccagcacgtaatggtggggactcgtgagagactgccggggtcaactcggaggaaggtggggatgacgtcaagtcatcatgccccttatgtccagggcttcacacatgctacaatggccggtacaaagggctgcgatgccgcgaggttaagcgaatccttaaaagccggtctcagttcggatcggggtctgcaactcgaccccgtgaagtcggagtcgctagtaatcgcagatcagcaacgctgcggtgaatacgttcccgggccttgtacacaccgcccgtcacgtcatgaaagtcggtaacacccgaagccagtggcctaaccctcgggagggagctgtcgaaggtgggatcggcgattgggacgaagtcgtaacaaggtagccgtaccggaaggtgcggctggatcacctcctttct'

    for gene_name in ['katG','rpoB','pncA','Rv2042c','rrs']:

        assert tb_reference.contains_gene(gene_name)

        gene=tb_reference.build_gene(gene_name)

        if gene_name=='rrs':
            assert not gene.codes_protein, gene_name
            assert gene.feature_type=='RNA', gene_name
        elif gene_name=='Rv2042c':
            assert gene.codes_protein, gene_name
            assert gene.feature_type=='LOCUS', gene_name
        else:
            assert gene.codes_protein, gene_name
            assert gene.feature_type=='GENE', gene_name

        if gene_name in ['katG','pncA','Rv2042c']:
            assert gene.reverse_complement, gene_name
        else:
            assert not gene.reverse_complement, gene_name

        if gene_name=='rrs':
            nuc_sequence=gene.nucleotide_sequence[gene.nucleotide_number>0]
            assert ''.join(i for i in nuc_sequence) == truth_gene_sequence[gene_name], gene_name+' nucleotide sequence incorrect!'
        else:
            assert ''.join(i for i in gene.amino_acid_sequence) == truth_gene_sequence[gene_name], gene_name+' amino acid sequence incorrect!'

def test_instanciate_genome_covid():

    reference = gumpy.Genome('config/NC_045512.2.gbk')

    assert len(reference)==29903
    assert reference.name=='NC_045512'
    assert reference.id=='NC_045512.2'

    # check to see if the stacking can cope with two genes overlapping
    assert set(reference.at_index(27756)) == {'ORF7a','ORF7b'}, 'not correctly detecting that ORF7a and ORF7b both include 22756-27759 incl. in SARS-CoV_2'
    assert set(reference.at_index(27757)) == {'ORF7a','ORF7b'}, 'not correctly detecting that ORF7a and ORF7b both include 22756-27759 incl. in SARS-CoV_2'
    assert set(reference.at_index(27758)) == {'ORF7a','ORF7b'}, 'not correctly detecting that ORF7a and ORF7b both include 22756-27759 incl. in SARS-CoV_2'
    assert set(reference.at_index(27759)) == {'ORF7a','ORF7b'}, 'not correctly detecting that ORF7a and ORF7b both include 22756-27759 incl. in SARS-CoV_2'

    #
    assert reference.contains_gene('ORF1ab')
    assert reference.contains_gene('S')
    assert reference.contains_gene('ORF7a')
    assert not reference.contains_gene('Rv2042c')

    # # check the first and last dozen bases of the whole sequence
    assert ''.join(i for i in reference.nucleotide_sequence[:12]) == 'attaaaggttta'
    assert ''.join(i for i in reference.nucleotide_sequence[-12:]) == 'aaaaaaaaaaaa'


def test_instanciate_genes_covid():

    reference = gumpy.Genome('config/NC_045512.2.gbk')

    assert reference.contains_gene('S')
    gene=reference.build_gene('S')
    assert not gene.reverse_complement
    assert gene.name=='S'
    assert gene.codes_protein
    assert gene.feature_type=='GENE'
    assert ''.join(i for i in gene.amino_acid_sequence) == "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT!"

    assert reference.contains_gene('ORF7a')
    gene=reference.build_gene('ORF7a')
    assert ''.join(i for i in gene.amino_acid_sequence) == 'MKIILFLALITLATCELYHYQECVRGTTVLLKEPCSSGTYEGNSPFHPLADNKFALTCFSTQFAFACPDGVKHVYQLRARSVSPKLFIRQEEVQELYSPIFLIVAAIVFITLCFTLKRKTE!'

    assert reference.contains_gene('ORF7b')
    gene=reference.build_gene('ORF7b')
    assert ''.join(i for i in gene.amino_acid_sequence) == 'MIELSLIDFYLCFLAFLLFLVLIMLIIFWFSLELQDHNETCHA!'

def test_instanciate_genes_covid_ORF1ab():

    reference = gumpy.Genome('config/NC_045512.2.gbk')

    assert reference.contains_gene('ORF1ab')

    gene=reference.build_gene('ORF1ab')

    assert ''.join(i for i in gene.amino_acid_sequence) == 'MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGGAYTRYVDNNFCGPDGYPLECIKDLLARAGKASCTLSEQLDFIDTKRGVYCCREHEHEIAWYTERSEKSYELQTPFEIKLAKKFDTFNGECPNFVFPLNSIIKTIQPRVEKKKLDGFMGRIRSVYPVASPNECNQMCLSTLMKCDHCGETSWQTGDFVKATCEFCGTENLTKEGATTCGYLPQNAVVKIYCPACHNSEVGPEHSLAEYHNESGLKTILRKGGRTIAFGGCVFSYVGCHNKCAYWVPRASANIGCNHTGVVGEGSEGLNDNLLEILQKEKVNINIVGDFKLNEEIAIILASFSASTSAFVETVKGLDYKAFKQIVESCGNFKVTKGKAKKGAWNIGEQKSILSPLYAFASEAARVVRSIFSRTLETAQNSVRVLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYITGGVVQLTSQWLTNIFGTVYEKLKPVLDWLEEKFKEGVEFLRDGWEIVKFISTCACEIVGGQIVTCAKEIKESVQTFFKLVNKFLALCADSIIIGGAKLKALNLGETFVTHSKGLYRKCVKSREETGLLMPLKAPKEIIFLEGETLPTEVLTEEVVLKTGDLQPLEQPTSEAVEAPLVGTPVCINGLMLLEIKDTEKYCALAPNMMVTNNTFTLKGGAPTKVTFGDDTVIEVQGYKSVNITFELDERIDKVLNEKCSAYTVELGTEVNEFACVVADAVIKTLQPVSELLTPLGIDLDEWSMATYYLFDESGEFKLASHMYCSFYPPDEDEEEGDCEEEEFEPSTQYEYGTEDDYQGKPLEFGATSAALQPEEEQEEDWLDDDSQQTVGQQDGSEDNQTTTIQTIVEVQPQLEMELTPVVQTIEVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLRVCVDTVRTNVYLAVFDKNLYDKLVSSFLEMKSEKQVEQKIAEIPKEEVKPFITESKPSVEQRKQDDKKIKACVEEVTTTLEETKFLTENLLLYIDINGNLHPDSATLVSDIDITFLKKDAPYIVGDVVQEGVLTAVVIPTKKAGGTTEMLAKALRKVPTDNYITTYPGQGLNGYTVEEAKTVLKKCKSAFYILPSIISNEKQEILGTVSWNLREMLAHAEETRKLMPVCVETKAIVSTIQRKYKGIKIQEGVVDYGARFYFYTSKTTVASLINTLNDLNETLVTMPLGYVTHGLNLEEAARYMRSLKVPATVSVSSPDAVTAYNGYLTSSSKTPEEHFIETISLAGSYKDWSYSGQSTQLGIEFLKRGDKSVYYTSNPTTFHLDGEVITFDNLKTLLSLREVRTIKVFTTVDNINLHTQVVDMSMTYGQQFGPTYLDGADVTKIKPHNSHEGKTFYVLPNDDTLRVEAFEYYHTTDPSFLGRYMSALNHTKKWKYPQVNGLTSIKWADNNCYLATALLTLQQIELKFNPPALQDAYYRARAGEAANFCALILAYCNKTVGELGDVRETMSYLFQHANLDSCKRVLNVVCKTCGQQQTTLKGVEAVMYMGTLSYEQFKKGVQIPCTCGKQATKYLVQQESPFVMMSAPPAQYELKHGTFTCASEYTGNYQCGHYKHITSKETLYCIDGALLTKSSEYKGPITDVFYKENSYTTTIKPVTYKLDGVVCTEIDPKLDNYYKKDNSYFTEQPIDLVPNQPYPNASFDNFKFVCDNIKFADDLNQLTGYKKPASRELKVTFFPDLNGDVVAIDYKHYTPSFKKGAKLLHKPIVWHVNNATNKATYKPNTWCIRCLWSTKPVETSNSFDVLKSEDAQGMDNLACEDLKPVSEEVVENPTIQKDVLECNVKTTEVVGDIILKPANNSLKITEEVGHTDLMAAYVDNSSLTIKKPNELSRVLGLKTLATHGLAAVNSVPWDTIANYAKPFLNKVVSTTTNIVTRCLNRVCTNYMPYFFTLLLQLCTFTRSTNSRIKASMPTTIAKNTVKSVGKFCLEASFNYLKSPNFSKLINIIIWFLLLSVCLGSLIYSTAALGVLMSNLGMPSYCTGYREGYLNSTNVTIATYCTGSIPCSVCLSGLDSLDTYPSLETIQITISSFKWDLTAFGLVAEWFLAYILFTRFFYVLGLAAIMQLFFSYFAVHFISNSWLMWLIINLVQMAPISAMVRMYIFFASFYYVWKSYVHVVDGCNSSTCMMCYKRNRATRVECTTIVNGVRRSFYVYANGGKGFCKLHNWNCVNCDTFCAGSTFISDEVARDLSLQFKRPINPTDQSSYIVDSVTVKNGSIHLYFDKAGQKTYERHSLSHFVNLDNLRANNTKGSLPINVIVFDGKSKCEESSAKSASVYYSQLMCQPILLLDQALVSDVGDSAEVAVKMFDAYVNTFSSTFNVPMEKLKTLVATAEAELAKNVSLDNVLSTFISAARQGFVDSDVETKDVVECLKLSHQSDIEVTGDSCNNYMLTYNKVENMTPRDLGACIDCSARHINAQVAKSHNIALIWNVKDFMSLSEQLRKQIRSAAKKNNLPFKLTCATTRQVVNVVTTKIALKGGKIVNNWLKQLIKVTLVFLFVAAIFYLITPVHVMSKHTDFSSEIIGYKAIDGGVTRDIASTDTCFANKHADFDTWFSQRGGSYTNDKACPLIAAVITREVGFVVPGLPGTILRTTNGDFLHFLPRVFSAVGNICYTPSKLIEYTDFATSACVLAAECTIFKDASGKPVPYCYDTNVLEGSVAYESLRPDTRYVLMDGSIIQFPNTYLEGSVRVVTTFDSEYCRHGTCERSEAGVCVSTSGRWVLNNDYYRSLPGVFCGVDAVNLLTNMFTPLIQPIGALDISASIVAGGIVAIVVTCLAYYFMRFRRAFGEYSHVVAFNTLLFLMSFTVLCLTPVYSFLPGVYSVIYLYLTFYLTNDVSFLAHIQWMVMFTPLVPFWITIAYIICISTKHFYWFFSNYLKRRVVFNGVSFSTFEEAALCTFLLNKEMYLKLRSDVLLPLTQYNRYLALYNKYKYFSGAMDTTSYREAACCHLAKALNDFSNSGSDVLYQPPQTSITSAVLQSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQSAVKRTIKGTHHWLLLTILTSLLVLVQSTQWSLFFFLYENAFLPFAMGIIAMSAFAMMFVKHKHAFLCLFLLPSLATVAYFNMVYMPASWVMRIMTWLDMVDTSLSGFKLKDCVMYASAVVLLILMTARTVYDDGARRVWTLMNVLTLVYKVYYGNALDQAISMWALIISVTSNYSGVVTTVMFLARGIVFMCVEYCPIFFITGNTLQCIMLVYCFLGYFCTCYFGLFCLLNRYFRLTLGVYDYLVSTQEFRYMNSQGLLPPKNSIDAFKLNIKLLGVGGKPCIKVATVQSKMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQGAVDINKLCEEMLDNRATLQAIASEFSSLPSYAAFATAQEAYEQAVANGDSEVVLKKLKKSLNVAKSEFDRDAAMQRKLEKMADQAMTQMYKQARSEDKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIVQLSEISMDNSPNLAWPLIVTALRANSAVKLQNNELSPVALRQMSCAAGTTQTACTDDNALAYYNTTKGGRFVLALLSDLQDLKWARFPKSDGTGTIYTELEPPCRFVTDTPKGPKVKYLYFIKGLNNLNRGMVLGSLAATVRLQAGNATEVPANSTVLSFCAFAVDAAKAYKDYLASGGQPITNCVKMLCTHTGTGQAITVTPEANMDQESFGGASCCLYCRCHIDHPNPKGFCDLKGKYVQIPTTCANDPVGFTLKNTVCTVCGMWKGYGCSCDQLREPMLQSADAQSFLNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAVAKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGIVGVLTLDNQDLNGNWYDFGDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQTYHPNCVNCLDDRCILHCANFNVLFSTVFPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSRLSFKELLVYAADPAMHAASGNLLLDKRTTCFSVAALTNNVAFQTVKPGNFNKDFYDFAVSKGFFKEGSSVELKHFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNNLDKSAGFPFNKWGKARLYYDSMSYEDQDALFAYTKRNVIPTITQMNLKYAISAKNRARTVAGVSICSTMTNRQFHQKLLKSIAATRGATVVIGTSKFYGGWHNMLKTVYSDVENPHLMGWDYPKCDRAMPNMLRIMASLVLARKHTTCCSLSHRFYRLANECAQVLSEMVMCGGSLYVKPGGTSSGDATTAYANSVFNICQAVTANVNALLSTDGNKIADKYVRNLQHRLYECLYRNRDVDTDFVNEFYAYLRKHFSMMILSDDAVVCFNSTYASQGLVASIKNFKSVLYYQNNVFMSEAKCWTETDLTKGPHEFCSQHTMLVKQGDDYVYLPYPDPSRILGAGCFVDDIVKTDGTLMIERFVSLAIDAYPLTKHPNQEYADVFHLYLQYIRKLHDELTGHMLDMYSVMLTNDNTSRYWEPEFYEAMYTPHTVLQAVGACVLCNSQTSLRCGACIRRPFLCCKCCYDHVISTSHKLVLSVNPYVCNAPGCDVTDVTQLYLGGMSYYCKSHKPPISFPLCANGQVFGLYKNTCVGSDNVTDFNAIATCDWTNAGDYILANTCTERLKLFAAETLKATEETFKLSYGIATVREVLSDRELHLSWEVGKPRPPLNRNYVFTGYRVTKNSKVQIGEYTFEKGDYGDAVVYRGTTTYKLNVGDYFVLTSHTVMPLSAPTLVPQEHYVRITGLYPTLNISDEFSSNVANYQKVGMQKYSTLQGPPGTGKSHFAIGLALYYPSARIVYTACSHAAVDALCEKALKYLPIDKCSRIIPARARVECFDKFKVNSTLEQYVFCTVNALPETTADIVVFDEISMATNYDLSVVNARLRAKHYVYIGDPAQLPAPRTLLTKGTLEPEYFNSVCRLMKTIGPDMFLGTCRRCPAEIVDTVSALVYDNKLKAHKDKSAQCFKMFYKGVITHDVSSAINRPQIGVVREFLTRNPAWRKAVFISPYNSQNAVASKILGLPTQTVDSSQGSEYDYVIFTQTTETAHSCNVNRFNVAITRAKVGILCIMSDRDLYDKLQFTSLEIPRRNVATLQAENVTGLFKDCSKVITGLHPTQAPTHLSVDTKFKTEGLCVDIPGIPKDMTYRRLISMMGFKMNYQVNGYPNMFITREEAIRHVRAWIGFDVEGCHATREAVGTNLPLQLGFSTGVNLVAVPTGYVDTPNNTDFSRVSAKPPPGDQFKHLIPLMYKGLPWNVVRIKIVQMLSDTLKNLSDRVVFVLWAHGFELTSMKYFVKIGPERTCCLCDRRATCFSTASDTYACWHHSIGFDYVYNPFMIDVQQWGFTGNLQSNHDLYCQVHGNAHVASCDAIMTRCLAVHECFVKRVDWTIEYPIIGDELKINAACRKVQHMVVKAALLADKFPVLHDIGNPKAIKCVPQADVEWKFYDAQPCSDKAYKIEELFYSYATHSDKFTDGVCLFWNCNVDRYPANSIVCRFDTRVLSNLNLPGCDGGSLYVNKHAFHTPAFDKSAFVNLKQLPFFYYSDSPCESHGKQVVSDIDYVPLKSATCITRCNLGGAVCRHHANEYRLYLDAYNMMISAGFSLWVYKQFDTYNLWNTFTRLQSLENVAFNVVNKGHFDGQQGEVPVSIINNTVYTKVDGVDVELFENKTTLPVNVAFELWAKRNIKPVPEVKILNNLGVDIAANTVIWDYKRDAPAHISTIGVCSMTDIAKKPTETICAPLTVFFDGRVDGQVDLFRNARNGVLITEGSVKGLQPSVGPKQASLNGVTLIGEAVKTQFNYYKKVDGVVQQLPETYFTQSRNLQEFKPRSQMEIDFLELAMDEFIERYKLEGYAFEHIVYGDFSHSQLGGLHLLIGLAKRFKESPFELEDFIPMDSTVKNYFITDAQTGSSKCVCSVIDLLLDDFVEIIKSQDLSVVSKVVKVTIDYTEISFMLWCKDGHVETFYPKLQSSQAWQPGVAMPNLYKMQRMLLEKCDLQNYGDSATLPKGIMMNVAKYTQLCQYLNTLTLAVPYNMRVIHFGAGSDKGVAPGTAVLRQWLPTGTLLVDSDLNDFVSDADSTLIGDCATVHTANKWDLIISDMYDPKTKNVTKENDSKEGFFTYICGFIQQKLALGGSVAIKITEHSWNADLYKLMGHFAWWTAFVTNVNASSSEAFLIGCNYLGKPREQIDGYVMHANYIFWRNTNPIQLSSYSLFDMSKFPLKLRGTAVMSLKEGQINDMILSLLSKGRLIIRENNRVVISSDVLVNN!'

def test_instanciate_genome_rna():
    genome = gumpy.Genome("config/TEST-RNA.gbk")
    g2 = gumpy.Genome("config/TEST-RNA.gbk", multithreaded=True)
    #Check that multithreading gives the same result.
    #If not testing on Linux the multithreaded code will never be touched...
    assert genome == g2

    #Testing generic attributes such as name and length
    assert len(genome) == 99
    assert genome.name == "TEST_RNA"
    assert genome.id == "TEST_RNA.1"
    assert list(genome.genes.keys()) == list("ABC")

    #Testing annotations parsed correctly
    assert genome.annotations['organism'] == "TEST_RNA_ORGANISM"
    assert genome.annotations['source'] == "TEST_RNA_SOURCE"
    assert genome.annotations['references'] == [{
                                                "location": [{
                                                     "_start": 0,
                                                     "_end": 99
                                                    }],
                                                "authors": "Test,1., Test,2., Test,3.",
                                                "consrtm": "",
                                                "title": "Test title for a reference for an RNA strand",
                                                "journal": "Test journal for an RNA strand",
                                                "medline_id": "",
                                                "pubmed_id": "1",
                                                "comment": ""
                                                 }]
    #Testing sequence was parsed correctly
    assert numpy.all(
        genome.nucleotide_sequence == numpy.array(
            list("aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc")
        ))
    assert numpy.all(
        genome.nucleotide_index == numpy.array(range(1,100)))

    #Testing all arrays were setup correctly
    #Stacked arrays are changed during promoter assignment, so should test promoter assignment too
    original_stacked_gene_name = numpy.array([
        ['', '', '', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', '', '', ''],
        ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
    ])
    #Grown out genes with promoters
    full_gene_name = numpy.array([
                ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'A', 'A', 'A'],
                ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ])
    assert numpy.all(
        genome.stacked_gene_name == full_gene_name)
    #No rev-comp as RNA only has a single strand
    assert numpy.all(
        genome.stacked_is_reverse_complement == numpy.array([[False for i in range(99)], [False for i in range(99)]])
    )
    assert numpy.all(
        genome.stacked_nucleotide_index == numpy.array([[i for i in range(1, 100)], [i for i in range(1, 100)]])
    )
    assert numpy.all(
        genome.stacked_is_promoter == numpy.array([
            [True, True, True]+[False for i in range(3, 60)]+[True for i in range(60, 90)]+[False, False, False, False, False, False, True, True, True],
            [False for i in range(99)]
        ])
    )
    assert numpy.all(
        genome.stacked_nucleotide_number == numpy.array([
            [-3, -2, -1] + [i for i in range(1, 28)] + [0 for i in range(28, 58)] + [i for i in range(-30, 0)] + [i for i in range(1, 7)] + [-6, -5, -4],
            [0 for i in range(27)] + [i for i in range(1, 34)] + [0 for i in range(39)]
        ])
    )
    assert numpy.all(
        genome.stacked_is_cds == (original_stacked_gene_name != "")
    )

    #Testing all genes instanciated (not correct genes)
    assert genome.contains_gene('A')
    gene=genome.build_gene('A')
    assert gene == gumpy.Gene(
        name="A",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[0] == "A"],
        index=genome.nucleotide_index[full_gene_name[0] == "A"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "A"],
        is_cds=genome.stacked_is_cds[full_gene_name == "A"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "A"],
        is_indel=genome.is_indel[full_gene_name[0] == "A"],
        indel_length=[0 for i in range(33)],
        feature_type="GENE"
    )
    assert genome.contains_gene('B')
    gene=genome.build_gene('B')
    assert gene == gumpy.Gene(
        name="B",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[1] == "B"],
        index=genome.nucleotide_index[full_gene_name[1] == "B"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "B"],
        is_cds=genome.stacked_is_cds[full_gene_name == "B"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "B"],
        is_indel=genome.is_indel[full_gene_name[1] == "B"],
        indel_length=[0 for i in range(33)],
        feature_type="GENE"
    )
    assert genome.contains_gene('C')
    gene=genome.build_gene('C')
    assert gene == gumpy.Gene(
        name="C",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[0] == "C"],
        index=genome.nucleotide_index[full_gene_name[0] == "C"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "C"],
        is_cds=genome.stacked_is_cds[full_gene_name == "C"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "C"],
        is_indel=genome.is_indel[full_gene_name[0] == "C"],
        indel_length=[0 for i in range(36)],
        feature_type="GENE"
    )

def test_instanciate_genes_rna():
    genome = gumpy.Genome("config/TEST-RNA.gbk")

    #If the previous checks pass, this should be equal to a freshly instanciated Gene
    gene = genome.build_gene("A")

    #Ground truth values for all genes - used to extract values as required
    #Grown out genes with promoters
    full_gene_name = numpy.array([
                ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'A', 'A', 'A'],
                ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ])
    #Nucleotides for the genome
    nucleotide_sequence = numpy.array(
            list("aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc")
        )
    #Indices for nucleotides
    nucleotide_index = numpy.array(range(1,100))

    #Test generic attributes which are passed into the constructor
    assert gene.name == "A"
    assert gene.feature_type == "GENE"
    assert numpy.all(gene.nucleotide_sequence == nucleotide_sequence[full_gene_name[0] == "A"])
    assert numpy.all(gene.nucleotide_index == nucleotide_index[full_gene_name[0] == "A"])
    assert numpy.all(gene.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1, 28))+[-6, -5, -4]))
    assert numpy.all(gene.is_cds == numpy.array([False, False, False]+[True for i in range(27)]+[False, False, False]))
    assert numpy.all(gene.is_promoter == numpy.array([True, True, True]+[False for i in range(27)]+[True, True, True]))
    assert numpy.all(gene.is_indel == numpy.array([False for i in range(33)]))
    assert numpy.all(gene.indel_length == numpy.array([0 for i in range(33)]))

    #Test attributes set during object instanciation
    assert numpy.all(gene.codon_number == numpy.array([math.ceil(i/3) for i in range(1, 28)]))
    assert numpy.all(gene.amino_acid_number == numpy.array(range(1,10)))
    assert numpy.all(gene.codons == numpy.array(["aaa", "aaa", "acc", "ccc", "ccc", "ccg", "ggg", "ggg", "ggg"]))
    assert numpy.all(gene.amino_acid_sequence == numpy.array(list("KKTPPPGGG")))

def test_instanciate_genome_dna():

    genome = gumpy.Genome("config/TEST-DNA.gbk")

    #Testing generic attributes such as name and length
    assert len(genome) == 99
    assert genome.name == "TEST_DNA"
    assert genome.id == "TEST_DNA.1"
    assert genome.description == "TEST_DNA, complete genome"
    assert not genome.is_reference
    assert genome.max_promoter_length == 100
    assert list(genome.genes.keys()) == list("ABC")
    assert genome.genes['C']['reverse_complement']

    #Testing annotations parsed correctly
    assert genome.annotations['organism'] == "TEST_DNA_ORGANISM"
    assert genome.annotations['source'] == "TEST_DNA_SOURCE"
    assert genome.annotations['references'] == [{
                                                "location": [{
                                                     "_start": 0,
                                                     "_end": 99
                                                    }],
                                                "authors": "Test,1., Test,2., Test,3.",
                                                "consrtm": "",
                                                "title": "Test title for a reference for an DNA strand",
                                                "journal": "Test journal for an DNA strand",
                                                "medline_id": "",
                                                "pubmed_id": "1",
                                                "comment": ""
                                                 }]
    #Testing sequence was parsed correctly
    assert numpy.all(
        genome.nucleotide_sequence == numpy.array(
            list("aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc")
        ))
    assert numpy.all(
        genome.nucleotide_index == numpy.array(range(1,100)))

    #Testing all arrays were setup correctly
    #Stacked arrays are changed during promoter assignment, so should test promoter assignment too
    original_stacked_gene_name = numpy.array([
        ['', '', '', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', '', '', ''],
        ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
    ])
    #Grown out genes with promoters
    full_gene_name = numpy.array([
                ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
                ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ])
    assert numpy.all(
        genome.stacked_gene_name == full_gene_name)
    assert numpy.all(
        genome.stacked_is_reverse_complement == numpy.array([[False for i in range(90)]+[True for i in range(9)], [False for i in range(99)]])
    )
    assert numpy.all(
        genome.stacked_nucleotide_index == numpy.array([[i for i in range(1, 100)], [i for i in range(1, 100)]])
    )
    assert numpy.all(
        genome.stacked_is_promoter == numpy.array([
            [True, True, True]+[False for i in range(3, 90)]+[False, False, False, False, False, False, True, True, True],
            [False for i in range(99)]
        ])
    )
    assert numpy.all(
        genome.stacked_nucleotide_number == numpy.array([
            [-3, -2, -1] + [i for i in range(1, 28)] + [0 for i in range(28, 88)] + [i for i in range(1, 7)][::-1] + [-1, -2, -3],
            [0 for i in range(27)] + [i for i in range(1, 34)] + [0 for i in range(39)]
        ])
    )
    assert numpy.all(
        genome.stacked_is_cds == (original_stacked_gene_name != "")
    )

    #Testing all genes instanciated (not correct genes)
    assert genome.contains_gene('A')
    gene=genome.build_gene('A')
    assert gene == gumpy.Gene(
        name="A",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[0] == "A"],
        index=genome.nucleotide_index[full_gene_name[0] == "A"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "A"],
        is_cds=genome.stacked_is_cds[full_gene_name == "A"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "A"],
        is_indel=genome.is_indel[full_gene_name[0] == "A"],
        indel_length=[0 for i in range(30)],
        feature_type="GENE"
    )
    assert genome.contains_gene('B')
    gene=genome.build_gene('B')
    assert gene == gumpy.Gene(
        name="B",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[1] == "B"],
        index=genome.nucleotide_index[full_gene_name[1] == "B"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "B"],
        is_cds=genome.stacked_is_cds[full_gene_name == "B"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "B"],
        is_indel=genome.is_indel[full_gene_name[1] == "B"],
        indel_length=[0 for i in range(33)],
        feature_type="GENE"
    )
    assert genome.contains_gene('C')
    gene=genome.build_gene('C')
    assert gene == gumpy.Gene(
        name="C",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[0] == "C"],
        index=genome.nucleotide_index[full_gene_name[0] == "C"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "C"],
        is_cds=genome.stacked_is_cds[full_gene_name == "C"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "C"],
        is_indel=genome.is_indel[full_gene_name[0] == "C"],
        indel_length=[0 for i in range(9)],
        feature_type="GENE",
        reverse_complement=True
    )

def test_instanciate_genes_dna():
    genome = gumpy.Genome("config/TEST-DNA.gbk")

    #If the previous checks pass, this should be equal to a freshly instanciated Gene
    assert genome.contains_gene('A')
    gene = genome.build_gene("A")

    #Ground truth values for all genes - used to extract values as required
    #Grown out genes with promoters
    full_gene_name = numpy.array([
                ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
                ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ])
    #Nucleotides for the genome
    nucleotide_sequence = numpy.array(
            list("aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc")
        )
    #Indices for nucleotides
    nucleotide_index = numpy.array(range(1,100))

    #Test generic attributes which are passed into the constructor
    assert gene.name == "A"
    assert gene.feature_type == "GENE"
    assert numpy.all(gene.nucleotide_sequence == nucleotide_sequence[full_gene_name[0] == "A"])
    assert numpy.all(gene.nucleotide_index == nucleotide_index[full_gene_name[0] == "A"])
    assert numpy.all(gene.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1, 28))))
    assert numpy.all(gene.is_cds == numpy.array([False, False, False]+[True for i in range(27)]))
    assert numpy.all(gene.is_promoter == numpy.array([True, True, True]+[False for i in range(27)]))
    assert numpy.all(gene.is_indel == numpy.array([False for i in range(30)]))
    assert numpy.all(gene.indel_length == numpy.array([0 for i in range(30)]))

    #Test attributes set during object instanciation
    assert numpy.all(gene.codon_number == numpy.array([math.ceil(i/3) for i in range(1, 28)]))
    assert numpy.all(gene.amino_acid_number == numpy.array(range(1,10)))
    assert numpy.all(gene.codons == numpy.array(["aaa", "aaa", "acc", "ccc", "ccc", "ccg", "ggg", "ggg", "ggg"]))
    assert numpy.all(gene.amino_acid_sequence == numpy.array(list("KKTPPPGGG")))

    #Checks for a reverse complement gene
    complementary_bases = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'x':'x', 'z':'z', 'o':'o', 'n':'n', 'r':'y', 'y':'r', 's':'w', 'w':'s'}
    assert genome.contains_gene('C')
    gene = genome.build_gene("C")

    #Test generic attributes which are passed into the constructor
    assert gene.name == "C"
    assert gene.feature_type == "GENE"
    assert gene.reverse_complement == True
    assert numpy.all(gene.nucleotide_sequence == [complementary_bases[n] for n in nucleotide_sequence[full_gene_name[0] == "C"]])
    assert numpy.all(gene.nucleotide_index == nucleotide_index[full_gene_name[0] == "C"][::-1])
    assert numpy.all(gene.nucleotide_number == numpy.array((list(range(1, 7))[::-1]+[-1, -2, -3])[::-1]))
    assert numpy.all(gene.is_cds == numpy.array([False, False, False]+[True for i in range(6)]))
    assert numpy.all(gene.is_promoter == numpy.array([True, True, True]+[False for i in range(6)]))
    assert numpy.all(gene.is_indel == numpy.array([False for i in range(9)]))
    assert numpy.all(gene.indel_length == numpy.array([0 for i in range(9)]))

    #Test attributes set during object instanciation
    assert numpy.all(gene.codon_number == numpy.array([math.ceil(i/3) for i in range(1, 7)]))
    assert numpy.all(gene.amino_acid_number == numpy.array(range(1,3)))
    assert numpy.all(gene.codons == numpy.array(["ggg", "ggg"]))
    assert numpy.all(gene.amino_acid_sequence == numpy.array(list("GG")))

def test_instanciate_vcf():
    vcf = gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf",
                            ignore_filter=True,format_fields_min_thresholds={'GT_CONF':0})

    #Testing some populated items for the object
    assert vcf.vcf_version == (4, 2)
    assert vcf.contig_lengths == {"TEST_DNA": 99}
    assert vcf.format_fields_metadata == {
        "COV": {
            "type":"Integer",
            "description":"Number of reads on ref and alt alleles",
            "id": 1
        },
        "GT": {
            "type": "String",
            "description": "Genotype",
            "id": 2
        },
        "DP": {
            "type": "Integer",
            "description": "total kmer depth from gramtools",
            "id": 3
        },
        "GT_CONF": {
            "type": "Float",
            "description": "Genotype confidence. Difference in log likelihood of most likely and next most likely genotype",
            "id": 4
        }
    }
    assert len(vcf.records) == 16

    #Due to the dict structure here, several asserts are required
    calls = {
        2: {'type': 'null','call': 'x','ref': 'a','pos': 0,'original_vcf_row': {'GT': (None, None),'DP': 2,'COV': (1,1),'GT_CONF': 2.05,'REF': 'a','ALTS': ('g',)}},

        4: {'type': 'null','call': 'x','ref': 'a','pos': 0,'original_vcf_row': {'GT': (None, None),'DP': 4,'COV': (1, 2, 2),'GT_CONF': 3.77,'REF': 'a','ALTS': ('g', 't')}},

        6: {'type': 'null','call': 'x','ref': 'a','pos': 0,'original_vcf_row': {'GT': (None, None),'DP': 4,'COV': (1, 1, 1,1),'GT_CONF': 2.76,'REF': 'aaa','ALTS': ('ggt', 'gta', 'ata')}},

        7: {'type': 'null','call': 'x','ref': 'a','pos': 1,'original_vcf_row': {'GT': (None, None),'DP': 4,'COV': (1, 1, 1,1),'GT_CONF': 2.76,'REF': 'aaa','ALTS': ('ggt', 'gta', 'ata')}},

        8: {'type': 'null','call': 'x','ref': 'a','pos': 2,'original_vcf_row': {'GT': (None, None),'DP': 4,'COV': (1, 1, 1,1),'GT_CONF': 2.76,'REF': 'aaa','ALTS': ('ggt', 'gta', 'ata')}},

        12: {'type': 'snp','call': 't','ref': 'c','pos': 0,'original_vcf_row': {'GT': (1, 1),'DP': 50,'COV': (0,50),'GT_CONF':200.58,'REF': 'c','ALTS': ('t',)}},

        14: {'type': 'snp','call': 'g','ref': 'c','pos': 0,'original_vcf_row': {'GT': (2, 2),'DP': 45,'COV': (0, 2,43),'GT_CONF':155.58,'REF': 'c','ALTS': ('t', 'g')}},

        16: {'type': 'snp','call': 't','ref': 'c','pos': 0,'original_vcf_row': {'GT': (1, 1),'DP': 70,'COV': (0, 68,8),'GT_CONF':300.25,'REF': 'ccc','ALTS': ('tgc', 'gtg')}},

        17: {'type': 'snp','call': 'g','ref': 'c','pos': 1,'original_vcf_row': {'GT': (1, 1),'DP': 70,'COV': (0, 68,8),'GT_CONF':300.25,'REF': 'ccc','ALTS': ('tgc', 'gtg')}},

        22: {'type': 'het','call': 'z','ref': 'g','pos': 0,'original_vcf_row': {'GT': (1, 2),'DP': 202,'COV': (1, 99, 100,2),'GT_CONF': 613.77,'REF': 'g','ALTS': ('t', 'c', 'a')}},

        24: {'type': 'het','call': 'z','ref': 'g','pos': 0,'original_vcf_row': {'GT': (0, 2),'DP': 202,'COV': (99, 1, 100,2),'GT_CONF': 613.77,'REF': 'g','ALTS': ('t', 'c', 'a')}},

        26: {'type': 'het','call': 'z','ref': 'g','pos': 0,'original_vcf_row': {'GT': (1, 2),'DP': 100,'COV': (0, 48, 50,2),'GT_CONF': 475.54,'REF': 'gg','ALTS': ('aa', 'ct', 'at')}},

        27: {'type': 'het','call': 'z','ref': 'g','pos': 1,'original_vcf_row': {'GT': (1, 2),'DP': 100,'COV': (0, 48, 50,2),'GT_CONF': 475.54,'REF': 'gg','ALTS': ('aa', 'ct', 'at')}},

        28: {'type': 'het','call': 'z','ref': 'g','pos': 0,'original_vcf_row': {'GT': (1, 3),'DP': 100,'COV': (0, 48, 2,50),'GT_CONF': 315.11,'REF': 'gg','ALTS': ('aa', 't', 'a')}},

        29: {'type': 'het','call': 'z','ref': 'g','pos': 1,'original_vcf_row': {'GT': (1, 3),'DP': 100,'COV': (0, 48, 2,50),'GT_CONF': 315.11,'REF': 'gg','ALTS': ('aa', 't', 'a')}},

        34: {'type': 'indel','call': ('ins', 'tt'),'ref': 't','pos': 1,'original_vcf_row': {'GT': (1, 1),'DP': 200,'COV': (1,199),'GT_CONF': 145.21,'REF': 't','ALTS': ('ttt',)}},

        37: {'type': 'indel','call': ('del', 't'),'ref': 't','pos': 1,'original_vcf_row': {'GT': (1, 1),'DP': 200,'COV': (1,199),'GT_CONF': 145.21,'REF': 'tt','ALTS': ('t',)}},

        40: {'type': 'indel','call': ('ins', 'g'),'ref': 't','pos': 1,'original_vcf_row': {'GT': (1, 1),'DP': 200,'COV': (1,199),'GT_CONF': 145.21,'REF': 'tt','ALTS': ('agt',)}},

        39: {'type': 'snp','call': 'a','ref': 't','pos': 0,'original_vcf_row': {'GT': (1, 1),'DP': 200,'COV': (1,199),'GT_CONF':145.21,'REF': 'tt','ALTS': ('agt',)}},

        65: {'type': 'indel','call': ('ins', 'ca'),'ref': 'g','pos': 0,'original_vcf_row': {'GT': (1, 1),'DP': 200,'COV': (1,199),'GT_CONF': 145.21,'REF': 'gg','ALTS': ('cagg',)}},

        74: {'type': 'indel','call': ('ins', 'a'),'ref': 't','pos': 1,'original_vcf_row': {'GT': (1, 1),'DP': 200,'COV': (1, 198,1), 'GT_CONF': 145.21,'REF': 't','ALTS': ('ta', 'at')}}
    }


    # could use assertDictEqual from unittest framework, but not using at present
    assert vcf.calls == calls
    #Testing record objects

    #Features common to all record objects:
    for record in vcf.records:
        assert record.chrom == "TEST_DNA"
    #Pos
    assert vcf.records[0].pos == 2
    assert vcf.records[1].pos == 4
    assert vcf.records[2].pos == 6
    assert vcf.records[3].pos == 12
    assert vcf.records[4].pos == 14
    assert vcf.records[5].pos == 16
    assert vcf.records[6].pos == 22
    assert vcf.records[7].pos == 24
    assert vcf.records[8].pos == 26
    assert vcf.records[9].pos == 28
    assert vcf.records[10].pos == 33
    assert vcf.records[11].pos == 36
    assert vcf.records[12].pos == 39
    assert vcf.records[13].pos == 65
    assert vcf.records[14].pos == 69
    assert vcf.records[15].pos == 73
    #Ref
    assert vcf.records[0].ref == "a"
    assert vcf.records[1].ref == "a"
    assert vcf.records[2].ref == "aaa"
    assert vcf.records[3].ref == "c"
    assert vcf.records[4].ref == "c"
    assert vcf.records[5].ref == "ccc"
    assert vcf.records[6].ref == "g"
    assert vcf.records[7].ref == "g"
    assert vcf.records[8].ref == "gg"
    assert vcf.records[9].ref == "gg"
    assert vcf.records[10].ref == "t"
    assert vcf.records[11].ref == "tt"
    assert vcf.records[12].ref == "tt"
    assert vcf.records[13].ref == "gg"
    assert vcf.records[14].ref == "gg"
    assert vcf.records[15].ref == "t"

    #Alt
    assert numpy.all(vcf.records[0].alts == ("g", ))
    assert numpy.all(vcf.records[1].alts == ("g",'t'))
    assert numpy.all(vcf.records[2].alts == ("ggt", "gta",'ata'))
    assert numpy.all(vcf.records[3].alts == ("t", ))
    assert numpy.all(vcf.records[4].alts == ("t", "g"))
    assert numpy.all(vcf.records[5].alts == ('tgc','gtg'))
    assert numpy.all(vcf.records[6].alts == ('t','c','a'))
    assert numpy.all(vcf.records[7].alts == ('t','c','a'))
    assert numpy.all(vcf.records[8].alts == ('aa','ct','at'))
    assert numpy.all(vcf.records[9].alts == ('aa','t','a'))
    assert numpy.all(vcf.records[10].alts == ('ttt',))
    assert numpy.all(vcf.records[11].alts == ('t',))
    assert numpy.all(vcf.records[12].alts == ('agt',))
    assert numpy.all(vcf.records[13].alts == ('cagg',))
    assert numpy.all(vcf.records[14].alts == ('gg',))
    assert numpy.all(vcf.records[15].alts == ('ta', 'at'))

    for i in range(15):
        assert vcf.records[i].qual == None
        assert vcf.records[i].info == {'KMER': 15}

    for i in range(15):
        if i == 1:
            assert vcf.records[i].filter == 'MIN_GT'
        else:
            assert vcf.records[i].filter == 'PASS'

    #GT
    gt = [record.values["GT"] for record in vcf.records]
    #None is given as a GT value for null values for alts
    assert numpy.all(gt == [(None, None),(None, None),(None, None),(1, 1),(2, 2),(1, 1),(1, 2),(0, 2),(1, 2),(1, 3),(1, 1),(1, 1), (1, 1), (1, 1), (0, 0), (1,1)])

    #DP
    dp = [record.values["DP"] for record in vcf.records]
    assert numpy.all(dp == [2, 4, 4, 50, 45, 70, 202, 202, 100, 100, 200, 200, 200, 200, 200, 200])
    #COV
    cov = [record.values["COV"] for record in vcf.records]
    assert numpy.all(cov == [(1, 1), (1, 2, 2), (1, 1, 1, 1), (0, 50), (0, 2, 43), (0, 68, 8), (1, 99, 100, 2), (99, 1, 100, 2), (0, 48, 50, 2), (0, 48, 2, 50), (1, 199), (1, 199), (1, 199), (1, 199), (1,199), (1, 198, 1)])

    #GT_CONF
    gt_conf = [record.values["GT_CONF"] for record in vcf.records]
    assert numpy.all(gt_conf == [2.05, 3.77, 2.76, 200.58, 155.58, 300.25, 613.77, 613.77, 475.54, 315.11, 145.21,
                                145.21, 145.21, 145.21, 145.21, 145.21])

    #Quick test for VCFRecord.__repr__()
    assert vcf.records[0].__repr__() == "TEST_DNA\t2\ta\t('g',)\tNone\tPASS\tGT:DP:COV:GT_CONF\t(None, None):2:(1, 1):2.05\n"

    #Quick test for vcf.__repr__()
    rep = [line.replace("\n","") for line in str(vcf).split("\n") if line.replace("\n","").strip() != ""]
    assert rep[0] == 'VCF variant file, version 4.2'
    assert rep[1] == 'tests/test-cases/TEST-DNA.vcf'
    assert rep[2] == "16 records"
    assert rep[3] == "FORMAT columns: COV, DP, GT, GT_CONF"
    assert rep[4] == str(vcf.records[0]).strip()
    assert rep[5] == str(vcf.records[1]).strip()
    assert rep[6] == str(vcf.records[2]).strip()
    assert rep[7] == "..."
    assert rep[8] == str(vcf.records[-1]).strip()
