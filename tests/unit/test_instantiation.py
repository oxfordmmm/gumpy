import numpy, gumpy, pytest, math


#As BioPython thinks that the locus line of the TEST-RNA.gbk is malformed, it gives a warning
#So ignore it to stop failing tests...
pytestmark = pytest.mark.filterwarnings("ignore")

def test_instanciate_genome_tb():

    reference = gumpy.Genome('config/NC_000962.3.gbk',gene_subset=['katG','rpoB','pncA','Rv2042c'])

    assert len(reference)==4411532
    assert reference.name=='NC_000962'
    assert reference.id=='NC_000962.3'

    # check to see if the stacking can cope with two genes overlapping
    mask=reference.stacked_nucleotide_index==2288681
    assert set(reference.stacked_gene_name[mask]) == {'pncA','Rv2042c'}, 'stacking failing at 2288681 in H37rV v3'

    # ..and the at_index function
    set(reference.at_index(2288681)) == {'pncA','Rv2042c'}, 'not correctly detecting that pncA and Rv2042c both include 2288681 in NC_000962.3'

    assert reference.contains_gene('katG')
    assert reference.contains_gene('rpoB')
    assert reference.contains_gene('pncA')
    assert reference.contains_gene('Rv2042c')
    assert ~reference.contains_gene('rpoC')

    # check the first and last dozen bases of the whole sequence
    assert ''.join(i for i in reference.nucleotide_sequence[:12]) == 'ttgaccgatgac'
    assert ''.join(i for i in reference.nucleotide_sequence[-12:]) == 'ggagatacgtcg'

def test_instanciate_genes_tb():

    reference = gumpy.Genome('config/NC_000962.3.gbk',gene_subset=['katG','rpoB','pncA','Rv2042c','rrs'])

    truth_gene_sequence={}
    truth_gene_sequence['katG']='VPEQHPPITETTTGAASNGCPVVGHMKYPVEGGGNQDWWPNRLNLKVLHQNPAVADPMGAAFDYAAEVATIDVDALTRDIEEVMTTSQPWWPADYGHYGPLFIRMAWHAAGTYRIHDGRGGAGGGMQRFAPLNSWPDNASLDKARRLLWPVKKKYGKKLSWADLIVFAGNCALESMGFKTFGFGFGRVDQWEPDEVYWGKEATWLGDERYSGKRDLENPLAAVQMGLIYVNPEGPNGNPDPMAAAVDIRETFRRMAMNDVETAALIVGGHTFGKTHGAGPADLVGPEPEAAPLEQMGLGWKSSYGTGTGKDAITSGIEVVWTNTPTKWDNSFLEILYGYEWELTKSPAGAWQYTAKDGAGAGTIPDPFGGPGRSPTMLATDLSLRVDPIYERITRRWLEHPEELADEFAKAWYKLIHRDMGPVARYLGPLVPKQTLLWQDPVPAVSHDLVGEAEIASLKSQIRASGLTVSQLVSTAWAAASSFRGSDKRGGANGGRIRLQPQVGWEVNDPDGDLRKVIRTLEEIQESFNSAAPGNIKVSFADLVVLGGCAAIEKAAKAAGHNITVPFTPGRTDASQEQTDVESFAVLEPKADGFRNYLGKGNPLPAEYMLLDKANLLTLSAPEMTVLVGGLRVLGANYKRLPLGVFTEASESLTNDFFVNLLDMGITWEPSPADDGTYQGKDGSGKVKWTGSRVDLVFGSNSELRALVEVYGADDAQPKFVQDFVAAWDKVMNLDRFDVR!'
    truth_gene_sequence['rpoB']='LADSRQSKTAASPSPSRPQSSSNNSVPGAPNRVSFAKLREPLEVPGLLDVQTDSFEWLIGSPRWRESAAERGDVNPVGGLEEVLYELSPIEDFSGSMSLSFSDPRFDDVKAPVDECKDKDMTYAAPLFVTAEFINNNTGEIKSQTVFMGDFPMMTEKGTFIINGTERVVVSQLVRSPGVYFDETIDKSTDKTLHSVKVIPSRGAWLEFDVDKRDTVGVRIDRKRRQPVTVLLKALGWTSEQIVERFGFSEIMRSTLEKDNTVGTDEALLDIYRKLRPGEPPTKESAQTLLENLFFKEKRYDLARVGRYKVNKKLGLHVGEPITSSTLTEEDVVATIEYLVRLHEGQTTMTVPGGVEVPVETDDIDHFGNRRLRTVGELIQNQIRVGMSRMERVVRERMTTQDVEAITPQTLINIRPVVAAIKEFFGTSQLSQFMDQNNPLSGLTHKRRLSALGPGGLSRERAGLEVRDVHPSHYGRMCPIETPEGPNIGLIGSLSVYARVNPFGFIETPYRKVVDGVVSDEIVYLTADEEDRHVVAQANSPIDADGRFVEPRVLVRRKAGEVEYVPSSEVDYMDVSPRQMVSVATAMIPFLEHDDANRALMGANMQRQAVPLVRSEAPLVGTGMELRAAIDAGDVVVAEESGVIEEVSADYITVMHDNGTRRTYRMRKFARSNHGTCANQCPIVDAGDRVEAGQVIADGPCTDDGEMALGKNLLVAIMPWEGHNYEDAIILSNRLVEEDVLTSIHIEEHEIDARDTKLGAEEITRDIPNISDEVLADLDERGIVRIGAEVRDGDILVGKVTPKGETELTPEERLLRAIFGEKAREVRDTSLKVPHGESGKVIGIRVFSREDEDELPAGVNELVRVYVAQKRKISDGDKLAGRHGNKGVIGKILPVEDMPFLADGTPVDIILNTHGVPRRMNIGQILETHLGWCAHSGWKVDAAKGVPDWAARLPDELLEAQPNAIVSTPVFDGAQEAELQGLLSCTLPNRDGDVLVDADGKAMLFDGRSGEPFPYPVTVGYMYIMKLHHLVDDKIHARSTGPYSMITQQPLGGKAQFGGQRFGEMECWAMQAYGAAYTLQELLTIKSDDTVGRVKVYEAIVKGENIPEPGIPESFKVLLKELQSLCLNVEVLSSDGAAIELREGEDEDLERAAANLGINLSRNESASVEDLA!'
    truth_gene_sequence['pncA']='MRALIIVDVQNDFCEGGSLAVTGGAALARAISDYLAEAADYHHVVATKDFHIDPGDHFSGTPDYSSSWPPHCVSGTPGADFHPSLDTSAIEAVFYKGAYTGAYSGFEGVDENGTPLLNWLRQRGVDEVDVVGIATDHCVRQTAEDAVRNGLATRVLVDLTAGVSADTTVAALEEMRTASVELVCSS!'
    truth_gene_sequence['Rv2042c']='MAPPNRDELLAAVERSPQAAAAHDRAGWVGLFTGDARVEDPVGSQPQVGHEAIGRFYDTFIGPRDITFHRDLDIVSGTVVLRDLELEVAMDSAVTVFIPAFLRYDLRPVTGEWQIAALRAYWELPAMMLQFLRTGSGATRPALQLSRALLGNQGLGGTAGFLTGFRRAGRRHKKLVETFLNAASRADKSAAYHALSRTATMTLGEDELLDIVELFEQLRGASWTKVTGAGSTVAVSLASDHRRGIMFADVPWRGNRINRIRYFPA!'
    truth_gene_sequence['rrs']='ttttgtttggagagtttgatcctggctcaggacgaacgctggcggcgtgcttaacacatgcaagtcgaacggaaaggtctcttcggagatactcgagtggcgaacgggtgagtaacacgtgggtgatctgccctgcacttcgggataagcctgggaaactgggtctaataccggataggaccacgggatgcatgtcttgtggtggaaagcgctttagcggtgtgggatgagcccgcggcctatcagcttgttggtggggtgacggcctaccaaggcgacgacgggtagccggcctgagagggtgtccggccacactgggactgagatacggcccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagcgacgccgcgtgggggatgacggccttcgggttgtaaacctctttcaccatcgacgaaggtccgggttctctcggattgacggtaggtggagaagaagcaccggccaactacgtgccagcagccgcggtaatacgtagggtgcgagcgttgtccggaattactgggcgtaaagagctcgtaggtggtttgtcgcgttgttcgtgaaatctcacggcttaactgtgagcgtgcgggcgatacgggcagactagagtactgcaggggagactggaattcctggtgtagcggtggaatgcgcagatatcaggaggaacaccggtggcgaaggcgggtctctgggcagtaactgacgctgaggagcgaaagcgtggggagcgaacaggattagataccctggtagtccacgccgtaaacggtgggtactaggtgtgggtttccttccttgggatccgtgccgtagctaacgcattaagtaccccgcctggggagtacggccgcaaggctaaaactcaaaggaattgacgggggcccgcacaagcggcggagcatgtggattaattcgatgcaacgcgaagaaccttacctgggtttgacatgcacaggacgcgtctagagataggcgttcccttgtggcctgtgtgcaggtggtgcatggctgtcgtcagctcgtgtcgtgagatgttgggttaagtcccgcaacgagcgcaacccttgtctcatgttgccagcacgtaatggtggggactcgtgagagactgccggggtcaactcggaggaaggtggggatgacgtcaagtcatcatgccccttatgtccagggcttcacacatgctacaatggccggtacaaagggctgcgatgccgcgaggttaagcgaatccttaaaagccggtctcagttcggatcggggtctgcaactcgaccccgtgaagtcggagtcgctagtaatcgcagatcagcaacgctgcggtgaatacgttcccgggccttgtacacaccgcccgtcacgtcatgaaagtcggtaacacccgaagccagtggcctaaccctcgggagggagctgtcgaaggtgggatcggcgattgggacgaagtcgtaacaaggtagccgtaccggaaggtgcggctggatcacctcctttct'

    for gene_name in ['katG','rpoB','pncA','Rv2042c','rrs']:

        assert reference.genes[gene_name].name==gene_name

        if gene_name=='rrs':
            assert ~reference.genes[gene_name].codes_protein, gene_name
            assert reference.genes[gene_name].feature_type=='RNA', gene_name
        elif gene_name=='Rv2042c':
            assert reference.genes[gene_name].codes_protein, gene_name
            assert reference.genes[gene_name].feature_type=='LOCUS', gene_name
        else:
            assert reference.genes[gene_name].codes_protein, gene_name
            assert reference.genes[gene_name].feature_type=='GENE', gene_name

        if gene_name in ['katG','pncA','Rv2042c']:
            assert reference.genes[gene_name].reverse_complement, gene_name
        else:
            assert ~reference.genes[gene_name].reverse_complement, gene_name

        if gene_name=='rrs':
            nuc_sequence=reference.genes[gene_name].nucleotide_sequence[reference.genes[gene_name].nucleotide_number>0]
            assert ''.join(i for i in nuc_sequence) == truth_gene_sequence[gene_name], gene_name+' nucleotide sequence incorrect!'
        else:
            assert ''.join(i for i in reference.genes[gene_name].amino_acid_sequence) == truth_gene_sequence[gene_name], gene_name+' amino acid sequence incorrect!'

def test_instanciate_genome_covid():

    reference = gumpy.Genome('config/NC_045512.2.gbk')

    assert len(reference)==29903
    assert reference.name=='NC_045512'
    assert reference.id=='NC_045512.2'

    # check to see if the stacking can cope with two genes overlapping
    set(reference.at_index(27756)) == {'ORF7a','ORF7b'}, 'not correctly detecting that ORF7a and ORF7b both include 22756-27759 incl. in SARS-CoV_2'
    set(reference.at_index(27757)) == {'ORF7a','ORF7b'}, 'not correctly detecting that ORF7a and ORF7b both include 22756-27759 incl. in SARS-CoV_2'
    set(reference.at_index(27758)) == {'ORF7a','ORF7b'}, 'not correctly detecting that ORF7a and ORF7b both include 22756-27759 incl. in SARS-CoV_2'
    set(reference.at_index(27759)) == {'ORF7a','ORF7b'}, 'not correctly detecting that ORF7a and ORF7b both include 22756-27759 incl. in SARS-CoV_2'

    #
    assert reference.contains_gene('ORF1ab')
    assert reference.contains_gene('S')
    assert reference.contains_gene('ORF7a')
    assert ~reference.contains_gene('Rv2042c')

    # # check the first and last dozen bases of the whole sequence
    assert ''.join(i for i in reference.nucleotide_sequence[:12]) == 'attaaaggttta'
    assert ''.join(i for i in reference.nucleotide_sequence[-12:]) == 'aaaaaaaaaaaa'


def test_instanciate_genes_covid():

    reference = gumpy.Genome('config/NC_045512.2.gbk')

    assert ~reference.genes['S'].reverse_complement
    assert reference.genes['S'].name=='S'
    assert reference.genes['S'].codes_protein
    assert reference.genes['S'].feature_type=='GENE'

    assert ''.join(i for i in reference.genes['ORF7a'].amino_acid_sequence) == 'MKIILFLALITLATCELYHYQECVRGTTVLLKEPCSSGTYEGNSPFHPLADNKFALTCFSTQFAFACPDGVKHVYQLRARSVSPKLFIRQEEVQELYSPIFLIVAAIVFITLCFTLKRKTE!'
    assert ''.join(i for i in reference.genes['ORF7b'].amino_acid_sequence) == 'MIELSLIDFYLCFLAFLLFLVLIMLIIFWFSLELQDHNETCHA!'
    assert ''.join(i for i in reference.genes['S'].amino_acid_sequence) == "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT!"

def test_instanciate_genes_covid_ORF1ab():

    reference = gumpy.Genome('config/NC_045512.2.gbk')

    assert ''.join(i for i in reference.genes['ORF1ab'].amino_acid_sequence) == 'MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGGAYTRYVDNNFCGPDGYPLECIKDLLARAGKASCTLSEQLDFIDTKRGVYCCREHEHEIAWYTERSEKSYELQTPFEIKLAKKFDTFNGECPNFVFPLNSIIKTIQPRVEKKKLDGFMGRIRSVYPVASPNECNQMCLSTLMKCDHCGETSWQTGDFVKATCEFCGTENLTKEGATTCGYLPQNAVVKIYCPACHNSEVGPEHSLAEYHNESGLKTILRKGGRTIAFGGCVFSYVGCHNKCAYWVPRASANIGCNHTGVVGEGSEGLNDNLLEILQKEKVNINIVGDFKLNEEIAIILASFSASTSAFVETVKGLDYKAFKQIVESCGNFKVTKGKAKKGAWNIGEQKSILSPLYAFASEAARVVRSIFSRTLETAQNSVRVLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYITGGVVQLTSQWLTNIFGTVYEKLKPVLDWLEEKFKEGVEFLRDGWEIVKFISTCACEIVGGQIVTCAKEIKESVQTFFKLVNKFLALCADSIIIGGAKLKALNLGETFVTHSKGLYRKCVKSREETGLLMPLKAPKEIIFLEGETLPTEVLTEEVVLKTGDLQPLEQPTSEAVEAPLVGTPVCINGLMLLEIKDTEKYCALAPNMMVTNNTFTLKGGAPTKVTFGDDTVIEVQGYKSVNITFELDERIDKVLNEKCSAYTVELGTEVNEFACVVADAVIKTLQPVSELLTPLGIDLDEWSMATYYLFDESGEFKLASHMYCSFYPPDEDEEEGDCEEEEFEPSTQYEYGTEDDYQGKPLEFGATSAALQPEEEQEEDWLDDDSQQTVGQQDGSEDNQTTTIQTIVEVQPQLEMELTPVVQTIEVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLRVCVDTVRTNVYLAVFDKNLYDKLVSSFLEMKSEKQVEQKIAEIPKEEVKPFITESKPSVEQRKQDDKKIKACVEEVTTTLEETKFLTENLLLYIDINGNLHPDSATLVSDIDITFLKKDAPYIVGDVVQEGVLTAVVIPTKKAGGTTEMLAKALRKVPTDNYITTYPGQGLNGYTVEEAKTVLKKCKSAFYILPSIISNEKQEILGTVSWNLREMLAHAEETRKLMPVCVETKAIVSTIQRKYKGIKIQEGVVDYGARFYFYTSKTTVASLINTLNDLNETLVTMPLGYVTHGLNLEEAARYMRSLKVPATVSVSSPDAVTAYNGYLTSSSKTPEEHFIETISLAGSYKDWSYSGQSTQLGIEFLKRGDKSVYYTSNPTTFHLDGEVITFDNLKTLLSLREVRTIKVFTTVDNINLHTQVVDMSMTYGQQFGPTYLDGADVTKIKPHNSHEGKTFYVLPNDDTLRVEAFEYYHTTDPSFLGRYMSALNHTKKWKYPQVNGLTSIKWADNNCYLATALLTLQQIELKFNPPALQDAYYRARAGEAANFCALILAYCNKTVGELGDVRETMSYLFQHANLDSCKRVLNVVCKTCGQQQTTLKGVEAVMYMGTLSYEQFKKGVQIPCTCGKQATKYLVQQESPFVMMSAPPAQYELKHGTFTCASEYTGNYQCGHYKHITSKETLYCIDGALLTKSSEYKGPITDVFYKENSYTTTIKPVTYKLDGVVCTEIDPKLDNYYKKDNSYFTEQPIDLVPNQPYPNASFDNFKFVCDNIKFADDLNQLTGYKKPASRELKVTFFPDLNGDVVAIDYKHYTPSFKKGAKLLHKPIVWHVNNATNKATYKPNTWCIRCLWSTKPVETSNSFDVLKSEDAQGMDNLACEDLKPVSEEVVENPTIQKDVLECNVKTTEVVGDIILKPANNSLKITEEVGHTDLMAAYVDNSSLTIKKPNELSRVLGLKTLATHGLAAVNSVPWDTIANYAKPFLNKVVSTTTNIVTRCLNRVCTNYMPYFFTLLLQLCTFTRSTNSRIKASMPTTIAKNTVKSVGKFCLEASFNYLKSPNFSKLINIIIWFLLLSVCLGSLIYSTAALGVLMSNLGMPSYCTGYREGYLNSTNVTIATYCTGSIPCSVCLSGLDSLDTYPSLETIQITISSFKWDLTAFGLVAEWFLAYILFTRFFYVLGLAAIMQLFFSYFAVHFISNSWLMWLIINLVQMAPISAMVRMYIFFASFYYVWKSYVHVVDGCNSSTCMMCYKRNRATRVECTTIVNGVRRSFYVYANGGKGFCKLHNWNCVNCDTFCAGSTFISDEVARDLSLQFKRPINPTDQSSYIVDSVTVKNGSIHLYFDKAGQKTYERHSLSHFVNLDNLRANNTKGSLPINVIVFDGKSKCEESSAKSASVYYSQLMCQPILLLDQALVSDVGDSAEVAVKMFDAYVNTFSSTFNVPMEKLKTLVATAEAELAKNVSLDNVLSTFISAARQGFVDSDVETKDVVECLKLSHQSDIEVTGDSCNNYMLTYNKVENMTPRDLGACIDCSARHINAQVAKSHNIALIWNVKDFMSLSEQLRKQIRSAAKKNNLPFKLTCATTRQVVNVVTTKIALKGGKIVNNWLKQLIKVTLVFLFVAAIFYLITPVHVMSKHTDFSSEIIGYKAIDGGVTRDIASTDTCFANKHADFDTWFSQRGGSYTNDKACPLIAAVITREVGFVVPGLPGTILRTTNGDFLHFLPRVFSAVGNICYTPSKLIEYTDFATSACVLAAECTIFKDASGKPVPYCYDTNVLEGSVAYESLRPDTRYVLMDGSIIQFPNTYLEGSVRVVTTFDSEYCRHGTCERSEAGVCVSTSGRWVLNNDYYRSLPGVFCGVDAVNLLTNMFTPLIQPIGALDISASIVAGGIVAIVVTCLAYYFMRFRRAFGEYSHVVAFNTLLFLMSFTVLCLTPVYSFLPGVYSVIYLYLTFYLTNDVSFLAHIQWMVMFTPLVPFWITIAYIICISTKHFYWFFSNYLKRRVVFNGVSFSTFEEAALCTFLLNKEMYLKLRSDVLLPLTQYNRYLALYNKYKYFSGAMDTTSYREAACCHLAKALNDFSNSGSDVLYQPPQTSITSAVLQSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQSAVKRTIKGTHHWLLLTILTSLLVLVQSTQWSLFFFLYENAFLPFAMGIIAMSAFAMMFVKHKHAFLCLFLLPSLATVAYFNMVYMPASWVMRIMTWLDMVDTSLSGFKLKDCVMYASAVVLLILMTARTVYDDGARRVWTLMNVLTLVYKVYYGNALDQAISMWALIISVTSNYSGVVTTVMFLARGIVFMCVEYCPIFFITGNTLQCIMLVYCFLGYFCTCYFGLFCLLNRYFRLTLGVYDYLVSTQEFRYMNSQGLLPPKNSIDAFKLNIKLLGVGGKPCIKVATVQSKMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQGAVDINKLCEEMLDNRATLQAIASEFSSLPSYAAFATAQEAYEQAVANGDSEVVLKKLKKSLNVAKSEFDRDAAMQRKLEKMADQAMTQMYKQARSEDKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIVQLSEISMDNSPNLAWPLIVTALRANSAVKLQNNELSPVALRQMSCAAGTTQTACTDDNALAYYNTTKGGRFVLALLSDLQDLKWARFPKSDGTGTIYTELEPPCRFVTDTPKGPKVKYLYFIKGLNNLNRGMVLGSLAATVRLQAGNATEVPANSTVLSFCAFAVDAAKAYKDYLASGGQPITNCVKMLCTHTGTGQAITVTPEANMDQESFGGASCCLYCRCHIDHPNPKGFCDLKGKYVQIPTTCANDPVGFTLKNTVCTVCGMWKGYGCSCDQLREPMLQSADAQSFLNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAVAKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGIVGVLTLDNQDLNGNWYDFGDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQTYHPNCVNCLDDRCILHCANFNVLFSTVFPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSRLSFKELLVYAADPAMHAASGNLLLDKRTTCFSVAALTNNVAFQTVKPGNFNKDFYDFAVSKGFFKEGSSVELKHFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNNLDKSAGFPFNKWGKARLYYDSMSYEDQDALFAYTKRNVIPTITQMNLKYAISAKNRARTVAGVSICSTMTNRQFHQKLLKSIAATRGATVVIGTSKFYGGWHNMLKTVYSDVENPHLMGWDYPKCDRAMPNMLRIMASLVLARKHTTCCSLSHRFYRLANECAQVLSEMVMCGGSLYVKPGGTSSGDATTAYANSVFNICQAVTANVNALLSTDGNKIADKYVRNLQHRLYECLYRNRDVDTDFVNEFYAYLRKHFSMMILSDDAVVCFNSTYASQGLVASIKNFKSVLYYQNNVFMSEAKCWTETDLTKGPHEFCSQHTMLVKQGDDYVYLPYPDPSRILGAGCFVDDIVKTDGTLMIERFVSLAIDAYPLTKHPNQEYADVFHLYLQYIRKLHDELTGHMLDMYSVMLTNDNTSRYWEPEFYEAMYTPHTVLQAVGACVLCNSQTSLRCGACIRRPFLCCKCCYDHVISTSHKLVLSVNPYVCNAPGCDVTDVTQLYLGGMSYYCKSHKPPISFPLCANGQVFGLYKNTCVGSDNVTDFNAIATCDWTNAGDYILANTCTERLKLFAAETLKATEETFKLSYGIATVREVLSDRELHLSWEVGKPRPPLNRNYVFTGYRVTKNSKVQIGEYTFEKGDYGDAVVYRGTTTYKLNVGDYFVLTSHTVMPLSAPTLVPQEHYVRITGLYPTLNISDEFSSNVANYQKVGMQKYSTLQGPPGTGKSHFAIGLALYYPSARIVYTACSHAAVDALCEKALKYLPIDKCSRIIPARARVECFDKFKVNSTLEQYVFCTVNALPETTADIVVFDEISMATNYDLSVVNARLRAKHYVYIGDPAQLPAPRTLLTKGTLEPEYFNSVCRLMKTIGPDMFLGTCRRCPAEIVDTVSALVYDNKLKAHKDKSAQCFKMFYKGVITHDVSSAINRPQIGVVREFLTRNPAWRKAVFISPYNSQNAVASKILGLPTQTVDSSQGSEYDYVIFTQTTETAHSCNVNRFNVAITRAKVGILCIMSDRDLYDKLQFTSLEIPRRNVATLQAENVTGLFKDCSKVITGLHPTQAPTHLSVDTKFKTEGLCVDIPGIPKDMTYRRLISMMGFKMNYQVNGYPNMFITREEAIRHVRAWIGFDVEGCHATREAVGTNLPLQLGFSTGVNLVAVPTGYVDTPNNTDFSRVSAKPPPGDQFKHLIPLMYKGLPWNVVRIKIVQMLSDTLKNLSDRVVFVLWAHGFELTSMKYFVKIGPERTCCLCDRRATCFSTASDTYACWHHSIGFDYVYNPFMIDVQQWGFTGNLQSNHDLYCQVHGNAHVASCDAIMTRCLAVHECFVKRVDWTIEYPIIGDELKINAACRKVQHMVVKAALLADKFPVLHDIGNPKAIKCVPQADVEWKFYDAQPCSDKAYKIEELFYSYATHSDKFTDGVCLFWNCNVDRYPANSIVCRFDTRVLSNLNLPGCDGGSLYVNKHAFHTPAFDKSAFVNLKQLPFFYYSDSPCESHGKQVVSDIDYVPLKSATCITRCNLGGAVCRHHANEYRLYLDAYNMMISAGFSLWVYKQFDTYNLWNTFTRLQSLENVAFNVVNKGHFDGQQGEVPVSIINNTVYTKVDGVDVELFENKTTLPVNVAFELWAKRNIKPVPEVKILNNLGVDIAANTVIWDYKRDAPAHISTIGVCSMTDIAKKPTETICAPLTVFFDGRVDGQVDLFRNARNGVLITEGSVKGLQPSVGPKQASLNGVTLIGEAVKTQFNYYKKVDGVVQQLPETYFTQSRNLQEFKPRSQMEIDFLELAMDEFIERYKLEGYAFEHIVYGDFSHSQLGGLHLLIGLAKRFKESPFELEDFIPMDSTVKNYFITDAQTGSSKCVCSVIDLLLDDFVEIIKSQDLSVVSKVVKVTIDYTEISFMLWCKDGHVETFYPKLQSSQAWQPGVAMPNLYKMQRMLLEKCDLQNYGDSATLPKGIMMNVAKYTQLCQYLNTLTLAVPYNMRVIHFGAGSDKGVAPGTAVLRQWLPTGTLLVDSDLNDFVSDADSTLIGDCATVHTANKWDLIISDMYDPKTKNVTKENDSKEGFFTYICGFIQQKLALGGSVAIKITEHSWNADLYKLMGHFAWWTAFVTNVNASSSEAFLIGCNYLGKPREQIDGYVMHANYIFWRNTNPIQLSSYSLFDMSKFPLKLRGTAVMSLKEGQINDMILSLLSKGRLIIRENNRVVISSDVLVNN!'

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
    assert list(genome.genes_lookup.keys()) == list("ABC")

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
    assert list(genome.genes.keys()) == ["A", "B", "C"]

    assert genome.genes["A"] == gumpy.Gene(
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
    assert genome.genes["B"] == gumpy.Gene(
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

    assert genome.genes["C"] == gumpy.Gene(
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
    gene = genome.genes["A"]

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
    assert numpy.all(gene.index == nucleotide_index[full_gene_name[0] == "A"])
    assert numpy.all(gene.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1, 28))+[-6, -5, -4]))
    assert numpy.all(gene.is_cds == numpy.array([False, False, False]+[True for i in range(27)]+[False, False, False]))
    assert numpy.all(gene.is_promoter == numpy.array([True, True, True]+[False for i in range(27)]+[True, True, True]))
    assert numpy.all(gene.is_indel == numpy.array([False for i in range(33)]))
    assert numpy.all(gene.indel_length == numpy.array([0 for i in range(33)]))

    #Test attributes set during object instanciation
    assert numpy.all(gene.triplet_number == numpy.array([math.ceil(i/3) for i in range(1, 28)]))
    assert numpy.all(gene.amino_acid_number == numpy.array(range(1,10)))
    assert numpy.all(gene.codons == numpy.array(["aaa", "aaa", "acc", "ccc", "ccc", "ccg", "ggg", "ggg", "ggg"]))
    assert numpy.all(gene.amino_acid_sequence == numpy.array(list("KKTPPPGGG")))

def test_instanciate_genome_dna():

    genome = gumpy.Genome("config/TEST-DNA.gbk")

    #Testing generic attributes such as name and length
    assert len(genome) == 99
    assert genome.name == "TEST_DNA"
    assert genome.id == "TEST_DNA.1"
    assert list(genome.genes_lookup.keys()) == list("ABC")

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
    #No rev-comp as DNA only has a single strand
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
    assert list(genome.genes.keys()) == ["A", "B", "C"]

    assert genome.genes["A"] == gumpy.Gene(
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
    assert genome.genes["B"] == gumpy.Gene(
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
    assert genome.genes["C"] == gumpy.Gene(
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
    gene = genome.genes["A"]

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
    assert numpy.all(gene.index == nucleotide_index[full_gene_name[0] == "A"])
    assert numpy.all(gene.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1, 28))))
    assert numpy.all(gene.is_cds == numpy.array([False, False, False]+[True for i in range(27)]))
    assert numpy.all(gene.is_promoter == numpy.array([True, True, True]+[False for i in range(27)]))
    assert numpy.all(gene.is_indel == numpy.array([False for i in range(30)]))
    assert numpy.all(gene.indel_length == numpy.array([0 for i in range(30)]))

    #Test attributes set during object instanciation
    assert numpy.all(gene.triplet_number == numpy.array([math.ceil(i/3) for i in range(1, 28)]))
    assert numpy.all(gene.amino_acid_number == numpy.array(range(1,10)))
    assert numpy.all(gene.codons == numpy.array(["aaa", "aaa", "acc", "ccc", "ccc", "ccg", "ggg", "ggg", "ggg"]))
    assert numpy.all(gene.amino_acid_sequence == numpy.array(list("KKTPPPGGG")))

    #Checks for a reverse complement gene
    complementary_bases = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'x':'x', 'z':'z', 'o':'o', 'n':'n', 'r':'y', 'y':'r', 's':'w', 'w':'s'}
    gene = genome.genes["C"]

    #Test generic attributes which are passed into the constructor
    assert gene.name == "C"
    assert gene.feature_type == "GENE"
    assert gene.reverse_complement == True
    assert numpy.all(gene.nucleotide_sequence == [complementary_bases[n] for n in nucleotide_sequence[full_gene_name[0] == "C"]])
    assert numpy.all(gene.index == nucleotide_index[full_gene_name[0] == "C"][::-1])
    assert numpy.all(gene.nucleotide_number == numpy.array((list(range(1, 7))[::-1]+[-1, -2, -3])[::-1]))
    assert numpy.all(gene.is_cds == numpy.array([False, False, False]+[True for i in range(6)]))
    assert numpy.all(gene.is_promoter == numpy.array([True, True, True]+[False for i in range(6)]))
    assert numpy.all(gene.is_indel == numpy.array([False for i in range(9)]))
    assert numpy.all(gene.indel_length == numpy.array([0 for i in range(9)]))

    #Test attributes set during object instanciation
    assert numpy.all(gene.triplet_number == numpy.array([math.ceil(i/3) for i in range(1, 7)]))
    assert numpy.all(gene.amino_acid_number == numpy.array(range(1,3)))
    assert numpy.all(gene.codons == numpy.array(["ggg", "ggg"]))
    assert numpy.all(gene.amino_acid_sequence == numpy.array(list("GG")))

def test_instanciate_vcf():
    vcf = gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf")
    #Testing some populated items for the object
    assert vcf.VCF_VERSION == (4, 2)
    assert vcf.contig_lengths == {"TEST_DNA": 99}
    assert vcf.formats == {
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
    assert len(vcf.records) == 6

    #Due to the dict structure here, several asserts are required
    changes = {
        1: ('g', [(0, '-'), (68, 'g')]),
        15: ('t', [(0, '-'), (68, 't')]),
        27: ('z', [(1, '-'), (99, 't'), (100, 'c')]),
        71: (numpy.array(['g', 'c', 'c']), [(0, '-'), (68, numpy.array(['g', 'c', 'c']))]),
        77: ('z', [(0, '-'), (48, numpy.array(['g', 't', 't'])), (20, 'g')]),
        89: ('x', [(0, '-'), (68, 'x')])
    }
    assert vcf.changes.keys() == changes.keys()
    for key in changes.keys():
        assert numpy.all(vcf.changes[key][0] == changes[key][0])
        for ((n_reads1, call1), (n_reads2, call2)) in zip(vcf.changes[key][1], changes[key][1]):
            assert n_reads1 == n_reads2
            assert numpy.all(call1 == call2)

    #Testing record objects

    #Features common to all record objects:
    for record in vcf.records:
        assert record.chrom == "TEST_DNA"
    #Pos
    assert vcf.records[0].pos == 2
    assert vcf.records[1].pos == 16
    assert vcf.records[2].pos == 28
    assert vcf.records[3].pos == 72
    assert vcf.records[4].pos == 78
    assert vcf.records[5].pos == 90

    #Ref
    assert vcf.records[0].ref == "A"
    assert vcf.records[1].ref == "C"
    assert vcf.records[2].ref == "G"
    assert vcf.records[3].ref == "T"
    assert vcf.records[4].ref == "T"
    assert vcf.records[5].ref == "A"

    #Alt
    assert numpy.all(vcf.records[0].alts == ("G", ))
    assert numpy.all(vcf.records[1].alts == ("T", ))
    assert numpy.all(vcf.records[2].alts == ("T", "C"))
    assert numpy.all(vcf.records[3].alts == ("GCC", ))
    assert numpy.all(vcf.records[4].alts == ("GTT", "G"))
    assert numpy.all(vcf.records[5].alts == None)

    #Qual
    assert vcf.records[0].qual == None
    assert vcf.records[1].qual == None
    assert vcf.records[2].qual == None
    assert vcf.records[3].qual == None
    assert vcf.records[4].qual == None
    assert vcf.records[5].qual == None

    #Filter
    assert vcf.records[0].filter == "PASS"
    assert vcf.records[1].filter == "."
    assert vcf.records[2].filter == "PASS"
    assert vcf.records[3].filter == "PASS"
    assert vcf.records[4].filter == "."
    assert vcf.records[5].filter == "."

    #Other fields
    assert vcf.records[0].info == {'KMER': 15}
    assert vcf.records[1].info == {'KMER': 15}
    assert vcf.records[2].info == {'KMER': 15}
    assert vcf.records[3].info == {'KMER': 15}
    assert vcf.records[4].info == {'KMER': 15}
    assert vcf.records[5].info == {'KMER': 15}

    print(vcf.records[0].values.keys())
    #GT
    gt = [record.values["GT"] for record in vcf.records]
    #None is given as a GT value for null values for alts
    assert numpy.all(gt == [(1, 1), (1, 1), (1, 2), (1, 1), (1, 1), (None, None)])
    #DP
    dp = [record.values["DP"] for record in vcf.records]
    assert numpy.all(dp == [68, 68, 200, 68, 68, 68])
    #COV
    cov = [record.values["COV"] for record in vcf.records]
    assert numpy.all(cov == [(0, 68), (0, 68), (1, 99, 100), (0, 68), (0, 48, 20), (0, 68)])
    #GT_CONF
    gt_conf = [record.values["GT_CONF"] for record in vcf.records]
    assert numpy.all(gt_conf == [613.77, 613.77, 613.77, 63.77, 63.77, 63.77])

    #Quick test for VCFRecord.__repr__()
    assert vcf.records[0].__repr__() == "TEST_DNA\t2\tA\t('G',)\tNone\tPASS\tGT:DP:COV:GT_CONF\t(1, 1):68:(0, 68):613.77\n"
