from physcraper import wrappers_numTax
import os



#################################
seqaln =  "/home/blubb/Documents/gitdata/physcraper/tiny_test_example/test.fas"
trfn= "/home/blubb/Documents/gitdata/physcraper/tiny_test_example/test.tre"
id_to_spn = r"/home/blubb/Documents/gitdata/physcraper/tiny_test_example/test_nicespl.csv"
workdir="tiny4Blast"
mattype="fasta"
schema_trf = "newick"
configfi = "example.config"
cwd = os.getcwd() 
treshold=4
selectby="blast" 
downto= None
add_local_seq = None
id_to_spn_addseq_json = None

otu_json = wrappers_numTax.OtuJsonDict(id_to_spn, configfi)



wrappers_numTax.own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 treshold,
                 selectby,
                 downto,
                 otu_json,
                 add_local_seq,
                 id_to_spn_addseq_json,
                 configfi)
