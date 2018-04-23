import sys
import os
import json
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax
#


#
#gene1
seqaln1= "tiny_test_example/test.fas"
mattype1="fasta"
trfn1= "tiny_test_example/test.tre"
schema_trf1 = "newick"
workdir1="test_ods_tiny"
id_to_spn1 = r"tiny_test_example/test_nicespl.csv"
otu_jsonfi1 = "{}/otu_dict.json".format(workdir)

#gene2
seqaln2= "tiny_test_example/test.fas"
mattype2="fasta"
trfn2= "tiny_test_example/test.tre"
schema_trf2 = "newick"
workdir2="test_ods_tiny"
id_to_spn2 = r"tiny_test_example/test_nicespl.csv"
otu_jsonfi2 = "{}/otu_dict.json".format(workdir)


configfi = "tests/data/localblast.config"

treshold=2
selectby="blast"
downtorank = "species"
add_local_seq = None
id_to_spn_addseq_json = None




if os.path.exists(otu_jsonfi1):
    otu_json1 = json.load(open(otu_jsonfi1))
if os.path.exists(otu_jsonfi2):
    otu_json1 = json.load(open(otu_jsonfi2))

if os.path.isfile("{}/scrape_checkpoint.p".format(workdir1)): 
    scraper1 = pickle.load(open("{}/scrape_checkpoint.p".format(workdir1),'rb'))
if os.path.isfile("{}/scrape_checkpoint.p".format(workdir2)): 
    scraper2 = pickle.load(open("{}/scrape_checkpoint.p".format(workdir2),'rb'))



wrappers.concat(scraper1, scraper2)