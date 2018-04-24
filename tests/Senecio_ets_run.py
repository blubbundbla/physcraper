import sys
import os
import json
from physcraper import wrappers, generate_ATT_from_files, AlignTreeTax
#


#
seqaln= "Senecio_ets/senecio_ets.fasta"
mattype="fasta"
trfn= "Senecio_ets/senecio_ets.tre"
schema_trf = "newick"
workdir="Senecio_ets_out"
id_to_spn = r"Senecio_ets/nicespl.csv"


configfi = "tests/data/localblast.config"
# configfi = "tests/data/aws.config"
otu_jsonfi = "{}/otu_dict.json".format(workdir)
treshold=4
selectby="blast"
downtorank = "species"
add_local_seq = None
id_to_spn_addseq_json = None

cwd = os.getcwd()  

if os.path.exists(otu_jsonfi):
    otu_json = json.load(open(otu_jsonfi))
else:
    otu_json = wrappers.OtuJsonDict(id_to_spn, configfi)
    if not os.path.exists(workdir):
       os.mkdir(workdir)
    json.dump(otu_json, open(otu_jsonfi,"w"))

wrappers.filter_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 treshold,
                 selectby,
                downtorank,
                 otu_jsonfi,
                  add_local_seq,
                 id_to_spn_addseq_json,
                 configfi)





