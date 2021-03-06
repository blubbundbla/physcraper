#!/usr/bin/env python
import sys
from dendropy import Tree, DnaCharacterMatrix
from peyotl import gen_otu_dict, iter_node
from peyotl.manip import iter_trees, iter_otus
from peyotl.api.phylesystem_api import PhylesystemAPI
from peyotl.sugar import tree_of_life, taxonomy
from peyotl.nexson_syntax import extract_tree_nexson, get_subtree_otus, PhyloSchema
from peyotl.api import APIWrapper
api_wrapper = APIWrapper()
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO, Entrez
import re
import os
import subprocess
import time
import datetime
import glob
import configparser
import json



'''
study_id=sys.argv[1]
tree_id=sys.argv[2]
seqaln=sys.argv[3]
mattype=sys.argv[4]
runname=sys.argv[5]
'''
def otu_to_ottid(otu): #returns none if tip does not have an ottid
    label = otu.get('^ot:ottId')
    if label is None:
            o = otu.get('^ot:originalLabel', '<unknown>')
            label = "'*tip not mapped to OTT. Original label - {o}'"
            label = label.format(o=o)





#blast_loc = config['blast']['location']
#E_VALUE_THRESH = config['blast']['e_value_thresh']
#Entrez.email = config['blast']['Entrez.email']
#ott_ncbi = config['ncbi.taxonomy']['ott_ncbi']
#get_ncbi_taxonomy = config['ncbi.taxonomy']['get_ncbi_taxonomy']
#ncbi_dmp = config['ncbi.taxonomy']['ncbi_dmp']
#MISSINGNESS_THRESH = 0.5



class physcraper_setup:
    """This needs to vet the inputs, standardize names, and prep for first run through of blasting. Does not blast itself!"""
    _found_mrca = 0
    _config_read = 0
    _id_dicts = 0
    _study_get = 0
    _phylesystem = 0
    def __init__(self, study_id, tree_id, seqaln, mattype, runname, configfi='/home/ejmctavish/projects/otapi/physcraper/config', run=1):
        """initialized object, most attributes generated through self._checkArgs using config file."""
        self.configfi = configfi
        self.study_id = study_id
        self.tree_id = tree_id
        self.seqaln = seqaln
        self.mattype = mattype
        self.runname = runname
    def _read_config(self):
        _config_read=1
        self.config = configparser.ConfigParser()
        self.config.read(self.configfi)
        if not (os.path.isfile("last_update") and os.path.isfile("{}_stream.tre".format(runname)) and os.path.isfile("{}_aln_ott.fas".format(runname))): 
            self.firstrun = 1
        else: #TODO move this shiz to make's problem
            self.firstrun = 0
    def _make_id_dicts(self):
        if not self._config_read:
            self._read_config()
        ott_to_ncbi = {}
        ncbi_to_ott = {}
        fi =open(ott_ncbi)
        for lin in fi:
            lii= lin.split(",")
            ott_to_ncbi[int(lii[0])]=int(lii[1])
            ncbi_to_ott[int(lii[1])]=int(lii[0])
        fi.close()
        self._id_dicts = 1
    def _get_study(self):
        if not self._phylesystem:
            self._phylesystem_setup()
        self.nexson = self.phy.return_study(study_id)[0]
        self._study_get = 1
    def _get_mrca(self):
        if not _study_get:
            _get_study()
        ott_ids = get_subtree_otus(self.nexson, tree_id=self.tree_id) #TODO are these ottids or OTUs?
        self.mrca_node = tree_of_life.mrca(ott_ids=ottids, wrap_response=True)
        sys.stdout.write("mrca_node found, {}\n".format(mrca_node.nearest_taxon.ott_id))
        self._found_mrca = 1
    def _phylesystem_setup(self):
        phylesystem_loc = self.config['phylesystem']['location']
        phylesystem_api_wrapper = PhylesystemAPI(get_from=phylesystem_loc)
        self.phy = phylesystem_api_wrapper.phylesystem_obj
        self._phylesystem = 1
    def _reconcile_names(self):
        d = DnaCharacterMatrix.get(path=self.seqaln, schema=self.mattype)
        d.taxon_namespace.is_mutable = True
        "so here I need to be getting the original names off of the "





'''
        tre = Tree.get(data=newick,
                        schema="newick",
                        taxon_namespace=d.taxon_namespace)

    # get all of the taxa associated with tips of the tree, and make sure that
    #   they include all of the members of the data's taxon_namespace...
        treed_taxa = [i.taxon for i in tre.leaf_nodes()]
        if len(treed_taxa) != len(d.taxon_namespace):
            missing = [i.label for i in d.taxon_namespace if i not in treed_taxa]
            emf = 'Some of the taxa in the alignment are not in the tree. Missing "{}"\n'
            em = emf.format('", "'.join(missing))
            raise ValueError(em)



    #Get original label to ott id mapping to match alignement to tree 
    ps = PhyloSchema('nexson', content='otumap', version='1.2.1')
    map_dict = map_dict=ps.convert(n)
    mapped_ids = set()
    for taxon in d.taxon_namespace:
        if taxon.label.replace("_"," ") in map_dict:
            if  map_dict[taxon.label.replace("_"," ")] == {}:
                taxon.label = taxon.label.replace("/","_") # it's legal nexus, but RAxML chokes. Probably other chars this is true of as well...
            else:
                if map_dict[taxon.label.replace("_"," ")]['^ot:ottId'] not in mapped_ids: #Can't have two tips with same name. Need better alternative tho!
                    mapped_ids.add(map_dict[taxon.label.replace("_"," ")]['^ot:ottId'])
                    taxon.label = str(map_dict[taxon.label.replace("_"," ")]['^ot:ottId'])
                else:
                    taxon.label = taxon.label.replace("/","_")
        else:
            sys.sterr.write("taxon label problem")

#This whole section is because de-concatenating the data leaves some taxa wholly missing (should move to preprocessing?)
    prune = []
    dp = {}
    for taxon, seq in d.items():
        if len(seq.symbols_as_string().translate(None, "-?")) == 0:
            prune.append(taxon.label)
        else:
            dp[taxon.label] = seq
  
    dna_prune = DnaCharacterMatrix.from_dict(dp)
    tre.prune_taxa_with_labels(prune)

    fi = open("pruned_taxa",'w')
    fi.write("taxa pruned from tree and alignemnt due to excessive missing data\n")
    for tax in prune:
        fi.write("\n".format(tax))
    fi.close()

    tre.resolve_polytomies()
    tre.write(path = "{}_random_resolve.tre".format(runname), schema = "newick", unquoted_underscores=True, suppress_edge_lengths=True)
    
    dna_prune.write(path="{}_aln_ott.phy".format(runname), schema="phylip")
    dna_prune.write(path="{}_aln_ott.fas".format(runname), schema="fasta")


    #This section grabs the MRCA node and blasts for seqs that are desc from that node
    mrca_node = tree_of_life.mrca(ott_ids=ottids, wrap_response=True)
    sys.stdout.write("mrca_node found, {}\n".format(mrca_node.nearest_taxon.ott_id))
    fi=open("last_update","w")
    today = datetime.date.today()
    fi.write("{}\n".format(str(today).replace("-","/")))
    fi.close()
#Below here only get runs on later iterations, but doesn't yet account for changes to the db, does full search again.
else:
    sys.stdout.write("updating tree\n")
    d = DnaCharacterMatrix.get(path="{}_aln_ott.fas".format(runname),
                               schema="fasta")
    d.taxon_namespace.is_mutable = False
    tre = Tree.get(path="{}_stream.tre".format(runname),
                    schema="newick",
                    taxon_namespace=d.taxon_namespace,
                    preserve_underscores =True)
    d.write(path="{}_aln_ott.phy".format(runname), schema="phylip")

    # get all of the taxa associated with tips of the tree, and make sure that
    #   they include all of the members of the data's taxon_namespace.
    # This shouldn't fail as we repeat through the loop!
    outgroup = open("outgroup.txt")
    outnames=[lin.split(',')[2] for lin in outgroup]
    for taxon in d.taxon_namespace:
        try:
            if int(taxon.label) in ott_to_ncbi.keys():
                if taxon.label not in outgroup:
                    ottids.append(taxon.label)
            else:
                pass
        except:
            pass
    mrca_node = tree_of_life.mrca(ott_ids=ottids, wrap_response=True)
    sys.stdout.write("mrca_node found, {}\n".format(mrca_node.nearest_taxon.ott_id))

assert(len(ottids) > 2)



if firstrun:
    equery = "txid{}[orgn]".format(ott_to_ncbi[mrca_node.nearest_taxon.ott_id])
    for i, record in enumerate(SeqIO.parse(seqaln, mattype)):
        record.seq._data = record.seq._data.replace("-","") # blast gets upset about too many gaps from aligned file
        record.seq._data = record.seq._data.replace("?","")
        sys.stdout.write("blasting seq {}\n".format(i))
        if len(record.seq._data) > 10:
            if not os.path.isfile("{}_{}.xml".format(runname,record.name)): 
                result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"),  entrez_query=equery)
                save_file = open("{}_{}.xml".format(runname,record.name), "w")
                save_file.write(result_handle.read())
                save_file.close()
                result_handle.close()
else:
    lastupdate = subprocess.check_output(['tail', '-1', 'last_update']).strip()
    today = datetime.date.today()
    fi=open("last_update","a")
    fi.write("{}\n".format(str(today).replace("-","/")))
    fi.close()
    equery = "txid{}[orgn] AND {}:{}[mdat]".format(ott_to_ncbi[mrca_node.nearest_taxon.ott_id], lastupdate, str(today).replace("-","/"))
    sys.stdout.write("searching with limit {}\n".format(equery))
    for i, record in enumerate(SeqIO.parse(seqaln, mattype)):
        record.seq._data = record.seq._data.replace("-","") # blast gets upset about too many gaps from aligned file
        sys.stdout.write("blasting seq {}\n".format(i))
        if not os.path.isfile("{}_{}_{}.xml".format(runname,record.name,today)): 
            result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"),  entrez_query=equery)
            save_file = open("{}_{}_{}.xml".format(runname,record.name,today), "w")
            save_file.write(result_handle.read())
            save_file.close()
            result_handle.close()
    


#else: #This doesn't work BC is not taxon limited...
#    for i, record in enumerate(SeqIO.parse(seqaln, mattype)):
#        output_handle = open("query.fasta", "w")
#        SeqIO.write(record, output_handle, "fasta")
#        output_handle.close()
#        blastx_cline =  NcbiblastxCommandline(cmd='blastn', out="{}_{}.xml".format(runname,i), outfmt=5, query="query.fasta", db='nt', evalue=E_VALUE_THRESH)
#        stdout, stderr = blastx_cline()

gi_ncbi_map = {}
if os.path.isfile("id_map.txt"):
    fi = open("id_map.txt")
    for lin in fi:
        gi_ncbi_map[int(lin.split(",")[0])]=lin.split(",")[1]

new_seqs={}

if firstrun:
  for i, record in enumerate(SeqIO.parse(seqaln, mattype)):
    if os.path.isfile("{}_{}.xml".format(runname,record.name)):
        result_handle = open("{}_{}.xml".format(runname,record.name))
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                       new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct

else:
  for i, record in enumerate(SeqIO.parse(seqaln, mattype)):
    result_handle = open("{}_{}_{}.xml".format(runname,record.name,today))
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH and int(alignment.title.split('|')[1]) not in gi_ncbi_map:
                   new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct

##SET MINIMUM SEQUENCE LENGTH
if len(new_seqs)==0:
    sys.stdout.write("No new sequences found.\n")
    sys.exit()



mapped_taxon_ids=open("id_map.txt","a")
sys.stdout.write("finding taxon ids\n")
for gi in new_seqs.keys():
    sys.stdout.write(".")
    if gi in gi_ncbi_map:
        sys.stdout.write("*")
        tax_id = int(gi_ncbi_map[gi])
    else:
        tax_id = int(subprocess.check_output(["bash", get_ncbi_taxonomy, "{}".format(gi), "{}".format(ncbi_dmp)]).split('\t')[1])
        mapped_taxon_ids.write("{}, {}\n".format(gi, tax_id))
        gi_ncbi_map[gi] = tax_id
sys.stdout.write("\n")

mapped_taxon_ids.close()

newseqs_file = "{}.fasta".format(runname)
fi = open(newseqs_file,'w')
sys.stdout.write("writing out sequences\n")
for gi in new_seqs.keys():
        try:
            ott_id = ncbi_to_ott[gi_ncbi_map[gi]]
        except:
            print("ncbi taxon ID {} not in ott".format(gi_ncbi_map[gi]))
            continue
        if ott_id not in ottids: # only adding seqs we don't have taxa for
                ottids.append(ott_id)
                fi.write(">{}_q\n".format(ncbi_to_ott[gi_ncbi_map[gi]]))
                fi.write("{}\n".format(new_seqs[gi]))
                print("success {}".format(ott_id))
        if ott_id in ottids: 
                print("ncbi taxon ID {} already in tree".format(gi_ncbi_map[gi]))



fi.close()

for fl in glob.glob("papara_*"):
    os.remove(fl)

#parallelelize across more threads!!
p1 = subprocess.call(["papara", "-t","{}_random_resolve.tre".format(runname), "-s", "{}_aln_ott.phy".format(runname), "-q",  newseqs, "-n", "extended"]) 



newaln = DnaCharacterMatrix.get(path="papara_alignment.extended",schema="phylip")

#prune out identical sequences
d = {}

starts = []
stops = []
for taxon, seq in newaln.items():
    if not taxon.label[-2:] == "_q":
        seqstr = seq.symbols_as_string()
        starts.append(min(seqstr.find('A'), seqstr.find('C'),seqstr.find('G'),seqstr.find('T')))
        stops.append(min(seqstr.rfind('A'), seqstr.rfind('C'),seqstr.rfind('G'),seqstr.rfind('T')))


start = sum(starts)/len(starts)
stop = sum(stops)/len(stops)

exclude = []
d = {}
for taxon, seq in newaln.items():
    if len(seq.symbols_as_string().translate(None, "-?")) < (stop-start)*MISSINGNESS_THRESH:
        d[taxon.label] = seq.values()[start:stop]
    else:
        exclude.append(taxon.label)

  
dna_cut = DnaCharacterMatrix.from_dict(d)
tre.prune_taxa_with_labels(exclude)

tre.write(path = "{}_cut.tre".format(runname), schema = "newick", unquoted_underscores=True, suppress_edge_lengths=True)

dna_cut.write(path="{}_aln_ott_cut.phy".format(runname), schema="phylip")
dna_cut.write(path="{}_aln_ott_cut.fas".format(runname), schema="fasta")


#realin the parts that are left
p1a = subprocess.call(["papara", "-t","{}_cut.tre".format(runname), "-s", "{}_aln_ott_cut.phy".format(runname), "-q",  newseqs, "-n", "cut"]) 



#run RAXML EPA on the alignments
for fl in glob.glob("RAxML_*"):
    os.remove(fl)
#placement

p2 = subprocess.call(["raxmlHPC", "-m", "GTRCAT", "-f", "v", "-s", "papara_alignment.cut", "-t","{}_cut.tre".format(runname), "-n", "{}_PLACE".format(runname)])

# this next line is on the assumption that you have ended up with some identical sequences. They get randomly pruned I think by raxml.
#p3 = subprocess.call(["raxmlHPC", "-m", "GTRCAT", "-f", "v", "-s", "papara_alignment.extended.reduced", "-t","{}_random_resolve.tre".format(runname), "-n", "{}_PLACE".format(runname)]) 


placetre = Tree.get(path="RAxML_labelledTree.{}_PLACE".format(runname),
                schema="newick",
                preserve_underscores=True)
placetre.resolve_polytomies()

for taxon in placetre.taxon_namespace:
    if taxon.label.startswith("QUERY"):
        taxon.label=taxon.label.replace("QUERY___","")


placetre.write(path = "{}_place_resolve.tre".format(runname), schema = "newick", unquoted_underscores=True)


#Full run with starting tree from placements

subprocess.call(["raxmlHPC", "-m", "GTRCAT", "-s", "papara_alignment.cut", "-t","{}_place_resolve.tre".format(runname), "-p", "1", "-n", "{}".format(runname)])


e = DnaCharacterMatrix.get(path="papara_alignment.extended",
                           schema="phylip")

e.taxon_namespace.is_mutable = False # Names should be the same at this point.

newtre = Tree.get(path="RAxML_bestTree.{}".format(runname),
                schema="newick",
                taxon_namespace=e.taxon_namespace,
                preserve_underscores =True)


newtre.write(path = "{}_stream.tre".format(runname), schema = "newick", unquoted_underscores=True)
e.write(path="{}_aln_ott.fas".format(runname), schema="fasta")

if tnrs_wrapper is None:
    from peyotl.sugar import tnrs
    tnrs_wrapper = tnrs

from peyotl.sugar import taxonomy

for taxon in newtre.taxon_namespace:
    if taxon.label.split("_")[0] in ottids:
        info = taxonomy.taxon(taxon.label.split("_")[0],
                                  include_lineage=False,
                                  list_terminal_descendants=True,
                                  wrap_response=True)
        taxon.label="{}{}".format(info.name,taxon.label.split("_")[1:])
        taxon.label = taxon.label.replace(" ","_")
newtre.write(path = "{}_stream_names.tre".format(runname), schema = "newick", unquoted_underscores=True)

'''
