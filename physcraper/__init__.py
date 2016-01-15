#!/usr/bin/env python
import sys
from dendropy import Tree, DnaCharacterMatrix
from peyotl import gen_otu_dict, iter_node
from peyotl.manip import iter_trees, iter_otus
from peyotl.api.phylesystem_api import PhylesystemAPI
from peyotl.sugar import tree_of_life, taxonomy
from peyotl.nexson_syntax import extract_tree, extract_tree_nexson, get_subtree_otus, extract_otu_nexson, PhyloSchema
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
import pickle
import unicodedata
#import dill



'''
study_id=sys.argv[1]
tree_id=sys.argv[2]
seqaln=sys.argv[3]
mattype=sys.argv[4]
runname=sys.argv[5]
'''


def seq_dict_build(seq, label, seq_dict): #Sequence needs to be passed in as string.
    new_seq = seq.replace("-", "")
    for tax in seq_dict.keys():
        inc_seq = seq_dict[tax].replace("-", "")
        if len(inc_seq) > len(new_seq):
            if inc_seq.find(new_seq) != -1:
                sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax))
                return
        else:
            if new_seq.find(inc_seq) != -1:
                del seq_dict[tax]
                seq_dict[label] = seq
                sys.stdout.write("seq {} is supersequence of {}, {} added and {} removed\n".format(label, tax, label, tax))
                return
    sys.stdout.write(".")
    seq_dict[label] = seq
    return



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
        self.lastupdate = '1900/01/01'
        self.PS_otu = 1
        self._read_config()
        self._get_study()
        try:
            self._reconcile_names()
        except:
            print "name reconcillation error"
        self.workdir = runname
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
    def _read_config(self):
        _config_read = 1
        self.config = configparser.ConfigParser()
        self.config.read(self.configfi)
        self.E_VALUE_THRESH = self.config['blast']['e_value_thresh']
        self.hitlist_size = int(self.config['blast']['hitlist_size'])
        self.seq_len_perc = float(self.config['physcraper']['seq_len_perc'])
        self.get_ncbi_taxonomy = self.config['ncbi.taxonomy']['get_ncbi_taxonomy']
        self.ncbi_dmp = self.config['ncbi.taxonomy']['ncbi_dmp']
    def _make_id_dicts(self):
        if not self._config_read:
            self._read_config()
        self.ott_to_ncbi = {}
        self.ncbi_to_ott = {}
        fi = open(self.config['ncbi.taxonomy']['ott_ncbi']) #TODO need to keep updated
        for lin in fi:
            lii = lin.split(",")
            self.ott_to_ncbi[int(lii[0])] = int(lii[1])
            self.ncbi_to_ott[int(lii[1])] = int(lii[0])
        fi.close()
        self.gi_ncbi_dict = {}
        if os.path.isfile("{}/id_map.txt".format(self.workdir)): #todo config?!
            fi = open("{}/id_map.txt".format(self.workdir))
            for lin in fi:
                self.gi_ncbi_dict[int(lin.split(",")[0])] = lin.split(",")[1]
        self._id_dicts = 1
    def _get_study(self):
        if not self._phylesystem:
            self._phylesystem_setup()
        self.nexson = self.phy.return_study(self.study_id)[0]
      #  self.newick = extract_tree(self.nexson, tree_id, PhyloSchema('newick', output_nexml2json = '1.2.1', content="tree", tip_label="ot:originalLabel"))
        self._study_get = 1
    def _get_mrca(self):
        if not self._study_get:
            self._get_study()
        self._make_id_dicts()
        self.orig_ott_ids = get_subtree_otus(self.nexson, tree_id=self.tree_id, subtree_id="ingroup", return_format="ottid")
        try:
            self.orig_ott_ids.remove(None)
        except:
            pass
        mrca_node = tree_of_life.mrca(ott_ids=list(self.orig_ott_ids), wrap_response=True)
        self.mrca_ott = mrca_node.nearest_taxon.ott_id
        self.mrca_ncbi = self.ott_to_ncbi[self.mrca_ott]
        self._found_mrca = 1
    def _phylesystem_setup(self):
        phylesystem_loc = self.config['phylesystem']['location']
        phylesystem_api_wrapper = PhylesystemAPI(get_from=phylesystem_loc)
        self.phy = phylesystem_api_wrapper.phylesystem_obj
        self._phylesystem = 1
    def _reconcile_names(self): #TODO here is where should build up the dicts with the info from open tree, and then mint new PS OTU_ids. SHould they be output with OTU ids?!
        otus = get_subtree_otus(self.nexson, tree_id=self.tree_id)
        self.treed_taxa = {}
        self.otu_dict = {}
        self.orig_lab_to_otu = {}
        for otu_id in otus:
            self.otu_dict[otu_id] = extract_otu_nexson(self.nexson, otu_id)[otu_id]
            self.otu_dict[otu_id]['physcraper:status'] = "original"
            orig = self.otu_dict[otu_id].get(u'^ot:originalLabel').replace(" ", "_")
            self.orig_lab_to_otu[orig] = otu_id
            self.treed_taxa[orig] = self.otu_dict[otu_id].get(u'^ot:ottId')
        self.aln = DnaCharacterMatrix.get(path=self.seqaln, schema=self.mattype)
        for tax in self.aln.taxon_namespace:
            tax.label = tax.label.replace(" ", "_")#Forcing all spaces to underscore
        missing = [i.label for i in self.aln.taxon_namespace if i.label not in self.treed_taxa.keys()] #some intense name forcing...
        if missing:
            emf = 'Some of the taxa in the alignment are not in the tree. Missing "{}"\n'
            em = emf.format('", "'.join(missing))
            raise ValueError(em)
#        for taxon in self.aln:
#            otu_id = self.orig_lab_to_otu[taxon.label]
#            taxon.label = otu_id
    def _prune(self):
        prune = []
        dp = {}
        self.orig_seqlen = []
        for taxon, seq in self.aln.items():
            if len(seq.symbols_as_string().translate(None, "-?")) == 0:
                prune.append(taxon.label)
            else:
                dp[taxon.label] = seq
                self.orig_seqlen.append(len(seq.symbols_as_string().translate(None, "-?")))
        self.newick = extract_tree(self.nexson, self.tree_id, PhyloSchema('newick', output_nexml2json='1.2.1', content="tree", tip_label="ot:originalLabel"))
        self.newick = self.newick.replace(" ", "_") #UGHHHH FIX ME (but actually works fine for now...)
        self.orig_aln = self.aln
        self.aln = DnaCharacterMatrix.from_dict(dp)
        self.tre = Tree.get(data=self.newick,
                    schema="newick",
                    preserve_underscores=True,
                    taxon_namespace=self.aln.taxon_namespace)
        self.tre.prune_taxa_with_labels(prune)
        if prune:
            fi = open("{}/pruned_taxa".format(self.workdir), 'w')
            fi.write("taxa pruned from tree and alignment due to excessive missing data\n")
            for tax in prune:
                fi.write("\n".format(tax))
            fi.close()
        for taxon in self.aln.taxon_namespace:
            otu_id = self.orig_lab_to_otu[taxon.label]
            taxon.label = otu_id.encode('ascii')
    def __getstate__(self):
        result = self.__dict__.copy()
        del result['config']
        del result['phy']
        return result
    def setup_physcraper(self):
        self._read_config()
        self._get_mrca()
        self._reconcile_names()
        self._prune()


class physcraper_scrape():
    def __init__(self, pickle_dump): #SO ideally this is getting loaded, inclueding the alignment, from a physcrpaer setup pickle file.... But would a different approach be better?
        self.Load(pickle_dump)
        self.today = str(datetime.date.today())
        self.new_seqs = {}
        self.otu_by_gi = {}
        self.ident_removed = 0
        self.blast_subdir = "{}/blast_run_{}".format(self.workdir, self.today)
        self._write_files()
    def Load(self, pickfi): # this is dangerous beacuse unp picklinga nd loading from a pickle file result in slightly different objest (i.e. date updating)
        f = open(pickfi, 'rb')
        tmp_dict = pickle.load(f)
        f.close()
        self.__dict__.update(tmp_dict.__dict__)
    def set_date(self, date):
        self.today = date
        self.blast_subdir = "{}/blast_run_{}".format(self.workdir, self.today)
    def _write_files(self):
        #First write rich annotation json file with everything needed for later?
        self.tre.resolve_polytomies()
        self.tre.write(path="{}/{}_random_resolve.tre".format(self.workdir, self.runname), schema="newick", unquoted_underscores=True, suppress_edge_lengths=True)
        self.aln.write(path="{}/{}_aln_ott.phy".format(self.workdir, self.runname), schema="phylip")
        self.aln.write(path="{}/{}_aln_ott.fas".format(self.workdir, self.runname), schema="fasta")
#        self.prun_aln = "{}/{}_aln_ott.fas".format(self.workdir, self.runname)
#        self.prun_mattype = "fasta"
        self.tre.write(path="{}/{}_random_resolve{}.tre".format(self.workdir, self.runname, self.today), schema="newick", unquoted_underscores=True, suppress_edge_lengths=True)
        self.aln.write(path="{}/{}_aln_ott{}.fas".format(self.workdir, self.runname, self.today), schema="fasta")
        pickle.dump(self, open('{}/{}_scrape{}.p'.format(self.workdir, self.runname, self.today), 'wb'))
    def run_blast(self): #TODO Should this be happening elsewhere?
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        equery = "txid{}[orgn] AND {}:{}[mdat]".format(self.mrca_ncbi, self.lastupdate, self.today.replace("-", "/"))
        #for i, record in enumerate(SeqIO.parse(self.prun_aln, self.prun_mattype)): #TODO, why not using the dendropy alignement?!! Switch to that.
        for taxon, seq in self.aln.items():
            query = seq.symbols_as_string().replace("-", "").replace("?", "") # blast gets upset about too many gaps from aligned file
            sys.stdout.write("blasting seq {}\n".format(taxon.label))
            xml_fi = "{}".format(self.gen_xml_name(taxon))
            if not os.path.isfile(xml_fi):
                result_handle = NCBIWWW.qblast("blastn", "nt", query, entrez_query=equery, hitlist_size=self.hitlist_size)
                save_file = open(xml_fi, "w")
                save_file.write(result_handle.read())
                save_file.close()
                result_handle.close()
        self.lastupdate = self.today #TODO - how to propagate this forward?????? Pickle whole class for re-use, or write to file?!
  #      self._blast_complete = 1 #TODO DAMNNNN is this the best way?!
    def gen_xml_name(self, taxon):
        return "{}/{}_{}.xml".format(self.blast_subdir, self.runname, taxon.label, self.today)#TODO pull the repeated runname
    def read_blast(self):
        self.gi_dict = {}
        self.run_blast()
        for taxon, seq in self.aln.items():
            xml_fi = "{}".format(self.gen_xml_name(taxon))
            if os.path.isfile(xml_fi):
                result_handle = open(xml_fi)
                blast_records = NCBIXML.parse(result_handle)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.expect < self.E_VALUE_THRESH:
                                self.new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct
                                self.gi_dict[int(alignment.title.split('|')[1])] = alignment.__dict__
    def seq_in_orig(self, new_seq, label):
        for taxon, seq in self.aln.items():
            seq = seq.symbols_as_string().replace("-", "").replace("?", "")
            new_seq = new_seq.replace("-", "").replace("?", "")
            if seq.find(new_seq) != -1:
                sys.stdout.write("seq gi{} is subsequence of original seq {}, not added\n".format(label, taxon))
                return 1
        return 0
    def remove_identical_seqs(self):
        '''goes through the new seqs pulled down, and removes ones that are shorter than LENGTH_THRESH
        percent of the orig seq lengths, and chooses the longer of two that are other wise identical,
        and puts them in a dict with new name as gi_ott_id. Does not test if they are identical to ones in the original alignment....'''
        d = {}
        avg_seqlen = sum(self.orig_seqlen)/len(self.orig_seqlen)
        seq_len_cutoff = avg_seqlen*self.seq_len_perc
        for gi, seq in self.new_seqs.items():
            if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                #first check if already in the original alignment
                if not self.seq_in_orig(seq, gi):
                    otu_id = self._add_otu(gi)
                    seq_dict_build(seq, otu_id, d)
        self.new_seqs_otu_id = d # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        self._ident_removed = 1 #TODO: now there can be ones that are identical to the original included sequences...
    def _add_otu(self, gi):
        otu_id = "otuPS{}".format(self.PS_otu)
        self.PS_otu += 1
        self.otu_dict[otu_id] = {}
        self.otu_dict[otu_id]['ncbi:gi'] = gi
        self.otu_dict[otu_id]['ncbi:accession'] = self.gi_dict[gi]['accession']
        self.otu_dict[otu_id]['ncbi:title'] = self.gi_dict[gi]['title']
        self.otu_dict[otu_id]['ncbi:taxon'] = self.map_gi_ncbi(gi)
        self.otu_dict[otu_id]['physcraper:status'] = "query"
        return otu_id
    def map_gi_ncbi(self, gi):
        mapped_taxon_ids = open("{}/id_map.txt".format(self.workdir), "a")
        if gi in self.gi_ncbi_dict:
            tax_id = int(self.gi_ncbi_dict[gi])
        else:
            tax_id = int(subprocess.check_output(["bash", self.get_ncbi_taxonomy, "{}".format(gi), "{}".format(self.ncbi_dmp)]).split('\t')[1])
            mapped_taxon_ids.write("{}, {}\n".format(gi, tax_id))
            self.gi_ncbi_dict[gi] = tax_id
            assert tax_id  #if this doesn't work then the gi to taxon mapping needs to be updated - shouldhappen anyhow perhaps?!
        return tax_id
        mapped_taxon_ids.close()
    def map_gi_to_ott_id(self, gi):
        if self._id_dicts != 1:
            self._make_id_dicts()
        try:
            ott_id = self.ncbi_to_ott[self.map_gi_ncbi(gi)]
            return ott_id
        except:
            return None
            sys.stderror.write("ncbi taxon id {} has no ottID".format())
    def write_query_seqs(self):
        self.newseqs_file = "{}/{}_{}.fasta".format(self.workdir, self.runname, self.today)
        fi = open(self.newseqs_file, 'w')
        sys.stdout.write("writing out sequences\n")
        if not self._ident_removed:
            self.remove_identical_seqs()
        for otu_id in self.new_seqs_otu_id.keys():
            fi.write(">{}_q\n".format(otu_id))
            fi.write("{}\n".format(self.new_seqs_otu_id[otu_id]))
        self._query_written = 1
    def align_query_seqs(self, papara_runname="extended"):
        os.chdir(self.workdir)#Clean up dir moving
        p1 = subprocess.call(["papara", "-t", "{}_random_resolve.tre".format(self.runname), "-s", "{}_aln_ott.phy".format(self.runname), "-q", "{}_{}.fasta".format(self.runname, self.today), "-n", papara_runname]) #FIx directpry ugliness
        os.chdir('..')
    def place_query_seqs(self):
        os.chdir(self.workdir)
        p2 = subprocess.call(["raxmlHPC", "-m", "GTRCAT", "-f", "v", "-s", "papara_alignment.extended", "-t", "{}_random_resolve.tre".format(self.runname), "-n", "{}_PLACE".format(self.runname)])
        placetre = Tree.get(path="RAxML_labelledTree.{}_PLACE".format(self.runname),
                schema="newick",
                preserve_underscores=True)
        placetre.resolve_polytomies()
        for taxon in placetre.taxon_namespace:
            if taxon.label.startswith("QUERY"):
                taxon.label = taxon.label.replace("QUERY___", "")
        placetre.write(path="{}_place_resolve.tre".format(self.runname), schema="newick", unquoted_underscores=True)
        os.chdir('..')
    def est_full_tree(self):
        os.chdir(self.workdir)
        subprocess.call(["raxmlHPC", "-m", "GTRCAT", "-s", "papara_alignment.extended", "-t", "{}_place_resolve.tre".format(self.runname), "-p", "1", "-n", "{}".format(self.today)])
        os.chdir('..')
    def generate_streamed_alignment(self):
        self.read_blast()
        self.remove_identical_seqs()
        self.write_query_seqs()
        self.align_query_seqs()
        self.place_query_seqs()
        self.est_full_tree()
        self.aln = DnaCharacterMatrix.get(path="{}/papara_alignment.extended".format(self.workdir), schema="phylip")
        self.aln.taxon_namespace.is_mutable = False #This should enforce name matching throughout...
        self.tre = Tree.get(path="{}/RAxML_bestTree.{}".format(self.workdir, self.today),
                schema="newick",
                preserve_underscores=True,
                taxon_namespace=self.aln.taxon_namespace)
        self._write_files()