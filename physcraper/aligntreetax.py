#!/usr/bin/env python
"""Physcraper module"""
import sys
import re
import os
import json
import pickle
from ete2 import NCBITaxa
from Bio import Entrez
from dendropy import Tree, DnaCharacterMatrix, datamodel

_DEBUG = 1


def otu_tiplabel(aln):
    """makes a dictionary with tiplabels of aln as key and
    removes n from tiplabels with start with a number """
    otu_taxonlabel_problem = {}
    for key in aln.taxon_namespace:
        match = re.match("'n[0-9]{1,3}", str(key))
        if match:
            newname = str(key)[2:]
            newname = newname[:-1]
        else:
            newname = str(key)
            newname = newname[:-1]
            newname = newname[1:]
        otu_taxonlabel_problem[key.label] = newname
    return otu_taxonlabel_problem

class AlignTreeTax(object):
    """wrap up the key parts together, requires OTT_id, and names must already match """
    def __init__(self, newick, otu_dict, alignment, ingroup_mrca, workdir, schema=None):
        #TODO add assertions that inputs are correct type!!!
        self.otu_taxonlabel_problem = {}
        self.aln = alignment
        self.otu_taxonlabel_problem = otu_tiplabel(alignment)
        if schema == None:
            self.tre = Tree.get(data=newick,
                                schema="newick",
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
        else:
            self.tre = Tree.get(data=newick,
                                schema=schema,
                                preserve_underscores=True,
                                taxon_namespace=self.aln.taxon_namespace)
        assert self.tre.taxon_namespace is self.aln.taxon_namespace
        assert isinstance(self.aln, datamodel.charmatrixmodel.DnaCharacterMatrix)
        self.tre = Tree.get(data=newick,
                            schema="newick",
                            preserve_underscores=True,
                            taxon_namespace=self.aln.taxon_namespace)
        assert isinstance(otu_dict, dict)
        self.otu_dict = otu_dict
        self.ps_otu = 1 #iterator for new otu IDs
        self._reconcile_names()
        self.workdir = workdir #TODO - is this where the workdir should live?
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        assert int(ingroup_mrca)
        self.ott_mrca = ingroup_mrca
        self.orig_seqlen = [] #FIXME
        self.gi_dict = {}
        self.orig_aln = alignment
        self.orig_newick = newick
        self._reconciled = 0

    def _reconcile_names(self):
        """This checks that the tree "original labels" from phylsystem
        align with those found in the taxonomy. Spaces vs underscores
        kept being an issue, so all spaces are coerced to underscores throughout!"""
        treed_taxa = set()
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon)
        aln_tax = set([tax for tax in self.aln.taxon_namespace])
        prune = treed_taxa ^ aln_tax
        missing = [i.label for i in prune]
        if missing:
            errmf = 'NAME RECONCILIATION Some of the taxa in the tree are not in the alignment or vice versa and will be pruned. Missing "{}"\n'
            errm = errmf.format('", "'.join(missing))
            sys.stderr.write(errm)
        self.aln.remove_sequences(prune)
        self.tre.prune_taxa(prune)
        for tax in prune:
            try:
                self.otu_dict[tax.label]['physcraper:status'] = "deleted in name reconciliation"
            except:
                self.otu_dict[self.otu_taxonlabel_problem[tax.label]]['physcraper:status'] = "deleted in name reconciliation"
            self.aln.taxon_namespace.remove_taxon(tax)
        assert self.aln.taxon_namespace == self.tre.taxon_namespace

    def prune_short(self, min_seqlen=0):
        """Sometimes in the de-concatenating of the original alignment
        taxa with no sequence are generated.
        This gets rid of those from both the tre and the alignement. MUTATOR"""
        print("prune short")
        prune = []
        for tax, seq in self.aln.items():
            if len(seq.symbols_as_string().translate(None, "-?")) <= min_seqlen:
                prune.append(tax)
        if prune:
            self.aln.remove_sequences(prune)
            self.tre.prune_taxa(prune)
            #self.aln.taxon_namespace.remove_taxon(tax)
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in prune short step due to sequence shorter than {}\n".format(min_seqlen))
            for tax in prune:
                fi.write("{}, {}\n".format(tax.label, self.otu_dict[tax.label].get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            self.otu_dict[tax.label]['physcraper:status'] = "deleted in prune short"
            self.aln.taxon_namespace.remove_taxon(tax)
        assert self.aln.taxon_namespace == self.tre.taxon_namespace
        self.orig_seqlen = [len(self.aln[tax].symbols_as_string().replace("-", "").replace("N", "")) for tax in self.aln]
        self.reconcile()

    def reconcile(self, seq_len_perc=0.75):
        """all missing data seqs are sneaking in, but from where?!"""
        #assert self.aln.taxon_namespace == self.tre.taxon_namespace
        print("reconcile")
        prune = []
        avg_seqlen = sum(self.orig_seqlen)/len(self.orig_seqlen)
        seq_len_cutoff = avg_seqlen*seq_len_perc
        for tax, seq in self.aln.items():
            # print(tax)
            if len(seq.symbols_as_string().translate(None, "-?")) < seq_len_cutoff:
                prune.append(tax)
        if prune:
            fi = open("{}/pruned_taxa".format(self.workdir), 'a')
            fi.write("Taxa pruned from tree and alignment in reconcilation step due to sequence shorter than {}\n".format(seq_len_cutoff))
            for tax in prune:
                fi.write("{}, {}\n".format(tax.label, self.otu_dict.get(tax.label).get('^ot:originalLabel')))
            fi.close()
        for tax in prune:
            self.otu_dict[tax.label]['physcraper:status'] = "deleted in reconcile"
            self.remove_taxa_aln_tre(tax.label)
        aln_ids = set()
        for tax in self.aln:
            aln_ids.add(tax.label)
        assert aln_ids.issubset(self.otu_dict.keys())
        treed_taxa = set()
        orphaned_leafs = set()
        #assert self.aln.taxon_namespace == self.tre.taxon_namespace
        ## here leaf_nodes have taxa that were dropped before. Why do we have this anyways?
        for leaf in self.tre.leaf_nodes():
            treed_taxa.add(leaf.taxon.label)
            if leaf.taxon.label not in aln_ids:
                self.otu_dict[leaf.taxon.label]['physcraper:status'] = "deleted due to presence in tree but not aln. ?!"
                orphaned_leafs.add(leaf)
                print(self.otu_taxonlabel_problem.keys())
                self.otu_dict[leaf.taxon.label]['physcraper:status'] = "deleted due to presence in tree but not aln. ?!"
                #orphaned_leafs.add(leaf)
            # else:
            #     treed_taxa.add(leaf.taxon.label)
        assert treed_taxa.issubset(aln_ids)
       # for key in  self.otu_dict.keys():
       #      if key not in aln_ids:
       #           sys.stderr.write("{} was in otu dict but not alignment. it should be in new seqs...\n".format(key)
        self.trim()
        self._reconciled = 1

    def trim(self, taxon_missingness=0.75):
        '''cuts off ends of alignment, maintaining similar to original seq len
        Important bc other while whole chromosomes get dragged in!'''
        seqlen = len(self.aln[0])
        for tax in self.aln:
            if len(self.aln[tax]) != seqlen:
                sys.stderr.write("can't trim un-aligned inputs, moving on")
                return
        start = 0
        stop = seqlen
        cutoff = len(self.aln) *  taxon_missingness
        for i in range(seqlen):
            counts = {'?':0, '-':0}
            for tax in self.aln:
                call = self.aln[tax][i].label
                if call in ['?', '-']:
                    counts[call] += 1
            if counts['?']+counts['-'] <= cutoff: #first ok column
                start = i
                break
        for i in range(seqlen-1, 0, -1):
            counts = {'?':0, '-':0}
            for tax in self.aln:
                call = self.aln[tax][i].label
                if call in ['?', '-']:
                    counts[call] += 1
            if counts['?']+counts['-'] <= cutoff:
                stop = i
                break
        for taxon in self.aln:
            self.aln[taxon] = self.aln[taxon][start:stop]
        sys.stdout.write("trimmed alignement ends to < {} missing taxa, start {}, stop {}\n".format(taxon_missingness, start, stop))
        return

    def add_otu(self, gi, ids_obj, email):
        """generates an otu_id for new sequences and adds them into the otu_dict.
        Needs to be passed an IdDict to do the mapping"""
        # otu_id = "otuPS{}".format(self.ps_otu)
        # self.ps_otu += 1
        ## this exists somewhere more or less identical, write function!!
        try:
            ncbi_id = int(ids_obj.map_gi_ncbi(gi))
            try:
                ott = int(ids_obj.ncbi_to_ott[ncbi_id])
            except:
                ott = "OTT_{}".format(self.ps_otu)
                self.ps_otu += 1
            spn = str(ids_obj.ott_to_name[ott]).replace(" ", "_")
        except:
            ncbi = NCBITaxa()
            Entrez.email = email
            tries = 5
            for i in range(tries):
                try:
                    handle = Entrez.efetch(db="nucleotide", id=gi, retmode="xml")
                except:
                    if i < tries - 1: # i is zero indexed
                        continue
                    else:
                        raise
                break
            read_handle = Entrez.read(handle)[0]
            spn = read_handle['GBSeq_feature-table'][0]['GBFeature_quals'][0]['GBQualifier_value'].replace(" ", "_")
            print(spn)
            try:
                ncbi_id = int(Entrez.read(Entrez.esearch(db="taxonomy", term=spn, RetMax=100))['IdList'][0])
            except:
                tax_info = ncbi.get_name_translator([spn])
                ncbi_id = int(tax_info.items()[0][1][0])
            try:
                ott = int(ids_obj.ncbi_to_ott[ncbi_id])
            except:
                ott = "OTT_{}".format(self.ps_otu)
                self.ps_otu += 1
        otu_id = "{}_{}".format(spn, gi)
        self.otu_dict[otu_id] = {}
        self.otu_dict[otu_id]['^ncbi:gi'] = gi
        self.otu_dict[otu_id]['^ncbi:accession'] = self.gi_dict[gi]['accession']
        self.otu_dict[otu_id]['^ncbi:title'] = self.gi_dict[gi]['title']
        # self.otu_dict[otu_id]['^ncbi:taxon'] = ids_obj.map_gi_ncbi(gi)
        # self.otu_dict[otu_id]['^ot:ottId'] = ids_obj.ncbi_to_ott.get(ids_obj.map_gi_ncbi(gi))
        # self.otu_dict[otu_id]['^physcraper:status'] = "query"
        # self.otu_dict[otu_id]['^ot:ottTaxonName'] = ids_obj.ott_to_name.get(self.otu_dict[otu_id]['^ot:ottId'])
        self.otu_dict[otu_id]['^ncbi:taxon'] = ncbi_id
        self.otu_dict[otu_id]['^ot:ottId'] = ott
        self.otu_dict[otu_id]['^physcraper:status'] = "query"
        self.otu_dict[otu_id]['^ot:ottTaxonName'] = spn
        self.otu_dict[otu_id]['^physcraper:last_blasted'] = "1800/01/01" #1800 = never blasted; 1900 = blasted 1x, not added, this century = blasted and added
        if _DEBUG:
            sys.stderr.write("gi:{} assigned new otu: {}\n".format(gi, otu_id))
        return otu_id

    def write_papara_files(self, treefilename="random_resolve.tre", alnfilename="aln_ott.phy"):
        """Papara is finicky about trees and needs phylip, this writes out needed files for papara
        (except query sequences)"""
        #CAN I even evaulte things in the function definitions?
        self.tre.resolve_polytomies()
        self.tre.deroot()
        tmptre = self.tre.as_string(schema="newick",
                                    unquoted_underscores=True,
                                    suppress_rooting=True)
        tmptre = tmptre.replace(":0.0;", ";")#Papara is diffffffficult about root
        fi = open("{}/{}".format(self.workdir, treefilename), "w")
        fi.write(tmptre)
        fi.close()
        self.aln.write(path="{}/{}".format(self.workdir, alnfilename), schema="phylip")

    def write_files(self, treepath="physcraper.tre", treeschema="newick", alnpath="physcraper.fas", alnschema="fasta"):
        """Outputs both the streaming files and a ditechecked"""
        #First write rich annotation json file with everything needed for later?
        self.tre.write(path="{}/{}".format(self.workdir, treepath),
                       schema=treeschema,
                       unquoted_underscores=True)
        self.aln.write(path="{}/{}".format(self.workdir, alnpath),
                       schema=alnschema)

    def write_otus(self, filename, schema='table'):
        """Writes out OTU dict as json"""
        assert schema in ['table', 'json']
        with open("{}/{}".format(self.workdir, filename), 'w') as outfile:
            json.dump(self.otu_dict, outfile)

    def remove_taxa_aln_tre(self, taxon_label):
        """removes taxa from aln, tre and otu_dict"""
        #self.aln.taxon_namespace=self.tre.taxon_namespace
        #assert self.aln.taxon_namespace == self.tre.taxon_namespace
        tax = self.aln.taxon_namespace.get_taxon(taxon_label)
        #assert(self.aln.taxon_namespace==self.tre.taxon_namespace)
        if tax:
            self.aln.remove_sequences([tax])
            self.aln.taxon_namespace.remove_taxon(tax)
            self.tre.prune_taxa([tax])
            self.otu_dict[tax.label]['physcraper:status'] = "deleted"
        else:
            self.otu_dict[taxon_label]['physcraper:status'] = "deleted, but it wasn't in teh alignemnet..."
        # print(self.tre)
        # print(self.tre.taxon_namespace)

    def write_labelled(self, label, treepath="labelled.tre", alnpath="labelled.fas"):
        """output tree and alignement with human readable labels
        Jumps through abunch of hoops to make labels unique.
        NOT MEMORY EFFICIENT AT ALL"""
        # print("write labelled files")
        assert label in ['^ot:ottTaxonName', 'user:TaxonName', "^ot:originalLabel", "^ot:ottId", "^ncbi:taxon"]
        # print(label)
        tmp_newick = self.tre.as_string(schema="newick")
        tmp_tre = Tree.get(data=tmp_newick,
                           schema="newick",
                           preserve_underscores=True)
        tmp_fasta = self.aln.as_string(schema="fasta")
        tmp_aln = DnaCharacterMatrix.get(data=tmp_fasta,
                                         schema="fasta")
        new_names = set()
        for taxon in tmp_tre.taxon_namespace:
            print(taxon)
            ### here the double names of the labelled tre files are generated.
            new_label = self.otu_dict[taxon.label].get(label)
            if new_label:
                if new_label in new_names:
                    new_label = " ".join([new_label, taxon.label])
            elif self.otu_dict[taxon.label].get("user:TaxonName"):
                new_label = self.otu_dict[taxon.label].get("user:TaxonName")
                if new_label in new_names:
                    new_label = " ".join([new_label, taxon.label])
            elif self.otu_dict[taxon.label].get('^ot:ottTaxonName'):
                new_label = self.otu_dict[taxon.label].get('^ot:ottTaxonName')
                if new_label in new_names:
                    new_label = " ".join([new_label, taxon.label])
            elif self.otu_dict[taxon.label].get("^ot:originalLabel"):
                new_label = self.otu_dict[taxon.label].get("^ot:originalLabel")
                if new_label in new_names:
                    new_label = " ".join([new_label, taxon.label])
            elif self.otu_dict[taxon.label].get("^ncbi:taxon"):
                new_label = " ".join(["ncbi", str(self.otu_dict[taxon.label].get("^ncbi:taxon"))])
                if new_label in new_names:
                    new_label = " ".join([new_label, taxon.label])
            new_names.add(new_label)
            taxon.label = new_label
            print(taxon.label)
        tmp_tre.write(path="{}/{}".format(self.workdir, treepath),
                      schema="newick",
                      unquoted_underscores=True,
                      suppress_edge_lengths=False)

        tmp_aln.write(path="{}/{}".format(self.workdir, alnpath),
                      schema="fasta")

    def dump(self, filename="att_checkpoint.p"):
        """w""rites object to file"""
#        frozen = jsonpickle.encode(self)
#        with open('{}/{}'.format(self.workdir, filename), 'w') as pjson:
#            pjson.write(frozen)
        pickle.dump(self, open("{}/{}".format(self.workdir, filename), "wb"))
        #TODO... write as proper nexml?!
