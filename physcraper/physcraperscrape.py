#!/usr/bin/env python
"""Physcraper module"""
import sys
import os
import subprocess
import datetime
import glob
import json
from copy import deepcopy
from urllib2 import URLError
import pickle
from Bio.Blast import NCBIWWW, NCBIXML
import physcraper.AWSWWW as AWSWWW
from ete2 import NCBITaxa
from dendropy import Tree, \
                     DnaCharacterMatrix

import numpy
import random

_DEBUG = 1

class PhyscraperScrape(object): #TODO do I wantto be able to instantiate this in a different way?!
    #set up needed variables as nones here?!
    #TODO better enforce ordering
    """This is the class that does the perpetual updating"""
    def __init__(self, data_obj, ids_obj, config_obj):
        #todo check input types assert()
        self.workdir = data_obj.workdir
        self.logfile = "{}/logfile".format(self.workdir)
        self.data = data_obj
        self.ids = ids_obj
        self.config = config_obj
        self.new_seqs = {}
        self.new_seqs_otu_id = {}
        self.otu_by_gi = {}
        self._to_be_pruned = []
        self.mrca_ncbi = ids_obj.ott_to_ncbi[data_obj.ott_mrca]
        self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir)
        self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.newseqs_file = "tmp.fasta"
        self.date = str(datetime.date.today()) #Date of the run - may lag behind real date!
        self.repeat = 1
        self.reset_markers()

 #TODO is this the right place for this?
    def reset_markers(self):
        """method to reset internal markers"""
        self._blasted = 0
        self._blast_read = 0
        self._identical_removed = 0
        self._query_seqs_written = 0
        self._query_seqs_aligned = 0
        self._query_seqs_placed = 0
        self._reconciled = 0
        self._full_tree_est = 0

    def run_blast(self): #TODO Should this be happening elsewhere?
        """generates the blast queries and saves them to xml"""
        print("run_blast")
        if not os.path.exists(self.blast_subdir):
            os.makedirs(self.blast_subdir)
        with open(self.logfile, "a") as log:
            log.write("Blast run {} \n".format(datetime.date.today()))
        for taxon, seq in self.data.aln.items():
            otu_id = taxon.label
            #TODO temp until I fix delete
            if otu_id in self.data.otu_dict:
                last_blast = self.data.otu_dict[otu_id]['^physcraper:last_blasted']
                today = str(datetime.date.today()).replace("-", "/")
                if abs((datetime.datetime.strptime(today, "%Y/%m/%d") - datetime.datetime.strptime(last_blast, "%Y/%m/%d")).days) > 14: #TODO make configurable
                    equery = "{}:{}[mdat]".format(self.mrca_ncbi,
                                                  last_blast,
                                                  today)
                    query = seq.symbols_as_string().replace("-", "").replace("?", "")
                    xml_fi = "{}/{}.xml".format(self.blast_subdir, taxon.label)
                    if not os.path.isfile(xml_fi):
                        sys.stdout.write("blasting seq {}\n".format(taxon.label))
                        if self.config.blast_loc == 'local':
                            fi_old = open("{}/tmp.fas".format(self.blast_subdir), 'w')
                            fi_old.write(">{}\n".format(taxon.label))
                            fi_old.write("{}\n".format(query))
                            fi_old.close()
                            blastcmd = "blastn -query " + \
                                       "{}/tmp.fas".format(self.blast_subdir) + \
                                       " -db {}nt -out ".format(self.config.blastdb) + \
                                       xml_fi + \
                                       " -outfmt 5 -num_threads {}".format(self.config.num_threads) + \
                                       " -max_target_seqs  {} -max_hsps {}".format(self.config.hitlist_size, self.config.hitlist_size) #TODO query via stdin
                            # print(blastcmd)
                            os.system(blastcmd)
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                        if self.config.blast_loc == 'remote':
                            if self.config.url_base:
                                result_handle = AWSWWW.qblast("blastn",
                                                             "nt",
                                                              query,
                                                              url_base=self.config.url_base,
                                                              entrez_query=equery,
                                                              hitlist_size=self.config.hitlist_size,
                                                              num_threads=self.config.num_threads)
                            else:
                                print("use BLAST websevice")
                                result_handle = AWSWWW.qblast("blastn",
                                                              "nt",
                                                              query,
                                                              # enabling entrez query makes the blast output file to be empty
                                                              # entrez_query=equery,
                                                              hitlist_size=self.config.hitlist_size)
                                # print(result_handle.read())
                            save_file = open(xml_fi, "w")
                            save_file.write(result_handle.read())
                            save_file.close()
                            self.data.otu_dict[otu_id]['^physcraper:last_blasted'] = today
                            result_handle.close()
                       # except (ValueError, URLError): TODO what to do when NCBI down?! how to handle error
                       #     sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(otu_id))
        self._blasted = 1
        return

    def add_local_seq(self, path_to_local_seq, id_to_spn_addseq_json):
        """add sequences from local database, without checking that it fits into the alignment, to your phylogeny, then do normal blast"""
        #### to DO make local blast query if new seq fit
        print("add_local_seq")
        self.sp_d = {}
        # get list of sequences,
        localfiles = os.listdir(path_to_local_seq)
        for index, item in enumerate(localfiles):
            print(index)
            print(item)
            item = str(item)
            print(item.startswith(".~"))
            if item.startswith(".~"):
                print("in if")
                localfiles[index] = None
        print(localfiles)
        localfiles = filter(None, localfiles)
        print(localfiles)
        gi_counter = 1
        ## add assert that tests that every file is a fasta file in the folder
        for file in localfiles:
            print(file)
            filepath = "{}/{}".format(path_to_local_seq, file)
            open_file = open(filepath)
            # for line in openFile:
            content = open_file.readlines()
            content = [x.strip() for x in content]
            content = filter(None, content) # fastest
            count = 0
            gi_list = content[::2]
            seq_list = content[1::2]
            print(gi_list)
            print(seq_list)
            for i in xrange(0, len(gi_list)):
                key = gi_list[i].replace(">", "")
                count = count+1
                seq = seq_list[i]
                print(key)
                print(seq)
                self.filtered_seq[key] = seq
                gi_counter += gi_counter
                self.data.gi_dict[key] = {'accession': "000000{}".format(gi_counter), 'title': "unpublished"} ## numbers startuibg with 0000 are unpublished data
                print("id_to_spn_addseq_json")
                print(id_to_spn_addseq_json)
                self.data.otu_dict[key] = {}
                self.data.otu_dict[key]['^ncbi:gi'] = self.data.gi_dict[key]['accession']
                self.data.otu_dict[key]['^ncbi:accession'] = self.data.gi_dict[key]['accession']
                self.data.otu_dict[key]['^ncbi:title'] = self.data.gi_dict[key]['title']
                # self.otu_dict[otu_id]['^ncbi:taxon'] = ids_obj.map_gi_ncbi(gi)
                # self.otu_dict[otu_id]['^ot:ottId'] = ids_obj.ncbi_to_ott.get(ids_obj.map_gi_ncbi(gi))
                # self.otu_dict[otu_id]['^physcraper:status'] = "query"
                self.data.otu_dict[key]['^ot:ottTaxonName'] = id_to_spn_addseq_json[key]['^ot:ottTaxonName']
                self.data.otu_dict[key]['^ncbi:taxon'] = id_to_spn_addseq_json[key]['^ncbiID']
                self.data.otu_dict[key]['^ot:ottId'] = id_to_spn_addseq_json[key]['ot:ottId']
                self.data.otu_dict[key]['^physcraper:status'] = "local seq"
                self.data.otu_dict[key]['^physcraper:last_blasted'] = "1800/01/01"
                # use those seq to make a blast search, that saves one round of calculating a tree.
                print("blast newly added seq")
                query = seq.replace("-", "").replace("?", "")
                xml_fi = "{}/{}.xml".format(self.blast_subdir, key)
                if not os.path.isfile(xml_fi):
                    sys.stdout.write("blasting seq {}\n".format(key))
                    try:
                        ###try to make a list of blast queries
                        print(datetime.datetime.now())
                        result_handle = NCBIWWW.qblast("blastn", "nt",
                                                       query,
                                                       url_base="http://ec2-18-144-9-156.us-west-1.compute.amazonaws.com/cgi-bin/blast.cgi",
                                                       # entrez_query=equery,
                                                       hitlist_size=self.config.hitlist_size)
                        print(datetime.datetime.now())
                        save_file = open(xml_fi, "w")
                        save_file.write(result_handle.read())
                        save_file.close()
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = str(datetime.date.today()).replace("-", "/")
                        result_handle.close()
                    except (ValueError, URLError):
                        sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(key))
        self.new_seqs = deepcopy(self.filtered_seq)
        self.new_seqs_otu_id = deepcopy(self.filtered_seq) ### !!! key is not exactly same format as before
        print("read blast results")
        # read blast results and add to new_seqs
        for key in self.filtered_seq:
            xml_fi = "{}/{}.xml".format(self.blast_subdir, key)
            if os.path.isfile(xml_fi):
                result_handle = open(xml_fi)
                try:
                    blast_records = NCBIXML.parse(result_handle)
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                if float(hsp.expect) < float(self.config.e_value_thresh):
                                    if int(alignment.title.split('|')[1]) not in self.data.gi_dict: #skip ones we already have (does it matter if these were delted? No...)
                                        self.new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct
                                        self.data.gi_dict[int(alignment.title.split('|')[1])] = alignment.__dict__
                except (ValueError, URLError):
                    sys.stderr.write("NCBIWWW error. Carrying on, but skipped {}\n".format(key))
        #set back to empty dict
        self.sp_d.clear()
        self.filtered_seq.clear()

    def read_blast(self, blast_dir=None):
        """reads in and prcesses the blast xml files"""
        print("read blast")
        if not blast_dir:
            blast_dir = self.blast_subdir
        if not self._blasted:
            self.run_blast()
        for taxon in self.data.aln:
            # xml_fi =  str(taxon.label) + ".xml"
            # xml_fi = "{}/{}.xml".format(blast_dir, taxon.label)
            xml_fi = "{}/{}.xml".format(self.blast_subdir, taxon.label)
            if os.path.isfile(xml_fi):
                result_handle = open(xml_fi)
                try:
                    blast_records = NCBIXML.parse(result_handle)
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                if float(hsp.expect) < float(self.config.e_value_thresh):
                                    giId = int(alignment.title.split('|')[1])
                                    if giId not in self.data.gi_dict: #skip ones we already have (does it matter if these were delted? No...)
                                        self.new_seqs[giId] = hsp.sbjct
                                        self.data.gi_dict[giId] = alignment.__dict__
                except ValueError:
                    sys.stderr.write("Problem reading {}, skipping\n".format(xml_fi))
        self.date = str(datetime.date.today())
        self._blast_read = 1

    # TODO this should go back in the class and should prune the tree
    def seq_dict_build(self, seq, label, seq_dict): #Sequence needs to be passed in as string.
        """takes a sequence, a label (the otu_id) and a dictionary and adds the
        sequence to the dict only if it is not a subsequence of a
        sequence already in the dict.
        If the new sequence is a super suquence of one in the dict, it
        removes that sequence and replaces it"""
        new_seq = seq.replace("-", "")
        tax_list = deepcopy(seq_dict.keys())
        i = 0
        for tax_lab in tax_list:
            i += 1
            inc_seq = seq_dict[tax_lab].replace("-", "")
            if len(inc_seq) >= len(new_seq):
                if inc_seq.find(new_seq) != -1:
                    sys.stdout.write("seq {} is subsequence of {}, not added\n".format(label, tax_lab))
                    self.data.otu_dict[tax_lab]['physcraper:status'] = "subsequence, not added"
                    # self.data.otu_dict[tax_lab]['^physcraper:last_blasted'] = "1900/01/01"

                    return
            else:
                if new_seq.find(inc_seq) != -1:#
                    if self.data.otu_dict[tax_lab].get('^physcraper:status') == "original":
                        sys.stdout.write("seq {} is supersequence of original seq {}, both kept in alignment\n".format(label, tax_lab))
                        # self.data.otu_dict[tax_lab]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[tax_lab]['physcraper:status'] = "new seq added"
                        seq_dict[label] = seq
                        return
                    else:
                        del seq_dict[tax_lab]
                        seq_dict[label] = seq
                        self.data.remove_taxa_aln_tre(tax_lab)
                        # self.data.otu_dict[tax_lab]['^physcraper:last_blasted'] = "1900/01/01"

                        sys.stdout.write("seq {} is supersequence of {}, {} added and {} removed\n".format(label, tax_lab, label, tax_lab))
                        self.data.otu_dict[tax_lab]['physcraper:status'] = "new seq added in place of {}".format(tax_lab)
                        return
        sys.stdout.write(".")
        if i%50 == 0:
            sys.stdout.write("\n")
        seq_dict[label] = seq
        return

    def remove_identical_seqs(self):
        """goes through the new seqs pulled down, and removes ones that are
        shorter than LENGTH_THRESH percent of the orig seq lengths, and chooses
        the longer of two that are other wise identical, and puts them in a dict
        with new name as gi_ott_id.
        Does not test if they are identical to ones in the original alignment...."""
        print("remove identical seqs")
        tmp_dict = dict((taxon.label, self.data.aln[taxon].symbols_as_string()) for taxon in self.data.aln)
        old_seqs = tmp_dict.keys()
        #Adding seqs that are different, but needs to be maintained as diff than aln that the tree has been run on
        avg_seqlen = sum(self.data.orig_seqlen)/len(self.data.orig_seqlen) #HMMMMMMMM
        seq_len_cutoff = avg_seqlen*self.config.seq_len_perc

        for gi, seq in self.new_seqs.items():
            if len(seq.replace("-", "").replace("N", "")) > seq_len_cutoff:
                otu_id = self.data.add_otu(gi, self.ids, self.config.email)
                self.seq_dict_build(seq, otu_id, tmp_dict)
        for tax in old_seqs:
            try:
                del tmp_dict[tax]
            except KeyError:
                pass
        self.new_seqs_otu_id = tmp_dict # renamed new seq to their otu_ids from GI's, but all info is in self.otu_dict
        with open(self.logfile, "a") as log:
            log.write("{} new sequences added from genbank, of {} before filtering\n".format(len(self.new_seqs_otu_id), len(self.new_seqs)))
        self.data.dump()

    def dump(self, filename="scrape_checkpoint.p"):
#        frozen = jsonpickle.encode(self)
#        with open('{}/{}'.format(self.workdir, filename), 'w') as pjson:
#            pjson.write(frozen)
        pickle.dump(self, open("{}/{}".format(self.workdir, filename), "wb"))
        #TODO... write as proper nexml?!



    def sp_dict(self, downtorank=None):
        """makes dict with Species name as key and the corresponding seq information from aln and blast seq"""
        print("make sp_dict")
        self.sp_d = {}
        for key in self.data.otu_dict:
            if downtorank != None:
                ncbi = NCBITaxa()
                print("downto is not None")
                if '^ncbi:taxon' in self.data.otu_dict[key]:
                    ncbiID = self.data.otu_dict[key]['^ncbi:taxon']
                elif '^ncbiID' in self.data.otu_dict[key]:
                    ncbiID = self.data.otu_dict[key]['^ncbiID']
                elif '^ot:ottTaxonName' in self.data.otu_dict[key]:
                    tax_name = self.data.otu_dict[key]['^ot:ottTaxonName']
                    ncbiID = ncbi.get_name_translator(tax_name)
                    ## does not always seem to work
                else:
                    print("no taxon name provided! It will fail")
                # print(ncbiID)
                lineage = ncbi.get_lineage(ncbiID)
                # names = ncbi.get_taxid_translator(lineage)
                lineage2ranks = ncbi.get_rank(lineage)
                for k, v in lineage2ranks.iteritems():
                    if v == downtorank:
                        tax_id = k
                        value_d = ncbi.get_taxid_translator([tax_id])
                        # print(value_d)
                        value = value_d[int(tax_id)]
            else:
                if '^user:TaxonName' in self.data.otu_dict[key].keys():
                    value = self.data.otu_dict[key]['^user:TaxonName']
                elif '^ot:ottTaxonName' in self.data.otu_dict[key].keys():
                    value = self.data.otu_dict[key]['^ot:ottTaxonName']
                    if value == None:
                        value = self.data.otu_dict[key]['^ncbi:taxon']
            # print("populate dict with values")
            value = str(value).replace(" ", "_")
            if value in self.sp_d:
                self.sp_d[value].append(self.data.otu_dict[key])
            else:
                self.sp_d[value] = [self.data.otu_dict[key]]
            # print(self.sp_d[value])
        print(self.sp_d)
        print(something_stupid)
        return

    def write_query_seqs(self):
        """writes out the query sequence file"""
        print("write query seq")
        if not self._blast_read:
            self.read_blast()
        self.newseqs_file = "{}.fasta".format(self.date)
        fi = open("{}/{}".format(self.workdir, self.newseqs_file), 'w')
        sys.stdout.write("writing out sequences\n")
        for otu_id in self.new_seqs_otu_id.keys():
            if otu_id not in self.data.aln: #new seqs only
                fi.write(">{}\n".format(otu_id))
                fi.write("{}\n".format(self.new_seqs_otu_id[otu_id]))
        self._query_seqs_written = 1

    def align_query_seqs(self, papara_runname="extended"):
        """runs papara on the tree, the alinment and the new query sequences"""
        if not self._query_seqs_written:
            self.write_query_seqs()
        for filename in glob.glob('{}/papara*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[1]))
        sys.stdout.write("aligning query sequences \n")
        self.data.write_papara_files()
        os.chdir(self.workdir)#Clean up dir moving
        # print("trying to call papara")
        try:
            print("I call papara")
            subprocess.call(["papara",
                             "-t", "random_resolve.tre",
                             "-s", "aln_ott.phy",
                             "-q", self.newseqs_file,
                             "-n", papara_runname]) #FIx directory ugliness
            sys.stdout.write("Papara done")
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("failed running papara. Is it installed?\n")
                sys.exit()
            # handle file not found error.
            else:
            # Something else went wrong while trying to run `wget`
                raise
        os.chdir('..')
        assert os.path.exists(path="{}/papara_alignment.{}".format(self.workdir, papara_runname))
        self.data.aln = DnaCharacterMatrix.get(path="{}/papara_alignment.{}".format(self.workdir, papara_runname), schema="phylip")
        self.data.aln.taxon_namespace.is_mutable = True #Was too strict...
        sys.stdout.write("Papara done")
        lfd = "{}/logfile".format(self.workdir)
        with open(lfd, "a") as log:
            log.write("Following papara alignment, aln has {} seqs \n".format(len(self.data.aln)))
        self.data.reconcile()
        self._query_seqs_aligned = 1

    def place_query_seqs(self):
        """runs raxml on the tree, and the combined alignment including the new quesry seqs
        Just for placement, to use as starting tree."""
        if os.path.exists("RAxML_labelledTree.PLACE"):
                os.rename("RAxML_labelledTree.PLACE", "RAxML_labelledTreePLACE.tmp")
        sys.stdout.write("placing query sequences \n")
        #print(os.getcwd())
        os.chdir(self.workdir)
        try:
            subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                             "-f", "v",
                             "-s", "papara_alignment.extended",
                             "-t", "random_resolve.tre",
                             "-n", "PLACE"])
            placetre = Tree.get(path="RAxML_labelledTree.PLACE",
                                schema="newick",
                                preserve_underscores=True)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("failed running raxmlHPC. Is it installed?")
                sys.exit()
            # handle file not
        # handle file not found error.
            else:
        # Something else went wrong while trying to run `wget`
                raise
        placetre.resolve_polytomies()
        for taxon in placetre.taxon_namespace:
            if taxon.label.startswith("QUERY"):
                taxon.label = taxon.label.replace("QUERY___", "")
        placetre.write(path="place_resolve.tre", schema="newick", unquoted_underscores=True)
        os.chdir('..')
        self._query_seqs_placed = 1

    def est_full_tree(self):
        """Full raxml run from the placement tree as starting tree"""
        os.chdir(self.workdir)
        for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
            os.rename(filename, "{}/{}_tmp".format(self.workdir, filename.split("/")[1]))
        subprocess.call(["raxmlHPC", "-m", "GTRCAT",
                              "-s", "papara_alignment.extended",
                              "-t", "place_resolve.tre",
                              "-p", "1",
                              "-n", "{}".format(self.date)])
        os.chdir('..')#TODO mordir not always one down!
        self._full_tree_est = 1

    def generate_streamed_alignment(self, treshold=None):
        """runs the key steps and then replaces the tree and alignme nt with the expanded ones"""
        ## first if should not be necessary, as this is in the wrapper
        if len(self.new_seqs) > 0:
            self.remove_identical_seqs()
            self.data.write_files() #should happen before aligning in case of pruning
            if len(self.new_seqs_otu_id) > 0:#TODO rename to something more intutitive
                self.write_query_seqs()
                self.align_query_seqs()
                self.data.reconcile()
                self.place_query_seqs()
                self.est_full_tree()
                self.data.tre = Tree.get(path="{}/RAxML_bestTree.{}".format(self.workdir, self.date),
                                         schema="newick",
                                         preserve_underscores=True,
                                         taxon_namespace=self.data.aln.taxon_namespace)
                self.data.write_files()
                if os.path.exists("{}/previous_run".format(self.workdir)):
                    prev_dir = "{}/previous_run{}".format(self.workdir, self.date)
                    i = 0
                    while os.path.exists(prev_dir):
                        i += 1
                        prev_dir = "{}/previous_run{}".format(self.workdir, self.date) + str(i)
                    os.rename("{}/previous_run".format(self.workdir), prev_dir)
                os.rename(self.blast_subdir, "{}/previous_run".format(self.workdir))
                if os.path.exists("{}/last_completed_update".format(self.workdir)):
                    os.rename(self.tmpfi, "{}/last_completed_update".format(self.workdir))
                for filename in glob.glob('{}/RAxML*'.format(self.workdir)):
                    os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[1]))
                for filename in glob.glob('{}/papara*'.format(self.workdir)):
                    os.rename(filename, "{}/previous_run/{}".format(self.workdir, filename.split("/")[1]))
                os.rename("{}/{}".format(self.workdir, self.newseqs_file), "{}/previous_run/newseqs.fasta".format(self.workdir))
                try:
                    self.data.write_labelled(label='^ot:ottTaxonName')
                except:
                    self.data.write_labelled(label='user:TaxonName')
                self.data.write_otus("otu_info", schema='table')
                self.new_seqs = {} #Wipe for next run
                self.new_seqs_otu_id = {}
                self.repeat = 1
            else:
                sys.stdout.write("No new sequences after filtering.\n")
                self.repeat = 0
        else:
            sys.stdout.write("No new sequences found.\n")
            self.repeat = 0
        self.reset_markers()
        self.data.dump()
#        frozen = jsonpickle.encode(self.data)
#        pjson = open('{}/att_checkpoint.json'.format(self.workdir), 'wb')
#        pjson.write(frozen)
        json.dump(self.data.otu_dict, open('{}/otu_dict.json'.format(self.workdir), 'wb'))

###############################

class FilterBlast(PhyscraperScrape):
    def __init__(self, data_obj, ids_obj, config_obj):
        super(FilterBlast, self).__init__(data_obj, ids_obj, config_obj)
        self.workdir = data_obj.workdir
        self.logfile = "{}/logfile".format(self.workdir)
        self.data = data_obj
        self.ids = ids_obj
        self.config = config_obj
        self.new_seqs = {}
        self.new_seqs_otu_id = {}
        self.otu_by_gi = {}
        self._to_be_pruned = []
        self.mrca_ncbi = ids_obj.ott_to_ncbi[data_obj.ott_mrca]
        self.tmpfi = "{}/physcraper_run_in_progress".format(self.workdir)
        self.blast_subdir = "{}/current_blast_run".format(self.workdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.newseqs_file = "tmp.fasta"
        self.date = str(datetime.date.today()) #Date of the run - may lag behind real date!
        self.repeat = 1
        self.reset_markers()

        self.sp_d = {}
        self.sp_seq_d = {}
        self.filtered_seq = {}

    def make_new_seq_dict(self, treshold, selectby):
        '''select how many sequences per species to keep in the alignment.'''
        self.sp_seq_d = {}
        self.filtered_seq = {}
        print("in make_new_seq_dict")
        
        for key in self.sp_d:
            """loop to populate dict. key1 = sp name, key2= gi number, value = seq, key2 amount determined by treshold and already present seq"""
            print(key)
            seq_d = {}
            tres_minimizer = 0
           
            for giID in self.sp_d[key]:
                # print(giID)
                # print(giID['^physcraper:last_blasted'] )





                if giID['^physcraper:last_blasted'] != '1800/01/01':
                    if '^user:TaxonName'  in giID:
                        """generate entry for already existing sp"""
                        tres_minimizer += 1
                        userName = giID['^user:TaxonName']
                        for userName_aln, seq in self.data.aln.items():
                            if '^user:TaxonName' in self.data.otu_dict[userName_aln.label]:
                                print("initial user input seq")
                                if userName == self.data.otu_dict[userName_aln.label]['^user:TaxonName']:
                                    # fn_oldseq = key +  "_tobeblasted"
                                    seq = seq.symbols_as_string().replace("-", "")
                                    seq = seq.replace("?", "")
                                    #seq_d[userName.label] = [seq, len(seq)]
                                    seq_d[userName_aln.label] = seq
     
                    elif '^ot:ottTaxonName'  in giID:
                        """generate entry for already existing sp"""
                        tres_minimizer += 1
                        print("added in earlier round, get aln directly")
                        userName = giID['^ot:ottTaxonName']
                        for userName_aln, seq in self.data.aln.items():
                            if '^ot:ottTaxonName' in self.data.otu_dict[userName_aln.label]:
                                if userName == self.data.otu_dict[userName_aln.label]['^ot:ottTaxonName']:
                                    seq = seq.symbols_as_string().replace("-", "")
                                    seq = seq.replace("?", "")
                                    #seq_d[userName.label] = [seq, len(seq)]
                                    seq_d[userName_aln.label] = seq
                    # print(self.sp_d[key])
                else:
                    print("make list of gi which are new")


                    # still problem with gi's. 
                    # self. gi_dict has every ever blasted seq
                    # self. otu_dict only stuff which has been processed by PhyscraperScrape

                    # look for things in gi_dict if they are really new!!



                    if "^ncbi:gi" in giID: # this line should not be necessary as all new blast seq have gi
                        #maybe this line?
                        print(giID)
                        # print(self.otu_dict.keys())
                        print( self.data.gi_dict.keys())
                        if giID["^ncbi:gi"] not in self.data.gi_dict.keys(): # if gi has never been queried before...
                            gi = int(giID['^ncbi:gi'])
                            print(gi)
                            
                            # print(self.data.otu_dict)
                            # fp = False
                            # for k ,v in self.data.otu_dict.iteritems():
                            #     if '^ncbi:gi' in v:
                            #         print(v['^ncbi:gi'])
                            #         if gi == v['^ncbi:gi']:
                            #             fp = True
                            # if fp == False:

                                # print(self.data.otu_dict.keys())
                                # # if gi in self.data.otu_dict.keys():
                                # print(self.new_seqs.keys())
                            if gi in self.new_seqs.keys():
                                seq = self.new_seqs[gi]
                                seq = seq.replace("-", "")
                                seq = seq.replace("?", "")
                                print(seq)
                                seq_d[gi] = seq
                                    # if selectby=="blast":
                                    #     # print("writing: {}".format(fi_new))
                                    #     fi_new = open(fn_newseq, 'a')
                                    #     fi_new.write(">{}\n".format(gi))
                                    #     fi_new.write("{}\n".format(seq))
                                    #     fi_new.close()
                self.sp_seq_d[key] = seq_d

        print(self.sp_seq_d)
        print(self.sp_seq_d.keys())
        print(something_stupid)

        return

    def select_seq_by_local_blast(self, seq_d, blast_seq, blast_db, treshold, count):
        """select sequences according to  local blast.  only species included which have a blast score of mean plus/minus sd """
        """input for now are the sequence names, better to generate them here and just give the sp_d/seq_d?"""
        print("select_seq_by_local_blast")
        general_wd = os.getcwd()
        os.chdir(os.path.join(self.workdir, "blast"))
        blast_seq = blast_seq.replace(" ", "_")

        fn = "blast/{}_tobeblasted".format(str(blast_seq))
        #fn = str(blast_seq) +  "_tobeblasted"
        blast_db_fn = "blast/{}"format(blast_db.replace(" ", "_"))
        cmd1 = "makeblastdb -in {}_db -dbtype nucl".format(str(blast_db_fn))
        print(cmd1)
        os.system(cmd1)
        cmd2 = "blastn -query {} -db {}_db -out output_{}.xml -outfmt 5".format(str(fn), str(blast_db_fn), str(fn))

        # cmd2 = "blastn -query " + str(fn) +" -db " + str(blast_db) +"_db -out output_" + str(fn) + ".xml -outfmt 5"
        print(cmd2)
        os.system(cmd2)
        output_blast = "output_" + str(fn) + ".xml"
        print(output_blast)
        xml_file = open(output_blast)
        os.chdir(general_wd)
        blast_out = NCBIXML.parse(xml_file)
        hsp_scores = {}
        for record in blast_out:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    # filter by e-value
                    ## !!! maybe don"t do that....
                    if hsp.expect < 0.001:
                        gi = alignment.title.split(" ")[1]
                        hsp_scores[gi] = {"hsp.bits" : hsp.bits, "hsp.score" : hsp.score, "alignment.length" : alignment.length, "hsp.expect" : hsp.expect}
                    else:
                        print("sequences highly different")
                        print(hsp)
        # make values to select for blast search
        total_seq = 0
        bit_sum = 0
        bit_l = []
        for gi in hsp_scores:
            total_seq += 1
            bit_sum += hsp_scores[gi]["hsp.bits"]
            bit_l.append(hsp_scores[gi]["hsp.bits"])
        bit_sd = numpy.std(bit_l)
        mean_hsp_bits = float(bit_sum/total_seq)
        seq_w_maxlen = {}
        for gi in hsp_scores:
            if  (hsp_scores[gi]["hsp.bits"] >= mean_hsp_bits-float(bit_sd)) & (hsp_scores[gi]["hsp.bits"] <= mean_hsp_bits+float(bit_sd)):
                # print(seq_d.keys())
                seq_w_maxlen[int(gi)] = seq_d[int(gi)]
        random_seq_ofsp = {}
        if (treshold - count) <= 0:
            print("already too many samples of sp in aln, skip adding more.")
        elif len(seq_w_maxlen.keys()) == (treshold - count):
            random_seq_ofsp = seq_w_maxlen
        elif len(seq_w_maxlen.keys()) > (treshold - count):
            random_seq_ofsp = random.sample(seq_w_maxlen.items(), (treshold - count))
            random_seq_ofsp = dict(random_seq_ofsp)
        elif len(seq_w_maxlen.keys()) < (treshold - count):
            random_seq_ofsp = seq_w_maxlen
        if len(random_seq_ofsp) > 0:
            for gi_n in random_seq_ofsp.keys():
                self.filtered_seq[gi_n] = random_seq_ofsp[gi_n]
        return(self.filtered_seq)

    def select_seq_by_length(self, treshold):
        """select new sequences by length"""
        print("do something")
        for key in self.sp_seq_d:
            count = 0
            if len(self.sp_seq_d[key]) > treshold:
                print(key)
                for sp_keys in self.sp_seq_d[key].keys():
                    if isinstance(sp_keys, str) == True:
                        count += 1
                ReqMaxLength = max(self.sp_seq_d[key].values())
                ### !!! sometimes the only seq in seq_w_maxlen is the original seq, then this is the one to be added, but it will be removed,
                # later as it is no new seq! thus no new seq for that species is added
                ##
                seq_w_maxlen = {}
                for k, v in self.sp_seq_d[key].iteritems():
                    if len(v) == len(ReqMaxLength):
                        seq_w_maxlen[k] = v
                if (treshold - count) <= 0:
                    print("already to many samples of sp in aln, skip adding more.")
                    random_seq_ofsp = None
                elif len(seq_w_maxlen) == (treshold - count):
                    random_seq_ofsp = seq_w_maxlen
                elif len(seq_w_maxlen) > (treshold - count):
                    random_seq_ofsp = random.sample(seq_w_maxlen.items(), (treshold - count))
                    random_seq_ofsp = dict(random_seq_ofsp)
                else:
                    toselect = range(len(seq_w_maxlen), (treshold - count))
                    keymax = seq_w_maxlen.keys()
                    subdict = {k:v for k, v in self.sp_seq_d[key].iteritems() if k not in keymax}
                    SencondLen = max(subdict.values())
                    seq2len = {}
                    for k, v in subdict.iteritems():
                        if len(v) == len(SencondLen):
                            seq2len[k] = v
                    random_seq_ofsp = random.sample(seq2len.items(), len(toselect))
                    random_seq_ofsp = dict(random_seq_ofsp)
                    random_seq_ofsp.update(seq_w_maxlen)
                if random_seq_ofsp != None:
                    for key in random_seq_ofsp.keys():
                        # print(key)
                        self.filtered_seq[key] = random_seq_ofsp[key]
        # print(self.filtered_seq)

    def how_many_sp_to_keep(self, treshold, selectby):
        """selects number of sequences according to treshhold"""
        print("how_many_sp_to_keep")
        print(self.sp_d)
        print(self.sp_d.keys())

        general_wd = os.getcwd()

        for key in self.sp_d:
            print("key")
            print(key)
            if len(self.sp_d[key]) <= treshold:
                print("add all seq to tre")
                ##add all stuff to self.filtered_seq[gi_n], because the sequence number is below threshold, thus everythng will be added
                for giID in self.sp_d[key]:
                    if giID['^physcraper:last_blasted'] != '1800/01/01':
                        """seq was already blasted"""
                        if '^user:TaxonName' in giID:
                            userName = giID['^user:TaxonName']                               
                            for userName_aln, seq in self.data.aln.items():
                                if '^user:TaxonName' in self.data.otu_dict[userName_aln.label]:
                                    if userName == self.data.otu_dict[userName_aln.label]['^user:TaxonName']:
                                        self.filtered_seq[userName_aln.label] = seq.symbols_as_string()                                           
                        elif '^ot:ottTaxonName' in giID:
                            userName = giID['^ot:ottTaxonName']
                            for userName_aln, seq in self.data.aln.items():
                                if '^ot:ottTaxonName' in self.data.otu_dict[userName_aln.label]:
                                    if userName == self.data.otu_dict[userName_aln.label]['^ot:ottTaxonName']:
                                        self.filtered_seq[userName] = seq.symbols_as_string()
                    elif giID['^physcraper:last_blasted'] == '1800/01/01':

                        gi = giID['^ncbi:gi']
                        if gi in self.new_seqs.keys():
                            print(self.new_seqs[gi])
                            seq = self.new_seqs[gi]
                            # seq = self.sp_seq_d[key][gi]
                            self.filtered_seq[gi] = seq

            elif len(self.sp_d[key]) > treshold:
                if selectby == "length":
                    self.select_seq_by_length(treshold)
                else:
                    # print("write files for blast")
                    for giID in self.sp_d[key]:
                        # print(giID)
                        # print(giID['^physcraper:last_blasted'])
                        if giID['^physcraper:last_blasted'] != '1800/01/01':
                            # print("genrate files used for blasting")
                            if '^user:TaxonName' in giID:
                                userName = giID['^user:TaxonName']
                                # print(userName)
                                for userName_aln, seq in self.data.aln.items():
                                    # print(userName_aln)
                                    if '^user:TaxonName' in self.data.otu_dict[userName_aln.label]:
                                        if userName == self.data.otu_dict[userName_aln.label]['^user:TaxonName']:
                                            if selectby=="blast":
                                                print("user")
                                                # print(userName, seq.symbols_as_string())
                                                self.write_blast_files(userName, seq)
                            elif '^ot:ottTaxonName' in giID:
                                # print("seq already in aln")
                                userName = giID['^ot:ottTaxonName']
                                for userName_aln, seq in self.data.aln.items():
                                    if '^ot:ottTaxonName' in self.data.otu_dict[userName_aln.label]:
                                        if userName == self.data.otu_dict[userName_aln.label]['^ot:ottTaxonName']:
                                            if selectby=="blast":
                                                print("ott")
                                                # print(userName, seq)

                                                self.write_blast_files(userName, seq)
                        else:
                            # print("make gilist as local blast database")
                            if "^ncbi:gi" in giID:
                                gi = int(giID['^ncbi:gi'])
                                if selectby=="blast":
                                    print("new")
                                    # print(gi, seq)





                                    fp = False
                                    for k ,v in self.data.otu_dict.iteritems():
                                        if '^ncbi:gi' in v:
                                            print(v['^ncbi:gi'])
                                            if gi == v['^ncbi:gi']:
                                                fp = True
                                    if fp == False:
                                        seq = self.sp_seq_d[key][gi]

                                        self.write_blast_files(gi, seq, db=True , fn=key)
                            
        """add sequences according to blast similarity. """
        for taxonID in self.sp_seq_d:
            print(taxonID)
            print(treshold)
            print(len(self.sp_seq_d[taxonID]))
            if len(self.sp_seq_d[taxonID]) > treshold:
                print("blast locally")
                print(taxonID)
                count= 0
                print(self.sp_seq_d[taxonID].keys())
                for sp_keys in self.sp_seq_d[taxonID].keys():
                    if isinstance(sp_keys, str) == True:
                        count += 1
                if count >= 1:
                    print("count>0")
                    """species is not new in alignment, make blast with existing seq"""
                    if taxonID in self.sp_d.keys():
                        for element in self.sp_d[taxonID]:
                            if '^ot:ottTaxonName' in element:
                                    blast_seq = "{}".format(element['^ot:ottTaxonName'])
                                    blast_db = "{}".format(element['^ot:ottTaxonName'])
                    self.select_seq_by_local_blast(self.sp_seq_d[taxonID], blast_seq, blast_db, treshold, count)



                else:
                    print("completely new taxon to blast")
                    """species is completely new in alignment, make blast with random species"""
                    for gi in self.sp_seq_d[taxonID]:
                        print(gi)
                        self.data.add_otu(gi, self.ids, self.config.email)
                    blast_seq = self.sp_seq_d[taxonID].keys()[0]
                    if type(blast_seq) == int:
                        str_db = taxonID
                    else:
                        str_db = str(blast_seq)

                    blast_db = self.sp_seq_d[taxonID].keys()[1:]
                    # write files for local blast first:
                    seq = self.sp_seq_d[taxonID][blast_seq]
                    self.write_blast_files(taxonID, seq) #blast qguy
                    str_toblast = str(taxonID)
                    for blast_key in blast_db:
                        seq = self.sp_seq_d[taxonID][blast_key]
                        self.write_blast_files(blast_key, seq, db=True, fn=taxonID) #local db
                    # make local blast of sequences
                    print("new")
                    self.select_seq_by_local_blast(self.sp_seq_d[taxonID], str_toblast, str_db, treshold, count)
        return


    def write_blast_files(self, userName, seq, db=False, fn=None):
        print("writing files")
        if not os.path.exists("{}/blast".format(self.data.workdir)):
            os.makedirs("{}/blast/".format(self.data.workdir))
        if db == True:
            fnw = "./{}/blast/{}_db".format(self.workdir, fn)
            fi_o = open(fnw, 'a')
        else:
            fnw = "./{}/blast/{}_tobeblasted".format(self.workdir, userName)
            fi_o = open(fnw, 'w')
        print(fnw)
        fi_o.write(">{}\n".format(userName))
        fi_o.write("{}\n".format(seq))
        fi_o.close()

    def replace_new_seq(self):
        print("replace new seq")
        print(self.filtered_seq)
        keylist = self.filtered_seq.keys()
        print(keylist)
        keylist = [x for x in keylist if type(x) == int]
        print(keylist)
        SeqNotAdded = self.new_seqs.keys()
        SeqNotAdded = [x for x  in SeqNotAdded if type(x) == int]
        for gi in SeqNotAdded:
            for key in self.data.otu_dict.keys():
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if  self.data.otu_dict[key]['^ncbi:gi'] == gi:
                        # if self.data.otu_dict[key]['^physcraper:last_blasted'] == "1800/01/01":
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'not added, because there are enough seq per sp in tre'
        for gi in keylist:
            for key in self.data.otu_dict.keys():
                # print(self.data.otu_dict[key])
                if '^ncbi:gi' in self.data.otu_dict[key]:
                    if  self.data.otu_dict[key]['^ncbi:gi'] == gi:
                        # if self.data.otu_dict[key]['^physcraper:last_blasted'] == "1800/01/01":
                        self.data.otu_dict[key]['^physcraper:last_blasted'] = "1900/01/01"
                        self.data.otu_dict[key]['^physcraper:status'] = 'added, as one of the representatives of the taxon'
#             self.data.otu_dict[otu_id]['^ncbi:gi'] = gi
        reduced_gi_dict = {k: self.data.gi_dict[k] for k in keylist}
        print(reduced_gi_dict)
        self.data.gi_dict.clear()
        self.data.gi_dict = reduced_gi_dict
        reduced_new_seqs_dic = {k: self.filtered_seq[k] for k in keylist}
        print(reduced_new_seqs_dic)
        self.new_seqs = deepcopy(reduced_new_seqs_dic)
        self.new_seqs_otu_id = deepcopy(reduced_new_seqs_dic) ### !!! key is not exactly same format as before
        #set back to empty dict
        self.sp_d.clear()
        self.filtered_seq.clear()
        # print(self.new_seqs)
        # print(something_stupid)
        return self.new_seqs
