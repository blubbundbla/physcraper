#!/usr/bin/env python
"""Physcraper module"""
import sys
import os
import subprocess
import pickle
from Bio import Entrez
from ete2 import NCBITaxa

_DEBUG = 1

class IdDicts(object):
    """Wraps up the annoying conversions"""
    #TODO - could - should be shared acrosss runs?! .... nooo.
    def __init__(self, config_obj, workdir):
        """Generates a series of name disambiguation dicts"""
        self.workdir = workdir
        self.config = config_obj
        self.ott_to_ncbi = {}
        self.ncbi_to_ott = {}
        self.ott_to_name = {}
        self.gi_ncbi_dict = {}
        fi = open(config_obj.ott_ncbi) #TODO need to keep updated
        for lin in fi:
            lii = lin.split(",")
            self.ott_to_ncbi[int(lii[0])] = int(lii[1])
            self.ncbi_to_ott[int(lii[1])] = int(lii[0])
            self.ott_to_name[int(lii[0])] = lii[2].strip()
            assert len(self.ott_to_ncbi) > 0
            assert len(self.ncbi_to_ott) > 0
            assert len(self.ott_to_name) > 0
        fi.close()
        if os.path.isfile("{}/id_map.txt".format(workdir)): #todo config?!
            fi = open("{}/id_map.txt".format(workdir))
            for lin in fi:
                self.gi_ncbi_dict[int(lin.split(",")[0])] = lin.split(",")[1]

    def map_gi_ncbi(self, gi):
        """get the ncbi taxon id's for a gi input"""
        if gi in self.gi_ncbi_dict:
            tax_id = int(self.gi_ncbi_dict[gi])
        else:
            try:
                Entrez.email = self.config.email
                handle = Entrez.efetch(db="nucleotide", id=gi, retmode="xml")
                read_handle = Entrez.read(handle)[0]
                tax_name = read_handle['GBSeq_feature-table'][0]['GBFeature_quals'][0]['GBQualifier_value']
                try:
                    tax_id = Entrez.read(Entrez.esearch(db="taxonomy", term=tax_name, RetMax=100))['IdList'][0]
                except:
                    ncbi = NCBITaxa()
                    tax_info = ncbi.get_name_translator([tax_name])
                    tax_id = tax_info.items()[0][1][0]
            except:

                if _DEBUG:
                    sys.stderr.write("Entrez fetch failed, for GI:{} running subprocess to get taxon id\n".format(gi))
                tax_id = int(subprocess.check_output(["bash", self.config.get_ncbi_taxonomy,
                                                      "{}".format(gi),
                                                      "{}".format(self.config.ncbi_dmp)]).split('\t')[1])
            self.gi_ncbi_dict[gi] = tax_id
        return tax_id

    def dump(self):
        """save object to file"""
        filename = "id_pickle.p"
        pickle.dump(IdDicts, open("{}/{}".format(self.workdir,filename), "wb"))
 