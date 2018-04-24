#!/usr/bin/env python
"""Physcraper module"""
import sys
import os
import configparser


_DEBUG = 1

def is_number(s):
    """test if string can be coerced to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False

class ConfigObj(object):
    """Pulls out the configuration information from
    the config file and makes it easier to pass
    around and store."""
    def __init__(self, configfi):
        if _DEBUG:
            sys.stderr.write("Building config object\n")
        assert os.path.isfile(configfi)
        config = configparser.ConfigParser()
        config.read(configfi)
        self.e_value_thresh = config['blast']['e_value_thresh']
        assert is_number(self.e_value_thresh)
        self.hitlist_size = int(config['blast']['hitlist_size'])
        self.seq_len_perc = float(config['physcraper']['seq_len_perc'])
        assert 0 < self.seq_len_perc < 1
        self.get_ncbi_taxonomy = config['taxonomy']['get_ncbi_taxonomy']
        assert os.path.isfile(self.get_ncbi_taxonomy)
        self.ncbi_dmp = config['taxonomy']['ncbi_dmp']
        if not os.path.isfile(self.ncbi_dmp):
            os.system("rsync -av ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz {}.gz".format(self.config.ncbi_dmp))
            os.system("tar -xzvf taxonomy/gi_taxid_nucl.dmp.gz")
            self.ncbi_dmp = "taxonomy/gi_taxid_nucl.dmp.gz"
        self.phylesystem_loc = config['phylesystem']['location']
        assert(self.phylesystem_loc in ['local', 'api'])
        self.ott_ncbi = config['taxonomy']['ott_ncbi']
        assert os.path.isfile(self.ott_ncbi)
        self.id_pickle = os.path.abspath(config['taxonomy']['id_pickle'])#TODO what is theis doing again?
        self.email = config['blast']['Entrez.email']
        self.blast_loc = config['blast']['location']
        self.num_threads = config['blast'].get('num_threads')
        assert self.blast_loc in ['local', 'remote']
        if self.blast_loc == 'local':
            self.blastdb = config['blast']['localblastdb']
        if self.blast_loc == 'remote':
            self.url_base = config['blast'].get('url_base')
        if _DEBUG:
            sys.stderr.write("{}\n".format(self.email))
            if self.blast_loc == 'remote':
                sys.stderr.write("url base = {}\n".format(self.url_base))
            sys.stderr.write("{}\n".format(self.blast_loc))
            if self.blast_loc == 'local':
                sys.stderr.write("local blast db {}\n".format(self.blastdb))