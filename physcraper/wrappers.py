#!/usr/bin/env python
import pickle
import sys
import os
import subprocess
import json
import csv
from ete2 import NCBITaxa
from physcraper import generate_ATT_from_phylesystem, generate_ATT_from_files, ConfigObj, IdDicts, PhyscraperScrape, FilterBlast, Concat
from dendropy import DnaCharacterMatrix


def sync_ncbi(configfi):
    conf = ConfigObj(configfi)
    subprocess.call(["rsync", "av", "ftp.ncbi.nih.gov::pub/taxonomy/gi_taxid_nucl.dmp.gz", "{}/gi_taxid_nucl.dmp.gz".format(conf.ncbi_dmp)])
    subprocess.call(["gunzip", "{}/gi_taxid_nucl.dmp.gz".format(dir)])


def sync_ott(configfi):
    conf = ConfigObj(configfi)
    subprocess.call(["process_ott.sh", "".format(conf.ott_ncbi)])

##generates IdDicts physcrapper class        
def get_ottid(configfi, cwd):
                     conf = ConfigObj(configfi)
                     ids = IdDicts(conf, cwd)  
                     return(ids)            



def standard_run(study_id,
                 tree_id,
                 seqaln,
                 mattype,
                 workdir,
                 configfi):
    '''looks for a json file to continue run, or builds and runs
    new analysis for as long as new seqs are found'''
    conf = ConfigObj(configfi)
#    if os.path.isfile("{}/att_checkpoint.json".format(workdir)):
#        sys.stdout.write("Reloading data object from json scrapefile\n")
#        thawed = open("{}/att_checkpoint.json".format(workdir), 'r').readlines()
#        data_obj = jsonpickle.decode(thawed)
#        scraper.repeat = 1
    if os.path.isfile("{}/att_checkpoint.p".format(workdir)):
        sys.stdout.write("Reloading data object from pickle file\n")
        data_obj = pickle.load( open( "{}/att_checkpoint.p".format(workdir), "rb" ) )
#        scraper.repeat = 1
    else:
#            sync_names()
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        #read the config file into a configuration object
        conf = ConfigObj(configfi)
        aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
        #Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_phylesystem(aln=aln,
                             workdir=workdir,
                             study_id = study_id,
                             tree_id = tree_id,
                             phylesystem_loc = conf.phylesystem_loc)
        #Mapping identifiers between OpenTree and NCBI requires and identifier dict object
        ids = IdDicts(conf, workdir="example")
        #Prune sequnces below a certain length threshold
        #This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()
        data_obj.write_files()
        data_obj.write_labelled(label='^ot:ottTaxonName')
        data_obj.write_otus('otu_info', schema='table')
        data_obj.dump()
        #Mapping identifiers between OpenTree and NCBI requires and identifier dict object
    if os.path.isfile(conf.id_pickle):
        sys.stdout.write("Reloading id dicts from {}\n".format(conf.id_pickle))
#        thawed_id = open(conf.id_json, 'r').readlines()
#        ids = jsonpickle.decode(thawed_id)
#        scraper.repeat = 1
        ids = pickle.load(open(conf.id_pickle, "rb" ))
    else:
        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()
        ids = IdDicts(conf, workdir=workdir)
        ids.dump()
#Now combine the data, the ids, and the configuration into a single physcraper scrape object
    scraper = PhyscraperScrape(data_obj, ids, conf)
    #run the ananlyses
    scraper.run_blast()
    scraper.read_blast()
    scraper.remove_identical_seqs()
    scraper.generate_streamed_alignment()
    while scraper.repeat == 1:
        scraper.data.write_labelled(label=  '^ot:ottTaxonName')
        scraper.data.write_otus("otu_info", schema='table')
        scraper.run_blast()
        scraper.read_blast()
        scraper.remove_identical_seqs()
        scraper.generate_streamed_alignment()


def OtuJsonDict(id_to_spn, configfi):
    """Make otu json dict, which is also produces within the openTreeLife-query"""
    cwd = os.getcwd()  
    ## reads input file into the var spInfo
    with open(id_to_spn, mode='r') as idtospn:
        reader = csv.reader(idtospn)
        spInfo = dict((rows[0], rows[1]) for rows in reader)
    #print(spInfoDict) 
 
    ###generate spinfodict
    
    ottdic = get_ottid(configfi, cwd) 
    print(ottdic)
    ncbi = NCBITaxa()    
    
    spInfoDict = {}
    for item in spInfo:
        spn = spInfo[item].replace("_", " ")
        name2taxid = ncbi.get_name_translator([spn])
       
        
        if len(name2taxid.items())>=1:
            ncbiid = name2taxid.items()[0][1][0]
            ott = ottdic.ncbi_to_ott[ncbiid]
            spn = ottdic.ott_to_name[ott]
            get_info = {'^ncbiID': ncbiid, '^ot:ottTaxonName': spn, '^ot:ottId': ott, '^user:TaxonName': spInfo[item],  '^physcraper:status': 'original','^physcraper:last_blasted' : "1900/01/01" }  
            spInfoDict[item] = get_info
        else:
            
            spInfoDict[item] = {'^user:TaxonName': spInfo[item],  '^physcraper:status': 'original','^physcraper:last_blasted' : "1900/01/01"}
    return  spInfoDict 



def own_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 treshold,
                 selectby,
                 downtorank,
                 spInfoDict,
                 add_local_seq,
                 id_to_spn_addseq_json,
                 configfi):
    '''looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found'''
    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        scraper = pickle.load(open("{}/scrape_checkpoint.p".format(workdir),'rb'))
        scraper.repeat = 1    
    else:   
#            sync_names()
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        #read the config file into a configuration object
        conf = ConfigObj(configfi)
        # print("config")
        print(dir(conf))
        print(conf.email)
        # print(seqaln, mattype)
        # #aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)

        #Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=seqaln, 
                             mattype=mattype, 
                             workdir=workdir,
                             treefile=trfn,
                             schema_trf = schema_trf,
                             otu_json=spInfoDict,
                             email = conf.email,
                             ingroup_mrca=None)

        #Prune sequnces below a certain length threshold
        #This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()

        data_obj.write_files()

        data_obj.write_labelled( label='user:TaxonName')
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()
        
        #ids = IdDicts(conf, workdir="example")
    #         if os.path.isfile("{}/id_pickle.p".format(workdir)): 

    #         #if os.path.isfile(conf.id_pickle):
    #             sys.stdout.write("Reloading id dicts from {}\n".format(conf.id_pickle))
    # #        thawed_id = open(conf.id_json, 'r').readlines()
    # #        ids = jsonpickle.decode(thawed_id)
    # #        scraper.repeat = 1
    #             ids = pickle.load(open("{}/id_pickle.p".format(workdir),'rb'))
    #         else:
        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()
        # if os.path.isfile("{}/id_pickle.p".format(workdir)): 
        #     sys.stdout.write("Reloading from pickled scrapefile: id\n")
        #     ids = pickle.load(open("{}/id_pickle.p".format(workdir),'rb'))

        # else:   
        ids = IdDicts(conf, workdir=workdir)
        ids.dump()

        # if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        #     sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        #     scraper = pickle.load(open("{}/scrape_checkpoint.p".format(workdir),'rb'))
        #     scraper.repeat = 1    
        # else:   
            #Now combine the data, the ids, and the configuration into a single physcraper scrape object
        scraper =  PhyscraperScrape(data_obj, ids, conf)
        if add_local_seq != None:
            print("will add local sequences now")
            scraper.add_local_seq(add_local_seq, id_to_spn_addseq_json)
            # scraper.replace_new_seq()
            scraper.remove_identical_seqs()
            scraper.generate_streamed_alignment(treshold)




        #run the ananlyses
        scraper.run_blast()
        scraper.read_blast()
        scraper.remove_identical_seqs()
        scraper.dump()
        print(treshold)
        if treshold != None:  

            scraper.sp_dict(downtorank)
            scraper.make_new_seq_dict(treshold=treshold, selectby=selectby)
            scraper.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
            scraper.replace_new_seq()
        print("from replace to streamed aln")
        scraper.generate_streamed_alignment(treshold)
    while scraper.repeat == 1: 
        scraper.data.write_labelled(label='user:TaxonName')
        scraper.data.write_otus("otu_info", schema='table')
        scraper.run_blast()
        scraper.read_blast()
        scraper.remove_identical_seqs()
#        scraper.how_many_sp_to_keep(treshold=treshhold)
        print("make sp_dict")    
        if treshold != None:  
            scraper.sp_dict(downtorank)
            filteredScrape = FilterBlast(scraper)
            scraper.make_new_seq_dict(treshold=treshold, selectby=selectby)
            scraper.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
            scraper.replace_new_seq()
        scraper.generate_streamed_alignment(treshold)

def concat(gene1, gene2, user_concat = None):
    scraper = PhyscraperScrape(data_obj, ids, conf)
    comb = scraper.combine(gene1, gene2)



def filter_data_run(seqaln,
                 mattype,
                 trfn,
                 schema_trf,
                 workdir,
                 treshold,
                 selectby,
                 downtorank,
                 spInfoDict,
                 add_local_seq,
                 id_to_spn_addseq_json,
                 configfi):
    '''looks for pickeled file to continue run, or builds and runs 
    new analysis for as long as new seqs are found'''
    if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        filteredScrape = pickle.load(open("{}/scrape_checkpoint.p".format(workdir),'rb'))
        filteredScrape.repeat = 1    
    else:   
#            sync_names()
        sys.stdout.write("setting up Data Object\n")
        sys.stdout.flush()
        #read the config file into a configuration object
        conf = ConfigObj(configfi)
        # print("config")
        print(dir(conf))
        print(conf.email)
        # print(seqaln, mattype)
        # #aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)

        #Generate an linked Alignment-Tree-Taxa object
        data_obj = generate_ATT_from_files(seqaln=seqaln, 
                             mattype=mattype, 
                             workdir=workdir,
                             treefile=trfn,
                             schema_trf = schema_trf,
                             otu_json=spInfoDict,
                             email = conf.email,
                             ingroup_mrca=None)

        #Prune sequnces below a certain length threshold
        #This is particularly important when using loci that have been de-concatenated, as some are 0 length which causes problems.
        data_obj.prune_short()

        data_obj.write_files()

        data_obj.write_labelled( label='user:TaxonName')
        data_obj.write_otus("otu_info", schema='table')
        data_obj.dump()
        
        #ids = IdDicts(conf, workdir="example")
    #         if os.path.isfile("{}/id_pickle.p".format(workdir)): 

    #         #if os.path.isfile(conf.id_pickle):
    #             sys.stdout.write("Reloading id dicts from {}\n".format(conf.id_pickle))
    # #        thawed_id = open(conf.id_json, 'r').readlines()
    # #        ids = jsonpickle.decode(thawed_id)
    # #        scraper.repeat = 1
    #             ids = pickle.load(open("{}/id_pickle.p".format(workdir),'rb'))
    #         else:
        sys.stdout.write("setting up id dictionaries\n")
        sys.stdout.flush()
        # if os.path.isfile("{}/id_pickle.p".format(workdir)): 
        #     sys.stdout.write("Reloading from pickled scrapefile: id\n")
        #     ids = pickle.load(open("{}/id_pickle.p".format(workdir),'rb'))

        # else:   
        ids = IdDicts(conf, workdir=workdir)
        ids.dump()

        # if os.path.isfile("{}/scrape_checkpoint.p".format(workdir)): 
        #     sys.stdout.write("Reloading from pickled scrapefile: scrape\n")
        #     scraper = pickle.load(open("{}/scrape_checkpoint.p".format(workdir),'rb'))
        #     scraper.repeat = 1    
        # else:   
            #Now combine the data, the ids, and the configuration into a single physcraper scrape object
        filteredScrape =  FilterBlast(data_obj, ids, conf)
        if add_local_seq != None:
            print("will add local sequences now")
            filteredScrape.add_local_seq(add_local_seq, id_to_spn_addseq_json)
            # scraper.replace_new_seq()
            filteredScrape.remove_identical_seqs()
            filteredScrape.generate_streamed_alignment(treshold)




        #run the ananlyses
        filteredScrape.run_blast()
        filteredScrape.read_blast()
        filteredScrape.remove_identical_seqs()
        filteredScrape.dump()
        print(treshold)
        if treshold != None:  

            filteredScrape.sp_dict(downtorank)
            filteredScrape.make_new_seq_dict(treshold=treshold, selectby=selectby)
            filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
            filteredScrape.replace_new_seq()
        print("from replace to streamed aln")
        filteredScrape.generate_streamed_alignment(treshold)
    while filteredScrape.repeat == 1: 
        filteredScrape.data.write_labelled(label='user:TaxonName')
        filteredScrape.data.write_otus("otu_info", schema='table')
        filteredScrape.run_blast()
        filteredScrape.read_blast()
        filteredScrape.remove_identical_seqs()
#        scraper.how_many_sp_to_keep(treshold=treshhold)
        print("make sp_dict")    
        if treshold != None:  
            filteredScrape.sp_dict(downtorank)
            filteredScrape.make_new_seq_dict(treshold=treshold, selectby=selectby)
            filteredScrape.how_many_sp_to_keep(treshold=treshold, selectby=selectby)
            filteredScrape.replace_new_seq()
        filteredScrape.generate_streamed_alignment(treshold)