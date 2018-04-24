#!/usr/bin/env python
"""Physcraper module"""
import sys
import json
import collections
from dendropy import Tree,\
                     DnaCharacterMatrix,\
                     DataSet,\
                     datamodel
from peyotl.api.phylesystem_api import PhylesystemAPI
from peyotl.sugar import tree_of_life
from peyotl.nexson_syntax import extract_tree,\
                                 get_subtree_otus,\
                                 extract_otu_nexson,\
                                 PhyloSchema

from configobj import ConfigObj
from iddicts import IdDicts
from aligntreetax import AlignTreeTax
from physcraperscrape import PhyscraperScrape, FilterBlast
from concatClass import Concat


_DEBUG = 1


#ATT is a dumb acronym for Alignment Tree Taxa object
def get_dataset_from_treebase(study_id,
                              phylesystem_loc='api'):
    try:
        nexson = get_nexson(study_id, phylesystem_loc)
    except:
        sys.stderr.write("couldn't find study id {} in phylesystem location {}\n".format(study_id, phylesystem_loc))
    treebase_url = nexson['nexml'][u'^ot:dataDeposit'][u'@href']
    if 'treebase' not in nexson['nexml'][u'^ot:dataDeposit'][u'@href']:
        sys.stderr.write("No treebase record associated with study ")
        sys.exit()
    else:
        tb_id = treebase_url.split(':S')[1]
        url = "https://treebase.org/treebase-web/search/downloadAStudy.html?id={}&format=nexus".format(tb_id)
#        url = "http://treebase.org/treebase-web/phylows/study/TB2:S{}?format=nexml".format(tb_id)
        if _DEBUG:
            sys.stderr.write(url +"\n")
        dna = DataSet.get(url=url,
                          schema="nexml")
        return dna

def generate_ATT_from_phylesystem(aln,
                                  workdir,
                                  study_id,
                                  tree_id,
                                  phylesystem_loc='api'):
    """gathers together tree, alignment, and study info - forces names to otu_ids.
    Outputs AlignTreeTax object.
    Input can be either a study ID and tree ID from OpenTree
    Alignemnt need to be a Dendropy DNA character matrix!"""

    #TODO CHECK ARGS
    assert isinstance(aln, datamodel.charmatrixmodel.DnaCharacterMatrix)
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_") #Forcing all spaces to underscore UGH
    nexson = get_nexson(study_id, phylesystem_loc)
    ott_ids = get_subtree_otus(nexson,
                               tree_id=tree_id,
                               subtree_id="ingroup",
                               return_format="ottid")
    ott_mrca = get_mrca_ott(ott_ids)
    newick = extract_tree(nexson,
                          tree_id,
                          PhyloSchema('newick',
                                      output_nexml2json='1.2.1',
                                      content="tree",
                                      tip_label="ot:originalLabel"))
    newick = newick.replace(" ", "_") #UGH Very heavy handed, need to make sure happens on alignement side as well.
    tre = Tree.get(data=newick,
                   schema="newick",
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)
    otus = get_subtree_otus(nexson, tree_id=tree_id)
    otu_dict = {}
    orig_lab_to_otu = {}
    treed_taxa = {}
    for otu_id in otus:
        otu_dict[otu_id] = extract_otu_nexson(nexson, otu_id)[otu_id]
        otu_dict[otu_id]['^physcraper:status'] = "original"
        otu_dict[otu_id]['^physcraper:last_blasted'] = "1800/01/01"
        orig = otu_dict[otu_id].get(u'^ot:originalLabel').replace(" ", "_")
        orig_lab_to_otu[orig] = otu_id
        treed_taxa[orig] = otu_dict[otu_id].get(u'^ot:ottId')
    for tax in aln.taxon_namespace:
        try:
            tax.label = orig_lab_to_otu[tax.label].encode('ascii')
        except KeyError:
            sys.stderr.write("{} doesn't have an otu id. It is being removed from the alignement. This may indicate a mismatch between tree and alignement\n".format(tax.label))
   #need to prune tree to seqs and seqs to tree...
    otu_newick = tre.as_string(schema="newick")
    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir)
    #newick should be bare, but alignement should be DNACharacterMatrix


def convert(data):
    """convert json 2.7 as string problem"""
    if isinstance(data, basestring):
        return str(data)
    elif isinstance(data, collections.Mapping):
        return dict(map(convert, data.iteritems()))
    elif isinstance(data, collections.Iterable):
        return type(data)(map(convert, data))
    else:
        return data

def generate_ATT_from_files(seqaln,
                            mattype,
                            workdir,
                            treefile,
                            otu_json,
                            email,
                            schema_trf="newick",
                            ingroup_mrca=None):
    """Build an ATT object without phylesystem.
    If no ingroup mrca ott_id is provided, will use all taxa in tree to calc mrca.
    otu_json should encode the taxon names for each tip"""
    aln = DnaCharacterMatrix.get(path=seqaln, schema=mattype)
    for tax in aln.taxon_namespace:
        tax.label = tax.label.replace(" ", "_") #Forcing all spaces to underscore UGH

    tre = Tree.get(path=treefile,
                   schema=schema_trf,
                   preserve_underscores=True,
                   taxon_namespace=aln.taxon_namespace)

    assert tre.taxon_namespace is aln.taxon_namespace

    otu_newick = tre.as_string(schema=schema_trf)

    if ingroup_mrca:
        ott_mrca = int(ingroup_mrca)
    else:
        otu_dict = json.load(open(otu_json, "r"))
        ott_ids = [otu_dict[otu].get(u'^ot:ottId',) for otu in otu_dict]
        ott_ids = filter(None, ott_ids)
        ott_ids = set(ott_ids)
        ott_mrca = get_mrca_ott(ott_ids)
        print(ott_mrca)

    return AlignTreeTax(otu_newick, otu_dict, aln, ingroup_mrca=ott_mrca, workdir=workdir, schema=schema_trf)


#####################################
def get_nexson(study_id, phylesystem_loc):
    """Grabs nexson from phylesystem"""
    phy = PhylesystemAPI(get_from=phylesystem_loc)
    nexson = phy.get_study(study_id)['data']
    return  nexson

def get_mrca_ott(ott_ids):
    """finds the mrca of the taxa in the ingroup of the original
    tree. The blast search later is limited to descendents of this
    mrca according to the ncbi taxonomy"""
    if None in ott_ids:
        ott_ids.remove(None)
    synth_tree_ott_ids = []
    ott_ids_not_in_synth = []
    for ott in ott_ids:
        try:
            tree_of_life.mrca(ott_ids=[ott], wrap_response=False)
            synth_tree_ott_ids.append(ott)
        except:
            ott_ids_not_in_synth.append(ott) 
    mrca_node = tree_of_life.mrca(ott_ids=synth_tree_ott_ids, wrap_response=False)# need to fix wrap eventually
    try:
        tax_id = mrca_node[u'mrca'][u'taxon'][u'ott_id']
        sys.stdout.write('MRCA of sampled taxa is {}\n'.format(mrca_node[u'mrca'][u'taxon'][u'name']))
    except: #Hackaround for V2 apis
        tax_id = mrca_node[u'nearest_taxon_mrca_ott_id']
        sys.stdout.write('MRCA of sampled taxa is {}\n'.format(mrca_node[u'nearest_taxon_mrca_unique_name']))
    return tax_id

def get_ott_ids_from_otu_dict(otu_dict): #TODO put into data obj?
    """Get the ott ids from an otu dict object"""
    ott_ids = []
    for otu in otu_dict:
        try:
            ott_ids.append(otu['^ot:ottId'])
        except KeyError:
            pass
