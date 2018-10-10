from physcraper import get_dataset_from_treebase, generate_ATT_from_files, generate_ATT_from_phylesystem, ConfigObj, IdDicts,  PhyscraperScrape
import pickle
import sys
import os


study_id = "pg_1799"
tree_id = "tree3637"

configfi = "tests/data/test.config"

conf = ConfigObj(configfi)

dataset = get_dataset_from_treebase(study_id,
                                phylesystem_loc='api')

aln = dataset[0]

ATTs = generate_ATTs_from_phylesystem(aln=aln,
                                  workdir='treebase',
                                  study_id=study_id,
                                  tree_id=tree_id,
                                  phylesystem_loc='local')