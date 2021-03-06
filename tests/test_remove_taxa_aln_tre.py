import os
import pickle
import sys
from physcraper import FilterBlast, ConfigObj, IdDicts

sys.stdout.write("\ntests remove_taxa_aln_tre\n")

# tests if I can remove sequences from aln and tre
workdir = "tests/output/test_remove_taxa_aln_tre"
configfi = "tests/data/test.config"
absworkdir = os.path.abspath(workdir)

try:
    conf = ConfigObj(configfi)
    data_obj = pickle.load(open("tests/data/precooked/tiny_dataobj.p", 'rb'))
    data_obj.workdir = absworkdir
    ids = IdDicts(conf, workdir=data_obj.workdir)
    ids.gi_ncbi_dict = pickle.load(open("tests/data/precooked/tiny_gi_map.p", "rb"))
except:
    sys.stdout.write("\n\nTest FAILED\n\n")
    sys.exit()
filteredScrape =  FilterBlast(data_obj, ids)

len_aln_before = len(filteredScrape.data.aln.as_string('phylip'))
len_tre_before = len(filteredScrape.data.tre.as_string(schema="newick"))
namespace_before = len(filteredScrape.data.aln.taxon_namespace)
namespace_tre_before = len(filteredScrape.data.tre.taxon_namespace)

for tax in filteredScrape.data.aln.taxon_namespace:
    filteredScrape.data.remove_taxa_aln_tre(tax.label)
    break

len_aln_after = len(filteredScrape.data.aln.as_string('phylip'))
len_tre_after = len(filteredScrape.data.tre.as_string(schema="newick"))
namespace_after = len(filteredScrape.data.aln.taxon_namespace)
namespace_tre_after = len(filteredScrape.data.tre.taxon_namespace)

try:
    assert len_aln_before != len_aln_after
    assert len_tre_before != len_tre_after
    assert namespace_before != namespace_after
    assert namespace_tre_before != namespace_tre_after
    sys.stdout.write("\nall subtests passed\n")
except:
    # print("test remove tax aln tre failed")
    sys.stderr.write("\ntest failed\n")

