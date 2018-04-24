import pickle
import os

class Concat():
    """combine several physcraper runs of the same lineage with
     different genes into one final concatenated aln and tre"""
    def __init__(self, workdir_comb):
        # super(PhyscraperScrape, self).__init__()
        self.workdir = workdir_comb
        self.gene_comb = {}
        self.sp_gi_comb = {}
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

    def combine(self, genelist):
    	count=1
    	for item in genelist:
    		genename = "gene{}".format(count)
        	# print(item.data.otu_dict)
        	
        	count+=1
        	self.gene_comb[genename] = item
        for genename in self.gene_comb:
            print(genename)

            # print(self.gene_comb[genename].data.otu_dict)
            for key2, value in self.gene_comb[genename].data.otu_dict.iteritems():

                # print(key2)
                # print(value['^ot:ottTaxonName'])
            #     if value['^ot:ottTaxonName'] not in self.sp_gi_comb:
            #         self.sp_gi_comb[value['^ot:ottTaxonName']] = {genename : value['^ncbi:gi']}
            #     elif value['^ot:ottTaxonName'] in self.sp_gi_comb:
            # #         # print(v["^ncbi:gi"])
            # #         self.gilist.append(v["^ncbi:gi"])
            #         self.sp_gi_comb[value['^ot:ottTaxonName']].append({genename : value['^ncbi:gi']})
                if value['^ot:ottTaxonName'] not in self.sp_gi_comb:
                    if '^ncbi:gi' in value:
                        self.sp_gi_comb[value['^ot:ottTaxonName']] = {genename : value['^ncbi:gi']}
                    else:
                        self.sp_gi_comb[value['^ot:ottTaxonName']] = {genename : key2}
                else:   
                    if '^ncbi:gi' in value:
                        self.sp_gi_comb[value['^ot:ottTaxonName']][genename] = value['^ncbi:gi']
                    else:
                        self.sp_gi_comb[value['^ot:ottTaxonName']][genename] = key2
     
        print(self.sp_gi_comb)
        return

    def dump(self, filename="concat_checkpoint.p"):
        pickle.dump(self, open("{}/{}".format(self.workdir, filename), "wb"))
