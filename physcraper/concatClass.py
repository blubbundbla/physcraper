class Concat():
    """combine several physcraper runs of the same lineage with
     different genes into one final concatenated aln and tre"""
    def __init__(self, workdir_comb):
        # super(PhyscraperScrape, self).__init__()
        self.gene_comb = {}
        self.workdir = workdir_comb
        
    def combine(self, genelist):
    	count=1
    	for item in genelist:
    		genename = "gene{}".format(count)
        	# print(item.data.otu_dict)
        	
        	count+=1
        	self.gene_comb[genename] = item
        for key in self.gene_comb:
            print(key)   
            print(self.gene_comb[key].data)
            for key2, value in self.gene_comb[key].iteritems():
                print(key2)
                print(value)
