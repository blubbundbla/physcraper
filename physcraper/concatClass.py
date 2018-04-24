class Concat():
    """combine several physcraper runs of the same lineage with
     different genes into one final concatenated aln and tre"""
    def __init__(self):
        # super(PhyscraperScrape, self).__init__()
        self.gene_comb = {}
        
    def combine(self, gene1, gene2):
        print(gene1)
        print(gene2)