from __init__ import PhyscraperScrape

class FilterBlast(PhyscraperScrape):
	def __init__(self):
		PhyscraperScrape.__init__(self)
		self.sp_d = {}
        self.sp_seq_d = {}
        self.filtered_seq = {}


