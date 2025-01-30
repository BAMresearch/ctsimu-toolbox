class run:
    def __init__ (self, metdatafile):
        self.runName = 0	
        self.MWert = 0
        #self.Temperature = 0
        self.alpha = 0
        self.Tmeas = 0
        self.Tcal = 0

    def TempKorr(MWert, Tmeas, Tcal, alpha):
	    return MWert*(1-(alpha*(Tmeas-Tcal)))