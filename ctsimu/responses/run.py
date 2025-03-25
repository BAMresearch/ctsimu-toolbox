class run:
    def __init__ (self):
        self.runName = 0	
        self.Mvalue = 0
        #self.Temperature = 0
        self.alpha = 0
        self.Tmeas = 0
        self.Tcal = 0

    def TempKorr(Mvalue, Tmeas, Tcal, alpha):
        return Mvalue*(1-(alpha*(Tmeas-Tcal)))
    
    def ScalCorr(Mvalue, skalcoeff):
        return Mvalue*skalcoeff
    