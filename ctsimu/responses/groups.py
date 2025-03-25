class Messgröße:	
    def __init__ (self, Name):		
        self.Name=Name
        self.Mwert = Mwert

    def TempKorr(self):
        return self.MWert*(1-alpha*DeltaT)

    def Diff_MW_KW(self):
        return self.MWert-KWert

    def ENWert(self):
        return