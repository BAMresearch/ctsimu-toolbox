import pandas as pd

class Measurands:
    def __init__ (self, RealValues, SimValues, data, file_names):
    #def __init__ (self, Messgröße, CalPath, Kalibrierunsicherheit, Messwerte, sep, decimal, *metadata):
        #self.Name = Messgröße
        self.CalPath = data['calibration_info']['calibration_csv_path'] #CalPath
        #self.KUnsi = Kalibrierunsicherheit
        #elf.RealValuese = Messwerte
        self.separator = data['run_info']['csv_sep'] #sep
        self.decimal = data['run_info']['csv_decimal'] #decimal
        self.data_Refwerte =()
        #self.CalValues = ()
        self.data_RefWerte = pd.read_csv(self.CalPath, sep = self.separator, decimal = self.decimal, header=0, index_col=0)
        self.CalValues = self.data_RefWerte['Kalibrierwert (mm)']
        self.RealValues = RealValues
        self.SimValues = SimValues
        self.file_names = file_names
        self.MUSim = ()
        self.MU = ()


    def Diff_MW_KW(self):
        if isinstance(self.RealValues, float):
            print(self.file_names)
            return self.RealValues-self.CalValues[self.file_names]
        elif isinstance(self.RealValues, dict):
            print(self.RealValues)
            return self.RealValues-self.CalValues[self.RealValues.name]
        #self.file_names
        #print(self.CalValues[self.RealValues.name])
        #length(Messgröße)
        else:
            print("E")
            return self.RealValues-self.CalValues[self.RealValues.name]

    def Mittelwert(self):
        return self.RealValues.mean()
    
    def Unsicherheit(self):
        return self.RealValues.std()
    
    def ENWert(self):
        
        self.SimValues_avg = self.SimValues.mean()
        self.RealValues_avg = self.RealValues.mean()

        self.MUSim = self.SimValues.std()
        self.MU = self.RealValues.std()
        
        self.EN = (self.SimValues_avg - self.CalValues) - (self.RealValues_avg - self.CalValues) / (self.MUSim.pow(2) + self.MU.pow(2)).pow(1/2)
        return
    #def ENWert(self, runTemp, …):
