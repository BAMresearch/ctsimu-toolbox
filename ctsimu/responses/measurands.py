import pandas as pd

class Measurands:
    def __init__ (self, Values, data, file_names):
    #def __init__ (self, Messgröße, CalPath, Kalibrierunsicherheit, Messwerte, sep, decimal, *metadata):
        #self.Name = Messgröße
        self.CalPath = data['calibration']['calibration_csv_path'] #CalPath
        #self.KUnsi = Kalibrierunsicherheit
        #elf.Valuese = Messwerte
        self.separator = data['real_ct']['csv_sep'] #sep
        self.decimal = data['real_ct']['csv_decimal'] #decimal
        self.data_Refwerte =()
        #self.CalValues = ()
        self.data_RefWerte = pd.read_csv(self.CalPath, sep = self.separator, decimal = self.decimal, header=0, index_col=0)
        self.CalValues = self.data_RefWerte['Kalibrierwert (mm)']
        self.Values = Values
        #self.SimValues = SimValues
        self.file_names = file_names
        self.MUSim = ()
        self.MU = ()


    def Diff_MW_KW(self):
        if isinstance(self.Values, float):
            #print(self.file_names)
            return self.Values-self.CalValues[self.file_names]
        elif isinstance(self.Values, dict):
            #print(self.Values)
            return self.Values-self.CalValues[self.Values.name]
        #self.file_names
        #print(self.CalValues[self.Values.name])
        #length(Messgröße)
        else: # or somehow is instance of a Series
            #print("ERROR")
            #print(type(self.Values))
            return self.Values-self.CalValues[self.Values.name]

    def Mittelwert(self):
        return self.Values.mean()
    
    def Unsicherheit(self):
        return self.Values.std()