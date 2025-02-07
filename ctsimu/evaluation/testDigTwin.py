import json
import pandas as pd
from ..test import *
from ..helpers import *
from ..responses.run import *
#from ..responses.measurands import *
from ..primitives import Vector, Polygon
from ..scenario import Scenario
import numpy as np
from ..responses.measurands import Measurands
#from ..responses import *


class TestDigTwin(generalTest):
    def __init__(self, metadatafile):
        #empty declarations
        self.Zusatz_krit_real = []
        self.Zusatz_krit_sim = []
        self.combined_df = []
        
        #metadatfile --> Json File
        self.metadatafile = metadatafile
        self.data = read_json_file(self.metadatafile)

        #general Information for the test
        self.measurements_of_interest = self.data['general']['measurements_of_interest']
        self.csv_test_output_path = self.data['general']['csv_test_output_path']

        # Real CT Scan information
        self.real_folder_path = self.data['real_ct']['csv_real_path']
        self.real_output_path = self.data['real_ct']['csv_real_output_path']
        self.sep = self.data['real_ct']['csv_sep']
        self.decimal = self.data['real_ct']['csv_decimal']
        self.measurement_name_col = self.data['real_ct']['real_name_column']
        self.measurement_value_col = self.data['real_ct']['real_value_column']
        self.meas_name_col_nr = self.data['real_ct']['real_name_column_nr']
        self.meas_value_col_nr = self.data['real_ct']['real_value_column_nr']
        self.temperature = self.data['real_ct']['temperature']

        #Simulation Scan information
        self.sim_folder_path = self.data['simulation_ct']['csv_simu_path']
        self.sim_output_path = self.data['simulation_ct']['csv_sim_output_path']
            
        #Calibration Information
        self.CalPath = self.data['calibration']['calibration_csv_path'] #CalPath
        self.cal_sep = self.data['calibration']['csv_sep'] #sep
        self.cal_decimal = self.data['calibration']['csv_decimal'] #sep
        self.alpha = self.data['calibration']['material_alpha']
        self.Tcal = self.data['calibration']['temperature_calib']
        #self.calibration_values = self.data['calibration']['calibration_csv_path']
        
        
        self.data_RefWerte = pd.read_csv(self.CalPath, sep = self.cal_sep, decimal = self.cal_decimal, header=0, index_col=0)
        
        #Calibration Values
        self.CalValues = self.data_RefWerte['Kalibrierwert (mm)']  # muss an json angepasst werden!
        self.CalUncertainty = self.data_RefWerte['Kalibrierunsicherheit'] # muss an json angepasst werden!

       


    def data_input(self):
     
        with open(self.metadatafile, 'r') as file:
            self.data = json.load(file)
        return self.data

    def get_temperature(self):
        self.data = read_json_file(self.metadatafile)
        #with open(metadata_file, 'r') as file:
            #pd.data=[]
            #print(json.load(file))
            #self.data = json.load(file)
            #print(data)

        # Lese den Wert des Schlüssels 'run_info''temperature'
        

        # Ausgabe des Wertes
        print(f"the temperature is {self.temperature}")

    def read_and_filter_csv_files(self, *folder_paths):
    # Open and read the JSON file
        #with open(json_file, 'r') as file:
        #    data = json.load(file)
        
        # Extract the folder path, column names, and measurements of interest
        


        #self.runNames =  self.data['calibration_info']['material_alpha']

         # Dictionary to store measurement values for each file
        #measurements_dict = {measurement: [] for measurement in self.measurements_of_interest}
        # List to store DataFrames
        #dataframes = []
        combined_df_test = pd.DataFrame()
         # List to store file names
        file_names = []
        for folder_path in folder_paths:
            for filename in os.listdir(folder_path):
                if filename.endswith('.csv'):
                    file_path = os.path.join(folder_path, filename)
                    df = pd.read_csv(file_path,sep=self.sep, decimal=self.decimal, usecols=[ int(self.meas_name_col_nr), int(self.meas_value_col_nr)])

                    # Filter the DataFrame to include only the columns of interest
                    df_filtered = df[[self.measurement_name_col, self.measurement_value_col]]
                    
                    #df_filtered[self.measurement_value_col] = df_filtered[self.measurement_value_col].str[:-3].apply(convert_str_to_float)
                    # Filter the DataFrame to include only the measurements of interest
                    df_filtered = df_filtered[df_filtered[self.measurement_name_col].isin(self.measurements_of_interest)]
                    
                    df_filtered[self.measurement_value_col] = df_filtered[self.measurement_value_col].apply(convert_str_to_float)
                    df_filtered[self.measurement_value_col].apply(run.TempKorr, args=(float(self.temperature), float(self.Tcal), float(self.alpha)))
                    #print(df_filtered)
                    combined_df_test = pd.concat([combined_df_test, df_filtered[self.measurement_value_col]], axis=1)
                  

                    # Append the file name to the list
                    file_names.append(filename)
                   

            # Create a DataFrame from the dictionary
            
            combined_df_test.index = self.measurements_of_interest
            combined_df_test.columns = file_names
                      

            combined_df_subtracted = combined_df_test.apply(lambda x: Measurands(x, self.data, file_names).Diff_MW_KW(), axis=1)
            
            return combined_df_subtracted

    def En_calc(self, SimValues, RealValues):

        #self.RealValues = self.read_and_filter_csv_files(self, self.real_folder_path)
        #self.SimValues = self.read_and_filter_csv_files(self, self.sim_folder_path)
        df = pd.DataFrame()
        #print(SimValues)

        self.SimValues_avg = SimValues.mean(axis=1)
        self.RealValues_avg = RealValues.mean(axis=1)
        self.MUSim = SimValues.std(axis=1)
        self.MU = RealValues.std(axis=1)

        df["RealValues_avg"] = RealValues.mean(axis=1)
        df["RealValues_U"] = RealValues.std(axis=1)
        df["SimValues_avg"] = SimValues.mean(axis=1)
        df["SimValues_U"] = SimValues.std(axis=1)

        
        #print(self.SimValues_avg)
        self.EN = ((self.SimValues_avg) - (self.RealValues_avg)) / ((self.MUSim.pow(2) + self.MU.pow(2))).pow(1/2)
        #print(((self.MUSim.pow(2) + self.MU.pow(2))).pow(1/2))
        #print(self.EN)

        self.Zusatz_krit_real = self.MU / self.CalUncertainty
        self.Zusatz_krit_sim = self.MUSim / self.MU
        df["EN-wert"] = self.EN
        df["Zusatz_krit_real"] = self.MU / self.CalUncertainty
        df["Zusatz_krit_sim"] = self.MUSim / self.MU

        df["Test_Result"] = self.EN
        df["Test_Result"] = np.where((df["EN-wert"] < 1) & (df["EN-wert"] > 1) & (df["Zusatz_krit_real"] < 1) & (df["Zusatz_krit_sim"] < 1), True, False)
        
        print(df)
        df.to_csv(self.csv_test_output_path, self.sep, decimal=self.decimal, index=True)

        return

    def TwinTest(self):
        
        
        
        return
    



    # Example usage
    #json_file = 'test.json'
    #combined_df = read_and_filter_csv_files(json_file)

    # Export the combined DataFrame to a CSV file
    #combined_df.to_csv(output_path, index=False)

    # Output the combined DataFrame
    #print(combined_df)