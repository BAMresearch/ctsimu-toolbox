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
        self.name = 0
        self.temperature = 0
        self.data = []
        self.folder_path = 0
        self.output_path = 0
        self.measurement_name_col = 0
        self.measurement_value_col = 0
        self.measurements_of_interest = 0
        self.sep = 0
        self.decimmal = 0
        self.meas_name_col_nr = 0
        self.meas_value_col_nr = 0
        self.system = 0
        self.alpha = 0
        self.Tcal =0
        self.runNames = 0
        self.calibration_values = []

        self.combined_df = []
        self.combined_df_test = pd.DataFrame()

        self.combined_df_mean = []
        self.combined_df_subtracted =[]

    def data_input(self, metadata_file):
     
        with open(metadata_file, 'r') as file:
            self.data = json.load(file)
        return self.data

    def get_temperature(self, metadata_file):
        self.data = read_json_file(metadata_file)
        #with open(metadata_file, 'r') as file:
            #pd.data=[]
            #print(json.load(file))
            #self.data = json.load(file)
            #print(data)

        # Lese den Wert des Schlüssels 'run_info''temperature'
        self.temperature = self.data['run_info']['temperature']

        # Ausgabe des Wertes
        print(f"the temperature is {self.temperature}")

    def read_and_filter_csv_files(self, metadata_file):
    # Open and read the JSON file
        #with open(json_file, 'r') as file:
        #    data = json.load(file)
        self.data = read_json_file(metadata_file)
        # Extract the folder path, column names, and measurements of interest
        self.folder_path = self.data['run_info']['csv_folder_path']
        self.output_path = self.data['run_info']['csv_output_path']
        self.measurement_name_col = self.data['run_info']['measurement_name_column']
        self.measurement_value_col = self.data['run_info']['measurement_value_column']
        self.measurements_of_interest = self.data['run_info']['measurements_of_interest']
        self.sep = self.data['run_info']['csv_sep']
        self.decimmal = self.data['run_info']['csv_decimal']
        self.meas_name_col_nr = self.data['run_info']['measurement_name_column_nr']
        self.meas_value_col_nr = self.data['run_info']['measurement_value_column_nr']
        self.system = self.data['run_info']['system']
        self.temperature = self.data['run_info']['temperature']
        self.alpha = self.data['calibration_info']['material_alpha']
        self.Tcal = self.data['calibration_info']['temperature_calib']
        self.calibration_values = self.data['calibration_info']['calibration_csv_path']


        #self.runNames =  self.data['calibration_info']['material_alpha']

         # Dictionary to store measurement values for each file
        measurements_dict = {measurement: [] for measurement in self.measurements_of_interest}
        # List to store DataFrames
        dataframes = []

         # List to store file names
        file_names = []

        for filename in os.listdir(self.folder_path):
            if filename.endswith('.csv'):
                file_path = os.path.join(self.folder_path, filename)
                df = pd.read_csv(file_path,sep=self.sep, decimal=self.decimmal, usecols=[ int(self.meas_name_col_nr), int(self.meas_value_col_nr)])

                # Filter the DataFrame to include only the columns of interest
                df_filtered = df[[self.measurement_name_col, self.measurement_value_col]]
                
                #df_filtered[self.measurement_value_col] = df_filtered[self.measurement_value_col].str[:-3].apply(convert_str_to_float)
                # Filter the DataFrame to include only the measurements of interest
                df_filtered = df_filtered[df_filtered[self.measurement_name_col].isin(self.measurements_of_interest)]
                
                df_filtered[self.measurement_value_col] = df_filtered[self.measurement_value_col].apply(convert_str_to_float)
                df_filtered[self.measurement_value_col].apply(run.TempKorr, args=(float(self.temperature), float(self.Tcal), float(self.alpha)))
                #print(df_filtered)
                self.combined_df_test = pd.concat([self.combined_df_test, df_filtered[self.measurement_value_col]], axis=1)

                # Append the measurement values to the dictionary
                for measurement in self.measurements_of_interest:
                    values = df_filtered[df_filtered[self.measurement_name_col] == measurement][self.measurement_value_col].values
                    
                    values = [run.TempKorr(value, int(self.temperature), int(self.Tcal), float(self.alpha)) for value in values]
                    #print(values)
                    measurements_dict[measurement].append(values[0] if len(values) > 0 else None)

                # Append the file name to the list
                file_names.append(filename)
                #self.combined_df_test = 
                #df_filtered['Normal z'] = df_filtered['Normal z'].str[:-3].apply(convert_str_to_float)
                #print(type(df_filtered[self.measurement_value_col]))
                #print(measurements_dict)

        # Create a DataFrame from the dictionary
        self.combined_df = pd.DataFrame(measurements_dict)
        self.combined_df_test.index = self.measurements_of_interest
        self.combined_df_test.columns = file_names
        #print(self.combined_df_test)
        # Transpose the DataFrame and add the file names as columns
        self.combined_df = self.combined_df.T
        self.combined_df.columns = file_names

        # Add a column for measurement names
        #combined_df.insert(0, 'Measurement Name', self.measurements_of_interest)
        print(self.combined_df)

        #dataevl = Measurands(self.combined_df_test, self.data, file_names)
        #dataevl.Diff_MW_KW()
        #self.combined_df_subtracted = self.combined_df_test.apply(dataevl.Diff_MW_KW, axis=1)
        

        #self.combined_df_subtracted = self.combined_df_test.apply(lambda x: Measurands(x, self.data, file_names).Diff_MW_KW(), axis=1)
        self.combined_df_mean = self.combined_df_test.apply(lambda x: Measurands(x,x, self.data, file_names).Mittelwert(), axis=1)
        print(type(self.combined_df_mean))
        self.combined_df_subtracted = self.combined_df_mean.apply(lambda x: Measurands(x,x, self.data, self.combined_df_mean[self.combined_df_mean == x].index[0]).Diff_MW_KW())
        #data_RefWerte = pd.read_csv(self.calibration_values, sep=";", decimal=",", header=0, index_col=0)
        #self.combined_df_subtracted = self.combined_df_test.apply(measurands.Diff_MW_KW, axis=1)
        #self.combined_df_test.iloc[1] = measurands.Diff_MW_KW(self.combined_df_test.iloc[1], data_RefWerte['Kalibrierwert (mm)'][1]) args=(data_RefWerte['Kalibrierwert (mm)'])
        print(self.combined_df_subtracted)
        #return combined_df

    #def Test_En()
    



    # Example usage
    #json_file = 'test.json'
    #combined_df = read_and_filter_csv_files(json_file)

    # Export the combined DataFrame to a CSV file
    #combined_df.to_csv(output_path, index=False)

    # Output the combined DataFrame
    #print(combined_df)