# -*- coding: UTF-8 -*-
"""# Test of digital twin

.. include:: ./testTwin.md
"""

import json
import pandas as pd
from ..test import *
from ..helpers import *
#from ..responses.run import *
#from ..responses.measurands import *
from ..primitives import Vector, Polygon
from ..scenario import Scenario
import numpy as np
from ..responses.measurands import Measurands
#from ..responses import *
from reportlab.pdfgen import canvas 
from reportlab.pdfbase.ttfonts import TTFont 
from reportlab.pdfbase import pdfmetrics 
from reportlab.lib import colors 

import matplotlib.pyplot as plt


class TestDigTwin(generalTest):
    """ CTSimU2 test of digital twin. """

    def __init__(self, metadatafile):
        #empty declarations
        self.Zusatz_krit_real = []
        self.Zusatz_krit_sim = []
        self.combined_df = []
        self.df = pd.DataFrame()
        #metadatfile --> Json File
        self.metadatafile = metadatafile
        self.data = read_json_file(self.metadatafile)
        print(self.data['real_ct']['csv_path'])

        #general Information for the test
        self.measurements_of_interest = self.data['general']['measurements_of_interest']
        self.csv_test_output_path = self.data['general']['csv_test_output_path']

        # Real CT Scan information
        self.real_folder_path = self.data['real_ct']['csv_path']
        #print(self.real_folder_path)
        self.output_path = self.data['real_ct']['csv_output_path']
        self.sep = self.data['real_ct']['csv_sep']
        self.decimal = self.data['real_ct']['csv_decimal']
        self.measurement_name_col = self.data['real_ct']['name_column']
        self.measurement_value_col = self.data['real_ct']['value_column']
        self.meas_name_col_nr = self.data['real_ct']['name_column_nr']
        self.meas_value_col_nr = self.data['real_ct']['value_column_nr']
        self.temperature = self.data['real_ct']['temperature']
        self.scaling_factor = self.data['real_ct']["scaling_factor"]
        

        #Simulation Scan information
        self.sim_folder_path = self.data['simulation_ct']['csv_path']
        
        self.sim_output_path = self.data['simulation_ct']['csv_output_path']
        print(self.sim_output_path)

        #Calibration Information
        self.CalPath = self.data['calibration']['calibration_csv_path'] #CalPath
        self.cal_sep = self.data['calibration']['csv_sep'] #sep
        self.cal_decimal = self.data['calibration']['csv_decimal'] #sep
        self.alpha = self.data['calibration']['material_alpha']
        self.Tcal = self.data['calibration']['temperature_calib']
        self.CalValCol = self.data['calibration']['values_column']
        self.uncertaintyCol = self.data['calibration']['uncertainty_columns']
        self.STLvalCol = self.data['calibration']['STL-value_columns_nr']
        #self.calibration_values = self.data['calibration']['calibration_csv_path']

        self.data_RefWerte = pd.read_csv(self.CalPath, sep = self.cal_sep, decimal = self.cal_decimal, header=0, index_col=0)

        #Calibration Values
        self.CalValues = self.data_RefWerte[self.CalValCol]  # muss an json angepasst werden!
        self.CalUncertainty = self.data_RefWerte[self.uncertaintyCol] # muss an json angepasst werden!

        self.fileName = 'sample.pdf'
        self.documentTitle = 'sample'
        self.title = 'Report Test Digital Twin'
        self.subTitle = 'CTSimU2'
        self.textLines = [ 
            'Text for the Report', 
            'Anything that matters', 
            ] 
        #self.image = 'image.jpg'


    def read_and_filter_csv_files(self, folder_paths, ct_type):
    # Variables needed for the function
        #self.real_folder_path = self.data[ct_type]['csv_path']
        self.real_output_path = self.data[ct_type]['csv_output_path']
        self.sep = self.data[ct_type]['csv_sep']
        self.decimal = self.data[ct_type]['csv_decimal']
        self.measurement_name_col = self.data[ct_type]['name_column']
        self.measurement_value_col = self.data[ct_type]['value_column']
        self.meas_name_col_nr = self.data[ct_type]['name_column_nr']

        self.output_path = self.data[ct_type]['csv_output_path']

        self.meas_value_col_nr = self.data[ct_type]['value_column_nr']
        if ct_type == "simulation_ct":
            self.temperature = 20
            self.scaling_factor = 1
            #self.data_RefWerte = pd.read_csv(self.CalPath, sep = self.cal_sep, decimal = self.cal_decimal, header=0, index_col=0)
            #Calibration Values
            self.CalValues = self.data_RefWerte[self.STLvalCol]  # muss an json angepasst werden!
            #self.CalUncertainty = self.data_RefWerte['Kalibrierunsicherheit'] # muss an json angepasst werden!

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
        for folder_path in [folder_paths]:
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
                    #print(df_filtered)
                    #if len(Tmeas) != len()
                    #df_filtered[self.measurement_value_col] = df_filtered[self.measurement_value_col].apply(run.TempKorr, args=(float(self.temperature), float(self.Tcal), float(self.alpha)))
                    df_filtered[self.measurement_value_col] = df_filtered[self.measurement_value_col].apply(run.TempKorr, args=(float(self.temperature), float(self.Tcal), float(self.alpha)))
                    # print(df_filtered)
                    combined_df_test = pd.concat([combined_df_test, df_filtered[self.measurement_value_col]], axis=1)

                    # Append the file name to the list
                    file_names.append(filename)
            #pd.set_option('display.max_rows', None)
            #pd.set_option('display.max_columns', None)
            #pd.set_option('display.width', None)
            #pd.set_option('display.max_colwidth', -1)
            #print(combined_df_test)

            # Create a DataFrame from the dictionary

            combined_df_test.index = self.measurements_of_interest
            combined_df_test.columns = file_names

            #combined_df_subtracted = combined_df_test.apply(lambda x: Measurands(x, self.data, file_names).Diff_MW_KW(), axis=1)
            #combined_df_subtracted = combined_df_test.apply(lambda x: Measurands(x, self.data, file_names).Diff_MW_KW(self.CalValues), axis=1)
            combined_df_subtracted = combined_df_test.apply(lambda x: self.Diff_MW_KW(x, file_names, self.CalValues), axis=1)

            print(combined_df_subtracted)
            combined_df_test.to_csv(self.output_path, self.sep, decimal=self.decimal, index=True)
            return combined_df_subtracted

    def TempKorr(self, Mvalue, Tmeas, Tcal, alpha):
        return Mvalue*(1-(alpha*(Tmeas-Tcal)))

    def Diff_MW_KW(self, values, file_names, CalValues):
        #print(values.name)
        #print(CalValues)
        #print(values-CalValues[values.name])
        return values-CalValues[values.name]


    def En_calc(self, RealValues, SimValues):

        #self.RealValues = self.read_and_filter_csv_files(self, self.real_folder_path)
        #self.SimValues = self.read_and_filter_csv_files(self, self.sim_folder_path)

        #print(SimValues)
        #print(RealValues)
        #print(SimValues)
        self.RealValues_avg = RealValues.mean(axis=1)
        self.SimValues_avg = SimValues.mean(axis=1)

        u_sim = SimValues.std(axis=1)
        u_p = RealValues.std(axis=1)
        u_drift = 0
        u_b = 0
        self.MU = 2*(RealValues.std(axis=1).pow(2) + self.CalUncertainty.pow(2)).pow(1/2)
        self.MUSim = 2*u_sim

        self.df["RealValues_avg"] = RealValues.mean(axis=1)
        self.df["RealValues_u"] = RealValues.std(axis=1)
        self.df["SimValues_avg"] = SimValues.mean(axis=1)
        self.df["SimValues_u"] = SimValues.std(axis=1)

        #print(self.SimValues_avg)
        #Berechnung EN-Wert, noch zu verändern.
        self.EN = ((self.SimValues_avg) - (self.RealValues_avg)) * ((self.MUSim.pow(2) + self.MU.pow(2)).pow(1/2)).pow(-1)
        #print(((self.MUSim.pow(2) + self.MU.pow(2))).pow(1/2))
        #print(self.EN)

        self.Zusatz_krit_real = self.MU / self.CalUncertainty
        self.Zusatz_krit_sim = self.MUSim / self.MU
        self.df["EN-wert"] = self.EN
        self.df["Zusatz_krit_real"] = self.MU / self.CalUncertainty
        self.df["Zusatz_krit_sim"] = self.MUSim / self.MU
        self.df["Zusatz_krit_summe"] = self.MUSim / self.MU + self.MU / self.CalUncertainty
        self.df["Test_Result"] = self.EN
        #self.df["Test_Result"] = np.where((self.df["EN-wert"] < 1) & (self.df["EN-wert"] > 1) & (self.df["Zusatz_krit_real"] < 1) & (self.df["Zusatz_krit_sim"] < 1), True, False)
        self.df["Test_Result"] = np.where((abs(self.df["EN-wert"])< 1) & (self.df["Zusatz_krit_real"] < 1) & (self.df["Zusatz_krit_sim"] < 1), True, False)
        print(self.df)

        self.df.to_csv(self.csv_test_output_path, self.sep, decimal=self.decimal, index=True)

        return


    def plotResults(self):
        plt.figure(figsize=(10, 6))
        plt.plot(self.df['EN-wert'], marker='o')

        # Add horizontal lines at -1 and +1
        plt.axhline(y=-1, color='r', linestyle='--')
        plt.axhline(y=1, color='g', linestyle='--')

        # Label the x-axis with names from self.measurements_of_interest
        plt.xticks(ticks=range(len(self.measurements_of_interest)), labels=self.measurements_of_interest, rotation=90)

        # Add labels and title
        plt.xlabel('Measurements')
        plt.ylabel('EN-wert')
        plt.title('EN-wert Series')

        # Show the plot
        plt.show()


    def TwinTest_report(self):
        # creating a pdf object
        pdf = canvas.Canvas(self.fileName)
        # setting the title of the document
        pdf.setTitle(self.documentTitle)

        # registering a external font in python 
        #pdfmetrics.registerFont( 
        #TTFont('abc', 'SakBunderan.ttf') 
        #) 

        # creating the title by setting it's font  
        # and putting it on the canvas 
        #pdf.setFont('abc', 36)
        pdf.setFont("Courier-Bold", 36)
        pdf.drawCentredString(300, 770, self.title)

        # creating the subtitle by setting it's font,  
        # colour and putting it on the canvas 
        pdf.setFillColorRGB(0, 0, 255) 
        pdf.setFont("Courier-Bold", 24) 
        pdf.drawCentredString(290, 720, self.subTitle) 

        # drawing a line
        pdf.line(30, 710, 550, 710)
        # creating a multiline text using
        # textline and for loop 
        text = pdf.beginText(40, 680)
        text.setFont("Courier", 18)
        text.setFillColor(colors.red)

        for line in self.textLines:
            text.textLine(line)
        pdf.drawText(text)
        # drawing a image at the  
        # specified (x.y) position 
        #pdf.drawInlineImage(image, 130, 400) 

        # saving the pdf 
        pdf.save()
        return