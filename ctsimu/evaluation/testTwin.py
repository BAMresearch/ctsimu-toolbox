# -*- coding: UTF-8 -*-
"""# Test of digital twin

.. include:: ./testTwin.md
"""

import io
import pandas as pd
pd.options.mode.copy_on_write = True
from ..test import *
from ..helpers import *
#from ..responses.run import *
#from ..responses.measurands import *
from ..primitives import Vector, Polygon
from ..scenario import Scenario
import numpy as np
#from ..responses.measurands import Measurands
#from ..responses import *
from reportlab.pdfgen import canvas 
from reportlab.pdfbase.ttfonts import TTFont 
from reportlab.pdfbase import pdfmetrics 
from reportlab.lib import colors 
from PIL import Image
import matplotlib.pyplot as plt


table_ids = ['values-real', 'values-sim', 'reference-real', 'reference-sim']

ctsimu_twin_test_supported_version = {
    "major": 0,
    "minor": 3
}


class testTwin():
    """ CTSimU2 test of digital twin. """

    def __init__(self, filename:str=None):
        # empty declarations
        self.Zusatz_krit_real = []
        self.Zusatz_krit_sim = []
        self.combined_df = []
        self.df = pd.DataFrame()
        self.img_buf = io.BytesIO()
        # metadatafile --> Json File
        #self.metadatafile = metadatafile
        #self.metadata = read_json_file(self.metadatafile)
        self.metadata = {}
        self.current_metadata_path = None
        self.current_metadata_file = None
        self.current_metadata_basename = None
        self.current_metadata_directory = None
        self.metadata_is_set = False

        if filename is not None:
            self.read_metadata(filename)


    def reset_metadata(self):
        """Reset scenario's metadata information."""
        self.current_metadata_path = None
        self.current_metadata_file = None
        self.current_metadata_basename = None
        self.current_metadata_directory = None
        self.metadata_is_set = False

        # Create new, empty metadata:
        self.metadata = {}


    def read_metadata(self, filename:str=None):
        """Import metadata from a CTSimU Twin Test file or a given
        metadata dictionary.

        Parameters
        ----------
        filename : str
            Path to a CTSimU Twin Test file.

            Default value: `None`
        """
        if filename is not None:
            json_dict = read_json_file(filename=filename)
            self.current_metadata_path = filename
            self.current_metadata_directory = os.path.dirname(filename)
            self.current_metadata_file = os.path.basename(filename)
            self.current_metadata_basename, extension = os.path.splitext(self.current_metadata_file)

        # If a file is read, we want to make sure that it is a valid
        # and supported metadata file:
        if isinstance(json_dict, dict):
            file_type = get_value(json_dict, ["file", "file_type"])
            if file_type != "CTSimU Twin Test":
                raise Exception(f"Invalid metadata structure: the string 'CTSimU Twin Test' was not found in 'file.file_type' in the metadata file {filename}.")

            fileformatversion = get_value(json_dict, ["file", "file_format_version"])
            if not is_version_supported(ctsimu_twin_test_supported_version, fileformatversion):
                raise Exception(f"Unsupported or invalid metadata version. Currently supported: up to {ctsimu_twin_test_supported_version['major']}.{ctsimu_twin_test_supported_version['minor']}.")
        else:
            raise Exception(f"Error when reading the metadata file: {filename}. read_json_file did not return a Python dictionary.")

        self.metadata = json_dict
        self.metadata_is_set = True

        #general Information for the test
        self.name = self.metadata['general']['name']
        self.measurements_of_interest = self.metadata['general']['measurements_of_interest']
        self.output_path = self.metadata['general']['output_path']
        self.alpha = self.metadata['general']['material_alpha']


    def prepare_data(self):
        for table_id in table_ids[::-1]:
            print(table_id)
            folder_path = self.metadata[table_id]['csv_path']
            sep = self.metadata[table_id]['csv_sep']
            decimal = self.metadata[table_id]['csv_decimal']
            name_col = self.metadata[table_id]['name_column']
            value_col = self.metadata[table_id]['value_column']

            if table_id == table_ids[2]: # 'reference-real'
                self.RealRefT = self.metadata[table_id]['temperature_calib']

                data = pd.read_csv(folder_path, sep=sep, decimal=decimal, header=0, index_col=0)
                self.RealValuesRef = data[value_col]
                uncertaintyCol = self.metadata[table_id]['uncertainty_column']
                self.RealRefUncertainty = data[uncertaintyCol]
            elif table_id == table_ids[3]: # 'reference-sim'
                self.SimRefT = self.metadata[table_id]['temperature_calib']

                data = pd.read_csv(folder_path, sep=sep, decimal=decimal, header=0, index_col=0)
                self.SimValuesRef = data[value_col]
                # uncertaintyCol = self.metadata[table_id]['uncertainty_column']
                # self.CalUncertainty = self.data[uncertaintyCol]
            elif table_id == table_ids[0]:
                self.RealValues = self.read_and_filter_csv_files(folder_path, table_id)
            elif table_id == table_ids[1]:
                self.SimValues = self.read_and_filter_csv_files(folder_path, table_id)


    def read_and_filter_csv_files(self, folder_paths, ct_type):
        sep = self.metadata[ct_type]['csv_sep']
        decimal = self.metadata[ct_type]['csv_decimal']
        header_row = self.metadata[ct_type]['header_row']
        name_col = self.metadata[ct_type]['name_column']
        value_col = self.metadata[ct_type]['value_column']
        self.temperature = self.metadata[ct_type]['temperature']
        self.scaling_factor = self.metadata[ct_type]["scaling_factor"]
        
        #Calibration Information
        if ct_type == "values-sim":
            ref_values = self.SimValuesRef
            ref_temperature = self.SimRefT
        else:
            ref_values = self.RealValuesRef
            ref_temperature = self.RealRefT

        combined_df_test = pd.DataFrame()
        file_names = []
        for folder_path in [folder_paths]:
            for filename in os.listdir(folder_path):
                if filename.endswith('.csv'):
                    file_path = os.path.join(folder_path, filename)
                    df = pd.read_csv(file_path, sep=sep, decimal=decimal, header=int(header_row), usecols=[name_col, value_col])
                    # df.to_csv(f'{file_path}.df', sep=sep, decimal=decimal, index=True)

                    # Filter the DataFrame to include only the measurements of interest
                    df_filtered = df[df[name_col].isin(self.measurements_of_interest)]

                    df_filtered[value_col] = df_filtered[value_col].apply(convert_str_to_float)
                    #print(df_filtered)
                    if not self.temperature == ref_temperature:
                        print(self.temperature, ref_temperature, self.alpha)
                        df_filtered[value_col] = df_filtered[value_col].apply(self.TempKorr, args=(float(self.temperature), float(ref_temperature), float(self.alpha)))
                    # print(df_filtered)
                    combined_df_test = pd.concat([combined_df_test, df_filtered[value_col]], axis=1)

                    # Append the file name to the list
                    file_names.append(filename)
            #pd.set_option('display.max_rows', None)
            #pd.set_option('display.max_columns', None)
            #pd.set_option('display.width', None)
            #pd.set_option('display.max_colwidth', -1)
            #print(combined_df_test)

            combined_df_test.index = self.measurements_of_interest
            combined_df_test.columns = file_names

            combined_df_subtracted = combined_df_test.apply(lambda x: self.Deviation(x, ref_values), axis=1)

            print(combined_df_subtracted)

            os.makedirs(self.output_path, exist_ok=True)
            combined_df_test.to_csv(f"{self.output_path}/{self.name}_{ct_type}.csv", sep=sep, decimal=decimal, index=True)

        return combined_df_subtracted


    def TempKorr(self, Mvalue, Tmeas, Tcal, alpha):
        return Mvalue*(1-(alpha*(Tmeas-Tcal)))


    def Deviation(self, values, ref_values):
        #print(values.name)
        #print(ref_values)
        #print(values-ref_values[values.name])
        return values-ref_values[values.name]


    def En_calc(self, RealValues, SimValues):

        #print(SimValues)
        #print(RealValues)
        #print(SimValues)
        self.RealValues_avg = RealValues.mean(axis=1)
        self.SimValues_avg = SimValues.mean(axis=1)

        u_sim = SimValues.std(axis=1)
        u_p = RealValues.std(axis=1)
        u_drift = 0
        u_b = 0
        self.MU = 2*(RealValues.std(axis=1).pow(2) + self.RealRefUncertainty.pow(2)).pow(1/2)
        self.MUSim = 2*u_sim

        self.df["RealValues_avg"] = RealValues.mean(axis=1)
        self.df["RealValues_u"] = self.MU
        self.df["SimValues_avg"] = SimValues.mean(axis=1)
        self.df["SimValues_u"] = self.MUSim

        #print(self.SimValues_avg)
        #Berechnung EN-Wert, noch zu verändern.
        self.EN = ((self.SimValues_avg) - (self.RealValues_avg)) * ((self.MUSim.pow(2) + self.MU.pow(2)).pow(1/2)).pow(-1)
        #print(((self.MUSim.pow(2) + self.MU.pow(2))).pow(1/2))
        #print(self.EN)

        self.Zusatz_krit_real = self.MU / self.RealRefUncertainty
        self.Zusatz_krit_sim = self.MUSim / self.MU
        self.df["E_DM-value"] = self.EN
        self.df["Zusatz_krit_real"] = self.Zusatz_krit_real
        self.df["Zusatz_krit_sim"] = self.Zusatz_krit_sim
        self.df["Zusatz_krit_summe"] = self.Zusatz_krit_real + self.Zusatz_krit_sim
        self.df["Test_Result"] = self.EN
        self.df["Test_Result"] = np.where((abs(self.df["E_DM-value"])< 1) & (self.df["Zusatz_krit_real"] < 1) & (self.df["Zusatz_krit_sim"] < 1), True, False)
        print(self.df)

        # Save results.
        sep = self.metadata['general']['csv_sep']
        decimal = self.metadata['general']['csv_decimal']
        self.df.to_csv(f"{self.output_path}/{self.name}.csv", sep=sep, decimal=decimal, index=True)


    def plotResults(self):
        plt.figure(figsize=(10, 6))
        plt.plot(self.df['E_DM-value'], marker='o')

        # Add horizontal lines at -1 and +1
        plt.axhline(y=-1, color='r', linestyle='--')
        plt.axhline(y=1, color='g', linestyle='--')

        # Label the x-axis with names from self.measurements_of_interest
        plt.xticks(ticks=range(len(self.measurements_of_interest)), labels=self.measurements_of_interest, rotation=90)

        # Add labels and title
        plt.xlabel('Measurand')
        plt.ylabel('E_DM value')
        plt.title('CTSimU2 Test Result')

        plt.savefig(self.img_buf, format='png')

        # Show the plot
        #plt.show()
        plt.savefig(f"{self.output_path}/{self.name}.png")


    def TwinTest_report(self):
        from datetime import datetime

        documentTitle = f"{self.output_path}/{self.name}.pdf"
        title = 'Digital Twin Test Report'
        subTitle = 'CTSimU2'
        textLines = [
            f'Metadata file: {self.current_metadata_file}', 
            f'Name: {self.name}', 
            f'Measurands: {", ".join(self.measurements_of_interest)}', 
            ] 
        
        # creating a pdf object
        pdf = canvas.Canvas(documentTitle)
        # setting the title of the document
        #pdf.setTitle(self.documentTitle)

        # registering a external font in python 
        #pdfmetrics.registerFont( 
        #TTFont('abc', 'SakBunderan.ttf') 
        #) 

        # creating the title by setting it's font  
        # and putting it on the canvas 
        pdf.setFont("Helvetica-Bold", 36)
        pdf.drawCentredString(300, 770, title)

        # creating the subtitle by setting it's font,  
        # colour and putting it on the canvas 
        pdf.setFillColorRGB(0, 0, 255) 
        pdf.setFont("Helvetica-Bold", 24) 
        pdf.drawCentredString(290, 720, subTitle) 

        # drawing a line
        pdf.line(30, 710, 550, 710)
        # creating a multiline text using
        # textline and for loop 
        text = pdf.beginText(40, 680)
        text.setFont("Helvetica", 12)
        text.setFillColor(colors.black)

        now = datetime.now()
        text.textLine(f"Date: {now}")
        for line in textLines:
            text.textLine(line)
        pdf.drawText(text)
        # drawing a image at the  
        # specified (x.y) position 
        pdf.drawInlineImage(Image.open(self.img_buf), 0, 0, 600, None, True) 

        # saving the pdf 
        try:
            pdf.save()
        except Exception as e:
            log(f"   Error: {str(e)}")

    def run(self):
        #self.RealValues = self.read_and_filter_csv_files(self.real_folder_path, "values-real")
        #print(self.RealValues)
        self.prepare_data()
            
        #self.SimValues = self.read_and_filter_csv_files(self.sim_folder_path, "values-sim")
        #print(self.SimValues)

        self.En_calc(self.RealValues, self.SimValues)
        
        self.plotResults()
        self.TwinTest_report()
