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
from ..version import *


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
        self.measurands = self.metadata['general']['measurands']
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
                # Expand list of measurands
                if len(self.measurands) == 0:
                    self.measurands = data.index
                elif len([m for m in self.measurands if '*' in m]) > 0:
                    mds = []
                    for m in data.index:
                        for pattern in [m for m in self.measurands if '*' in m]:
                            if m.startswith(pattern.replace('*', "")):
                                mds.append(m)
                    self.measurands = mds
                print(self.measurands)
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
        temperature = self.metadata[ct_type]['temperature']
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
                    df_filtered = df[df[name_col].isin(self.measurands)]

                    df_filtered[value_col] = df_filtered[value_col].apply(convert_str_to_float)
                    #print(df_filtered)
                    if not temperature == ref_temperature:
                        print(temperature, ref_temperature, self.alpha)
                        df_filtered[value_col] = df_filtered[value_col].apply(self.TempKorr, args=(float(temperature), float(ref_temperature), float(self.alpha)))
                    # print(df_filtered)
                    combined_df_test = pd.concat([combined_df_test, df_filtered[value_col]], axis=1)

                    # Append the file name to the list
                    file_names.append(filename)
            #pd.set_option('display.max_rows', None)
            #pd.set_option('display.max_columns', None)
            #pd.set_option('display.width', None)
            #pd.set_option('display.max_colwidth', -1)
            #print(combined_df_test)

            combined_df_test.index = self.measurands
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
        return values - ref_values[values.name]


    def Edm_calc(self, RealValues, SimValues):
        T_real = float(self.metadata['values-real']['temperature'])
        T_ref = float(self.RealRefT)

        #print(SimValues)
        #print(RealValues)
        #print(SimValues)
        self.RealValues_avg = RealValues.mean(axis=1)
        self.SimValues_avg = SimValues.mean(axis=1)

        u_psim = SimValues.std(axis=1)
        k_sim = 2.0
        u_cal = self.RealRefUncertainty
        u_p = RealValues.std(axis=1)
        u_ab = 0.2 * self.alpha
        u_b = (T_real - T_ref) * u_ab * self.RealValues_avg
        print('u_b: ',u_b)
        k_real = 2.0
        self.U_real = k_real * (u_cal.pow(2) + u_p.pow(2) + u_b.pow(2)).pow(1/2)
        self.U_sim = k_sim * u_psim
        self.E_DM = ((self.SimValues_avg) - (self.RealValues_avg)) * (self.U_sim.pow(2) + self.U_real.pow(2)).pow(-0.5)
        #print(self.E_DM)

        # self.Zusatz_krit_real = self.U_real / self.RealRefUncertainty
        # self.Zusatz_krit_sim = self.U_sim / self.U_real
        self.Zusatz_krit_real = u_cal / u_p
        self.Zusatz_krit_sim = u_psim / u_p
        # self.Result = np.where((abs(self.E_DM) < 1) & (self.Zusatz_krit_real < 1) & (self.Zusatz_krit_sim < 1), True, False)
        self.Result = np.where((abs(self.E_DM) < 1) & (self.Zusatz_krit_real + self.Zusatz_krit_sim < 2), True, False)

        self.df["RealValues_avg"] = RealValues.mean(axis=1)
        self.df["RealValues_U"] = self.U_real
        self.df["SimValues_avg"] = SimValues.mean(axis=1)
        self.df["SimValues_U"] = self.U_sim
        self.df["E_DM-value"] = self.E_DM
        self.df["Zusatz_krit_real"] = self.Zusatz_krit_real
        self.df["Zusatz_krit_sim"] = self.Zusatz_krit_sim
        self.df["Zusatz_krit_summe"] = self.Zusatz_krit_real + self.Zusatz_krit_sim
        # self.df["Test_Result"] = np.where((abs(self.df["E_DM-value"])< 1) & (self.df["Zusatz_krit_real"] < 1) & (self.df["Zusatz_krit_sim"] < 1), True, False)
        self.df["Test_Result"] = np.where((abs(self.df["E_DM-value"])< 1) & (self.df["Zusatz_krit_real"] + self.df["Zusatz_krit_sim"] < 2), True, False)
        self.df.index.name = 'Masse'
        # self.df.reset_index()
        print(self.df)

        # Save results.
        sep = self.metadata['general']['csv_sep']
        decimal = self.metadata['general']['csv_decimal']
        self.df.to_csv(f"{self.output_path}/{self.name}.csv", sep=sep, decimal=decimal, index=True, index_label='Masse')


    def plotDeviations(self):
        import matplotlib.transforms
        from matplotlib.transforms import offset_copy

        fig = plt.figure(figsize=(7, 5))

        # Add horizontal lines at 0
        plt.axhline(y=0, color='lightgray')

        # plt.plot(self.df['RealValues_avg'], marker='o')
        plt.errorbar(self.df.index, self.df['RealValues_avg'], yerr=self.df['RealValues_U'], fmt='o')
        plt.errorbar([x+.15 for x in range(len(self.df.index))], self.df['SimValues_avg'], yerr=self.df['SimValues_U'], fmt='o')

        # Label the x-axis with names from self.measurands
        xtick_labels = [s[:7] for s in self.measurands]
        plt.xticks(ticks=range(len(xtick_labels)), labels=xtick_labels, rotation=90)

        # Add labels and title
        # plt.xlabel('Measurand')
        plt.ylabel('Deviation / mm')
        # plt.title('CTSimU2 Test Result')

        # plt.savefig(self.img_buf, format='png')

        # Show the plot
        #plt.show()
        plt.savefig(f"{self.output_path}/{self.name}_2.png")


    def plotResults(self):

        # assign categories and colormap
        cat = np.where(self.df['Zusatz_krit_summe'] < 2.0, 0, 1)
        col_map = np.array(['C0', 'C3'])
        # print(col_map[cat])

        plt.figure(figsize=(10, 6))

        # Add horizontal lines at -1 and +1
        plt.axhline(y=-1, color='lightgray')
        plt.axhline(y=1, color='lightgray')

        plt.scatter(self.df.index, self.df['E_DM-value'], c=col_map[cat])

        # Label the x-axis with names from self.measurands
        xtick_labels = [s[:7] for s in self.measurands]
        plt.xticks(ticks=range(len(xtick_labels)), labels=xtick_labels, rotation=90)
        # print(plt.xticks())

        # Add labels and title
        # plt.xlabel('Measurand')
        plt.ylabel('E_DM value')
        # plt.title('CTSimU2 Test Result')

        # plt.savefig(self.img_buf, format='png')

        # Show the plot
        #plt.show()
        plt.savefig(f"{self.output_path}/{self.name}_3.png")


    def plotResult(self):
        fig, (ax1, ax2) = plt.subplots(2, height_ratios=[2.25,1], sharex=True, dpi=120)
        ax1.set_title(self.name)
        # title = fig.suptitle(self.name)
        # title.set_position([.5,.9])
        fig.subplots_adjust(hspace = 0.05)

        # Deviation plot
        #
        # Add horizontal lines at 0
        ax1.axhline(y=0, color='lightgray')
        # Add the plot
        ax1.errorbar([s[:7] for s in self.df.index], self.df['RealValues_avg']*1e3, yerr=self.df['RealValues_U']*1e3, fmt='o', label='Real')
        ax1.errorbar([x+.15 for x in range(len(self.df.index))], self.df['SimValues_avg']*1e3, yerr=self.df['SimValues_U']*1e3, fmt='o', label='Sim')
        # ax1.set_ylim(-10,10)
        # Add labels and title
        ax1.set_axisbelow(False)
        ax1.tick_params(which='both', direction='in', top=True, right=True)
        # ax1.set_xlabel('Measurand')
        # ax1.title('CTSimU2 Test Result')
        ax1.set_ylabel('Deviation / µm')
        ax1.legend()

        # E_DM plot
        #
        # Add horizontal lines at -1 and +1
        ax2.axhline(y=-1, color='lightgray', zorder=0)
        ax2.axhline(y=1, color='lightgray', zorder=0)
        # Add the plot
        passed = np.where((self.df["Zusatz_krit_real"] + self.df["Zusatz_krit_sim"] < 2), self.df["E_DM-value"], None)
        missed = np.where((self.df["Zusatz_krit_real"] + self.df["Zusatz_krit_sim"] < 2), None, self.df["E_DM-value"])
        ax2.scatter(range(len(passed)), passed, label='passed exta criteria')
        ax2.scatter(range(len(missed)), missed, c='C3', label='missed exta criteria')
        ax2.set_ylim(-3,3)
        # Add labels and title
        ax2.set_axisbelow(False)
        ax2.tick_params(which='both', direction='in', top=True, right=True)
        # ax1.xlabel('Measurand')
        # ax1.title('CTSimU2 Test Result')
        ax2.set_ylabel('$E_{DM}$')
        # xtick_labels = [s[:7] for s in self.measurands]
        # ax2.xticks(ticks=range(len(xtick_labels)), labels=xtick_labels, rotation=90)
        ax2.legend()

        # X ticks and labels
        mul = 1 + int(len(self.df.index) / 16)
        # print(mul)
        from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
        ax2.xaxis.set_major_locator(MultipleLocator(mul))
        ax2.xaxis.set_minor_locator(MultipleLocator(1))
        plt.xticks(ha='right', rotation=30)

        fig.savefig(self.img_buf, format='png')

        # Show the plot
        #plt.show()
        fig.savefig(f"{self.output_path}/{self.name}.png")


    def TwinTest_report(self):
        from datetime import datetime

        documentTitle = f"{self.output_path}/{self.name}.pdf"
        title = 'Digital Twin Test Report'
        subTitle = 'CTSimU2'
        textLines = [
            f'Metadata file: {self.current_metadata_file}', 
            f'Name: {self.name}', 
            f'Measurands: {", ".join(self.measurands)}', 
            ] 
        
        # creating a pdf object
        pdf = canvas.Canvas(documentTitle)
        pdf.setTitle('Digital Twin Test Report')
        pdf.setCreator('CTSimU Toolbox - https://github.com/BAMresearch/ctsimu-toolbox')

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
        text.textLine(f"Date: {now:%Y-%m-%d %H:%M}")
        text.textLine(f"CTSimU Toolbox Version: {get_version()}")
        for line in textLines:
            text.textLine(line)
        pdf.drawText(text)
        # drawing a image at the  
        # specified (x.y) position 
        pdf.drawInlineImage(Image.open(self.img_buf), 30, 0, 550, None, True) 

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

        self.Edm_calc(self.RealValues, self.SimValues)

        self.plotDeviations()
        self.plotResults()
        self.plotResult()
        self.TwinTest_report()
