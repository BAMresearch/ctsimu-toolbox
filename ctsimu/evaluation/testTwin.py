# -*- coding: UTF-8 -*-
"""# Test of digital model

.. include:: ./testTwin.md
"""

import io
import pandas as pd
pd.options.mode.copy_on_write = True
from ..test import *
from ..helpers import *
import numpy as np
from reportlab.lib import colors 
import matplotlib.pyplot as plt
from ..version import *

table_ids = ['values-real', 'values-sim', 'reference-real', 'reference-sim']

ctsimu_twin_test_supported_version = {
    "major": 0,
    "minor": 4
}


class testTwin():
    """ CTSimU2 test of digital twin. """

    def __init__(self, filename:str=None):
        # empty declarations
        self.Criterion_real = []
        self.Criterion_sim = []
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
            self.name = self.current_metadata_basename

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
        self.measurands = self.metadata['measurement']['measurands']
        self.output_path = self.metadata['measurement']['output_path']
        self.sep = self.metadata['measurement']['csv_sep']
        self.decimal = self.metadata['measurement']['csv_decimal']
        self.alpha = self.metadata['measurement']['material_alpha']


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
                if len(self.measurands) == 0: # take all calibrated quantities of the test piece
                    self.measurands = data.index
                elif len([m for m in self.measurands if '*' in m]) > 0: # wildcard patterns only
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
        scaling_factor = self.metadata[ct_type]["scaling_factor"]

        # Choose calibration
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

                    # Filter the DataFrame to include only the measurands of interest
                    df_filtered = df[df[name_col].isin(self.measurands)]

                    df_filtered[value_col] = df_filtered[value_col].apply(convert_str_to_float)
                    #print(df_filtered)
                    if not scaling_factor == 1:
                        print(f'  scaling: {scaling_factor}')
                        df_filtered[value_col] = df_filtered[value_col] * scaling_factor
                    if not temperature == ref_temperature:
                        print(f'  temp. compensation: {temperature}, {ref_temperature}, {self.alpha}')
                        df_filtered[value_col] = df_filtered[value_col].apply(self.TemperatureCorrection, args=(float(temperature), float(ref_temperature), float(self.alpha)))
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

            if ct_type == "values-sim":
                self.files_sim = file_names
            else:
                self.files_real = file_names

            combined_df_subtracted = combined_df_test.apply(lambda x: self.Deviation(x, ref_values), axis=1)

            print(combined_df_subtracted)

            os.makedirs(self.output_path, exist_ok=True)
            combined_df_test.to_csv(f"{self.output_path}/{self.name}_{ct_type}.csv", sep=self.sep, decimal=self.decimal, index=True)

        return combined_df_subtracted


    def TemperatureCorrection(self, Mvalue, Tmeas, Tcal, alpha):
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
        self.Real_avg = RealValues.mean(axis=1)
        self.Sim_avg = SimValues.mean(axis=1)

        u_psim = SimValues.std(axis=1)
        k_sim = get_value(self.metadata, ["values-sim", "k_value"], 2.0)
        u_cal = self.RealRefUncertainty
        u_p = RealValues.std(axis=1)
        u_ab = 0.2 * self.alpha
        u_b = (T_real - T_ref) * u_ab * self.Real_avg
        #print('u_b: ',u_b)
        k_real = get_value(self.metadata, ["values-real", "k_value"], 2.0)
        self.U_real = k_real * (u_cal.pow(2) + u_p.pow(2) + u_b.pow(2)).pow(1/2)
        self.U_sim = k_sim * u_psim
        self.E_DM = ((self.Sim_avg) - (self.Real_avg)) * (self.U_sim.pow(2) + self.U_real.pow(2)).pow(-0.5)
        #print(self.E_DM)

        self.Criterion_real = u_cal / u_p
        self.Criterion_sim = u_psim / u_p
        # self.Result = np.where((abs(self.E_DM) < 1) & (self.Criterion_real < 1) & (self.Criterion_sim < 1), True, False)
        self.Result = np.where((abs(self.E_DM) < 1) & (self.Criterion_real + self.Criterion_sim < 2), True, False)

        self.df["Real_avg"] = self.Real_avg
        self.df["u_cal"] = u_cal
        self.df["u_p"] = u_p
        self.df["u_b"] = u_b
        self.df["Real_U"] = self.U_real
        self.df["Sim_avg"] = self.Sim_avg
        self.df["Sim_U"] = self.U_sim
        self.df["E_DM"] = self.E_DM
        self.df["Criterion_real"] = self.Criterion_real
        self.df["Criterion_sim"] = self.Criterion_sim
        self.df["Criterion_sum"] = self.Criterion_real + self.Criterion_sim
        # self.df["Test Result"] = np.where((abs(self.df["E_DM"])< 1) & (self.df["Criterion_real"] < 1) & (self.df["Criterion_sim"] < 1), True, False)
        self.df["Test Result"] = np.where((abs(self.df["E_DM"])< 1) & (self.df["Criterion_real"] + self.df["Criterion_sim"] < 2), True, False)
        self.df.index.name = 'Measurand'
        # self.df.reset_index()
        print(self.df)

        # Save results.
        self.df.to_csv(f"{self.output_path}/{self.name}.csv", sep=self.sep, decimal=self.decimal, index=True, index_label='Measurand')


    def plotDeviations(self):
        import matplotlib.transforms
        from matplotlib.transforms import offset_copy

        fig = plt.figure(figsize=(7, 5))

        # Add horizontal lines at 0
        plt.axhline(y=0, color='lightgray')

        # plt.plot(self.df['Real_avg'], marker='o')
        plt.errorbar(self.df.index, self.df['Real_avg'], yerr=self.df['Real_U'], fmt='o')
        plt.errorbar([x+.15 for x in range(len(self.df.index))], self.df['Sim_avg'], yerr=self.df['Sim_U'], fmt='o')

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
        plt.savefig(f"{self.output_path}/{self.name}_deviations.png")


    def plotResults(self):

        # assign categories and colormap
        cat = np.where(self.df['Criterion_sum'] < 2.0, 0, 1)
        col_map = np.array(['C0', 'C3'])
        # print(col_map[cat])

        plt.figure(figsize=(10, 6))

        # Add horizontal lines at -1 and +1
        plt.axhline(y=-1, color='lightgray')
        plt.axhline(y=1, color='lightgray')

        plt.scatter(self.df.index, self.df["E_DM"], c=col_map[cat])

        # Label the x-axis with names from self.measurands
        xtick_labels = [s[:7] for s in self.measurands]
        plt.xticks(ticks=range(len(xtick_labels)), labels=xtick_labels, rotation=90)
        # print(plt.xticks())

        # Add labels and title
        # plt.xlabel('Measurand')
        plt.ylabel('$E_{DM}$ value')
        # plt.title('CTSimU2 Test Result')

        # plt.savefig(self.img_buf, format='png')

        # Show the plot
        #plt.show()
        plt.savefig(f"{self.output_path}/{self.name}_E_DM.png")


    def plotResult(self, wide = False):
        fig, (ax1, ax2) = plt.subplots(2, height_ratios=[2.25,1], sharex=True, dpi=120)
        ax1.set_title(self.name)
        # title = fig.suptitle(self.name)
        # title.set_position([.5,.9])
        fig.subplots_adjust(hspace = 0.05)
        if wide:
            fig.set_size_inches(18, 6)
            max = 50
            add = '_wide'
            fts = 'smaller'
        else:
            fig.set_size_inches(10, 6)
            max = 30
            add = ''
            fts = None

        # Deviation plot
        #
        # Add horizontal lines at 0
        ax1.axhline(y=0, color='lightgray')
        # Add the plot
        ax1.errorbar([s[:7] for s in self.df.index], self.df['Real_avg']*1e3, yerr=self.df['Real_U']*1e3, fmt='o', label='Real')
        ax1.errorbar([x+.15 for x in range(len(self.df.index))], self.df['Sim_avg']*1e3, yerr=self.df['Sim_U']*1e3, fmt='o', label='Sim')
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
        passed = np.where((self.df["Criterion_real"] + self.df["Criterion_sim"] < 2), self.df["E_DM"], None)
        missed = np.where((self.df["Criterion_real"] + self.df["Criterion_sim"] < 2), None, self.df["E_DM"])
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
        mul = 1 + int(len(self.df.index) / max)
        # print(mul)
        from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
        ax2.xaxis.set_major_locator(MultipleLocator(mul))
        ax2.xaxis.set_minor_locator(MultipleLocator(1))
        plt.xticks(fontsize=fts, ha='right', rotation=30)

        fig.savefig(self.img_buf, format='png')

        # Show the plot
        #plt.show()
        fig.savefig(f"{self.output_path}/{self.name}{add}.png")

    def TwinTest_report(self):
        from datetime import datetime
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import mm
        from reportlab.platypus import Paragraph, PageBreak, Spacer, Image, Preformatted
        from ..report import Report, NumberedCanvas

        now = datetime.now()

        documentTitle = f"{self.output_path}/{self.name}.pdf"
        
        # Create a style sheet for the document
        styles = getSampleStyleSheet()

        # Define custom styles
        title_style = ParagraphStyle('Title', parent=styles['Heading1'], alignment=1, spaceBefore=10, spaceAfter=10)
        subtitle_style = ParagraphStyle('Subtitle', parent=styles['Heading2'], alignment=1, spaceBefore=6, spaceAfter=6)
        body_style = ParagraphStyle('Body', parent=styles['Normal'], fontSize=12, spaceBefore=0, spaceAfter=3)
        body_center_style = ParagraphStyle('Body', parent=styles['Normal'], alignment=1, fontSize=12, spaceBefore=0, spaceAfter=3)
        #print(body_style.parent.listAttrs())

        # Define the content for the document
        title = Paragraph('Digital Model Test Report', title_style)
        #subtitle = Paragraph(self.name, subtitle_style)

        # creating a pdf object
        pdf = Report(documentTitle,
                pagesize = A4,
                showBoundary = 0,
                leftMargin = 27*mm,
                rightMargin = 20*mm,
                topMargin = 25*mm,
                bottomMargin = 25*mm,
                allowSplitting = 1,
                title = 'Digital Twin Test Report',
                creator = 'CTSimU Toolbox - https://github.com/BAMresearch/ctsimu-toolbox')
        pdf.setHeader('CTSimU Toolbox', 'Digital Model Test Report', f'{now:%Y-%m-%d}')
        pdf.setFooter(f'{self.name}.pdf', left2 = f'generated by CTSimU Toolbox {get_version()}', left3 = 'https://github.com/BAMresearch/ctsimu-toolbox')

        # Prepare the data for the table
        df = self.df.get(['Real_avg', 'Real_U', 'Sim_avg', 'Sim_U', "E_DM", 'Criterion_real', 'Criterion_sim', 'Test Result']).reset_index()
        df["Real_avg"] = df['Real_avg'].apply(lambda x: round(x, 4))
        df["Real_U"] = df['Real_U'].apply(lambda x: round(x, 5))
        df["Sim_avg"] = df['Sim_avg'].apply(lambda x: round(x, 4))
        df["Sim_U"] = df['Sim_U'].apply(lambda x: round(x, 5))
        df["E_DM"] = df["E_DM"].apply(lambda x: round(x, 2))
        df['Criterion_real'] = df['Criterion_real'].apply(lambda x: round(x, 2))
        df['Criterion_sim'] = df['Criterion_sim'].apply(lambda x: round(x, 2))

        # Add the content to the PDF document
        elements = [title, 
                    # Paragraph("on", body_center_style), subtitle, 
                    Paragraph(f"generated at {now:%Y-%m-%d %H:%M}", body_center_style),
                    Spacer(1, 24, True),
                    # Paragraph(f"CTSimU Toolbox Version: {get_version()}", body_style), 
                    Paragraph(f'Metadata file: {self.current_metadata_file}', body_style), 
                    Paragraph(f"Description: {get_value(self.metadata, ['measurement','description'],'')}", body_style, bulletText='*'), 
                    Paragraph(f"Contact: {self.metadata['measurement']['contact']}", body_style, bulletText='*'), 
                    Paragraph(f"Experimental CTs: ", body_style, bulletText='*'), 
                    Paragraph(f"Description: {get_value(self.metadata, ['values-real','description'],'')}", body_style, bulletText='    -'), 
                    Paragraph(f"File path: {self.metadata['values-real']['csv_path']}", body_style, bulletText='    -'), 
                    Paragraph(f"Number of scans: {len(self.files_real)}", body_style, bulletText='    -'), 
                    Paragraph(f"Temperature: {self.metadata['values-real']['temperature']}", body_style, bulletText='    -'), 
                    Paragraph(f"Scaling: {self.metadata['values-real']['scaling_factor']}", body_style, bulletText='    -'), 
                    Paragraph(f"Reference: {self.metadata['reference-real']['description']}", body_style, bulletText='    -'), 
                    Paragraph(f"Simulated CTs: ", body_style, bulletText='*'), 
                    Paragraph(f"Description: {get_value(self.metadata, ['values-sim','description'],'')}", body_style, bulletText='    -'), 
                    Paragraph(f"File path: {self.metadata['values-sim']['csv_path']}", body_style, bulletText='    -'), 
                    Paragraph(f"Number of scans: {len(self.files_sim)}", body_style, bulletText='    -'), 
                    #Paragraph(f"Temperature: {self.metadata['values-sim']['temperature']}", body_style, bulletText='    -'), 
                    #Paragraph(f"Scaling: {self.metadata['values-sim']['scaling_factor']}", body_style, bulletText='    -'), 
                    Paragraph(f"Reference: {self.metadata['reference-sim']['description']}", body_style, bulletText='    -'), 
                    Image(f"{self.output_path}/{self.name}.png", 450, 270),
                    # PageBreak(), 
                    # Paragraph('This is some content for the PDF document Page 2. ', body_style),
                    pdf.df2table(df),
                    Spacer(1, 24, True),
                    Paragraph("Full imput metadata", subtitle_style),
                    Paragraph(f'{self.current_metadata_file}:', body_style), 
                    Preformatted(json.dumps(self.metadata, indent=2), styles['Code']),
                    Spacer(1, 24, True),
                    Paragraph("List of output files", subtitle_style),
                    Paragraph(f'{self.name}.csv', body_style), 
                    Paragraph(f'{self.name}.pdf (this document)', body_style), 
                    Paragraph(f'{self.name}.png', body_style), 
                    Paragraph(f'{self.name}_deviations.png', body_style), 
                    Paragraph(f'{self.name}_E_DM.png', body_style), 
                    Paragraph(f'{self.name}_values-real.csv', body_style), 
                    Paragraph(f'{self.name}_values-sim.csv', body_style), 
                    ]
        pdf.build(elements, canvasmaker=NumberedCanvas)


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
        self.plotResult(wide = True)
        self.TwinTest_report()
