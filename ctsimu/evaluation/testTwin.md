To run this test evaluation, you can simply pass the name of the test and the metadata file which describes the data files that you want to evaluate to the `ctsimu.toolbox.Toolbox`:

```python
from ctsimu.toolbox import Toolbox
Toolbox("twin", "TwinTest_metadata.json")
```

The metafile for the CTSimU Twin Test has the following structure:

```json
{
    "file": {
        "file_type": "CTSimU Twin Test",                   # do not change
        "file_format_version": {"major": 0, "minor": 4}    # do not change
    },

    "measurement": {
        "description": "DigitalTwinTest on sample data using Lochplatte.",  # in report
        "contact": "Jane Doe",                             # in report

        "sample_material": "Al",                           # not used, for information
        "material_alpha": 0.0000234,                       # required, thermal expansion coefficient

        "measurands": ["dia01", "dia02", "dia03", "dia04", "dia05"],
                                                           # required, column names
        "output_path": "./results",                        # required
        "csv_sep": ";",                                    # required, column separator
        "csv_decimal": ","                                 # required, decimal mark
    },

    "values-real": {
        "description": "CT XYZ",                           # in report
        "nr_of_runs":"",                                   # not used
        "files": "",                                       # not used, instead take all CSVs at 'csv_path'
        "csv_path": "Data/Measurements/",                  # required
        "csv_sep": ";",                                    # required, column separator
        "csv_decimal": ",",                                # required, decimal mark
        "header_row": 103,                                 # required
        "name_column": "Name",                             # required
        "value_column": "Ist [mm/°]",                      # required
        "k_value": 2,                                      # required
        "temperature": 21.3,                               # required
        "scaling_factor": 1.0                              # required
    },

    "values-sim":{
        "description": "simulation XYZ",                   # in report
        "nr_of_runs":"",                                   # not used
        "files": "",                                       # not used, instead take all CSVs at 'csv_path'
        "csv_path": "Data/Simulations/",                   # required
        "csv_sep": ";",                                    # required, column separator
        "csv_decimal": ",",                                # required, decimal mark
        "header_row": 103,                                 # required
        "name_column": "Name",                             # required
        "value_column": "Ist [mm/°]",                      # required
        "temperature": 20,                                 # required
        "scaling_factor": 1                                # required
    },

    "reference-real":{
        "description": "Calibration XYZ",                  # in report
        "temperature_calib": 20,                           # required
        "csv_path": "Daten/Calibration.csv",               # required
        "csv_sep": ";",                                    # required, column separator
        "csv_decimal": ",",                                # required, decimal mark
        "name_column": "",                                 # required
        "value_column": "Kalibrierwert (mm)",              # required
        "uncertainty_column": "Kalibrierunsicherheit"      # required
    },

    "reference-sim":{
      "description": "Calibration XYZ",                    # in report
      "temperature_calib": 20,                             # required
      "csv_path": "Daten/Calibration.csv",                 # required
      "csv_sep": ";",                                      # required, column separator
      "csv_decimal": ",",                                  # required, decimal mark
      "name_column": "",                                   # required
      "value_column": "Nominal STL",                       # required
      "uncertainty_column": ""                             # required
    }
}
```

> The "measurands" are defined as a list of identifiers. An entry can end with a placeholder '*', which activates all identifiers that begin with the character string before the placeholder. If there is an entry with a placeholder in the list, entries without a placeholder are ignored.