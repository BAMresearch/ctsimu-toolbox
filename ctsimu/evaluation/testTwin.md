To run this test evaluation, you can simply pass the name of the test and the metadata file which describes the data files that you want to evaluate to the `ctsimu.toolbox.Toolbox`:

```python
from ctsimu.toolbox import Toolbox
Toolbox("twin", "TwinTest_metadata.json")
```

The metafile for the CTSimU Twin Test has the following structure:

```json
{
    "file": {
        "file_type": "CTSimU Twin Test",                    # required
        "file_format_version": {"major": 0, "minor": 3}     # required
    },

    "general": {
        "name": "SampleData",                               # required
        "description": "DigitalTwinTest on sample data.",
        "contact": "Jane Doe",
        "date_created": "2025-02-04",
        "date_changed": "2025-02-20",

        "sample_material": "Al",
        "material_alpha": 0.0000234,                        # required

        "measurands": ["dia01", "dia02", "dia03", "dia04", "dia05"],
                                                            # required
        "output_path": "./results",
        "csv_sep": ";",
        "csv_decimal": ","
    },

    "values-real": {
        "name": "CT XYZ",
        "nr_of_runs":"",
        "date_measure": "",
        "temperature": "21.3258",                           # required
        "files": "",
        "csv_path": "Daten/Messung175kV/",
        "csv_sep": ";",
        "csv_decimal": ",",
        "header_row": 103,
        "name_column": "Name",
        "value_column": "Ist [mm/°]",
        "scaling_factor": "1"
    },

    "values-sim":{
        "name": "simulation XYZ",
        "nr_of_runs":"",
        "date": "",
        "temperature": "20",
        "files": "",
        "csv_path": "Daten/Simulation175kV/",
        "csv_sep": ";",
        "csv_decimal": ",",
        "header_row": 103,
        "name_column": "Name",
        "value_column": "Ist [mm/°]",
        "scaling_factor": "1"
    },

    "reference-real":{
        "name": "Calibration XYZ",
        "date_calib": "",
        "temperature_calib": "20",
        "csv_path": "Daten/Kalibrierwerte.csv",
        "csv_sep": ";",
        "csv_decimal": ",",
        "header_row": 0,
        "name_column": "",
        "value_column": "Kalibrierwert (mm)",
        "uncertainty_column": "Kalibrierunsicherheit"
    },

    "reference-sim":{
        "name": "Nominal STL",
        "date_calib": "",
        "temperature_calib": "20",
        "csv_path": "Daten/Kalibrierwerte.csv",
        "csv_sep": ";",
        "csv_decimal": ",",
        "header_row": 0,
        "name_column": "",
        "value_column": "Nominal STL",
        "uncertainty_column": ""
    }
}
```
