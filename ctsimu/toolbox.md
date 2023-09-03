[CTSimU metadata files] contain information about the projection images of a
CT scan and about the tomogram. They reference a [CTSimU scenario file] which
describes the acquisition parameters and scanner geometry, and they may
also contain information about a simulation (if the projections originate
from a CT simulation).

Metadata files can be passed to the toolbox along with single commands, such
as the instruction to run a flat-field correction on all projection images or
to create reconstruction configuration files for OpenCT or SIEMENS CERA.

The following listing shows an example for a metadata file. The meanings of
the parameters are described in the documentation for [CTSimU metadata files].

[CTSimU metadata files]: https://bamresearch.github.io/ctsimu-scenarios/metadata.html
[CTSimU scenario file]: https://bamresearch.github.io/ctsimu-scenarios/index.html

```json
{
    "file":                 {
        "name":                "example",
        "description":         "Tetrahedron in a rigid frame.",
        "contact":             "Jane Doe",
        "date_created":        "2023-09-03",
        "date_changed":        "2023-09-03",
        "file_type":           "CTSimU Metadata",
        "file_format_version": {"major": 1, "minor": 2}
    },
    "output":               {
        "system":        "",
        "date_measured": "2023-09-03",
        "projections":   {
            "filename":      "../projections/example_%04d.raw",
            "number":        300,
            "frame_average": 3,
            "max_intensity": 42000,
            "datatype":      "uint16",
            "byteorder":     "little",
            "headersize":    {"file":  0, "image": 0},
            "dimensions":    {
                "x": {"value": 200, "unit":  "px"},
                "y": {"value": 150, "unit":  "px"}
            },
            "pixelsize":     {
                "x": {"value": 1.3, "unit":  "mm"},
                "y": {"value": 1.3, "unit":  "mm"}
            },
            "dark_field":    {
                "number":        1,
                "frame_average": 20,
                "filename":      "../projections/example_dark.raw",
                "projections_corrected": false
            },
            "flat_field":    {
                "number":        3,
                "frame_average": 20,
                "filename":      "../projections/example_flat_%04d.raw",
                "projections_corrected": false
            },
            "bad_pixel_map": {
                "filename": null,
                "projections_corrected": false
            }
        },
        "tomogram":      {
            "filename":   "example_reconstruction.raw",
            "datatype":   "float32",
            "byteorder":  "little",
            "headersize": {"file":  0, "image": 0},
            "dimensions": {
                "x": {"value": 200, "unit":  "px"},
                "y": {"value": 200, "unit":  "px"},
                "z": {"value": 150, "unit":  "px"}
            },
            "voxelsize":  {
                "x": {"value": 0.9, "unit":  "mm"},
                "y": {"value": 0.9, "unit":  "mm"},
                "z": {"value": 0.9, "unit":  "mm"}
            }
        }
    },
    "acquisition_geometry": {
        "path_to_CTSimU_JSON": "../../example.json"
    }
}
```

# Flat-field correction

Metadata files can reference dark-field images and flat-field (bright,
free-beam) images. These images can be used by the toolbox to run a dark-field and
flat-field correction on the projection images. The corrected images
will be placed in a sub-folder called `corrected` at the location of
the uncorrected projection images.

To run corrections, pass the command `correction` to the constructor of
a new toolbox object, followed by one or more metadata files:

```python
from ctsimu.toolbox import Toolbox
Toolbox("correction", "example_metadata.json")
```

After the flat-field division, the projection images will be rescaled to the
value of `"max_intensity"` defined in the metadata file. If this value is not
defined, a default value of `60000` is assumed. Alternatively, the rescale
factor can be set manually, as well as an additional offset that is added
after the correction and gray value rescaling:

```python
from ctsimu.toolbox import Toolbox
Toolbox("correction",
    "example_metadata.json",
    rescaleFactor=42000,
    offsetAfterRescale=1500
)
```

Note that it is possible to run the flat-field correction for multiple
metadata files, for example located in sub-directories. The keyword argument
`overwrite` determines if existing corrected projection files will be overwritten
(default: `True`).

```python
from ctsimu.toolbox import Toolbox
Toolbox("correction",
    "run001/example_metadata_run001.json",
    "run002/example_metadata_run002.json",
    "run003/example_metadata_run003.json",
    overwrite=False
)
```

For more information about the keyword arguments, see the documentation
of the function `ctsimu.toolbox.Toolbox.correction()`

# Reconstruction configurations

To create reconstruction configuration files for OpenCT and SIEMENS CERA,
the command `"recon_config"` along with one or more reconstruction metadata
files can be passed to the toolbox.

By default, the Toolbox will calculate a projection matrix for each frame. These
matrices will be included in the reconstruction configuration files. This means
that, in principle, the reconstruction of free-trajectory scans is possible.

The reconstruction metadata file must reference a valid [CTSimU scenario file]
which describes the full acquisition geometry. It should also contain information
about the tomogram to be created (number of voxels, voxel size). If tomogram
information is missing, the tomogram size will be determined automatically
from the scenario's detector size and magnification.

```python
from ctsimu.toolbox import Toolbox
Toolbox("recon_config", "reconstruction/metadata.json")
```

Existing configuration files will not be overwritten unless this is enforced:

```python
from ctsimu.toolbox import Toolbox
Toolbox("recon_config", "reconstruction/metadata.json", overwrite=True)
```

To select which files should be created, the keyword arguments `openct`
and `cera` can either be set to `True` or `False`:

```python
from ctsimu.toolbox import Toolbox
Toolbox("recon_config", "reconstruction/metadata.json", openct=True, cera=True)
```

By default, VGI files are created which reference the volume output file from
the reconstruction software and which can be opened in VGSTUDIO. To deactivate
the creation of VGI files, the keyword argument `create_vgi` can be set to
`False`.

```python
from ctsimu.toolbox import Toolbox
Toolbox("recon_config", "reconstruction/metadata.json", create_vgi=False)
```

Under some conditions, OpenCT files must include absolute paths to projection
files (instead of relative paths). This can be activated with the keyword argument
`openct_abspaths`. To create an OpenCT file of the circular trajectory variant
(instead of the default free trajectory), `openct_variant` can be set to `"circular"`.
The circular variant will not include any projection matrices. It is only
meant for ideal circular trajectory scans with the rotation axis in the
center of the detector and no tilts.

```python
from ctsimu.toolbox import Toolbox
Toolbox("recon_config",
    "reconstruction/metadata.json",
    openct_abspaths=True,
    openct_variant="circular"
)
```

For more information about the keyword arguments, see the documentation
of the function `ctsimu.toolbox.Toolbox.recon_config()`


# Scenario standardization

Scenario files can be updated to the latest version of the [CTSimU file format]
using the `standardize` command. The function accepts one or more CTSimU scenario
description JSON files:

```python
from ctsimu.toolbox import Toolbox
Toolbox("standardize", "scenario.json")
```

The previous scenario files will be backed up (i.e., renamed to `{name}.json.old`).
If a backup file already exists, the standardization will be skipped. It is
therefore not possible to overwrite backed up scenario files.

[CTSimU file format]: https://bamresearch.github.io/ctsimu-scenarios/index.html

The file format version that is currently supported by the toolbox can be
found in the `ctsimu.helpers` module in the dictionaries
`ctsimu_supported_scenario_version` and `ctsimu_supported_metadata_version`:

```python
# Find supported file format versions
# -------------------------------------
import ctsimu.helpers

# Scenarios:
print("Supported scenario file format:")
print(ctsimu.helpers.ctsimu_supported_scenario_version)

# Metadata files:
print("Supported metadata file format:")
print(ctsimu.helpers.ctsimu_supported_metadata_version)
```

# Recursive post-processing

The toolbox can recursively scan complete folders for scenario and metadata
files to run flat-field correction, standardize scenario files and create
reconstruction config files automatically for a collection of folders.

The command `"post-processing"` is used for recursive post-processing, followed
by one or more folder paths which should be scanned and processed.

The following example script runs recursive post-processing in the script's
current working directory (`"."`): projection images will be corrected and
reconstruction config files will be created, existing files will not be
overwritten. Scenario file standardization is turned off.

```python
from ctsimu.toolbox import Toolbox
Toolbox("post-processing",
    ".",
    correction=True,
    recon_config=True,
        cera=True,
        openct=True,
            openct_variant='free',
            openct_abspaths=False,
    standardize=False,
    overwrite=False
)
```

More information about the keyword arguments and more available options
for flat-field correction can be found in the documentation of the function
`ctsimu.toolbox.Toolbox.post_processing()`.