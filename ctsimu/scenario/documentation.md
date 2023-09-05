# Reading and writing scenario files

To read a scenario file, you can pass it to the constructor of
a `ctsimu.scenario.Scenario` object:

```python
.. include:: ../../examples/scenario/01_read_scenario.py
```

You can also employ the function `ctsimu.scenario.Scenario.read()`:

```python
from ctsimu.scenario import Scenario
s = Scenario()
s.read("example.json")
```

To write a scenario file, use the `ctsimu.scenario.Scenario.write()` function.

```python
from ctsimu.scenario import Scenario
s = Scenario()
s.read("example.json")
s.write("example_new.json")
```

Scenarios are always written in the currently supported file format version
of the [CTSimU scenario file format]. Any additional information from
the originally imported scenario file that was stored in non-standard
parameters, such as comments, is lost. However, the `"simulation"` section
of the scenario file, which stores proprietary, simulation-software specific
parameters, is completely preserved.

[CTSimU scenario file format]: https://bamresearch.github.io/ctsimu-scenarios/

# Manipulating scenarios

## Parameter getters and setters

A `ctsimu.scenario.Scenario` object provides the same sub-object
structure as the JSON structure of a [CTSimU scenario].
[CTSimU scenario]: https://bamresearch.github.io/ctsimu-scenarios/

In the toolbox, each parameter of the scenario is an object of the class
`ctsimu.scenario.parameter.Parameter`. It provides a function
`ctsimu.scenario.parameter.Parameter.get()` to get its current value and a
function `ctsimu.scenario.parameter.Parameter.set()` to set its standard and
current value.

The **current value** is the value that the parameter has in the frame that
is currently set (see below). Because parameter values can drift throughout
a CT scan, the current value of a parameter depends on the current frame
number. The **standard value** is the parameter value without any drifts
considered.

In the following example, we get the value of the environment temperature and
then set it to a new value.

```python
.. include:: ../../examples/scenario/02_getters_and_setters.py
```

## Units

Within the toolbox, parameter values are always in their **native units,**
which means the getter functions return values in native units and the
setter functions expect a value in the parameter's native unit. The following
table gives an overview of the native units used in the toolbox.

| Quantity            | Native unit                                          |
| :------------------ | :--------------------------------------------------- |
| dimensionless       | `None`                                               |
| Length              | `"mm"`                                               |
| Angle               | `"rad"` by default,                                  |
|                     | `"deg"` for `start_angle` and `stop_angle` of stage  |
| Time                | `"s"`                                                |
| Current             | `"mA`                                                |
| Voltage             | `"kV`                                                |
| Mass density        | `"g/cm^3"`                                           |
| Spatial frequency   | `"lp/mm"`                                            |
| Angular velocity    | `"deg/s"`                                            |
| Temperature         | `"C"`                                                |
| Boolean             | `"bool"`                                             |
| String              | `"string"`                                           |

There are also further *dummy units* which the toolbox accepts in scenario
descriptions, but handles them like no unit at all (dimensionless):
`"px"`, `"1/J"`, `"relative"`.

Additionally, each parameter has a **preferred unit.** This is the parameter's
unit used in the JSON scenario file. When writing a scenario file, the
value is converted from the internal native unit to the preferred unit.

# Frames

## Setting the frame

A scenario defines a CT trajectory, typically a rotation of the sample stage.
While the motion takes place, projection images are taken: each projection
image represents a scenario **frame** with its specific CT geometry.

To set the frame number for the scenario, use the function
`ctsimu.scenario.Scenario.set_frame()`. The function accepts the frame number
(starting at zero) and as second argument a Boolean value that decides if the
geometry should be set up as seen by the reconstruction software (`True`):
any drifts and deviations that are not `known_to_reconstruction` will be ignored
in this case.

```python
.. include:: ../../examples/scenario/03_set_frame.py
```

## Geometry of the current frame

After the scenario frame is set, the current `ctsimu.geometry.Geometry` can
be accessed through the function `ctsimu.scenario.Scenario.current_geometry()`.
This can be helpful if you want to calculate projection matrices for each frame
or use the coordinate systems of the individual parts for something else.
The function `ctsimu.scenario.Scenario.n_frames()` returns
the number of frames defined in the scenario.

The following example shows how to iterate through the frames of a scenario,
calculate a projection matrix for each frame, and access the coordinate
systems of the scenario's parts (X-ray source, stage, detector and samples).

It also shows how to transform the sample's coordinate system into world
coordinates. Sample coordinates are given in terms of the stage coordinate
system if the sample is attached to the stage (the usual behaviour). Here,
the method `ctsimu.geometry.CoordinateSystem.change_reference_frame()` from
the geometry sub-module is used to transform the sample coordinate system.

```python
.. include:: ../../examples/scenario/04_frame_geometry.py
```

# Reconstruction config files

## OpenCT & CERA

The toolbox can generate reconstruction configuration files for SIEMENS CERA
and in the OpenCT format (used by VGSTUDIO) directly from a scenario file.
The simplest way is to use the functions
`ctsimu.scenario.Scenario.write_CERA_config()` and
`ctsimu.scenario.Scenario.write_OpenCT_config()`:

```python
.. include:: ../../examples/scenario/05_reconstruction_configs.py
```

For a detailed explanation of the parameters, please see the documentation
of the two functions.

## Projection and volume parameters

Reconstructions need to read in projection images and will create a tomogram
volume as output. In the example above, information about the projection images
and tomogram volume is generated automatically. For example, the number of
voxels in the tomogram and the voxel size are calculated automatically from the
magnification and detector geometry.

If you want to change these parameters, you must change the scenario's
**metadata** before you create the config files. If you have a [metadata file]
that describes the projections and tomogram volume, you can import this file
directly after reading the scenario:

```python
from ctsimu.scenario import Scenario
s = Scenario("example.json")
s.read_metadata("recon_metadata.json")
```

It is also possible to set up the scenario's metadata manually by following
the same structure as defined in the [metadata file]. In the following example,
the projection image files are assumed to be sequentially numbered (four digits),
starting at zero.

```python
.. include:: ../../examples/scenario/06_reconstruction_metadata.py
```

[metadata file]: https://bamresearch.github.io/ctsimu-scenarios/metadata.html