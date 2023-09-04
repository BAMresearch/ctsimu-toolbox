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

To write a scenario file, use the `ctsimu.scenario.Scenario.write()` function:

```python
from ctsimu.scenario import Scenario
s = Scenario()
s.read("example.json")
s.write("example_new.json")
```

# Manipulating scenarios

## Parameter getters and setters

A `ctsimu.scenario.Scenario` object provides the same sub-object
structure as the JSON structure of a [CTSimU scenario].
[CTSimU scenario]: https://bamresearch.github.io/ctsimu-scenarios

In the toolbox, each parameter of the scenario is an object of the class
`ctsimu.scenario.parameter.Parameter`. It provides a function
`ctsimu.scenario.parameter.Parameter.get()` to get its current value and a
function `ctsimu.scenario.parameter.Parameter.set()` to set its standard and
current value.

In the following example, we get the value of the environment temperature and
then set it to a new value.

```python
.. include:: ../../examples/scenario/02_getters_and_setters.py
```

## Units

Within the toolbox, parameter values are always in their **native units,**
which means the getter functions return values in native units, and the
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
descriptions, but handles them like no unit at all: `"px"`, `"1/J"`, `"relative"`.

Additionally, each parameter has a **preferred unit.** This is the parameter's
unit used in the JSON scenario file. When writing a scenario file, the
value is converted from the internal native unit to the preferred unit.


# Setting the frame

A scenario defines a CT trajectory, typically a rotation of the sample stage.
While the motion takes place, projection images are taken: each projection
image represents a scenario frame, with its specific CT geometry.

To set the frame number for the scenario, use the function
`ctsimu.scenario.Scenario.set_frame()`.