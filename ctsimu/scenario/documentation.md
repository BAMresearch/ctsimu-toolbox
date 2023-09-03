# Reading scenario files

To read a scenario file, you can pass it to the constructor of
a `ctsimu.scenario.Scenario` object:

```python
from ctsimu.scenario import Scenario
s = Scenario("examples/scenario/example.json")
```

You can also employ the function `ctsimu.scenario.Scenario.read()`:

```python
from ctsimu.scenario import Scenario
s = Scenario()
s.read("examples/scenario/example.json")
```

# Manipulating scenarios

A `ctsimu.scenario.Scenario` object provides the same sub-object
structure as the JSON structure of a [CTSimU scenario].
[CTSimU scenario]: https://bamresearch.github.io/ctsimu-scenarios

Each parameter provides a function `ctsimu.scenario.parameter.Parameter.get()`
to get its current value and a function
`ctsimu.scenario.parameter.Parameter.set()` to set its standard and
current value.

In the following example, we get the value of the environment temperature and
then set it to a new value.

```python
from ctsimu.scenario import Scenario
s = Scenario("examples/scenario/example.json")

# Get current temperature:
T = s.environment.temperature.get()
print(f"Currently, it's {T}°C in the room.")
# Currently, it's 20°C in the room.

# Set new temperature to 23°C:
s.environment.temperature.set(23)

# Get new temperature:
T = s.environment.temperature.get()
print(f"Now, it's {T}°C in the room.")
# Now, it's 23°C in the room.
```