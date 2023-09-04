# -*- coding: UTF-8 -*-
from ctsimu.scenario import Scenario
s = Scenario("example.json")

# Get current temperature:
T = s.environment.temperature.get()
print(f"Currently, it's {T}°C in the room.")
# Currently, it's 20°C in the room.

# Set new temperature to 23°C:
s.environment.temperature.set(23)

# Get new temperature:
T = s.environment.temperature.get()
print(f"Now, it's {T}°C in the room.")
# Currently, it's 23°C in the room.

# Native unit and preferred unit:

print(s.environment.temperature.native_unit)
print(s.environment.temperature.preferred_unit)