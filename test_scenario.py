# -*- coding: UTF-8 -*-
from ctsimu.scenario import Scenario

s = Scenario()
s.read("examples/scenario/example.json")
print(s.json_dict())
s.write("test3.json")