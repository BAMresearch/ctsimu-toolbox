# -*- coding: UTF-8 -*-
# File: examples/scenario/03_set_frame.py

from ctsimu.scenario import Scenario
s = Scenario("example.json")

s.set_frame(frame=10, reconstruction=False)
geo = s.current_geometry()

print("Stage coordinate system for frame 10:")
print(geo.stage)

"""
Stage coordinate system for frame 10:
Center: [275.   0.   0.]
u:      [-0.93896467  0.34045905 -0.04932528]
v:      [-0.34274809 -0.93813153  0.04932528]
w:      [-0.02948036  0.06322084  0.99756405]
"""