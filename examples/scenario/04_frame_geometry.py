# -*- coding: UTF-8 -*-
# File: examples/scenario/04_frame_geometry.py

from ctsimu.scenario import Scenario
import ctsimu.geometry # for ctsimu_world
s = Scenario("example.json")

for frame_number in range(s.n_frames()):
    # Set frame number:
    s.set_frame(frame=frame_number, reconstruction=False)

    # Calculate projection matrix:
    geo = s.current_geometry()
    pmatrix = geo.projection_matrix(mode="OpenCT")

    # Individual coordinate systems for this frame:
    cs_source   = s.source.coordinate_system
    cs_stage    = s.stage.coordinate_system
    cs_detector = s.detector.coordinate_system

    # First sample in stage coordinate system:
    sample    = s.samples[0]  # get first sample
    cs_sample = sample.coordinate_system.get_copy()

    # Transform sample coordinate system from
    # stage to world coordinates if it is
    # attached to the stage:
    if sample.is_attached_to_stage():
        cs_sample.change_reference_frame(
            cs_from=cs_stage,
            cs_to=ctsimu.geometry.ctsimu_world
        )

    print(f"Frame {frame_number}")
    print("================")
    print("Stage coordinates in world:")
    print(cs_stage)

    print(f"Sample coordinates in world:")
    print(cs_sample)
