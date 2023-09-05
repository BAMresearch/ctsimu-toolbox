# -*- coding: UTF-8 -*-
# File: examples/scenario/05_reconstruction_configs.py

from ctsimu.scenario import Scenario
s = Scenario("example.json")

s.write_CERA_config(
    save_dir="cera_recon",
    basename="example",
    create_vgi=True
)

s.write_OpenCT_config(
    save_dir="openct_recon",
    basename="example",
    create_vgi=True,
    variant="free",
    abspaths=True
)