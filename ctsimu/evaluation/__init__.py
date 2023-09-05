# -*- coding: UTF-8 -*-
"""# Evaluations for CTSimU 2D tests

While this sub-module contains the 2D test implementations,
the actual tests are usually run using `ctsimu.toolbox.Toolbox`
commands. This works by passing the test name as a command
to a new `ctsimu.toolbox.Toolbox` object. Some tests, like
the noise test in the example below, specify keywords that
you must use to identify projection images.

```python
from ctsimu.toolbox import Toolbox
Toolbox("2D-FB-1",
  SNR100="2D-FB-1_Detektor1_SNR-100_metadata.json",
  SNR250="2D-FB-1_Detektor2_SNR-250_metadata.json"
)
```

As can be seen from the example, projection images are not passed directly
to the toolbox for test evaluation. Instead, [CTSimU metadata files] must be used.
These files allow the complete description of projection images and
the flat and dark field images that belong to them. For many tests,
the toolbox needs to run a flat-field correction before the actual
test is done. This is the reason why all tests work with metadata
files.

[CTSimU metadata files]: https://bamresearch.github.io/ctsimu-scenarios/metadata.html

For details on how to run the test scenarios using metadata
files for the given projection images to be tested, see the
following documentation for each individual test.
"""

from . import *