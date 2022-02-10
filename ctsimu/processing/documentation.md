Pipeline Concept
================

![Example for a processing pipeline.](../pictures/pipeline.png "Example for a processing pipeline.")

The required minimum set of parameters are an input file pattern and an output file pattern. As illustrated in the following picture, a file pattern is a string that specifies a single file name (for a single 2D image or a raw chunk of multiple images), or a stack of 2D images using the `%[N]d` placeholder for a number that contains `[N]` digits, e.g. `%4d` for the numbers from 0000 to 9999.

![File stacks and their file patterns.](../pictures/filestacks.png "File stacks and their file patterns.")
