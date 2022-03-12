#!/bin/bash
# Builds the documentation as a PDF file,
# using pdoc3 to create a complete Markdown file,
# which is then passed to pandoc to create a LaTeX
# file and a PDF using XeLaTeX.
pdoc --force --pdf ctsimu > docs/ctsimu-toolbox.md

pandoc --metadata=title:"CTSimU Toolbox Documentation" --from=markdown+abbreviations+tex_math_single_backslash --pdf-engine=xelatex --variable=mainfont:"DejaVu Sans" --toc --toc-depth=1 --output=docs/ctsimu-toolbox.pdf docs/ctsimu-toolbox.md

rm docs/ctsimu-toolbox.md