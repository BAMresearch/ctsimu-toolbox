#!/bin/bash
# Builds the HTML documentation using pdoc3
# and moves it to the "docs" folder (instead of docs/ctsimu).
pdoc -c latex_math=True --force --html --template-dir "docs/templates" -o "docs/" ctsimu
cp -R docs/ctsimu/* docs/
rm -R docs/ctsimu