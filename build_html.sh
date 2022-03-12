#!/bin/bash
# Builds the HTML documentation using pdoc3
# and moves it to the "docs" folder (instead of docs/ctsimu).
pdoc --force --html --template-dir "docs/templates" -o "docs/" ctsimu
cp -R docs/ctsimu/* docs/
rm -R docs/ctsimu