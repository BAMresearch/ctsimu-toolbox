!#/bin/bash
pdoc --force --html --template-dir "docs/templates" -o "docs/" ctsimu
cp -R docs/ctsimu/* docs/
rm -R docs/ctsimu