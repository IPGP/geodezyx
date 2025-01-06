#!/bin/bash

## go to the doc source and run the script in this folder
# cd docsrc

## create the auto documentation based on the doc
sphinx-apidoc --module-first -f -o source/ ../geodezyx

## build the HTML doc
sphinx-build -b html source build

## compile (is this necessary ???? it create the HTML sub-folder we copy
## after but seems redundant with the build one)
make html

## copy the compiled version to the frontend doc folder
cp -v -p -r build/html/* ../docs

## do not delete the docs folder (like a from scratch update),
## it contains more files than docsrc like the nojekyll (cp is injective-only)
