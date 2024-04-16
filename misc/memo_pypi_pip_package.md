### Build
python3 setup.py sdist bdist_wheel

### Check
twine check dist/*

### Upload
python3 -m twine upload  --verbose   dist/<package_new_version>
