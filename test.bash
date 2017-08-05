#!/bin/bash
# 
# Runs tests
python3 setup.py build_ext --inplace
python3 setup.py install --user
pytest
