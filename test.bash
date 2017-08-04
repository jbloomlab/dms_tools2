#!/bin/bash
# 
# Runs tests
python setup.py build_ext --inplace
pytest
