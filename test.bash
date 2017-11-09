#!/bin/bash
# 
# Runs tests
pip install -e .[rplot] --user
pytest
rm _*.png # remove images created by some docstrings
