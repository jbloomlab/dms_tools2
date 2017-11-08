#!/bin/bash
# 
# Runs tests
pip install -e .[rplot] --user
pytest
