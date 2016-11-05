#!/bin/bash
set -euo pipefail
pip install -r requirements_dev.txt
python setup.py install
