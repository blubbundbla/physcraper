#!/usr/bin/env bash
# [wf] execute example_run stage

virtualenv venv2 
source venv2/bin/activate

pip install -r ./../../requirements.txt
python ./../../setup.py develop

python ./../../docs/example.py
