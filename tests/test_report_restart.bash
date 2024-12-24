#!/usr/bin/env bash

rm -R report_tmp
rm civic.json
rm pandrugs2.json

pytest test_report.py -s -v
