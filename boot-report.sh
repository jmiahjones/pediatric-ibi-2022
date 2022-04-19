#!/bin/env bash
nice -n 15 Rscript ./R/report_pipeline.R > logs/report.log 2>&1
tail logs/report.log
nice -n 15 Rscript ./R/other_output.R > logs/other.log 2>&1
tail logs/other.log
echo "Reporting completed. Make sure to copy the tmp.csv file into xlsx file!"
exit 0
