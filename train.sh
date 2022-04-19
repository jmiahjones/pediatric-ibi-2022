#!/bin/env bash
nice -n 15 Rscript ./R/train_pipeline_ibi.R > logs/train.log 2>&1
tail logs/train.log
exit 0
