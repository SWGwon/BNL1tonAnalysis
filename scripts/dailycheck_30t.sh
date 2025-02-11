#!/bin/bash

EXECUTABLE="/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/DailyCheck30T"
INPUT_FILE="majority_test_250211T0722_1.root"
echo "Processing file: $INPUT_FILE"

$EXECUTABLE -i "/Users/gwon/WbLS/30ton_analysis/data/${INPUT_FILE}" -o "/Users/gwon/WbLS/30ton_analysis/data/" -n 10000 -t 2

if [ $? -ne 0 ]; then
    echo "Error processing file: $INPUT_FILE"
fi
