#!/bin/bash

echo "input file path"
read -e INPUT_FILE_PATH
#INPUT_FILE_NAME=$1

#removing '_subrun.root' 
INPUT_FILE_NAME="${INPUT_FILE_PATH%_[0-9][0-9][0-9][0-9].root}"
INPUT_FILE_NAME="${INPUT_FILE_NAME%_[0-9][0-9][0-9].root}"
INPUT_FILE_NAME="${INPUT_FILE_NAME%_[0-9][0-9].root}"
INPUT_FILE_NAME="${INPUT_FILE_NAME%_[0-9].root}"

echo "output file path"
read -e OUTPUT_DIR

echo "what trigger? (0: tp, 1: alpha, 2: majority)"
read T

echo "# events to process"
read EVENT_NUMBER

CALIBRATION_CSV="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/bnl1t_spe_fit_results_241111.csv"

EXECUTABLE="/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/ProcessRawRoot"

for i in {0..100}; do
    INPUT_FILE="${INPUT_FILE_NAME}_${i}.root"
    echo "Processing file: $INPUT_FILE"
    
    $EXECUTABLE -i "$INPUT_FILE" -o "$OUTPUT_DIR/" -n $EVENT_NUMBER -t $T -c $CALIBRATION_CSV

    if [ $? -ne 0 ]; then
        echo "Error processing file: $INPUT_FILE"
    fi
done
