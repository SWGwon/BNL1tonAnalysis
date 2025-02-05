#!/bin/bash

LIQUID_CONFIG="07"
INPUT_DIR="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_${LIQUID_CONFIG}0_WbLS"

# 실행 파일 경로
EXECUTABLE="/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/DailyCheck"
TIMESTAMP=$1

FILENAME="muon_wbls${LIQUID_CONFIG}_${TIMESTAMP}_1"

OUTPUT_DIR="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/${FILENAME}"
mkdir ${OUTPUT_DIR}

INPUT_FILE="${INPUT_DIR}/${FILENAME}.root"
echo "Processing file: $INPUT_FILE"

$EXECUTABLE -i "$INPUT_FILE" -o "$OUTPUT_DIR/" -n 10000 -t 2
#$EXECUTABLE -i "/Users/gwon/WbLS/30ton_analysis/data/majority_test_250202T0733_0.root" -o "/Users/gwon/WbLS/30ton_analysis/data/" -n 10000 -t 2

if [ $? -ne 0 ]; then
    echo "Error processing file: $INPUT_FILE"
fi
