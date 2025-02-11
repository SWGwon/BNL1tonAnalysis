#!/bin/bash

EXCUTABLE_DIR="/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build"

#$EXCUTABLE_DIR/SelectCrossingMuon -L $1 -T $2
#$EXCUTABLE_DIR/MakeCrossingMuonRootFile -L $1 -T $2
# python /Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/src/ProcessCrosingMuon.py --Dir /Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_$1/ --TimeStamp $2

CALIBRATION_CSV="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/bnl1t_spe_fit_results_241111.csv"

LIQUIDCONFIG=$1
TIMESTAMP=$2

if [ -z "$LIQUIDCONFIG" ] || [ -z "$TIMESTAMP" ]; then
    echo "use: ./process_crossing_muon \$LIQUIDCONFIG(070_WbLS) \$TIMESTAMP(YYMMDD)"
  exit 1
fi

shopt -s nullglob

#INPUT_DIR="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_${LIQUIDCONFIG}/${PATTERN}"

for i in {0..1}; do
    PATTERN="*${TIMESTAMP}*_${i}.root"
    INPUT_DIR="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_${LIQUIDCONFIG}/${PATTERN}"
    echo $INPUT_DIR

    for INPUT_FILE in $INPUT_DIR; do
        echo "Processing file: ${INPUT_FILE}"
        BASE_NAME=$(basename "${INPUT_FILE}" .root)
        OUTPUT_DIR="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/${BASE_NAME}"
        mkdir -p "${OUTPUT_DIR}"

        $EXCUTABLE_DIR/ProcessRawRoot -i "${INPUT_FILE}" -o "${OUTPUT_DIR}/" -n 10000 -t 3 -c "${CALIBRATION_CSV}"

        if [ $? -ne 0 ]; then
            echo "Error processing file: ${INPUT_FILE}"
        fi
    done
done
#rm ~/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_$1/*$2*root
