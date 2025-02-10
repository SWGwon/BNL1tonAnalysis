#!/bin/bash
source /Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/setup.sh
source /Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/venv/bin/activate


EXCUTABLE_DIR="/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build"

#$EXCUTABLE_DIR/SelectCrossingMuon -L $1 -T $2
#$EXCUTABLE_DIR/MakeCrossingMuonRootFile -L $1 -T $2
# python /Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/src/ProcessCrosingMuon.py --Dir /Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_$1/ --TimeStamp $2

CALIBRATION_CSV="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/bnl1t_spe_fit_results_241111.csv"

shopt -s nullglob

for i in {0..200}; do
    PATTERN="*${2}*_${i}.root"
    echo "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_${1}/${PATTERN}"

    for INPUT_FILE in /Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_${1}/${PATTERN}; do
        echo "Processing file: ${INPUT_FILE}"
        BASE_NAME=$(basename "${INPUT_FILE}" .root)
        OUTPUT_DIR="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/${BASE_NAME}"
        mkdir -p "${OUTPUT_DIR}"

        $EXCUTABLE_DIR/ProcessRawRoot -i "${INPUT_FILE}" -o "${OUTPUT_DIR}/" -n 10000 -t 1 -c "${CALIBRATION_CSV}"

        if [ $? -ne 0 ]; then
            echo "Error processing file: ${INPUT_FILE}"
        fi
    done
done
#rm ~/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_$1/*$2*root
