#!/bin/bash
# Enable nullglob so that if no files match the pattern, an empty array is returned.
shopt -s nullglob

LIQUID_CONFIG="07"
INPUT_DIR="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_${LIQUID_CONFIG}0_WbLS"

# Path to the executable
EXECUTABLE="/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/DailyCheck"

# Get TIMESTAMP from the first script argument
TIMESTAMP=$1

# Construct the search pattern for files that contain TIMESTAMP and end with _1.root
pattern="${INPUT_DIR}/*${TIMESTAMP}*_1.root"

# Expand the pattern into an array of matching files
files=($pattern)

# Check if any matching file is found
if [ ${#files[@]} -eq 0 ]; then
    echo "No file found matching the pattern: ${pattern}"
    exit 1
fi

# Use the first matching file
INPUT_FILE="${files[0]}"
# Extract the filename without the directory and without the .root extension
FILENAME=$(basename "$INPUT_FILE" .root)

# Create the output directory based on the found filename
OUTPUT_DIR="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/${FILENAME}"
mkdir -p "${OUTPUT_DIR}"

echo "Processing file: $INPUT_FILE"

# Execute the processing command
$EXECUTABLE -i "$INPUT_FILE" -o "$OUTPUT_DIR/" -n 10000 -t 2 -c ""

# Check for errors in execution
if [ $? -ne 0 ]; then
    echo "Error processing file: $INPUT_FILE"
fi
