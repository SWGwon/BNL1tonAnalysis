#!/bin/bash
# Enable nullglob: if no files match the glob pattern, an empty array is returned.
shopt -s nullglob

# -------------------------------
# Function: process_root_file
# Input: a .root file path
# Process: extracts the file name, creates an output directory, and calls the executable.
# -------------------------------
process_root_file() {
    local file_path="$1"
    # Remove the .root extension to extract the file name (e.g., muon_wbls07_250204T0936_1)
    local file_name
    file_name=$(basename "$file_path" .root)

    echo "Processing file: $file_path"

    # Define the output directory (adjust the path as needed)
    local output_dir="$OUTPUT_FILE_DIR/${file_name}/"
    mkdir -p "$output_dir"

    # Call the executable with the required parameters
    $EXECUTABLE -i "$file_path" -o "$output_dir" -n "$EVENT_NUMBER" -t "$T" -c "$CALIBRATION_CSV"

    if [ $? -ne 0 ]; then
        echo "Error processing file: $file_path"
    fi
}

###############################################
# Main Script
###############################################

echo "Enter the directory containing .root files:"
read -e INPUT_FILE_DIR
INPUT_FILE_DIR=$(eval echo "$INPUT_FILE_DIR")


# Prompt the user for the TIMESTAMP (e.g., "250204")
echo "Enter TIMESTAMP to search for in file names: (YYMMDD)"
read TIMESTAMP

OUTPUT_FILE_DIR="/Users/gwon/WbLS/1ton_analysis/output/$TIMESTAMP"

# Prompt for trigger and number of events to process
echo "Enter trigger (-1: all, 0: tp, 1: alpha, 2: majority, 3: crossing_muon):"
read T

echo "Enter number of events to process:"
read EVENT_NUMBER

# Set constant paths (adjust as needed)
CALIBRATION_CSV="/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/bnl1t_spe_fit_results_241111.csv"
EXECUTABLE="/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/ProcessRawRoot"

# Construct the glob pattern.
# Remove any trailing "/" from INPUT_FILE_DIR and append the pattern for files containing TIMESTAMP.
pattern="${INPUT_FILE_DIR%/}/*${TIMESTAMP}*.root"

# Expand the glob pattern and store matching files in an array.
found_files=($pattern)

# Exit if no matching files are found.
if [ ${#found_files[@]} -eq 0 ]; then
    echo "No .root files found in '$INPUT_FILE_DIR' containing '$TIMESTAMP'."
    exit 1
fi

# Process each found file using the process_root_file function.
for file in "${found_files[@]}"; do
    process_root_file "$file"
done
