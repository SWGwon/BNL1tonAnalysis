set -e
source /Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/setup.sh
source /Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/venv/bin/activate
##-W 6 : 025 WbLS
##-W 7 : 04 WbLS
./SelectCrossingMuon -L $1 -T $2
./MakeCrossingMuonRootFile -L $1 -T $2
python /Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/src/ProcessCrosingMuon.py --Dir /Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_$1/ --TimeStamp $2

rm ~/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_$1/*$2*root
