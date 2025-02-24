"""
Load package. 

Try run this script in analysis directory.
"""

import uproot
import re
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.fft import fft, fftfreq, rfft, rfftfreq
import os
import sys
import glob
from datetime import datetime
from numpy import array, isscalar, uint16, uint32
import pandas as pd
from iminuit import Minuit
from numpy import sqrt
from scipy.stats import norm
import statistics
import uproot
import types
from mpl_toolkits.mplot3d import Axes3D
from concurrent.futures import ThreadPoolExecutor
from functools import partial
import threading

import sys

sys.path.append("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/lib")
sys.path.append("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/lib")
sys.path.append("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/src")
sys.path.append("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/tools")

# define enviromental variables
os.environ["SOURCE_DIR"] = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/src"
os.environ["YAML_DIR"] = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml"
from event_display import EventDisplay

import argparse

# 인자 파서 생성
parser = argparse.ArgumentParser(description="Process some arguments.")

# 인자 추가
parser.add_argument("--TimeStamp", type=str, required=False, help="Your name")
# parser.add_argument('--WaterLevel', type=int, required=True, help='Your name')
parser.add_argument("--Dir", type=str, required=True, help="Your name")
# 인자 파싱
args = parser.parse_args()

# 인자 사용
print(f"TimeStamp: {args.TimeStamp}")
# print(f"WaterLevel: {args.WaterLevel}")
print(f"file directory: {args.Dir}")

# 06/23, change event_ttt to evnent_ttt_1 (run_drop.py, waveform.py)
# from 07/04, change event_ttt to evnent_ttt_1 (run_drop.py, waveform.py)
# config = 1: event_ttt_1              , 0 : event_ttt
config = 1
time_config = 1
no_thread = 0

# 0: empty,
# 1: quarter,
# 2: half,
# 3: full
crossingMuonDir = args.Dir

# 검색하려는 디렉토리 경로 설정
# if waterLevel == 1:
#    crossingMuonDir = '/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_empty/'
#    crossingMuonDir = '/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_water_quality_test/'
# if waterLevel == 2:
#    crossingMuonDir = '/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_quarter/'
# if waterLevel == 3:
#    crossingMuonDir = '/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_half/'
# if waterLevel == 4:
#    if config == 0:
#        crossingMuonDir = '/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_full_0/'
#    if config == 1:
#        crossingMuonDir = '/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_full_1/'
# if waterLevel == 5:
#    crossingMuonDir = '/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_injection_test/'
# if waterLevel == 6:
#    crossingMuonDir = '/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_025_WbLS/'

# 특정 패턴에 맞는 파일 목록을 스캔하여 리스트에 저장
fileName = glob.glob(os.path.join(crossingMuonDir, f"*{args.TimeStamp}*.root"))

# 확장자를 제외한 파일명만 추출
fileName = [os.path.splitext(os.path.basename(file))[0] for file in fileName]

print(fileName)

for name in fileName:
    os.makedirs(
        "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/" + name, exist_ok=True
    )

bot_pmt_channels = [
    "adc_b1_ch1",
    "adc_b1_ch2",
    "adc_b1_ch3",
    "adc_b1_ch4",
    "adc_b1_ch5",
    "adc_b1_ch6",
    "adc_b1_ch7",
    "adc_b1_ch8",
    "adc_b1_ch9",
    "adc_b1_ch10",
    "adc_b1_ch11",
    "adc_b1_ch12",
    "adc_b1_ch13",
    "adc_b1_ch14",
    "adc_b1_ch15",
    "adc_b2_ch0",
    "adc_b2_ch1",
    "adc_b2_ch2",
    "adc_b2_ch3",
    "adc_b2_ch4",
    "adc_b2_ch5",
    "adc_b2_ch6",
    "adc_b2_ch7",
    "adc_b2_ch8",
    "adc_b2_ch9",
    "adc_b2_ch10",
    "adc_b2_ch11",
    "adc_b2_ch12",
    "adc_b2_ch13",
    "adc_b2_ch14",
]

side_pmt_channels = [
    "adc_b3_ch0",
    "adc_b3_ch1",
    "adc_b3_ch2",
    "adc_b3_ch3",
    "adc_b3_ch4",
    "adc_b3_ch5",
    "adc_b3_ch6",
    "adc_b3_ch7",
    "adc_b3_ch8",
    "adc_b3_ch9",
    "adc_b3_ch10",
    "adc_b3_ch11",
    "adc_b3_ch12",
    "adc_b3_ch13",
    "adc_b3_ch14",
    "adc_b3_ch15",
    #'adc_b4_ch0', 'adc_b4_ch1','adc_b4_ch2', 'adc_b4_ch3', 'adc_b4_ch4', 'adc_b4_ch5', 'adc_b4_ch6', 'adc_b4_ch7','adc_b4_ch8', 'adc_b4_ch9', 'adc_b4_ch10', 'adc_b4_ch11',
]

side_additional_channels = [
    "adc_b4_ch0",
    "adc_b4_ch1",
    "adc_b4_ch2",
    "adc_b4_ch3",
    "adc_b4_ch4",
    "adc_b4_ch5",
    "adc_b4_ch6",
    "adc_b4_ch7",
    "adc_b4_ch8",
    "adc_b4_ch9",
    "adc_b4_ch10",
    "adc_b4_ch11",
]

bot_and_side_channels = [
    "adc_b1_ch1",
    "adc_b1_ch2",
    "adc_b1_ch3",
    "adc_b1_ch4",
    "adc_b1_ch5",
    "adc_b1_ch6",
    "adc_b1_ch7",
    "adc_b1_ch8",
    "adc_b1_ch9",
    "adc_b1_ch10",
    "adc_b1_ch11",
    "adc_b1_ch12",
    "adc_b1_ch13",
    "adc_b1_ch14",
    "adc_b1_ch15",
    "adc_b2_ch0",
    "adc_b2_ch1",
    "adc_b2_ch2",
    "adc_b2_ch3",
    "adc_b2_ch4",
    "adc_b2_ch5",
    "adc_b2_ch6",
    "adc_b2_ch7",
    "adc_b2_ch8",
    "adc_b2_ch9",
    "adc_b2_ch10",
    "adc_b2_ch11",
    "adc_b2_ch12",
    "adc_b2_ch13",
    "adc_b2_ch14",
    "adc_b3_ch0",
    "adc_b3_ch1",
    "adc_b3_ch2",
    "adc_b3_ch3",
    "adc_b3_ch4",
    "adc_b3_ch5",
    "adc_b3_ch6",
    "adc_b3_ch7",
    "adc_b3_ch8",
    "adc_b3_ch9",
    "adc_b3_ch10",
    "adc_b3_ch11",
    "adc_b3_ch12",
    "adc_b3_ch13",
    "adc_b3_ch14",
    "adc_b3_ch15",
    #'adc_b4_ch0', 'adc_b4_ch1', 'adc_b4_ch2', 'adc_b4_ch3', 'adc_b4_ch4', 'adc_b4_ch5', 'adc_b4_ch6', 'adc_b4_ch7','adc_b4_ch8', 'adc_b4_ch9', 'adc_b4_ch10', 'adc_b4_ch11',
]

lock = threading.Lock()

import ast

good_eventid = []
time_loc = []


def select_muon(wfm, ievt, npe_channel):
    if wfm == None:
        return

    t = wfm.time_axis_ns

    mask = (t >= 0) & (t < 2000)
    x_data = t[mask]
    y_data = wfm.amp_pe["sum"][mask]

    max_y = np.max(y_data)
    max_index = np.where(y_data == max_y)[0][0]
    time_loc.append(max_index)
    print("event ", ievt)
    print("max x and y (time) locations:", max_index, max_y)

    if time_config == 0:
        if max_index < 300 or max_index > 450:
            return  # to select muon peak
    if time_config == 1:
        if max_index < 300 or max_index > 450:
            return  # to select muon peak

    mask = (t >= (max_index - 20) * 2) & (t < (max_index + 40) * 2)  # muon peak

    print("mask", (max_index - 20) * 2, "~", (max_index + 40) * 2)
    for chn in bot_and_side_channels:
        npe_channel[chn].append(np.sum(wfm.amp_pe[chn][mask]) * 2)
        # npe_channel[chn].append(np.sum(wfm.amp_pe[chn])*2)

    for chn in side_additional_channels:
        npe_channel[chn].append(np.sum(wfm.amp_pe[chn][mask]) * 2)
        # npe_channel[chn].append(np.sum(wfm.amp_pe[chn])*2)

    print(ievt, "done")


thread_local = threading.local()


def get_thread_local_dp(name):
    if not hasattr(thread_local, "dp"):
        print("making dp instance")
        file = crossingMuonDir + name + ".root"
        print(file)
        date_str = name.split("_")[3][:6]
        if date_str >= "241111":
            date_str = "241111"
        yaml = f"/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/config_1740_{date_str}.yaml"
        # yaml = '/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/config_1740.yaml'
        print(yaml)
        thread_local.dp = EventDisplay(file, yaml)
    return thread_local.dp


import sys
import io
import time


def fetch_waveform(ievt, name, npe_channel):
    print(ievt)
    dp = get_thread_local_dp(name)
    wfm = dp.get_all_waveform(ievt)

    if wfm == None:
        print("wfm == None")
        return None

    with lock:
        select_muon(wfm, ievt, npe_channel)
    del thread_local.dp


fileIndex = 2

if no_thread == 1:
    npe_channel = {
        channel: [] for channel in bot_and_side_channels + side_additional_channels
    }

    # 파일 경로 설정
    fpath = os.path.join(
        crossingMuonDir, fileName[fileIndex] + "_good_eventId_with_paddle_CXX.txt"
    )

    # 파일이 존재하는지 확인
    if not os.path.exists(fpath):
        print(f"File not found: {fpath}")

    # 파일 열기 및 내용 읽기
    # with open(fpath, 'r') as f:
    #   event_ids = ast.literal_eval(f.read())

    with open(fpath, "r") as file:
        event_ids = [int(line.strip()) for line in file if line.strip().isdigit()]

    print(event_ids)
    for a in event_ids:
        fetch_waveform(a, fileName[fileIndex], npe_channel)

    for chn in bot_and_side_channels + side_additional_channels:
        with open(
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/"
            + fileName[fileIndex]
            + "/"
            + f"npe_channel_{chn}.pkl",
            "wb",
        ) as f:
            pickle.dump(npe_channel[chn], f)


max_workers = os.cpu_count()  # 논리적 CPU 코어 수
if no_thread == 0:
    for name in fileName:
        npe_channel = {
            channel: [] for channel in bot_and_side_channels + side_additional_channels
        }
        try:
            # 파일 경로 설정
            fpath = os.path.join(
                crossingMuonDir, name + "_good_eventId_with_paddle_CXX.txt"
            )

            # 파일이 존재하는지 확인
            if not os.path.exists(fpath):
                print(f"File not found: {fpath}")
                continue

            # 파일 열기 및 내용 읽기
            # with open(fpath, 'r') as f:
            #   event_ids = ast.literal_eval(f.read())
            with open(fpath, "r") as file:
                event_ids = [
                    int(line.strip()) for line in file if line.strip().isdigit()
                ]

        except Exception as e:
            # 다른 오류가 발생할 경우 에러 메시지 출력
            print(f"Error processing file {fpath}: {e}")
            continue

        print(name)
        print("muon events:", event_ids)

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            partial_fetch_waveform = partial(
                fetch_waveform, name=name, npe_channel=npe_channel
            )
            executor.map(partial_fetch_waveform, event_ids)

        os.makedirs(
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/" + name,
            exist_ok=True,
        )
        # 파일로 저장하는 코드
        for chn in bot_and_side_channels + side_additional_channels:
            with open(
                "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/"
                + name
                + "/"
                + f"npe_channel_{chn}.pkl",
                "wb",
            ) as f:
                pickle.dump(npe_channel[chn], f)

        for chn in bot_and_side_channels + side_additional_channels:
            file_path = f"/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/{name}/npe_channel_{chn}.txt"
            with open(file_path, "w") as f:
                for value in npe_channel[chn]:
                    f.write(f"{value}\n")


print("all done")
