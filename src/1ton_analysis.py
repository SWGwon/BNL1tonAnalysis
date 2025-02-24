import os
import sys
import numpy as np
import multiprocessing
import pickle
from functools import partial
from numba import njit
import glob

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

side_pmt_row1 = [
    "adc_b3_ch0",
    "adc_b3_ch4",
    "adc_b3_ch8",
    "adc_b3_ch12",
]
side_pmt_row2 = [
    "adc_b3_ch1",
    "adc_b3_ch5",
    "adc_b3_ch9",
    "adc_b3_ch13",
]
side_pmt_row3 = [
    "adc_b3_ch2",
    "adc_b3_ch6",
    "adc_b3_ch10",
    "adc_b3_ch14",
]
side_pmt_row4 = [
    "adc_b3_ch3",
    "adc_b3_ch7",
    "adc_b3_ch11",
    "adc_b3_ch15",
]


# mean과 std를 계산하는 함수 (low_limit과 high_limit을 사용)
def calculate_stats(data, low_limit, high_limit):
    filtered_data = [x for x in data if low_limit <= x < high_limit]
    mean = np.mean(filtered_data)
    std = np.std(filtered_data)
    return mean, std


# mean과 std를 기반으로 새로운 범위 내의 mean과 std를 계산하는 함수
def refined_stats(data, low_limit, high_limit):
    # 첫 번째 단계: 전체 범위에서 mean과 std 계산
    mean, std = calculate_stats(data, low_limit, high_limit)

    # 두 번째 단계: mean ± 2 * std 범위 내의 데이터로 다시 mean과 std 계산
    refined_low = mean - 2 * std
    refined_high = mean + 2 * std

    refined_mean, refined_std = calculate_stats(data, refined_low, refined_high)

    return refined_mean, refined_std, refined_low, refined_high


def compute_results(
    npe_channel,
    oldData=False,
    npe_channel_old=None,
    efficiency_ratio=None,
    bot_pmt_channels=None,
    side_pmt_channels=None,
):

    # Ensure required parameters are provided
    if bot_pmt_channels is None or side_pmt_channels is None:
        raise ValueError("bot_pmt_channels and side_pmt_channels must be provided.")

    if efficiency_ratio is None:
        raise ValueError("efficiency_ratio must be provided.")

    if oldData and npe_channel_old is None:
        raise ValueError("npe_channel_old must be provided when oldData is True.")

    # Combine bot and side channels
    bot_and_side_channels = bot_pmt_channels + side_pmt_channels

    # Initialize weighted NPE channel dictionary
    weighted_npe_channel = {}

    # Compute weighted NPE values for each channel
    for chn in bot_and_side_channels:
        weighted_npe_channel[chn] = [
            value / efficiency_ratio[chn] for value in npe_channel[chn]
        ]

    # Initialize arrays for bot, side, and total channels
    arrays_bot = []
    arrays_bot_weight = []
    arrays_side = []
    arrays_side_weight = []
    arrays_tot = []
    arrays_tot_weight = []

    if oldData:
        arrays_bot_old = []
        arrays_side_old = []
        arrays_tot_old = []

    # Process bottom PMT channels
    for chn in bot_pmt_channels:
        if oldData:
            arrays_bot_old.append(npe_channel_old[chn])
            arrays_tot_old.append(npe_channel_old[chn])
        arrays_bot.append(npe_channel[chn])
        arrays_bot_weight.append(weighted_npe_channel[chn])
        arrays_tot.append(npe_channel[chn])
        arrays_tot_weight.append(weighted_npe_channel[chn])

    # Process side PMT channels
    for chn in side_pmt_channels:
        if oldData:
            arrays_side_old.append(npe_channel_old[chn])
            arrays_tot_old.append(npe_channel_old[chn])
        arrays_side.append(npe_channel[chn])
        arrays_side_weight.append(weighted_npe_channel[chn])
        arrays_tot.append(npe_channel[chn])
        arrays_tot_weight.append(weighted_npe_channel[chn])

    # Initialize result lists
    result_bot = []
    result_bot_weight = []
    result_side = []
    result_side_weight = []
    result_tot = []
    result_tot_weight = []

    if oldData:
        result_bot_old = []
        result_side_old = []
        result_tot_old = []

    # Sum NPE values across channels for bot, side, and total
    for elements in zip(*arrays_bot):
        result_bot.append(sum(elements))

    for elements in zip(*arrays_bot_weight):
        result_bot_weight.append(sum(elements))

    for elements in zip(*arrays_side):
        result_side.append(sum(elements))

    for elements in zip(*arrays_side_weight):
        result_side_weight.append(sum(elements))

    for elements in zip(*arrays_tot):
        result_tot.append(sum(elements))

    for elements in zip(*arrays_tot_weight):
        result_tot_weight.append(sum(elements))

    if oldData:
        for elements in zip(*arrays_bot_old):
            result_bot_old.append(sum(elements))

        for elements in zip(*arrays_side_old):
            result_side_old.append(sum(elements))

        for elements in zip(*arrays_tot_old):
            result_tot_old.append(sum(elements))

    # Prepare output results
    results = {
        "result_bot": result_bot,
        "result_bot_weight": result_bot_weight,
        "result_side": result_side,
        "result_side_weight": result_side_weight,
        "result_tot": result_tot,
        "result_tot_weight": result_tot_weight,
    }

    if oldData:
        results["result_bot_old"] = result_bot_old
        results["result_side_old"] = result_side_old
        results["result_tot_old"] = result_tot_old

    return results


import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 파일에서 불러오는 코드

xlim = 50


def loadNPEchannel(waterLevel, xlim, timestamp=None):
    bot_and_side_channels = bot_pmt_channels + side_pmt_channels
    if waterLevel == -2:
        npe_channel = {channel: [] for channel in bot_and_side_channels}
    else:
        npe_channel = {
            channel: [] for channel in bot_and_side_channels + side_additional_channels
        }

    # 검색하려는 디렉토리 경로 설정
    if waterLevel == -2:
        directories = ["/Users/gwon/WbLS/1ton_analysis/oldfile/analysis/"]
    elif waterLevel == -1:
        directories = [
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_water_quality_test/"
        ]
    elif waterLevel == 0:
        directories = [
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_empty/"
        ]
    elif waterLevel == 1:
        directories = [
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_quarter/"
        ]
    elif waterLevel == 2:
        directories = [
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_half/"
        ]
    elif waterLevel == 3:
        directories = [
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_full_0/",
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_full_1/",
        ]
    elif waterLevel == 4:
        directories = [
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_injection_test/"
        ]
    elif waterLevel == 5:
        directories = [
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_025_WbLS/"
        ]
    else:
        directories = []

    print("Directories: ", directories)

    # 파일명 리스트 초기화
    fileName = []

    # timestamp가 None이거나 문자열이면 리스트로 변환
    if timestamp is None:
        timestamp_list = [""]
    elif isinstance(timestamp, str):
        timestamp_list = [timestamp]
    else:
        timestamp_list = timestamp  # 이미 리스트인 경우

    # 각 디렉토리에서 파일 목록 읽어서 fileName에 추가
    for dir_path in directories:
        print("Scanning directory: ", dir_path)
        for ts in timestamp_list:
            pattern = os.path.join(dir_path, f"*{ts}*.root")
            files = glob.glob(pattern)
            file_names = [os.path.splitext(os.path.basename(file))[0] for file in files]
            fileName.extend(file_names)

    # 파일명 출력
    # 특정 패턴에 맞는 파일 목록을 스캔하여 리스트에 저장
    # fileName = glob.glob(os.path.join(directory, '*.root'))

    # 확장자를 제외한 파일명만 추출
    # fileName = [os.path.splitext(os.path.basename(file))[0] for file in fileName]
    print(fileName)

    for name in fileName:
        # print(name)
        for chn in bot_and_side_channels + side_additional_channels:
            file_path = f"/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/{name}/npe_channel_{chn}.pkl"

            # 파일이 존재하는지 확인
            if os.path.exists(file_path):
                # 파일이 비어 있지 않은지 확인
                if os.path.getsize(file_path) > 0:
                    with open(file_path, "rb") as f:
                        loaded_data = pickle.load(f)
                        filtered_data = [
                            value for value in loaded_data if 0 <= value <= xlim
                        ]
                        npe_channel[chn].extend(loaded_data)

    print("total events:", len(npe_channel["adc_b1_ch1"]))

    n_rows = 8  # Adjust as needed
    n_cols = 8  # Adjust as needed
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 20))
    fig.tight_layout(pad=4.0)

    for i, chn in enumerate(bot_and_side_channels):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]

        ax.hist(
            [x for x in npe_channel[chn] if x < xlim and x > 0],
            bins=25,
            density=True,
            color="blue",
            alpha=0.5,
        )
        ax.set_title(chn)
        ax.set_xlim(0, xlim)

    # Hide any empty subplots
    for j in range(i + 1, n_rows * n_cols):
        fig.delaxes(axes.flatten()[j])

    plt.show()

    return npe_channel


npe_channel_old = loadNPEchannel(waterLevel=5, xlim=50)
