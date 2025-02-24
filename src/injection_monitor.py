import os
import sys
import numpy as np
import multiprocessing
import pickle
from functools import partial
from numba import njit

WbLS_config = "07"
injection_dir = (
    f"/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_{WbLS_config}0_injection"
)
timestamp = "250122T1421"

# 필요한 모듈 임포트
sys.path.append("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/lib")
sys.path.append("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/lib")
sys.path.append("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/src")
sys.path.append("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/tools")

from event_display import EventDisplay

# 환경 변수 설정
os.environ["SOURCE_DIR"] = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/src"
os.environ["YAML_DIR"] = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml"
yaml = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/config_1740.yaml"

# 채널 리스트
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



# Numba를 사용하여 계산 함수 컴파일
@njit
def compute_pe_sum(amp_pe_array):
    return np.sum(amp_pe_array) * 2


def fill_injection_test_hists(wfm, ievt):
    t = wfm.time_axis_ns
    mask = (t >= 0) & (t < 2000)
    y_data = wfm.amp_pe["sum"][mask]

    max_y = np.max(y_data)
    max_index = np.where(y_data == max_y)[0][0]
    if max_index > 200:
        return None

    temp_bot_npe = 0.0
    temp_side_npe = 0.0
    npe_channel = {channel: [] for channel in bot_pmt_channels + side_pmt_channels}

    for chn in bot_pmt_channels:
        amp_pe_array = wfm.amp_pe[chn][mask]
        PE = compute_pe_sum(amp_pe_array)
        temp_bot_npe += PE
        npe_channel[chn].append(PE)

    for chn in side_pmt_channels:
        amp_pe_array = wfm.amp_pe[chn][mask]
        PE = compute_pe_sum(amp_pe_array)
        temp_side_npe += PE
        npe_channel[chn].append(PE)

    print(f"Event {ievt} done")
    return temp_bot_npe, temp_side_npe, npe_channel


def process_event(args):
    ievt, fpath = args
    dp = EventDisplay(fpath, yaml)
    wfm = dp.get_all_waveform(ievt)
    if wfm is None:
        print(f"wfm == None, ievt: {ievt}")
        return None

    result = fill_injection_test_hists(wfm, ievt)
    return result


def process_file(name):
    fpath = os.path.join(injection_dir, f"{name}.root")
    print(f"Processing file: {fpath}")

    dp = EventDisplay(fpath, yaml)
    dp.grab_events(-1)
    interval = 10000
    rangeStart = dp.min_event_id
    rangeEnd = -999
    if dp.max_event_id > rangeStart + interval:
        rangeEnd = rangeStart + interval
    else:
        rangeEnd = dp.max_event_id
    print("range: " + str(rangeStart) + "~" + str(rangeEnd))
    event_ids = range(rangeStart, rangeEnd)

    # 멀티프로세싱을 위한 인자 준비
    args_list = [(ievt, fpath) for ievt in event_ids]

    # 멀티프로세싱 Pool 사용
    cpu_count = multiprocessing.cpu_count()
    # cpu_count = 5
    with multiprocessing.Pool(processes=cpu_count) as pool:
        results = pool.map(process_event, args_list)

    # 결과 집계
    bot_npe = []
    side_npe = []
    npe_channel = {channel: [] for channel in bot_pmt_channels + side_pmt_channels}

    for result in results:
        if result is not None:
            temp_bot_npe, temp_side_npe, event_npe_channel = result
            bot_npe.append(temp_bot_npe)
            side_npe.append(temp_side_npe)
            for chn in npe_channel.keys():
                npe_channel[chn].extend(event_npe_channel[chn])

    # 데이터 저장
    with open(f"{injection_dir}/bot_npe_{name}.pkl", "wb") as f:
        pickle.dump(bot_npe, f)
    with open(f"{injection_dir}/side_npe_{name}.pkl", "wb") as f:
        pickle.dump(side_npe, f)
    for chn in npe_channel.keys():
        with open(f"{injection_dir}/npe_channel_{name}_{chn}.pkl", "wb") as f:
            pickle.dump(npe_channel[chn], f)


if __name__ == "__main__":
    args = sys.argv[1:]
    print(f"Received arguments: {args}")
    start = int(sys.argv[1])
    for i in range(start, 200):
        name = f"injection_{WbLS_config}_{timestamp}_{i}"
        try:
            process_file(name)
        except:
            continue

    print("Done")
