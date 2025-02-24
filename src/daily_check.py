"""
Load package. 

Try run this script in analysis directory.
"""

import uproot
import re
import numpy as np
import time

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
import types
from mpl_toolkits.mplot3d import Axes3D

# define constants
ADC_TO_MV = 2000 / (2**14 - 1)
SAMPLE_TO_NS = 2

interval = 2000
fileName="muon_wbls07_250201T0947_1"


outdir = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/analysis/" + fileName
os.makedirs(outdir, exist_ok=True)
# fpath='/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_full_1/' + fileName + '.root'
fpath = (
    "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_070_WbLS/"
    + fileName
    + ".root"
)

# import EventDisplay
import sys

sys.path.append("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/lib")
sys.path.append("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/lib")
sys.path.append("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/src")
sys.path.append("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/tools")

os.environ["SOURCE_DIR"] = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/src"
os.environ["YAML_DIR"] = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml"
yaml = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/config_1740.yaml"
from event_display import EventDisplay

# veiw waveform using EventDisplay
date_str = fileName.split("_")[3][:6]
# if date_str >= "240817":
#   date_str = "240816"
print("date: ", date_str)
# 추출한 날짜를 사용해 yaml 파일 경로 생성
yaml = (
    f"/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/config_1740_{date_str}.yaml"
)
yaml = f"/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/config_1740.yaml"

dp = EventDisplay(fpath, yaml)
asdf = dp.grab_events(-1)
rangeStart = dp.min_event_id
rangeEnd = -999
if dp.max_event_id > rangeStart + interval:
    rangeEnd = rangeStart + interval
else:
    rangeEnd = dp.max_event_id
print("range: " + str(rangeStart) + "~" + str(rangeEnd))
event_ids = range(rangeStart, rangeEnd)


from concurrent.futures import ThreadPoolExecutor
from functools import partial
import threading

lock = threading.Lock()

dp.fig_height = 10
dp.fig_width = 20
schn = [[] * i for i in range(97)]
h_hodo = [[] * i for i in range(4)]
time_loc = []
ichn = 0
hevt = 0
print("starting")
hodo5 = 0
hodo6 = 0
hodo7 = 0
top_fire = np.zeros((10, 10))
bot_fire = np.zeros((10, 10))
time_loc = []
chg_wtime = []
good_eventID = []


def process_waveform(wfm, ievt):
    t = wfm.time_axis_ns
    mask = (t >= 0) & (t < 2000)
    x_data = t[mask]
    y_data = wfm.amp_pe["sum"][mask]

    max_y = np.max(y_data)
    max_index = np.where(y_data == max_y)[0][0]

    time_loc.append(max_index)
    chg_wtime.append(np.sum(y_data))

    hodo1 = 0
    hodo2 = 0
    hodo3 = 0
    hodo4 = 0
    top1 = []
    top2 = []
    bot1 = []
    bot2 = []
    ichn = 0
    iichn = 0

    for chn in [
        "adc_b1_ch1",
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
        "adc_b4_ch11",
        "adc_b4_ch13",
        "adc_b4_ch14",
        "adc_b4_ch15",
        "adc_b5_ch1",
        "adc_b5_ch2",
        "adc_b5_ch3",
        "adc_b5_ch4",
        "adc_b5_ch5",
        "adc_b5_ch6",
        "adc_b5_ch7",
        "adc_b5_ch8",
        "adc_b5_ch9",
        "adc_b5_ch10",
        "adc_b5_ch11",
        "adc_b5_ch12",
        "adc_b5_ch13",
        "adc_b5_ch14",
        "adc_b5_ch15",
        "adc_b5_ch16",
        "adc_b5_ch17",
        "adc_b5_ch18",
        "adc_b5_ch19",
        "adc_b5_ch20",
        "adc_b5_ch21",
        "adc_b5_ch22",
        "adc_b5_ch23",
        "adc_b5_ch24",
        "adc_b5_ch25",
        "adc_b5_ch26",
        "adc_b5_ch27",
        "adc_b5_ch28",
        "adc_b5_ch29",
        "adc_b5_ch30",
        "adc_b5_ch31",
        "adc_b5_ch32",
    ]:
        schn[ichn].append(np.sum(wfm.amp_pe[chn]))
        ichn += 1
        if ichn == 1:
            continue
        if ichn == 32:
            continue

    for ih in range(8):
        if schn[64 + ih][-1] > 0.5:
            hodo1 += 1
            top1.append(ih)
        if schn[64 + 8 + ih][-1] > 0.5:
            hodo2 += 1
            top2.append(ih)
        if schn[64 + 16 + ih][-1] > 0.5:
            hodo3 += 1
            bot1.append(ih)
        if schn[64 + 24 + ih][-1] > 0.5:
            hodo4 += 1
            bot2.append(ih)

    if len(top1) != 0 or len(top2) != 0:
        # print ('top1 ',top1)
        # print ('top2 ',top2)
        if len(top1) == 0 and len(top2) != 0:
            for jj in range(len(top2)):
                top_fire[9][jj] += 1
        if len(top2) == 0 and len(top1) != 0:
            for ii in range(len(top1)):
                top_fire[ii][9] += 1
        if len(top2) != 0 and len(top1) != 0:
            for ii in range(len(top1)):
                for jj in range(len(top2)):
                    # print ('?? ', ii,' ',jj, ' ',top1[ii],' ', top2[jj])
                    top_fire[top1[ii]][top2[jj]] = top_fire[top1[ii]][top2[jj]] + 1

    # print (top_fire)
    if len(bot1) != 0 or len(bot2) != 0:
        # print ('bot1 ',bot1)
        # print ('bot2 ',bot2)
        if len(bot1) == 0 and len(bot2) != 0:
            for jj in range(len(bot2)):
                bot_fire[9][jj] += 1
        if len(bot2) == 0 and len(bot1) != 0:
            for ii in range(len(bot1)):
                bot_fire[ii][9] += 1
        if len(bot1) != 0 and len(bot2) != 0:
            for ii in range(len(bot1)):
                for jj in range(len(bot2)):
                    bot_fire[bot1[ii]][bot2[jj]] = bot_fire[bot1[ii]][bot2[jj]] + 1

    if hodo1 == 1 and hodo2 == 1 and hodo3 == 1 and hodo4 == 1:
        # print ("::::::::::: good event ::::: ",ievt)
        good_eventID.append(ievt)
    print(
        "file: ",
        fileName,
        ", event: ",
        ievt,
        ", samples: ",
        t.size,
        ", max x and y (time) locations:",
        max_index,
    )


thread_local = threading.local()


def get_thread_local_dp():
    if not hasattr(thread_local, "dp"):
        # print('making dp instance')
        thread_local.dp = EventDisplay(
            fpath,
            "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/config_1740.yaml",
        )
    return thread_local.dp


def fetch_waveform(ievt):
    dp = get_thread_local_dp()
    wfm = dp.get_all_waveform(ievt)

    if wfm == None:
        print("wfm == None")
        return None

    with lock:
        process_waveform(wfm, ievt)


start_time = time.time()

# max_workers = os.cpu_count()  # 논리적 CPU 코어 수
max_workers = 5  # 논리적 CPU 코어 수

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    # with ThreadPoolExecutor() as executor:
    executor.map(fetch_waveform, event_ids)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"running time: {elapsed_time} sec")

##############################################################################


def saveCSV(data, dataName):
    df = pd.DataFrame(data)
    df.to_csv(outdir + "/" + fileName + "_" + dataName + ".csv", index=False)


plt.hist(time_loc, bins=70, range=(60, 500))
plt.ylabel("counts")
plt.xlabel("relative time (x 2ns)")
plt.yscale("log")
plt.savefig(outdir + "/" + fileName + "_double_peak.png")
plt.show()
saveCSV(time_loc, "time_loc")


print(top_fire)
print(bot_fire)
plt.imshow(top_fire, norm=colors.LogNorm())
plt.colorbar()
plt.savefig(outdir + "/" + fileName + "_top_fire.png")
plt.show()
saveCSV(top_fire, "top_fire")

plt.imshow(bot_fire, norm=colors.LogNorm())
plt.colorbar()
plt.savefig(outdir + "/" + fileName + "_bot_fire.png")
plt.show()
saveCSV(bot_fire, "bot_fire")


# print (schn)
# plt.figure()

f, a = plt.subplots(4, 4)
a = a.ravel()
for idx, ax in enumerate(a):
    if idx != 0:
        ax.hist(
            schn[idx + 0],
            bins=50,
            range=[0, 20],
            label="chn {}".format(idx + 1),
            histtype="step",
            linewidth=2,
            color="c",
            hatch="",
            edgecolor="b",
            fill=False,
        )
        ax.set_xlabel("mV ns")
        ax.legend()
        ax.set_yscale("log")
plt.tight_layout()
plt.savefig(outdir + "/" + fileName + "_b1.png")
saveCSV(schn, "schn")

f, a = plt.subplots(4, 4)
a = a.ravel()
for idx, ax in enumerate(a):
    if idx + 16 != 31:
        ax.hist(
            schn[idx + 16],
            bins=50,
            range=[0, 20],
            label="chn {}".format(idx + 1 + 16),
            histtype="step",
            linewidth=2,
            color="c",
            hatch="",
            edgecolor="b",
            fill=False,
        )
        ax.set_xlabel("mV ns")
        ax.legend()
        ax.set_yscale("log")
plt.tight_layout()
plt.savefig(outdir + "/" + fileName + "_b2.png")

f, a = plt.subplots(4, 4)
a = a.ravel()
for idx, ax in enumerate(a):
    if idx + 32 == 39:
        ax.hist(
            schn[idx + 32],
            bins=50,
            range=[0, 20],
            label="chn {}".format(idx + 32),
            histtype="step",
            linewidth=2,
            color="c",
            hatch="",
            edgecolor="b",
            fill=False,
        )
    else:
        ax.hist(
            schn[idx + 32],
            bins=50,
            range=[0, 20],
            label="chn {}".format(idx + 32),
            histtype="step",
            linewidth=2,
            color="c",
            hatch="",
            edgecolor="b",
            fill=False,
        )
    ax.set_xlabel("mV ns")
    ax.legend()
    ax.set_yscale("log")
plt.tight_layout()
plt.savefig(outdir + "/" + fileName + "_b3.png")

f, a = plt.subplots(4, 4)
a = a.ravel()
for idx, ax in enumerate(a):
    if idx + 48 != 60:
        if idx + 48 == 61 or idx + 48 == 62:
            ax.hist(
                schn[idx + 48],
                bins=50,
                range=[0, 1000],
                label="chn {}".format(idx + 48),
                histtype="step",
                linewidth=2,
                color="c",
                hatch="",
                edgecolor="b",
                fill=False,
            )
        else:
            ax.hist(
                schn[idx + 48],
                bins=50,
                range=[0, 20],
                label="chn {}".format(idx + 48),
                histtype="step",
                linewidth=2,
                color="c",
                hatch="",
                edgecolor="b",
                fill=False,
            )
        ax.set_xlabel("mV ns")
        ax.legend()
        ax.set_yscale("log")
plt.tight_layout()
plt.savefig(outdir + "/" + fileName + "_b4.png")

f, a = plt.subplots(4, 2)
a = a.ravel()
for idx, ax in enumerate(a):
    if idx == 6:
        ax.hist(
            schn[idx + 64],
            bins=50,
            range=[0, 1],
            label="chn {}".format(idx + 1),
            histtype="step",
            linewidth=2,
            color="c",
            hatch="",
            edgecolor="b",
            fill=False,
        )
    else:
        ax.hist(
            schn[idx + 64],
            bins=50,
            range=[0, 4],
            label="chn {}".format(idx + 1),
            histtype="step",
            linewidth=2,
            color="c",
            hatch="",
            edgecolor="b",
            fill=False,
        )
    ax.set_xlabel("mV ns")
    ax.legend()
    ax.set_yscale("log")
plt.tight_layout()
plt.savefig(outdir + "/" + fileName + "_b5_1.png")

f, a = plt.subplots(4, 2)
a = a.ravel()
for idx, ax in enumerate(a):
    ax.hist(
        schn[idx + 72],
        bins=50,
        range=[0, 4],
        label="chn {}".format(idx + 9),
        histtype="step",
        linewidth=2,
        color="c",
        hatch="",
        edgecolor="b",
        fill=False,
    )
    ax.set_xlabel("mV ns")
    ax.legend()
    ax.set_yscale("log")
plt.tight_layout()
plt.savefig(outdir + "/" + fileName + "_b5_2.png")

f, a = plt.subplots(4, 2)
a = a.ravel()
for idx, ax in enumerate(a):
    ax.hist(
        schn[idx + 80],
        bins=50,
        range=[0, 1],
        label="chn {}".format(idx + 17),
        histtype="step",
        linewidth=2,
        color="c",
        hatch="",
        edgecolor="b",
        fill=False,
    )
    ax.set_xlabel("mV ns")
    ax.legend()
    ax.set_yscale("log")
plt.tight_layout()
plt.savefig(outdir + "/" + fileName + "_b5_3.png")

f, a = plt.subplots(4, 2)
a = a.ravel()
for idx, ax in enumerate(a):
    if idx + 25 == 33:
        break
    ax.hist(
        schn[idx + 88],
        bins=50,
        range=[0, 1],
        label="chn {}".format(idx + 25),
        histtype="step",
        linewidth=2,
        color="c",
        hatch="",
        edgecolor="b",
        fill=False,
    )
    ax.set_xlabel("mV ns")
    ax.legend()
    ax.set_yscale("log")
plt.tight_layout()
plt.savefig(outdir + "/" + fileName + "_b5_4.png")
