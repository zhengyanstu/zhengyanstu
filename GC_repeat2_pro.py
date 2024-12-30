# -*- coding:utf-8 -*-
import os
import numpy as np
import random
import argparse
from scipy import stats
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from pandas import DataFrame


# 定义命令行参数
def parse_arguments():
    parser = argparse.ArgumentParser(description="GC含量与重复序列分析")

    # 文件路径相关的参数
    parser.add_argument('--repeat_file', type=str, required=True, help="重复序列文件的路径")
    parser.add_argument('--fasta_file', type=str, required=True, help="FASTA文件的路径")

    # 其他参数
    parser.add_argument('--flank_length', type=int, default=500, help=" flank的长度，默认值为500")
    parser.add_argument('--window_size', type=int, default=50, help="计算GC含量的窗口大小，默认值为50")
    parser.add_argument('--step_size', type=int, default=30, help="滑动窗口的步长，默认值为30")
    parser.add_argument('--length_min', type=int, default=30, help="最小重复序列长度，默认值为30")
    parser.add_argument('--length_max', type=int, default=9999, help="最大重复序列长度，默认值为9999")

    return parser.parse_args()


# F1 - 读取重复序列位置
def read_repeat(repeat, length_max, length_min):
    with open(repeat, "r") as r:
        repeat_pos = []
        data = r.readlines()
        for line in data:
            if line.startswith("Chr"):
                info_list = line.split()
                self_chr = info_list[0].strip()
                self_start = info_list[1].strip()
                self_end = info_list[2].strip()
                mate_chr = info_list[3].strip()
                mate_start = info_list[4].strip()
                mate_end = info_list[5].strip()
                length = int(info_list[7].strip())
                # 过滤
                if int(length_max) >= length >= int(length_min):
                    repeat_pos.append(self_chr + "_" + self_start + "_" + self_end)
                    repeat_pos.append(mate_chr + "_" + mate_start + "_" + mate_end)
    return repeat_pos


# 读取fasta文件
def read_fasta(fasta):
    fasta_dict = {}
    with open(fasta, "r") as f:
        data = f.readlines()
        for line in data:
            if line.startswith(">"):
                seq_name = line.strip()[1:]
                fasta_dict[seq_name] = []
            else:
                fasta_dict[seq_name].append(line.strip())
    for k, v in fasta_dict.items():
        fasta_dict[k] = "".join(v)
    return fasta_dict


# 计算GC含量
def GC_cal(fasta):
    A = fasta.count("A")
    T = fasta.count("T")
    G = fasta.count("G")
    C = fasta.count("C")
    gc_content = round(float(G + C) / float(A + T + G + C) * 100, 2)
    result = float(gc_content)
    return result


# F2 - 计算窗口GC含量
def window_gc_dict(fasta, window_size, step_size):
    window_dict = {}
    for chr_num in fasta.keys():
        window_dict[chr_num] = []
        for i in range(0, len(fasta[chr_num]), int(step_size)):
            if i + int(window_size) < len(fasta[chr_num]):
                window_start = i
                window_end = i + int(window_size)
                window_gc = GC_cal(fasta[chr_num][window_start:(window_end + 1)])
                win_info = [str(window_start), str(window_end), str(window_gc)]
                window_dict[chr_num].append("_".join(win_info))
            elif i + int(window_size) >= len(fasta[chr_num]):
                window_start = i
                window_end = len(fasta[chr_num])
                window_gc = GC_cal(fasta[chr_num][window_start:(window_end + 1)])
                win_info = [str(window_start), str(window_end), str(window_gc)]
                window_dict[chr_num].append("_".join(win_info))
                if window_end - window_start + 1 >= 25:
                    window_dict[chr_num].append("_".join(win_info))
                    break
    return window_dict


# F3 - 获取窗口位置的GC含量列表
def pos_gc_list(repeat_pos, fasta, window_dict, flank):
    flank_window_list = []
    for r in repeat_pos:
        chr_num = str(r).split("_")[0]
        start = str(r).split("_")[1]
        end = str(r).split("_")[2]
        # left flank
        fend = int(flank) - int(start)
        # right flank
        fstart = (int(end) + int(flank)) - len(fasta[chr_num])
        for w in window_dict[chr_num]:
            # left flank
            window_start = int(w.split("_")[0])
            if int(start) >= window_start >= int(start) - int(flank):
                flank_window_list.append(chr_num + "_" + w)
            elif int(start) >= window_start >= len(fasta[chr_num]) - fend:
                flank_window_list.append(chr_num + "_" + w)
            # right flank
            window_end = int(w.split("_")[1])
            if int(end) <= window_end <= int(end) + int(flank):
                flank_window_list.append(chr_num + "_" + w)
            elif int(end) <= window_end <= fstart:
                flank_window_list.append(chr_num + "_" + w)

    return list(set(flank_window_list))


def get_gc_from_list(alist):
    gc_list = []
    for i in alist:
        gc_list.append(float(i.split("_")[3]))
    return gc_list


# 主程序
def main():
    # 解析命令行参数
    args = parse_arguments()

    # 读取数据
    active_repeat_pos = read_repeat(args.repeat_file, args.length_max, args.length_min)
    unactive_repeat_pos = read_repeat(args.repeat_file, 99, 30)
    my_fasta = read_fasta(args.fasta_file)
    all_windows = window_gc_dict(my_fasta, args.window_size, args.step_size)

    # 获取GC含量列表
    active_repeat_gc_list = pos_gc_list(active_repeat_pos, my_fasta, all_windows, args.flank_length)
    unactive_repeat_gc_list = pos_gc_list(unactive_repeat_pos, my_fasta, all_windows, args.flank_length)

    # 删除重复项
    for i in active_repeat_gc_list:
        if i in unactive_repeat_gc_list:
            unactive_repeat_gc_list.remove(i)

    all_windows_gc_list = []
    for k in all_windows.keys():
        for i in all_windows[k]:
            all_windows_gc_list.append(k + "_" + i)
    all_rest_windows_gc_list = []
    for i in all_windows_gc_list:
        if i not in active_repeat_gc_list and i not in unactive_repeat_gc_list:
            all_rest_windows_gc_list.append(i)

    # 获取GC值
    actvie_gc_list = get_gc_from_list(active_repeat_gc_list)
    unactvie_gc_list = get_gc_from_list(unactive_repeat_gc_list)
    rest_gc_list = get_gc_from_list(all_rest_windows_gc_list)
    all_gc_list = rest_gc_list + unactvie_gc_list + actvie_gc_list

    # 随机GC均值
    random_gc_mean_list = []
    for i in range(10000):
        random_gc_list = random.sample(all_gc_list, int(len(unactvie_gc_list)))
        random_gc_mean_list.append(np.mean(random_gc_list))

    # 比较结果
    extreme_count = 0
    for i in random_gc_mean_list:
        if float(i) < float(np.mean(unactvie_gc_list)):
            extreme_count = extreme_count + 1
    print("极端情况计数：", extreme_count)


if __name__ == "__main__":
    main()
