# -*- coding:utf-8 -*-
import os
import random
import numpy as np


os.chdir("/Users/joe/Desktop/向日葵列当线粒体重排/GC含量/断点附近GC含量和随机抽/两两Mauve求边界")

boundary_file = "A_v_B_A_Boundary"
fasta_file = "A_mt.fasta"
chr_list = "A_chr_list.txt"
proximity = 50
window = 100
step = 50


def read_boundary(boundary):
    temp_list = []
    with open(boundary, "r") as b:
        data = b.readlines()
        line_count = 0
        for line in data:
            if line_count > 0:
                chr = line.split()[1].strip()
                start = line.split()[2].strip()
                end = line.split()[3].strip()
                temp_list.append(chr + "_" + start)
                temp_list.append(chr + "_" + end)
            line_count = line_count + 1
    # remove adjacent point (shorter than twice the testing window size:100)
    boundary_pos = []
    key_list = []
    for chr_pos in temp_list:
        chr = chr_pos.split("_")[0]
        pos = chr_pos.split("_")[1]
        key = int(pos) // 100
        if key not in key_list:
            boundary_pos.append(chr + "_" + pos)
            key_list.append(key)
    return boundary_pos


# read fasta
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


# compute GC
def GC_cal(fasta):
    A = fasta.count("A")
    T = fasta.count("T")
    G = fasta.count("G")
    C = fasta.count("C")
    gc_content = round(float(G + C) / float(A + T + G + C) * 100, 2)
    result = float(gc_content)
    return result


# concatenate fasta dict
def concatenate_fasta(fa_dict):
    concatenate_list = []
    for v in fa_dict.values():
        concatenate_list.append(v)
    return "".join(concatenate_list)


def median_quarter(fasta, window_size, step):
    window_gc_list = []
    for i in range(0, len(fasta), int(step)):
        if i + int(window_size) < len(fasta):
            window_fasta = fasta[i:i + int(window_size)]
            window_gc_list.append(GC_cal(window_fasta))
        elif i + int(window_size) >= len(fasta):
            window_fasta = fasta[i:]
            window_gc_list.append(GC_cal(window_fasta))
            break
    gc_median = np.median(window_gc_list)
    lower_q = np.quantile(window_gc_list, 0.01, interpolation='lower')
    higher_q = np.quantile(window_gc_list, 0.01, interpolation='higher')
    return [lower_q, gc_median, higher_q]


# random boundary
def random_boundary(chr_list, boundary_list):
    rand_list = []
    chr_dict = {}
    with open(chr_list, "r") as c:
        chrs = c.readlines()
        for chr in chrs:
            chr_name = chr.split()[0].strip()
            chr_length = chr.split()[1].strip()
            chr_dict[chr_name] = int(chr_length)
        # random chr
        for i in range(len(boundary_list)):
            chr = random.sample(chr_dict.keys(), 1)
            # random pos
            pos = random.randint(0, chr_dict[chr[0]])
            rand_list.append(chr[0] + "_" + str(pos))
    return rand_list


# get pos have low gc
def gc_pos(fa_dict, window_size, step, lower_q):
    gc_pos_list = []
    for chr, seq in fa_dict.items():
        for i in range(0, len(seq), int(step)):
            if i + int(window_size) < len(seq):
                window_fasta = seq[i:i + int(window_size)]
                window_gc = GC_cal(window_fasta)
                if float(window_gc) <= float(lower_q):
                    gc_pos_list.append(chr + "_" + str(i) + "_" + str(i + int(window_size)))
            elif i + int(window_size) >= len(seq):
                window_fasta = seq[i:]
                window_gc = GC_cal(window_fasta)
                if float(window_gc) <= float(lower_q):
                    gc_pos_list.append(chr + "_" + str(i) + "_" + str(i + int(window_size)))
                break
    return gc_pos_list


# def CoOccurrence(gc_pos_list, boundary_list, proximity):
#     vicinity_count = 0
#     for low_gc in gc_pos_list:
#         low_gc_chr = low_gc.split("_")[0].strip()
#         low_gc_start = low_gc.split("_")[1].strip()
#         low_gc_end = low_gc.split("_")[2].strip()
#         for boundary in boundary_list:
#             boundary_chr = boundary.split("_")[0].strip()
#             boundary_pos = boundary.split("_")[1].strip()
#             if boundary_chr == low_gc_chr:
#                 # left
#                 if abs(int(low_gc_start) - int(boundary_pos)) < proximity:
#                     vicinity_count = vicinity_count + 1
#                 # right
#                 if abs(int(low_gc_end) - int(boundary_pos)) < proximity:
#                     vicinity_count = vicinity_count + 1
#     return float(vicinity_count / float(len(gc_pos_list)))

def CoOccurrence(gc_pos_list, boundary_list, proximity):
    vicinity_count = 0
    for boundary in boundary_list:
        boundary_chr = boundary.split("_")[0].strip()
        boundary_pos = boundary.split("_")[1].strip()
        for low_gc in gc_pos_list:
            low_gc_chr = low_gc.split("_")[0].strip()
            low_gc_start = low_gc.split("_")[1].strip()
            low_gc_end = low_gc.split("_")[2].strip()
            if boundary_chr == low_gc_chr:
                # left
                if abs(int(low_gc_start) - int(boundary_pos)) < proximity:
                    vicinity_count = vicinity_count + 1
                # right
                if abs(int(low_gc_end) - int(boundary_pos)) < proximity:
                    vicinity_count = vicinity_count + 1
    return float(vicinity_count / float(len(gc_pos_list)))


# the position of low gc
fasta = read_fasta(fasta_file)
con_fasta = concatenate_fasta(fasta)
lower_q = median_quarter(con_fasta, window, step)[0]
lower_q_pos = gc_pos(fasta, window, step, lower_q)
print(len(lower_q_pos))
# the position of my boundary
my_boundary_list = read_boundary(boundary_file)
my_low_gc_pro = CoOccurrence(lower_q_pos, my_boundary_list, proximity)
print(my_low_gc_pro)


# sample 1000 replicates with the same size
replicate_list = []
for i in range(10000):
    random_boundary_list = random_boundary(chr_list, my_boundary_list)
    #print(CoOccurrence(lower_q_pos, random_boundary_list, proximity))
    replicate_list.append(CoOccurrence(lower_q_pos, random_boundary_list, proximity))
# compare
extreme_count = 0
for i in replicate_list:
    if float(i) > float(my_low_gc_pro):
            extreme_count = extreme_count + 1
print(extreme_count)
