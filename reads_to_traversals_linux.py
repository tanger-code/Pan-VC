"""
    >> 生成reads_name_path字典之后，先存入文件，然后从文件读取到内存，能《减少内存的使用》
    根据new_traversals文件，将比对到snarl path的reads名称加到对应路径之后
    >需要改动的参数
        >142：执行vg-chunk命令时，所选的step大小（10w）
        >193：log文件地址
        >200-205：各种文件地址
        >220：从graphaligner_aln_sort.json文件中提取read_name和read_path的chunk大小
        >224：220的pool数量
        >245：均匀打乱traversals文件的step大小
        >252：traversals文件的chunk大小
        >255：252的pool数量
"""
import pandas as pd
import json
import re
import time
import multiprocessing as mp
from functools import partial
from collections import defaultdict
import os
import logging
import argparse


# ------------------------------------------------判断snarl_path和read_path的重合----------------------------------------------
def path_overlap2(l_read, l_snarl):
    # 判断snarl_path第一个节点在reads_path中的位置
    idx = l_read.index(l_snarl[0])
    length = len(l_read) - idx
    if l_read[idx:] == l_snarl[:length]:
        return True
    if l_read[idx::-1] == l_snarl[:idx + 1]:
        return True
    return False


def path_overlap(snarl_path, reads_path, l_snarl, l_read):
    # snarl_path和reads_path都是字符串如['43632066', '43632064']，直接按字符串查找即可

    l_snarl_len = len(l_snarl)
    l_read_len = len(l_read)
    l_snarl_r = l_snarl[::-1]
    # reads的比对路径有正有反，都需要进行判定；找snarl和reads中较短的path进行正反判定
    if (l_snarl_len < l_read_len):
        # snarl短，则对snarl_path取正反
        snarl_path1 = snarl_path[1:-1]  # 将两端的[ ]去掉
        snarl_path2 = str(l_snarl_r)[1:-1]
        if (snarl_path1 in reads_path) or (snarl_path2 in reads_path):
            return True
    else:
        # read短，则对reads_path取正反
        l_read_r = l_read[::-1]
        read_path1 = reads_path[1:-1]
        read_path2 = str(l_read_r)[1:-1]
        if (read_path1 in snarl_path) or (read_path2 in snarl_path):
            return True

    # 左右都没有重叠，则返回false
    if (l_read[0] not in l_snarl) and (l_read[-1] not in l_snarl):
        return False
    else:
        # 如果snarl_path第一个节点或最后一个节点在reads_path中，进行进一步的判定
        if l_snarl[0] in l_read:
            return path_overlap2(l_read, l_snarl)
        elif l_snarl[-1] in l_read:
            return path_overlap2(l_read, l_snarl_r)
        # 如果都不在，则snarl_path与reads_path不重叠
        else:
            return False


def path_overlap_all(snarl_path, reads_path, l_snarl, l_read):
    # 对于小变异，只考虑snarl_path完全在reads_path中或者reads_path完全在snarl_path中
    # reads的比对路径有正有反，都需要进行判定
    # snarl_path和reads_path都是字符串如['43632066', '43632064']，直接按字符串查找即可

    if len(l_snarl) < len(l_read):
        # 计算snarl_path的正反
        l_snarl_r = l_snarl[::-1]
        snarl_path1 = snarl_path[1:-1]  # 将两端的[ ]去掉
        snarl_path2 = str(l_snarl_r)[1:-1]
        # 查找snarl路径的字符串是否包含在reads比对路径的字符串中
        if (snarl_path1 in reads_path) or (snarl_path2 in reads_path):
            return True
        else:
            return False
    else:
        if (l_read[0] not in l_snarl) or (l_read[-1] not in l_snarl):
            return False
        l_read_r = l_read[::-1]
        read_path1 = reads_path[1:-1]
        read_path2 = str(l_read_r)[1:-1]
        if (read_path1 in snarl_path) or (read_path2 in snarl_path):
            return True
        else:
            return False


def chunk_reads_extract_chunk_work(folder_path, file, chromosome):
    logging.info(f"---chunk_reads_extract---")

    # 提取文件中包含的所有reads
    chunk_reads = []
    with open(folder_path + file, 'r') as f:
        for line in f:
            line = line[:100]
            match = re.search(r'"name":\s*"([^"]+)"', line)
            # 找到chunk中的reads名称
            read_name = match.group(1)
            chunk_reads.append(read_name)

    # file形式：chunk_101_GRCh38#0#chr21_10100000_10199999.json
    match = re.search(chromosome + r'_(\d+_\d+)', file)
    # 10100000_10199999
    chunk_range = match.group(1)
    parts = chunk_range.split("_")
    parts[1] = str(int(parts[1]) + 1)
    # 10100000_10200000
    chunk_range = "_".join(parts)

    return [chunk_range, chunk_reads]


def chunk_reads_dict_file_write(chunk_info, chunk_reads_dict):
    chunk_range = chunk_info[0]
    chunk_reads = chunk_info[1]
    chunk_reads_dict[chunk_range] = chunk_reads


# -----------------------------将chr21.new_traversals.csv文件分成chunk，每个chunk的工作-----------------------------------------
def chunk_work(df_traver_chunk, ref_nodes_pos_df, read_name_path_dict, chunk_reads_dict):
    print('!!-------------------------------------------------')
    df_traver_chunk['support_reads'] = ''
    for idx, row in df_traver_chunk.iterrows():
        logging.info(f"{idx}")
        # 1、提取traversals每一行中的snarl节点
        snarl_name = row['name']
        snarl_path = row['path']
        snarl_name_nodes = re.findall(r"\d+", snarl_name)
        l_snarl_path = re.findall(r"\d+", snarl_path)

        # 2、判断snarl节点在参考基因组中的位置
        snarl_start_pos = ref_nodes_pos_df[ref_nodes_pos_df['node'] == snarl_name_nodes[0]].iloc[0, 1]
        # snarl_end_pos = ref_nodes_pos_df[ref_nodes_pos_df['node'] == snarl_name_nodes[-1]].iloc[0, 1]

        # 3、根据这个位置将snarl路径定位到对应的gam的chunk中
        # gam分割时，每step分割一次
        step = 100000
        start_chunk = snarl_start_pos // step
        # 需要判断是否是最后一个chunk，最后一个chunk范围不是10w的倍数，需要截断处理
        end_pos = (start_chunk + 1) * step
        if end_pos >= chromosome_len:
            start_chunk_range = f"{start_chunk * step}_{chromosome_len}"
        else:
            start_chunk_range = f"{start_chunk * step}_{end_pos}"
        # start_chunk_file_name = f"../data_HG002/chunks/chunk_0_GRCh38#0#chr21_{start_chunk * step}_{(start_chunk + 1) * step}.json"

        # 4、将snarl路径和gam的chunk中的reads路径进行比对
        # 如果snarl涉及到多个chunk（最多2个），则每个chunk都含有比对到snarl上的reads，大多数是重复的；所以只考虑第一个chunk即可
        # 遍历chunk文件，将与snarl path重合的read path对应的reads名称加到对应snarl path后
        reads_list = []
        chunk_reads = chunk_reads_dict[start_chunk_range]
        for read_name in chunk_reads:
            # 然后在dict中找到其对应的路径
            read_paths = read_name_path_dict[read_name]
            # 最后判断snarl路径和read路径是否重叠
            for read_path in read_paths:
                # 忽略只有一个节点的比对路径
                if len(read_path) == 1:
                    continue
                # 找到就可以停止
                read_path_str = str(read_path)
                if path_overlap(snarl_path, read_path_str, l_snarl_path, read_path):
                    reads_list.append(read_name)
                    break
        df_traver_chunk.loc[idx, 'support_reads'] = str(reads_list)
    return df_traver_chunk


# -----------------------------将gam文件分成chunk提取read_name和read_path，每个chunk的工作---------------------------------------
def read_name_path_chunk_work(gam_chunk):
    read_names = []
    read_paths = []
    # 提取reads名称和对应的比对路径，存储到字典中
    for line in gam_chunk:
        logging.info("1")
        line = json.loads(line)
        read_name = line["name"]
        read_path = []
        for item in line["path"]["mapping"]:
            read_path.append(item["position"]["node_id"])
        read_names.append(read_name)
        read_paths.append(read_path)
    return [read_names, read_paths]


def file_write_callback(read_names_paths, read_name_path_dict):
    read_names = read_names_paths[0]
    read_paths = read_names_paths[1]
    for i in range(len(read_names)):
        read_name_path_dict[read_names[i]].append(read_paths[i])


log_file = "./job_store/log_reads_to_traversals.txt"  # log文件地址
logging.basicConfig(filename=log_file, level=logging.INFO)
parser = argparse.ArgumentParser(description="ref_node_pos")
parser.add_argument('-c', '--chromosome', type=str, default="chr")
parser.add_argument('-l', '--chromosome_len', type=str, default="1")
args = parser.parse_args()
chromosome = args.chromosome
chromosome_len = int(args.chromosome_len)

if __name__ == '__main__':
    logging.info("start")
    print(time.localtime())
    time1 = time.time()
    # -----------------------------------------------------文件处理-----------------------------------------------------------
    gam_file_path = "./job_store/graphaligner_aln_sort.json"
    traversals_file = "./job_store/" + chromosome + ".new_traversals.csv"
    ref_nodes_pos_file = "./job_store/ref_nodes_pos.csv"
    folder_path = "./chunks/"  # 存放gam的chunk文件的文件夹地址，注意最后有个/，之后拼接路径用的
    traversals_shuffle_file = "./job_store/" + chromosome + "_traversals_shuffle.csv"  # 将traversals文件进行重新均匀排序，防止最后的snarl聚集在比较大的gam的chunk文件中，拖慢并行效果
    traversals_with_reads_file = "./job_store/traversals_with_reads_10w.csv"  # 最终生成的文件地址
    read_name_path_dict_file = "./job_store/read_name_path_dict.txt"  # 将read_name_path_dict存入文件，然后再读入，释放不必要的内存
    chunk_reads_dict_file = "./job_store/chunk_reads_dict.txt"
    
    # 从ref_nodes_pos.csv中加载每个ref_node的POS位置
    logging.info("ref_nodes_pos_df")
    ref_nodes_pos_df = pd.read_csv(ref_nodes_pos_file, sep='\t', header=None, names=['node', 'pos'])
    ref_nodes_pos_df['pos'] = ref_nodes_pos_df['pos'] + 1  # 准换为1-base的坐标
    ref_nodes_pos_df['node'] = ref_nodes_pos_df['node'].astype(str)  # 将node列的数据类型转换为str

    # ----------多进程提取read_name、read_path存储到read_name_path_dict中--------------------------------
    logging.info("read_name_path_dict")
    # 提取reads名称和对应的比对路径，存储到字典中
    read_name_path_dict_pre = defaultdict(list)
    # 将dataframe进行分块
    with open(gam_file_path, 'r') as f:
        lines_gam = f.readlines()
        chunk_step = 10000
        chunks_gam = [lines_gam[i:i + chunk_step] for i in range(0, len(lines_gam), chunk_step)]
    # 每个chunk一个进程
    new_file_write_callback = partial(file_write_callback, read_name_path_dict=read_name_path_dict_pre)
    with mp.Pool(processes=40) as pool:
        for chunk_gam in chunks_gam:
            pool.apply_async(read_name_path_chunk_work, (chunk_gam,), callback=new_file_write_callback)
        pool.close()
        pool.join()
    # 将read_name_path_dict_pre存入文件，然后再读入内存（能减少内存使用，python机制未知）
    with open(read_name_path_dict_file, 'w') as f:
        for k, v in read_name_path_dict_pre.items():
            f.write(f"{k}\t{v}\n")
    read_name_path_dict = dict()
    with open(read_name_path_dict_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            read_name_path_dict[line[0]] = eval(line[1])
    # -----------------------------------------------------------------------------------------------

    # 读取每个gam的chunk文件，并将chunk所包含的reads名称全提取出来，存入字典
    # 这样，之后就不用每个snarl_path都打开一次文件，然后提取reads，有很多冗余操作
    logging.info("chunk_reads_dict")
    chunk_reads_dict_pre = dict()
    new_chunk_reads_dict_file_write = partial(chunk_reads_dict_file_write, chunk_reads_dict=chunk_reads_dict_pre)
    json_files = [f for f in os.listdir(folder_path) if f.endswith('.json')]
    with mp.Pool(processes=40) as pool:
        for json_file in json_files:
            pool.apply_async(chunk_reads_extract_chunk_work, (folder_path, json_file, chromosome),
                             callback=new_chunk_reads_dict_file_write)
        pool.close()
        pool.join()
    # 将字典保存到文件中
    with open(chunk_reads_dict_file, 'w') as f:
        for k, v in chunk_reads_dict_pre.items():
            f.write(f"{k}\t{v}\n")
    chunk_reads_dict = dict()
    with open(chunk_reads_dict_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            chunk_reads_dict[line[0]] = eval(line[1])
    # ----------------------------------------------------------------------------------------------------------------------

    logging.info("traversals_with_reads")
    # 将traversals文件进行重新均匀排序，防止最后的snarl聚集在比较大的gam的chunk文件中，拖慢并行效果
    with open(traversals_file, 'r') as f:
        lines = f.readlines()
        header = lines[0]
        with open(traversals_shuffle_file, "w") as f_w:
            f_w.write(header)
            step = 1000
            total_len = len(lines)
            for i in range(1, step + 1):
                tmp_lines = lines[i:total_len:step]
                for line in tmp_lines:
                    f_w.write(line)
    # 将dataframe进行分块
    chunks = pd.read_csv(traversals_shuffle_file, sep='\t', chunksize=4000)
    new_chunk_work = partial(chunk_work, ref_nodes_pos_df=ref_nodes_pos_df, read_name_path_dict=read_name_path_dict,
                             chunk_reads_dict=chunk_reads_dict)
    with mp.Pool(processes=60) as pool:
        res_dfs = pool.map(new_chunk_work, chunks)
        result_df = pd.concat(res_dfs)
        result_df.to_csv(traversals_with_reads_file, sep='\t', index=False)

    time2 = time.time()
    print(time.localtime())
    print(time2 - time1)
