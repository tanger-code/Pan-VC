"""
将REF_GRCh38_chr21.txt中的参考node提取并写成一行一行的形式
之后进行vg find -N ref_nodes.txt -x chr21.xg -P GRCh38#chr21 > ref_nodes_pos.csv
将每个node的POS位置写入文件中
"""

import pandas as pd
import re
import time
import argparse

parser = argparse.ArgumentParser(description="ref_node_pos")
parser.add_argument('-c', '--chromosome', type=str, default="chr")
args = parser.parse_args()
chromosome = args.chromosome

ref_path_file = "./job_store/REF_GRCh38_" + chromosome + ".txt"
ref_nodes_file = "./job_store/ref_nodes.txt"

df = pd.read_csv(ref_path_file, sep='\t', header=None, names=['0', '1', '2', '3'])
# 提取nodes列
ref_nodes_str = df.iloc[0, 2]
# 保存所有结点，去掉+
nodes = re.findall(r"\d+", ref_nodes_str)
# 按行写入文件
with open(ref_nodes_file, 'w') as f:
    for node in nodes:
        f.write(f"{node}\n")
