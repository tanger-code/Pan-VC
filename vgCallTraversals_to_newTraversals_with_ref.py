"""
从vgCallTraversals文件中提取snarl名称、snarl路径、snarl方向、snarl路径的平均支持度、是否是snarl的参考路径
根据参考路径对提取的snarl名称、snarl路径进行矫正
"""
import json
import csv
import re
import argparse

parser = argparse.ArgumentParser(description="vgCallTraversals_to_newTraversals")
parser.add_argument('-c', '--chromosome', type=str, default="chr")
args = parser.parse_args()
chromosome = args.chromosome

vgCallTraversals_file = "./job_store/" + chromosome + ".vgCallTraversals.gaf"
ref_file = "./job_store/REF_GRCh38_" + chromosome + ".txt"
new_traversals_file = "./job_store/" + chromosome + ".new_traversals.csv"

# 从参考路径文件中提取参考路径
with open(ref_file, 'r') as f:
    # 文件中只有一行
    # P	GRCh38#chr21	42815844+,42815845+,42815846+,42815847+    *
    line = f.readline()
    line = line.strip().split()
    # line[2]是具体的path，根据"+,"进行分割后，得到path的node列表，但是最后一个node还有一个+号
    ref_path = line[2].split('+,')  # 42815844,42815845,42815846,42815847+
    # 去掉最后的+号
    ref_path[-1] = ref_path[-1][:-1]  # 42815844,42815845,42815846,42815847

ref_dict = {value: index for index, value in enumerate(ref_path)}

# 从vgCallTraversals文件中提取snarl名称、snarl路径、snarl路径的平均支持度
# 根据ref参考路径对snarl和snarl路径进行修正：[<67<64][67, 66, 64][<, >, <]==>[>64>67][64, 66, 67][>, <, >]
# 将有用信息存入chr21.new_traversals.csv
with open(vgCallTraversals_file, 'r') as f_read:
    with open(new_traversals_file, 'w', newline='') as f_write:
        # 写入header信息
        header = ['name', 'path', 'dir', 'depth', 'ref']
        csv_writer = csv.writer(f_write, delimiter='\t')
        csv_writer.writerow(header)

        # 遍历vgCallTraversals文件，提取有用信息，存入新文件
        for line in f_read:
            line = line.strip().split()
            # 第0部分是traversals的name，分割后snarl_name是第4个
            travel_name = line[0].split('#')
            snarl_name = travel_name[4]
            # 第5部分是traversal的path信息：>n1 >n2 >n3
            path_str = line[5]
            path = re.findall(r'\d+', path_str)
            dir = re.findall('>|<', path_str)
            # 平均reads支持度
            avg_depth = line[12].split(':')[-1]

            # 根据ref路径矫正信息
            # 提取snarl_name的两个node：snarl_orig_nodes和其方向snarl_orig_dir
            snarl_orig_nodes = re.findall(r'\d+', snarl_name)
            snarl_orig_dir = re.findall(r'<|>', snarl_name)

            # 在ref路径中找这两个node，比较索引的大小
            node1_idx = ref_dict[snarl_orig_nodes[0]]
            node2_idx = ref_dict[snarl_orig_nodes[1]]

            # 如果node1的索引>node2的索引，则应该矫正
            if node1_idx > node2_idx:
                # snarl_name的矫正
                snarl_now_nodes = snarl_orig_nodes[::-1]
                snarl_now_dir = ['<' if ch == '>' else '>' for ch in snarl_orig_dir[::-1]]
                new_snarl_name = snarl_now_dir[0] + snarl_now_nodes[0] + snarl_now_dir[1] + snarl_now_nodes[1]
                # snarl_path和dir的矫正
                new_path = path[::-1]
                new_dir = ['<' if ch == '>' else '>' for ch in dir[::-1]]
                # 判断snarl的这条路径是否是ref
                snarl_ref_path = ref_path[node2_idx:node1_idx + 1]
                if snarl_ref_path == new_path:
                    ref_flag = 0
                else:
                    ref_flag = 1
                csv_writer.writerow([new_snarl_name, new_path, new_dir, avg_depth, ref_flag])

            else:
                # 如果node1的索引<node2的索引，不用矫正
                snarl_ref_path = ref_path[node1_idx:node2_idx + 1]
                if snarl_ref_path == path:
                    ref_flag = 0
                else:
                    ref_flag = 1
                csv_writer.writerow([snarl_name, path, dir, avg_depth, ref_flag])
