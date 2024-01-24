"""
将traversals_with_reads_agg.csv中的每行按照行首node在REF路径中出现的顺序或者在ref_node_pos中出现的顺序进行排序
这样，之后进行snarl的拼接时，才能正确拼接序列
"""
import pandas as pd
import re

traversals_with_reads_file = "./job_store/traversals_with_reads_10w.csv"
traversals_with_reads_agg_file = "./job_store/traversals_with_reads_agg.csv"
traversals_with_reads_agg_sort_by_pos_file = "./job_store/traversals_with_reads_agg_sort_by_pos.csv"
ref_nodes_pos_file = "./job_store/ref_nodes_pos.csv"

#  -----------------------------traversals_with_reads文件按照snarl name进行聚合--------------------------------------------
# 将snarl_name相同的记录进行合并：
traversals_df = pd.read_csv(traversals_with_reads_file, sep='\t')
# groupby函数按照name进行分组，agg函数进行聚合，将snarl_name相同的记录进行合并
traversals_df = traversals_df.groupby('name', sort=False).agg(lambda x: ', '.join(x.astype(str))).reset_index()
traversals_df.to_csv(traversals_with_reads_agg_file, index=False, sep='\t')


# ----------------------------将traversals聚合之后的文件按照snarl在REF中的顺序进行排序-------------------------------------
traversals_df = pd.read_csv(traversals_with_reads_agg_file, sep='\t')

# 按照位置排序的所有的ref_nodes
ref_nodes_all = []
with open(ref_nodes_pos_file, 'r') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip().split('\t')
        node = line[0]
        ref_nodes_all.append(node)

# traversals中存在的ref_nodes
ref_node_in_traversals = []
with open(traversals_with_reads_agg_file, 'r') as f:
    # 跳过第一行的header
    f.readline()
    lines = f.readlines()
    for line in lines:
        line = line.strip().split('\t')
        temp_node = re.findall(r"\d+", line[0])[0]
        ref_node_in_traversals.append(temp_node)
ref_node_in_traversals_set = set(ref_node_in_traversals)

# 从ref_nodes_all中提取出现在ref_node_in_traversals中的node，保持原有顺序
ref_node_in_traversals_ordinal = []
for node in ref_nodes_all:
    if node in ref_node_in_traversals_set:
        ref_node_in_traversals_ordinal.append(node)

# 将snarl中的>node1>node2扩展为两列node1和node2
traversals_df[['', 'node1', 'node2']] = traversals_df['name'].str.split('>', expand=True)
# 上一步会多一列空行，这里删除它
traversals_df = traversals_df.drop(columns=[''])
# 以node1列为行索引
traversals_df['node1'] = traversals_df['node1'].astype(str)
traversals_df = traversals_df.set_index('node1')
# 按照node1在将ref_node_in_traversals_ordinal中的顺序进行取行排序，然后释放索引
traversals_df = traversals_df.loc[ref_node_in_traversals_ordinal].reset_index()
# 删除多余的列
traversals_df = traversals_df.drop(columns=['node1', 'node2'])

# 将排序好的记录写入文件
traversals_df.to_csv(traversals_with_reads_agg_sort_by_pos_file, index=False, sep='\t')
