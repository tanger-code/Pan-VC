"""
· 基因分型--总，即SNP和INDEL一起进行分型
· 用的是带有HG002的数据
· 主要改变：
    >>将通过node进行检测变异和分型--->改为通过比对直接进行分型
    >>如果两端ref节点有一个小于20bp，则向前或者向后进行扩展，将多个snarl进行拼接，联合进行检测和分型
    >去除了重复卫星序列处的变异，因为这些区域reads测序和比对质量都一般
    >不处理1/2的情况，而是将1和2都和ref进行比较
"""
import multiprocessing as mp
import pandas as pd
import re
import math
import json
import time
import logging
from vcf_format import VcfRecord, vcf_header
from alignment import variant_call
import argparse


# # 根据snarl所有路径所对应的reads，去重合并，得到此snarl处的所有reads
def extract_snarl_all_reads(support_reads):
    # 根据snarl所有路径所对应的reads，去重合并，得到此snarl处的所有reads
    # support_reads格式str：[read1, read2, read3], [read2, read3, read5], [...], [...]
    # 将support_reads转换为列表形式
    support_reads_list = eval('[' + support_reads + ']')
    snarl_reads = []
    for path_reads in support_reads_list:
        for read in path_reads:
            # 一个read基本上只会比对到同一snarl的一个路径上，不会出现同时比对到同一snarl的多个路径上，所以多个snarl_path的support_reads并不会有重叠
            if read not in snarl_reads:
                snarl_reads.append(read)
    return snarl_reads


# 计算基因型Gxy的似然
def likelihood_calculate(allele_x_reads, allele_y_reads, snarl_all_reads):
    # 判断每一个read对Gxy的支持度L(r | Gxy)，并计算最终似然L(R | Gxy)
    # 注，这里的likelihood记录的是指数部分，即（1/2）的次幂，这个值越大，总体似然值越小，最后返回-likelihood即可
    likelihood = 0
    for read in snarl_all_reads:
        if read in allele_x_reads:
            support_x = True
        else:
            support_x = False
        if read in allele_y_reads:
            support_y = True
        else:
            support_y = False

        if support_x and support_y:
            likelihood = likelihood
        elif support_x or support_y:
            likelihood = likelihood + 1
        else:
            likelihood = likelihood + 5

    return -likelihood


# 将图的GFA格式文件以pandas读取，并提取'S'--结点的DataFrame
def gfa_transform(filePath):
    # 读取gfa文件，记录以S开头的行数，即node行数
    node_rows = 0
    with open(filePath, 'r') as f:
        f.readline()
        for line in f:
            if line.startswith('S'):
                node_rows += 1

    # 提取node所有行，保存为DataFrame返回
    names = ['0', '1', '2']
    node_df = pd.read_csv(filePath, sep='\t', header=None, names=names, skiprows=1, nrows=node_rows)
    node_df.columns = ['flag', 'node', 'seq']

    return node_df


# 计算序列的反补序列
def dna_reverse_complement(sequence):
    trantab = str.maketrans('ACGT', 'TGCA')  # 制作翻译表
    rc_seq = sequence[::-1].translate(trantab)  # 转换字符
    return rc_seq


# 判断变异处是否是卫星序列，重复序列>10次认为是卫星序列
def whether_satellite_sequence(temp_pos, temp_alt):
    # 记录是否是卫星序列，如果是，记录变异起始位置和变异结束位置
    flag = 0
    start_pos = 0
    end_pos = 0

    # 注，有时可能有一个很大的INS或者DEL（可能上万），这时不去考虑它是串联重复序列；如果考虑的话，会出现很多需要考虑的子序列
    if len(temp_alt) > 500:
        return flag, start_pos, end_pos

    # 定位变异在fa文件中的位置：行索引+offset
    line_num = temp_pos // 60
    line_offset = temp_pos % 60 - 1

    # 识别重复模块
    require_seq = temp_alt[1:]
    pattern = re.compile(r'(?=(\w+)(?=\1))')
    matches = pattern.findall(require_seq)  # 只识别除了第一个碱基之外的碱基序列
    matches_truth = set([x for x in matches if require_seq.startswith(x)])
    matches_truth.add(require_seq)  # 添加temp_alt[1:]，增加完备性
    matches_truth.add(temp_alt[0])  # 添加temp_alt[0]，增加完备性

    # 判断此处是不是卫星序列
    parent_seq_behind = chr21_fa_lines[line_num] + chr21_fa_lines[line_num + 1]
    parent_seq_front = chr21_fa_lines[line_num - 1] + chr21_fa_lines[line_num]
    for dup_seq in matches_truth:
        # 从fa序列中分别向前后延申一行（60bp，可以检测6bp大小的重复，可以继续延申），看是否是真的重复序列
        dup_seq_len = len(dup_seq)
        # 向后延申
        behind_num = 1
        for i in range(line_offset + 1, len(parent_seq_behind), dup_seq_len):
            temp_seq = parent_seq_behind[i:i + dup_seq_len]
            if temp_seq == dup_seq:
                behind_num += 1
            else:
                end_pos = temp_pos + (behind_num + 1) * dup_seq_len
                break
        # 如果最后end_pos为0，说明往后取的一行序列全是参考序列的重复，没有进入else的逻辑里
        # 此时，应该继续往后取行；或者将pos位置往后多移几十个bp
        if end_pos == 0:
            end_pos = temp_pos + (behind_num + 1) * dup_seq_len + 50
        # 向前延申
        front_num = 0
        for i in range(60 + line_offset + 1, 0, -dup_seq_len):
            temp_seq = parent_seq_front[i - dup_seq_len: i]
            if temp_seq == dup_seq:
                front_num += 1
            else:
                start_pos = temp_pos - (front_num + 1) * dup_seq_len
                break
        if start_pos == 0:
            start_pos = temp_pos - (front_num + 1) * dup_seq_len - 50
        # 如果重复次数>10，则认为是重复卫星序列
        if behind_num + front_num >= 5:
            flag = 1
            break

    return flag, start_pos, end_pos


# -----------------------------------------------------------------------------------------------------------------------
# 进行具体的比对检测变异，返回相关POS, REF, ALT, INFO信息
def genotyping_all_normal(snarl_ref_path, snarl_alt_path, snarl_alt_path_dir, AD_ref, AD_alt):
    bed_records = []
    # 需要记录的变异信息
    REF = []
    ALT = []
    POS = []
    INFO = ""

    # snarl开始的位置，即最左边的位置
    snarl_start_pos = ref_nodes_pos_df[ref_nodes_pos_df['node'] == snarl_ref_path[0]].iloc[0, 1]
    # 参考路径序列和变异序列
    ref_seq = ""
    for node in snarl_ref_path:
        ref_seq += node_df[node_df['node'] == int(node)].iloc[0, 2]
    alt_seq = ""
    for node, node_dir in zip(snarl_alt_path, snarl_alt_path_dir):
        if node_dir == '>':
            alt_seq += node_df[node_df['node'] == int(node)].iloc[0, 2]
        else:
            alt_seq += dna_reverse_complement(node_df[node_df['node'] == int(node)].iloc[0, 2])

    # 进行比对检测变异和分型,得到变异的参考序列和变异序列
    ref_variants, alt_variants, temp_pos = variant_call(ref_seq, alt_seq)
    # 上面的temp_pos为snarl内部的相对位置，需要加上snarl_start_pos
    pos_variants = [x + snarl_start_pos - 1 for x in temp_pos]

    # *****************************************************************************************
    # 对于ATTATATA--ATATA这种变异，比对的结果中变异在开头即：'ATT'--''，会出现空字符的情况
    # 这时需要在前面加上一个碱基；且这种情况只会在开头第一个变异出现
    if len(ref_variants) > 0:
        if ref_variants[0] == '' or alt_variants[0] == '':
            # 详见variant_call函数定义处，开头是这种indel时，pos需要加上1
            pos_variants[0] += 1
            # 定位变异在fa文件中的位置：行索引+offset
            line_num = pos_variants[0] // 60
            line_offset = pos_variants[0] % 60 - 1
            # 这里往前取一行，防止变异发生行首位置（此时应该加上上一行最后一个碱基）
            parent_seq_front = chr21_fa_lines[line_num - 1] + chr21_fa_lines[line_num]
            # 59 + line_offset 即 60 + line_offset -1，减去1代表往前取一个碱基
            add_one_seq = parent_seq_front[59 + line_offset]

            ref_variants[0] = add_one_seq + ref_variants[0]
            alt_variants[0] = add_one_seq + alt_variants[0]
            pos_variants[0] -= 1
    # *****************************************************************************************

    # **************************对其中的每个变异进行判断是否是微卫星序列**************************
    for ref_v, alt_v, pos_v in zip(ref_variants, alt_variants, pos_variants):
        satellite_flag = 0
        start_pos = 0
        end_pos = 0
        # 只对indel进行卫星序列的判断,snp不用
        if len(ref_v) == 1 and len(alt_v) > 1:
            satellite_flag, start_pos, end_pos = whether_satellite_sequence(pos_v, alt_v)
        if len(ref_v) > 1 and len(alt_v) == 1:
            satellite_flag, start_pos, end_pos = whether_satellite_sequence(pos_v, ref_v)

        # 如果是卫星序列，则跳过这个变异；记录pos位置，用于生成bed文件，去除真值集中这些区域的变异
        if satellite_flag == 1 and (AD_ref >= 2 * AD_alt):
            temp_bed_record = f"{chromosome}\t{start_pos}\t{end_pos}\n"
            bed_records.append(temp_bed_record)
            continue
        # 否则，将其加入最后的保留信息中
        REF.append(ref_v)
        ALT.append(alt_v)
        POS.append(pos_v)
    # *************************************************************************************

    # INFO区域
    INFO = "AT="
    ref_info = ""
    alt_info = ""
    for node in snarl_ref_path:
        ref_info = ref_info + f">{node}"
    for node, node_dir in zip(snarl_alt_path, snarl_alt_path_dir):
        alt_info = alt_info + node_dir + f"{node}"
    INFO = INFO + ref_info + "," + alt_info

    return POS, REF, ALT, INFO, bed_records


# 正常的基因分型过程，宏观，调用genotyping_all_normal
def normal_genotyping(snarl_paths, snarl_paths_dirs, snarl_refs, paths_support_reads, snarl_all_reads):
    """前五个参数分别是（此次进行分型的traversals行）所对应的-所有路径、路径方向、是否参考、支持的reads、snarl处的所有reads；
    f_vcf是vcf写入文件的句柄；
    idx_traversal是这个traversal在traversals_df中的索引，便于找到相邻的traversals，进而判断此处的变异是否是微卫星序列"""
    vcf_records = []
    bed_records = []
    # 对snarl_paths中的路径进行两两组合Gxy，对snarl处的所有reads，判断每一个read对Gxy的支持度L(r|Gxy)，并计算最终似然L(R|Gxy)
    # max_Lxy记录最大的似然值，max_x和max_y分别记录最大似然所对应的两个allele在snarl_paths中的位置
    max_Lxy = -math.inf
    max_x = 0
    max_y = 0
    # AD_ref AD_x AD_y分别记录snarl_ref_path、allele1、allele2的depth
    AD_ref = 0
    AD_x = 0
    AD_y = 0
    for idx_x, allele_x in enumerate(snarl_paths):
        for idx_y, allele_y in enumerate(snarl_paths):
            if idx_y < idx_x:
                continue
            else:
                # 如果两个allele的reads支持度都是0，则跳过
                if len(paths_support_reads[idx_x]) == 0 and len(paths_support_reads[idx_y]) == 0:
                    continue
                Lxy = likelihood_calculate(paths_support_reads[idx_x], paths_support_reads[idx_y], snarl_all_reads)
                if Lxy > max_Lxy:
                    max_Lxy = Lxy
                    # max_x记录支持度最高的allele，max_y记录支持度次高的allele
                    if len(paths_support_reads[idx_x]) > len(paths_support_reads[idx_y]):
                        max_x = idx_x
                        max_y = idx_y
                        AD_x = len(paths_support_reads[idx_x])
                        AD_y = len(paths_support_reads[idx_y])
                    else:
                        max_x = idx_y
                        max_y = idx_x
                        AD_x = len(paths_support_reads[idx_y])
                        AD_y = len(paths_support_reads[idx_x])

    # print(max_x,max_y)
    # print(snarl_paths[max_x])
    # print(snarl_paths[max_y])
    # 选取似然最高的基因型Gxy当作此snarl处的genotype.
    # 注：在snarl处基因型全是ref即0/0的，不写入文件中
    if snarl_refs[max_x] == '0' and snarl_refs[max_y] == '0':
        return [vcf_records, bed_records]

    # 基因型似然最高的两个allele,及其每个节点的方向
    allele_1 = snarl_paths[max_x]
    allele_2 = snarl_paths[max_y]
    allele_1_dir = snarl_paths_dirs[max_x]
    allele_2_dir = snarl_paths_dirs[max_y]

    # 找snarl参考路径的索引，据此找到snarl参考路径
    snarl_ref_idx = snarl_refs.index('0')
    snarl_ref_path = snarl_paths[snarl_ref_idx]
    AD_ref = len(paths_support_reads[snarl_ref_idx])

    # -------------------------------------------------------------------------------------------------------
    flag = 0  # 记录是不是1/2的情况
    # allele_1是snarl参考路径: 0/1
    if max_x == snarl_ref_idx:
        POS, REF, ALT, INFO, temp_bed_records = genotyping_all_normal(allele_1, allele_2, allele_2_dir, AD_ref, AD_y)
        GT = "0/1"
        DP = AD_x + AD_y
        AD = f"{AD_ref}" + "," + f"{AD_y}"
    # allele_2是snarl参考路径: 1/0
    elif max_y == snarl_ref_idx:
        POS, REF, ALT, INFO, temp_bed_records = genotyping_all_normal(allele_2, allele_1, allele_1_dir, AD_ref, AD_x)
        GT = "1/0"
        DP = AD_x + AD_y
        AD = f"{AD_ref}" + "," + f"{AD_x}"
    # allele_1和allele_2一样，但都不是snarl参考路径: 1/1
    elif max_x == max_y:
        POS, REF, ALT, INFO, temp_bed_records = genotyping_all_normal(snarl_ref_path, allele_1, allele_1_dir, AD_ref,
                                                                      AD_x)
        GT = "1/1"
        DP = AD_ref + AD_x
        AD = f"{AD_ref}" + "," + f"{AD_y}"
    # allele_1和allele_2不一样，但都不是snarl参考路径: 1/2.
    # 这时可以选择支持度最高的allele1与ref进行比较分型，也可以两个都和ref进行比较（两个都比较可能好一些，因为真值集合中有一些1|2类型的变异，拆开成两个变异）
    else:
        flag = 1
        POS1, REF1, ALT1, INFO1, temp_bed_records = genotyping_all_normal(snarl_ref_path, allele_1, allele_1_dir,
                                                                          AD_ref, AD_x)
        GT1 = "1/0"
        DP1 = AD_ref + AD_x
        AD1 = f"{AD_ref}" + "," + f"{AD_x}"

        POS2, REF2, ALT2, INFO2, temp_bed_records = genotyping_all_normal(snarl_ref_path, allele_2, allele_2_dir,
                                                                          AD_ref, AD_y)
        GT2 = "1/0"
        DP2 = AD_ref + AD_y
        AD2 = f"{AD_ref}" + "," + f"{AD_y}"

    # ===============================================================================================================================================
    for temp_bed_record in temp_bed_records:
        bed_records.append(temp_bed_record)
    temp_id = f">{snarl_ref_path[0]}>{snarl_ref_path[-1]}"
    if flag == 0:
        if len(POS) == 0:
            return [vcf_records, bed_records]
        else:
            for i in range(len(POS)):
                # 建立一个vcf行记录
                vcf_record = VcfRecord()
                vcf_record.CHROM = chromosome
                vcf_record.ID = temp_id
                vcf_record.POS = POS[i]
                vcf_record.REF = REF[i]
                vcf_record.ALT = ALT[i]
                vcf_record.INFO = INFO
                vcf_record.SAMPLE = GT + ":" + f"{DP}" + ":" + AD + ":" + f"{max_Lxy}"
                vcf_records.append(vcf_record.to_string())
    else:
        if len(POS1) == 0 and len(POS2) == 0:
            return [vcf_records, bed_records]
        else:
            for i in range(len(POS1)):
                # 建立一个vcf行记录
                vcf_record = VcfRecord()
                vcf_record.CHROM = chromosome
                vcf_record.ID = temp_id
                vcf_record.POS = POS1[i]
                vcf_record.REF = REF1[i]
                vcf_record.ALT = ALT1[i]
                vcf_record.INFO = INFO1
                vcf_record.SAMPLE = GT1 + ":" + f"{DP1}" + ":" + AD1 + ":" + f"{max_Lxy}"
                vcf_records.append(vcf_record.to_string())
            for i in range(len(POS2)):
                # 建立一个vcf行记录
                vcf_record = VcfRecord()
                vcf_record.CHROM = chromosome
                vcf_record.ID = temp_id
                vcf_record.POS = POS2[i]
                vcf_record.REF = REF2[i]
                vcf_record.ALT = ALT2[i]
                vcf_record.INFO = INFO2
                vcf_record.SAMPLE = GT2 + ":" + f"{DP2}" + ":" + AD2 + ":" + f"{max_Lxy}"
                vcf_records.append(vcf_record.to_string())
    return [vcf_records, bed_records]


# -----------------------------------------------------------------------------------------------------------------------
# 进行具体的变异检测过程，返回相关信息
def genotyping_all_complex(snarl_ref_path, snarl_alt_path, snarl_alt_path_dir):
    # 用来记录需要返回的VCF记录
    vcf_records = []
    bed_records = []
    # 将参考路径和替代路径进行拼接，用于ref nodes的计数
    total_nodes = snarl_ref_path + snarl_alt_path
    # 记录同时出现在snarl_ref_path和snarl_alt_path中的ref node的索引
    common_nodes = []
    common_nodes_index = []
    for idx, node in enumerate(snarl_ref_path):
        node_count = total_nodes.count(node)
        if node_count == 2:
            common_nodes.append(node)
            common_nodes_index.append(idx)

    # 记录节点名顺序相反的节点（有的区域节点顺序相反，但是序列可能差不多），放入set中，涉及这些节点的snarl跳过不考虑
    # 记录common节点在ref和alt中的对应顺序
    ref_common_nodes = common_nodes
    alt_common_nodes = []
    common_nodes_set = set(common_nodes)
    for node in snarl_alt_path:
        if node in common_nodes_set:
            alt_common_nodes.append(node)
    # 将逆序节点加入集合
    inv_nodes = set()
    for i in range(len(ref_common_nodes) - 1):
        for j in range(i + 1, len(ref_common_nodes)):
            if ref_common_nodes[i:j + 1] == alt_common_nodes[j:i - 1:-1]:
                for inv_node in ref_common_nodes[i:j + 1]:
                    inv_nodes.add(inv_node)

    # 从第二个共有node开始，找snarl_ref_path和snarl_alt_path中，与前一个共有node不都相邻的位置
    for idx, node in enumerate(common_nodes[1:]):
        idx = idx + 1
        alt_idx1 = snarl_alt_path.index(common_nodes[idx - 1])
        alt_idx2 = snarl_alt_path.index(node)
        ref_gap = common_nodes_index[idx] - common_nodes_index[idx - 1]
        alt_gap = alt_idx2 - alt_idx1
        # 如果ref和alt路径中，两个共有节点都是彼此邻接的，说明此处没有变异
        if ref_gap == 1 and alt_gap == 1:
            continue
        # 否则，有变异
        else:
            small_snarl_ref_path = snarl_ref_path[common_nodes_index[idx - 1]: common_nodes_index[idx] + 1]
            small_snarl_alt_path = snarl_alt_path[alt_idx1:alt_idx2 + 1]
            small_snarl_alt_path_dir = snarl_alt_path_dir[alt_idx1:alt_idx2 + 1]

            # 如果小snarl中涉及到inv_nodes，则跳过
            if (small_snarl_ref_path[0] in inv_nodes) or (small_snarl_ref_path[-1] in inv_nodes):
                continue

            # 利用df_edge_pack信息，计算相应变异的覆盖度，推断基因型信息
            # ref路径的覆盖度
            small_snarl_ref_pack = []
            for temp_idx in range(len(small_snarl_ref_path) - 1):
                temp_node1 = small_snarl_ref_path[temp_idx]
                temp_node2 = small_snarl_ref_path[temp_idx + 1]
                mask = ((df_edge_pack['from.id'] == int(temp_node1)) & (df_edge_pack['to.id'] == int(temp_node2))) | (
                        (df_edge_pack['from.id'] == int(temp_node2)) & (df_edge_pack['to.id'] == int(temp_node1)))
                temp_pack = df_edge_pack[mask].iloc[0, 4]
                small_snarl_ref_pack.append(temp_pack)
            # alt路径的覆盖度
            small_snarl_alt_pack = []
            for temp_idx in range(len(small_snarl_alt_path) - 1):
                temp_node1 = small_snarl_alt_path[temp_idx]
                temp_node2 = small_snarl_alt_path[temp_idx + 1]
                mask = ((df_edge_pack['from.id'] == int(temp_node1)) & (df_edge_pack['to.id'] == int(temp_node2))) | (
                        (df_edge_pack['from.id'] == int(temp_node2)) & (df_edge_pack['to.id'] == int(temp_node1)))
                temp_pack = df_edge_pack[mask].iloc[0, 4]
                small_snarl_alt_pack.append(temp_pack)
            # 用路径上覆盖度最低的那个值当作路径的覆盖度（短板效应）
            if len(small_snarl_ref_pack) > 0:
                ref_min_pack = min(small_snarl_ref_pack)
            else:
                ref_min_pack = 0
            if len(small_snarl_alt_pack) > 0:
                alt_min_pack = min(small_snarl_alt_pack)
            else:
                alt_min_pack = 0

            # 只有alt的depth大于1，才认为可能有变异
            if alt_min_pack <= 1:
                continue
            # 如果alt的depth太小，说明此处变异不存在
            if ref_min_pack > 5 * alt_min_pack:
                continue
            if alt_min_pack > 5 * ref_min_pack:
                temp_GT = "1/1"
            else:
                temp_GT = "0/1"
            temp_DP = ref_min_pack + alt_min_pack
            temp_AD = f"{ref_min_pack}" + "," + f"{alt_min_pack}"

            # 调用正常分型的具体操作
            POS, REF, ALT, INFO, temp_bed_records = genotyping_all_normal(small_snarl_ref_path, small_snarl_alt_path,
                                                                          small_snarl_alt_path_dir, ref_min_pack,
                                                                          alt_min_pack)
            for temp_bed_record in temp_bed_records:
                bed_records.append(temp_bed_record)
            if len(POS) == 0:
                continue
            else:
                for i in range(len(POS)):
                    # 建立一个vcf行记录
                    vcf_record = VcfRecord()
                    vcf_record.CHROM = chromosome
                    temp_id = f">{small_snarl_ref_path[0]}>{small_snarl_ref_path[-1]}"
                    vcf_record.ID = temp_id
                    vcf_record.POS = POS[i]
                    vcf_record.REF = REF[i]
                    vcf_record.ALT = ALT[i]
                    # 在同一个小snarl中出现的变异，INFO、GT、DP、AD都一样；似然不好算，用DP代替
                    vcf_record.INFO = INFO
                    vcf_record.SAMPLE = temp_GT + ":" + f"{temp_DP}" + ":" + temp_AD + ":" + f"{temp_DP}"
                    vcf_records.append(vcf_record.to_string())
    return [vcf_records, bed_records]


# 复杂区域的基因分型过程
def complex_region_genotyping(snarl_paths, snarl_paths_dirs, snarl_refs, paths_support_reads):
    vcf_records = []
    bed_records = []
    # 找snarl参考路径的索引，据此找到snarl参考路径
    snarl_ref_idx = snarl_refs.index('0')
    snarl_ref_path = snarl_paths[snarl_ref_idx]

    # 将每条alt路径都和ref进行比较
    for idx, allele in enumerate(snarl_paths):
        if idx == snarl_ref_idx:
            continue
        allele_dir = snarl_paths_dirs[idx]
        temp_vcf_records, temp_bed_records = genotyping_all_complex(snarl_ref_path, allele, allele_dir)
        for temp_record in temp_vcf_records:
            vcf_records.append(temp_record)
        for temp_bed_record in temp_bed_records:
            bed_records.append(temp_bed_record)
    return [vcf_records, bed_records]


# -----------------------------------------------------------------------------------------------------------------------
# 对单个snarl进行分型，找到单个snarl的参考和变异的序列、覆盖度、基因型，将结果用于多snarl联合分型中
def single_snarl_genotyping(snarl_paths, snarl_paths_dirs, snarl_refs, paths_support_reads, snarl_all_reads):
    # 对snarl_paths中的路径进行两两组合Gxy，对snarl处的所有reads，判断每一个read对Gxy的支持度L(r|Gxy)，并计算最终似然L(R|Gxy)
    # max_Lxy记录最大的似然值，max_x和max_y分别记录最大似然所对应的两个allele在snarl_paths中的位置
    max_Lxy = -math.inf
    max_x = 0
    max_y = 0
    # AD_ref AD_x AD_y分别记录snarl_ref_path、allele1、allele2的depth
    AD_ref = 0
    AD_x = 0
    AD_y = 0
    for idx_x, allele_x in enumerate(snarl_paths):
        for idx_y, allele_y in enumerate(snarl_paths):
            if idx_y < idx_x:
                continue
            else:
                # 如果两个allele的reads支持度都是0，则跳过
                if len(paths_support_reads[idx_x]) == 0 and len(paths_support_reads[idx_y]) == 0:
                    continue
                Lxy = likelihood_calculate(paths_support_reads[idx_x], paths_support_reads[idx_y], snarl_all_reads)
                if Lxy > max_Lxy:
                    max_Lxy = Lxy
                    # max_x记录支持度最高的allele，max_y记录支持度次高的allele
                    if len(paths_support_reads[idx_x]) > len(paths_support_reads[idx_y]):
                        max_x = idx_x
                        max_y = idx_y
                        AD_x = len(paths_support_reads[idx_x])
                        AD_y = len(paths_support_reads[idx_y])
                    else:
                        max_x = idx_y
                        max_y = idx_x
                        AD_x = len(paths_support_reads[idx_y])
                        AD_y = len(paths_support_reads[idx_x])

    # 基因型似然最高的两个allele,及其每个节点的方向
    allele_1 = snarl_paths[max_x]
    allele_2 = snarl_paths[max_y]
    allele_1_dir = snarl_paths_dirs[max_x]
    allele_2_dir = snarl_paths_dirs[max_y]

    # 找snarl参考路径的索引，据此找到snarl参考路径
    snarl_ref_idx = snarl_refs.index('0')
    snarl_ref_path = snarl_paths[snarl_ref_idx]
    snarl_ref_path_dir = snarl_paths_dirs[snarl_ref_idx]
    AD_ref = len(paths_support_reads[snarl_ref_idx])

    # 获取参考序列长度，每个snarl最后一个ref_node没有算
    ref_seq = ""
    for node in snarl_ref_path[:-1]:
        ref_seq += node_df[node_df['node'] == int(node)].iloc[0, 2]
    ref_seq_len = len(ref_seq)

    # 根据不同情况，返回ref和alt对应的信息
    # 0/0的情况
    if snarl_refs[max_x] == '0' and snarl_refs[max_y] == '0':
        return [snarl_ref_path], [snarl_ref_path_dir], [AD_ref], [snarl_ref_idx], max_Lxy, "0/0", ref_seq_len
    # allele_1是snarl参考路径: 0/1
    elif max_x == snarl_ref_idx:
        return [snarl_ref_path, snarl_ref_path, allele_2], [snarl_ref_path_dir, snarl_ref_path_dir, allele_2_dir], \
               [AD_ref, AD_ref, AD_y], [snarl_ref_idx, snarl_ref_idx, max_y], max_Lxy, "0/1", ref_seq_len
    # allele_2是snarl参考路径: 1/0
    elif max_y == snarl_ref_idx:
        return [snarl_ref_path, snarl_ref_path, allele_1], [snarl_ref_path_dir, snarl_ref_path_dir, allele_1_dir], \
               [AD_ref, AD_ref, AD_x], [snarl_ref_idx, snarl_ref_idx, max_x], max_Lxy, "1/0", ref_seq_len
    # allele_1和allele_2一样，但都不是snarl参考路径: 1/1
    elif max_x == max_y:
        return [snarl_ref_path, allele_1], [snarl_ref_path_dir, allele_1_dir], [AD_ref, AD_x], \
               [snarl_ref_idx, max_x], max_Lxy, "1/1", ref_seq_len
    # allele_1和allele_2不一样，但都不是snarl参考路径: 1/2.
    else:
        return [snarl_ref_path, allele_1, allele_2], [snarl_ref_path_dir, allele_1_dir, allele_2_dir], \
               [AD_ref, AD_x, AD_y], [snarl_ref_idx, max_x, max_y], max_Lxy, "1/2", ref_seq_len


# 多个snarl联合的基因分型过程
def joint_genotyping(snarl_paths_list, snarl_paths_dirs_list, snarl_refs_list, paths_support_reads_list,
                     snarl_all_reads_list):
    vcf_records = []
    bed_records = []
    # 记录ref和alt对应的路径
    ref_path, alt1_path, alt2_path = [], [], []
    # 记录ref和alt对应路径的方向
    ref_path_dir, alt1_path_dir, alt2_path_dir = [], [], []
    # 对应的depth
    ref_depth, alt1_depth, alt2_depth = [], [], []
    # 用于区分alt1和alt2路径的snarl索引
    diff_index = -1
    # 其他信息
    Variant_indexes = []
    Likelihoods = []
    Genotypes = []
    temp_info = ""
    # 记录每个snarl的参考seq序列长度，便于之后从变异反推所属的snarl
    seq_len_list = []

    # 保存每个snarl的信息
    for i in range(len(snarl_refs_list)):
        # 1.找到每个snarl处的ref序列和alt序列，选择的路径在snarl_paths_list[i]中的索引，及对应的覆盖度、似然、基因型信息
        Paths, Paths_dirs, Depths, vr_idx, Lxy, gt, ref_seq_len = \
            single_snarl_genotyping(snarl_paths_list[i], snarl_paths_dirs_list[i], snarl_refs_list[i],
                                    paths_support_reads_list[i], snarl_all_reads_list[i])

        # 此处没有变异：[snarl_ref_path], [snarl_ref_path_dir], [AD_ref], [snarl_ref_idx], max_Lxy, "0/0", ref_seq_len
        if gt == "0/0":
            ref_path += Paths[0]
            alt1_path += Paths[0]
            alt2_path += Paths[0]
            ref_path_dir += Paths_dirs[0]
            alt1_path_dir += Paths_dirs[0]
            alt2_path_dir += Paths_dirs[0]
            ref_depth.append(Depths[0])
            alt1_depth.append(Depths[0])
            alt2_depth.append(Depths[0])

        # 此处为纯和变异：[snarl_ref_path, allele_1], [snarl_ref_path_dir, allele_1_dir], [AD_ref, AD_x], [snarl_ref_idx, max_x], max_Lxy, "1/1", ref_seq_len
        elif gt == "1/1":
            ref_path += Paths[0]
            alt1_path += Paths[1]
            alt2_path += Paths[1]
            ref_path_dir += Paths_dirs[0]
            alt1_path_dir += Paths_dirs[1]
            alt2_path_dir += Paths_dirs[1]
            ref_depth.append(Depths[0])
            alt1_depth.append(Depths[1])
            alt2_depth.append(Depths[1])

        # 此处为杂合变异：[snarl_ref_path, snarl_ref_path, allele_1], [snarl_ref_path_dir, snarl_ref_path_dir, allele_1_dir], \
        #                [AD_ref, AD_ref, AD_x], [snarl_ref_idx, snarl_ref_idx, max_x], max_Lxy, "1/0", ref_seq_len
        elif gt == "0/1" or gt == "1/0":
            # alt1和alt2在之前都一样
            if diff_index == -1:
                ref_path += Paths[0]
                alt1_path += Paths[1]
                alt2_path += Paths[2]
                ref_path_dir += Paths_dirs[0]
                alt1_path_dir += Paths_dirs[1]
                alt2_path_dir += Paths_dirs[2]
                ref_depth.append(Depths[0])
                alt1_depth.append(Depths[1])
                alt2_depth.append(Depths[2])
            # alt1和alt2在之前不一样
            else:
                # 能区分alt1和alt2的snarl处对应的路径索引（可能是0/1，也可能是1/2）
                temp_vr_idx = Variant_indexes[diff_index]
                read_find_flag = 0
                for read in paths_support_reads_list[i][vr_idx[1]]:
                    if read in paths_support_reads_list[diff_index][temp_vr_idx[1]]:
                        # 如果前一个alt1和这里的alt1对应，则拼接即可
                        ref_path += Paths[0]
                        alt1_path += Paths[1]
                        alt2_path += Paths[2]
                        ref_path_dir += Paths_dirs[0]
                        alt1_path_dir += Paths_dirs[1]
                        alt2_path_dir += Paths_dirs[2]
                        ref_depth.append(Depths[0])
                        alt1_depth.append(Depths[1])
                        alt2_depth.append(Depths[2])
                        read_find_flag = 1
                        break
                    elif read in paths_support_reads_list[diff_index][temp_vr_idx[2]]:
                        # 如果前一个alt1和这里的alt2对应，则拼接之后需要将alt1_path和alt2_path即对应的depth进行调换，以保证alt1总是靠前的那一条
                        ref_path += Paths[0]
                        alt1_path += Paths[2]
                        alt2_path += Paths[1]
                        temp_path = alt1_path
                        alt1_path = alt2_path
                        alt2_path = temp_path
                        ref_path_dir += Paths_dirs[0]
                        alt1_path_dir += Paths_dirs[2]
                        alt2_path_dir += Paths_dirs[1]
                        temp_path_dir = alt1_path_dir
                        alt1_path_dir = alt2_path_dir
                        alt2_path_dir = temp_path_dir
                        ref_depth.append(Depths[0])
                        alt1_depth.append(Depths[2])
                        alt2_depth.append(Depths[1])
                        temp_depth = alt1_depth
                        alt1_depth = alt2_depth
                        alt2_depth = temp_depth
                        read_find_flag = 1
                        break
                if read_find_flag == 0:
                    # 这里是为了增加完备性，防止有些reads有问题
                    ref_path += Paths[0]
                    alt1_path += Paths[1]
                    alt2_path += Paths[2]
                    ref_path_dir += Paths_dirs[0]
                    alt1_path_dir += Paths_dirs[1]
                    alt2_path_dir += Paths_dirs[2]
                    ref_depth.append(Depths[0])
                    alt1_depth.append(Depths[1])
                    alt2_depth.append(Depths[2])
            # 注意这一行放的位置，别移动，对上面的if和else都有效
            diff_index = i

        # 此处为1/2变异：[snarl_ref_path, allele_1, allele_2], [snarl_ref_path_dir, allele_1_dir, allele_2_dir], \
        #                [AD_ref, AD_x, AD_y], [snarl_ref_idx, max_x, max_y], max_Lxy, "1/2", ref_seq_len
        else:
            # alt1和alt2在之前都一样
            if diff_index == -1:
                ref_path += Paths[0]
                alt1_path += Paths[1]
                alt2_path += Paths[2]
                ref_path_dir += Paths_dirs[0]
                alt1_path_dir += Paths_dirs[1]
                alt2_path_dir += Paths_dirs[2]
                ref_depth.append(Depths[0])
                alt1_depth.append(Depths[1])
                alt2_depth.append(Depths[2])
            # alt1和alt2在之前不一样
            else:
                # 能区分alt1和alt2的snarl处对应的路径索引（可能是0/1，也可能是1/2）
                temp_vr_idx = Variant_indexes[diff_index]
                read_find_flag = 0
                for read in paths_support_reads_list[i][vr_idx[1]]:
                    if read in paths_support_reads_list[diff_index][temp_vr_idx[1]]:
                        # 如果前一个alt1和这里的alt1对应，则拼接即可
                        ref_path += Paths[0]
                        alt1_path += Paths[1]
                        alt2_path += Paths[2]
                        ref_path_dir += Paths_dirs[0]
                        alt1_path_dir += Paths_dirs[1]
                        alt2_path_dir += Paths_dirs[2]
                        ref_depth.append(Depths[0])
                        alt1_depth.append(Depths[1])
                        alt2_depth.append(Depths[2])
                        read_find_flag = 1
                        break
                    elif read in paths_support_reads_list[diff_index][temp_vr_idx[2]]:
                        # 如果前一个alt1和这里的alt2对应，则拼接之后需要将alt1_path和alt2_path即对应的depth进行调换，以保证alt1总是靠前的那一条
                        ref_path += Paths[0]
                        alt1_path += Paths[2]
                        alt2_path += Paths[1]
                        temp_path = alt1_path
                        alt1_path = alt2_path
                        alt2_path = temp_path
                        ref_path_dir += Paths_dirs[0]
                        alt1_path_dir += Paths_dirs[2]
                        alt2_path_dir += Paths_dirs[1]
                        temp_path_dir = alt1_path_dir
                        alt1_path_dir = alt2_path_dir
                        alt2_path_dir = temp_path_dir
                        ref_depth.append(Depths[0])
                        alt1_depth.append(Depths[2])
                        alt2_depth.append(Depths[1])
                        temp_depth = alt1_depth
                        alt1_depth = alt2_depth
                        alt2_depth = temp_depth
                        read_find_flag = 1
                        break
                if read_find_flag == 0:
                    # 这里是为了增加完备性，防止有些reads有问题
                    ref_path += Paths[0]
                    alt1_path += Paths[1]
                    alt2_path += Paths[2]
                    ref_path_dir += Paths_dirs[0]
                    alt1_path_dir += Paths_dirs[1]
                    alt2_path_dir += Paths_dirs[2]
                    ref_depth.append(Depths[0])
                    alt1_depth.append(Depths[1])
                    alt2_depth.append(Depths[2])
            # 注意这一行放的位置，别移动，对上面的if和else都有效
            diff_index = i

        # 记录info信息
        temp_info += f">{snarl_paths_list[i][0][0]}>{snarl_paths_list[i][0][-1]}|"
        # 记录max_Lxy似然信息
        Likelihoods.append(Lxy)
        # 记录基因型信息
        Genotypes.append(gt)
        # 记录每个snarl处选择路径的索引
        Variant_indexes.append(vr_idx)
        # 记录snarl的ref序列长度
        seq_len_list.append(ref_seq_len)

    # 前面已经找到了ref_path，alt1_path，alt2_path；现在需要判断alt1_path、alt2_path、ref_path是否相等
    if alt1_path == ref_path and alt2_path == ref_path:
        # 三者都一样，没有变异
        return [vcf_records, bed_records]
    elif alt1_path == ref_path:
        # alt1_path == ref_path， 只需要比较alt2_path和ref_path即可
        # ref序列
        ref_seq = ""
        last_node = ref_path[-1]
        for node in ref_path:
            if node == last_node:
                continue
            ref_seq += node_df[node_df['node'] == int(node)].iloc[0, 2]
            last_node = node
        # alt2序列
        alt_seq_2 = ""
        last_node = alt2_path[-1]
        for node, node_dir in zip(alt2_path, alt2_path_dir):
            if node == last_node:
                continue
            if node_dir == '>':
                alt_seq_2 += node_df[node_df['node'] == int(node)].iloc[0, 2]
            else:
                alt_seq_2 += dna_reverse_complement(node_df[node_df['node'] == int(node)].iloc[0, 2])
            last_node = node
        # 进行比对，检测变异并分型
        temp_vcf_records, temp_bed_records = aln_and_call(ref_seq, alt_seq_2, ref_depth, alt2_depth, snarl_paths_list,
                                                          temp_info, Genotypes, Likelihoods,
                                                          seq_len_list)
        for temp_record in temp_vcf_records:
            vcf_records.append(temp_record)
        for temp_bed_record in temp_bed_records:
            bed_records.append(temp_bed_record)

    elif alt2_path == ref_path:
        # alt2_path == ref_path， 只需要比较alt1_path和ref_path即可
        # ref序列
        ref_seq = ""
        last_node = ref_path[-1]
        for node in ref_path:
            if node == last_node:
                continue
            ref_seq += node_df[node_df['node'] == int(node)].iloc[0, 2]
            last_node = node
        # alt1序列
        alt_seq_1 = ""
        last_node = alt1_path[-1]
        for node, node_dir in zip(alt1_path, alt1_path_dir):
            if node == last_node:
                continue
            if node_dir == '>':
                alt_seq_1 += node_df[node_df['node'] == int(node)].iloc[0, 2]
            else:
                alt_seq_1 += dna_reverse_complement(node_df[node_df['node'] == int(node)].iloc[0, 2])
            last_node = node
        # 进行比对，检测变异并分型
        temp_vcf_records, temp_bed_records = aln_and_call(ref_seq, alt_seq_1, ref_depth, alt1_depth, snarl_paths_list,
                                                          temp_info, Genotypes, Likelihoods,
                                                          seq_len_list)
        for temp_record in temp_vcf_records:
            vcf_records.append(temp_record)
        for temp_bed_record in temp_bed_records:
            bed_records.append(temp_bed_record)
    else:
        # alt1_path、alt2_path和ref_path都不一样，都需要进行比较
        # ref序列
        ref_seq = ""
        last_node = ref_path[-1]
        for node in ref_path:
            if node == last_node:
                continue
            ref_seq += node_df[node_df['node'] == int(node)].iloc[0, 2]
            last_node = node
        # alt1序列
        alt_seq_1 = ""
        last_node = alt1_path[-1]
        for node, node_dir in zip(alt1_path, alt1_path_dir):
            if node == last_node:
                continue
            if node_dir == '>':
                alt_seq_1 += node_df[node_df['node'] == int(node)].iloc[0, 2]
            else:
                alt_seq_1 += dna_reverse_complement(node_df[node_df['node'] == int(node)].iloc[0, 2])
            last_node = node
        # alt2序列
        alt_seq_2 = ""
        last_node = alt2_path[-1]
        for node, node_dir in zip(alt2_path, alt2_path_dir):
            if node == last_node:
                continue
            if node_dir == '>':
                alt_seq_2 += node_df[node_df['node'] == int(node)].iloc[0, 2]
            else:
                alt_seq_2 += dna_reverse_complement(node_df[node_df['node'] == int(node)].iloc[0, 2])
            last_node = node

        # 进行比对，检测变异并分型
        temp_vcf_records, temp_bed_records = aln_and_call(ref_seq, alt_seq_1, ref_depth, alt1_depth, snarl_paths_list,
                                                          temp_info, Genotypes, Likelihoods,
                                                          seq_len_list)
        for temp_record in temp_vcf_records:
            vcf_records.append(temp_record)
        for temp_bed_record in temp_bed_records:
            bed_records.append(temp_bed_record)
        temp_vcf_records, temp_bed_records = aln_and_call(ref_seq, alt_seq_2, ref_depth, alt2_depth, snarl_paths_list,
                                                          temp_info, Genotypes, Likelihoods,
                                                          seq_len_list)
        for temp_record in temp_vcf_records:
            vcf_records.append(temp_record)
        for temp_bed_record in temp_bed_records:
            bed_records.append(temp_bed_record)

    return [vcf_records, bed_records]


def aln_and_call(ref_seq, alt_seq, ref_depth, alt_depth, snarl_paths_list, temp_info, Genotypes, Likelihoods,
                 seq_len_list):
    # 将对应的记录返回
    vcf_records = []
    bed_records = []
    # 3.进行比对分型
    # 进行比对检测变异和分型,得到变异的参考序列和变异序列
    ref_variants, alt_variants, temp_pos = variant_call(ref_seq, alt_seq)
    # 上面的temp_pos为snarl内部的相对位置，需要加上snarl_start_pos
    snarl_start_node = snarl_paths_list[0][0][0]
    snarl_start_pos = ref_nodes_pos_df[ref_nodes_pos_df['node'] == snarl_start_node].iloc[0, 1]
    pos_variants = [x + snarl_start_pos - 1 for x in temp_pos]

    # *****************************************************************************************
    # 对于ATTATATA--ATATA这种变异，比对的结果中变异在开头即：'ATT'--''，会出现空字符的情况
    # 这时需要在前面加上一个碱基；且这种情况只会在开头第一个变异出现
    if len(ref_variants) > 0:
        if ref_variants[0] == '' or alt_variants[0] == '':
            # 详见variant_call函数定义处，开头是这种indel时，pos需要加上1
            pos_variants[0] += 1
            # 定位变异在fa文件中的位置：行索引+offset
            line_num = pos_variants[0] // 60
            line_offset = pos_variants[0] % 60 - 1
            # 这里往前取一行，防止变异发生行首位置（此时应该加上上一行最后一个碱基）
            parent_seq_front = chr21_fa_lines[line_num - 1] + chr21_fa_lines[line_num]
            # 59 + line_offset 即 60 + line_offset -1，减去1代表往前取一个碱基
            add_one_seq = parent_seq_front[59 + line_offset]

            ref_variants[0] = add_one_seq + ref_variants[0]
            alt_variants[0] = add_one_seq + alt_variants[0]
            pos_variants[0] -= 1
    # *****************************************************************************************

    # 4.将其中的非卫星序列变异写入文件
    # INFO区域
    INFO = "AT="
    INFO += temp_info
    for ref_v, alt_v, pos_v in zip(ref_variants, alt_variants, pos_variants):
        vcf_record = VcfRecord()
        vcf_record.CHROM = chromosome
        vcf_record.ID = f">{snarl_paths_list[0][0][0]}>{snarl_paths_list[-1][0][-1]}"
        vcf_record.POS = pos_v
        vcf_record.REF = ref_v
        vcf_record.ALT = alt_v
        vcf_record.INFO = INFO
        # **************************对于GT这个重要指标，需要将变异再定位到snarl中进行判定************************
        # 重新计算变异在这个snarl序列群中的相对位置
        opposite_pos_v = pos_v - snarl_start_pos + 1
        accu_len = 0  # 累加序列长度
        # 找到第一个累加序列长度 >= 变异相对位置的索引i；则这个变异对应第i个snarl；再去snarl中找相应的信息
        snarl_idx = 0
        snarl_idx_flag = 0
        for i in range(len(seq_len_list)):
            accu_len += seq_len_list[i]
            if accu_len >= opposite_pos_v:
                snarl_idx = i
                snarl_idx_flag = 1
                break
        # 如果这个flag为0，说明循环中没有找到accu_len >= opposite_pos_v条件；此时变异发生在最后一个snarl的最后一个节点上，因为这个节点的长度没有加到seq_len_list中
        if snarl_idx_flag == 0:
            snarl_idx = len(seq_len_list) - 1
        GT = Genotypes[snarl_idx]
        # 防止变异位置移动（如左对齐右对齐）造成对齐到没有变异的snarl中
        if GT == '0/0':
            GT = '0/1'
        if GT == '1/2':
            GT = '0/1'
        max_Lxy = Likelihoods[snarl_idx]
        DP = ref_depth[snarl_idx] + alt_depth[snarl_idx]
        AD = f"{ref_depth[snarl_idx]},{alt_depth[snarl_idx]}"
        # **************************对其中的每个变异进行判断是否是微卫星序列**************************
        satellite_flag = 0
        start_pos = 0
        end_pos = 0
        # 只对indel进行卫星序列的判断,snp不用
        if len(ref_v) == 1 and len(alt_v) > 1:
            satellite_flag, start_pos, end_pos = whether_satellite_sequence(pos_v, alt_v)
        if len(ref_v) > 1 and len(alt_v) == 1:
            satellite_flag, start_pos, end_pos = whether_satellite_sequence(pos_v, ref_v)
        # 如果是卫星序列，则跳过这个变异；记录pos位置，用于生成bed文件，去除真值集中这些区域的变异
        if satellite_flag == 1 and (ref_depth[snarl_idx] >= 2 * alt_depth[snarl_idx]):
            temp_bed_record = f"{chromosome}\t{start_pos}\t{end_pos}\n"
            bed_records.append(temp_bed_record)
            continue
        # *************************************************************************************
        # 否则，将其写入文件中
        vcf_record.SAMPLE = GT + ":" + f"{DP}" + ":" + AD + ":" + f"{max_Lxy}"
        vcf_records.append(vcf_record.to_string())
    return [vcf_records, bed_records]


# ----------------------------------------------------并行任务，每个chunk的工作---------------------------------------------
def chunk_work(traversals_df):
    start_index=traversals_df.index[0]
    traversals_df.reset_index(drop=True, inplace=True)
    # 用来记录每个分型函数返回的vcf记录，便于return后调用回调函数
    vcf_records = []
    bed_records = []
    # 记录总共的traversals数量，便于判断结尾
    traversals_df_rows = traversals_df.shape[0]

    # -----------------------------------------------Genotyping-----------------------------------------------------------------
    # 当snarl发生联合时，记录联合的数量，也即下次traversals_df需要跳过的行数
    skip_nums = -1
    # traversals_df每一行都是一个潜在变异点
    for idx in traversals_df.index:
        # print(idx, "------------")
        logging.info(f"{start_index+idx}")
        # 跳过之前已经联合的snarl
        if skip_nums > 0:
            skip_nums -= 1
            continue

        # 定义包含多个snarl信息的列表
        snarl_paths_list = []
        snarl_paths_dirs_list = []
        snarl_refs_list = []
        paths_support_reads_list = []
        snarl_all_reads_list = []
        right_node_len = 0  # 开始设置为0以将第一个snarl放入列表中，防止空列表

        # 当snarl发生联合时，记录联合的数量，也即下次traversals_df需要跳过的行数
        skip_nums = -1

        # 记录是不是第一个snarl，如果是第一个snarl，则上一个snarl最后一个node和当前snarl第一个node就不需要比较
        first_snarl_flag = 1
        last_snarl_last_node = ""

        # 由于每次都是取到right_node的长度大于20bp，下一次的snarl范围中left_node（也即上一个snarl范围的right_node）序列长度也一定大于20bp
        while right_node_len < 50 and idx < traversals_df_rows:
            skip_nums += 1
            # 将数据转换为列表形式
            row1 = traversals_df.loc[idx]
            snarl_paths = eval('[' + row1['path'] + ']')
            snarl_paths_dirs = eval('[' + row1['dir'] + ']')
            snarl_refs = re.findall(r'\d+', row1['ref'])
            paths_support_reads = eval('[' + row1['support_reads'] + ']')
            # 将每个snarl_path处的reads进行合并、去重，得到比对到snarl处的所有reads. 注:main函数里边的变量，在其他函数中可以用
            snarl_all_reads = extract_snarl_all_reads(row1['support_reads'])

            # *******************如果是复杂snarl，则联合中断***************************
            # 找snarl参考路径的索引，据此找到snarl参考路径
            snarl_ref_idx = snarl_refs.index('0')
            snarl_ref_path = snarl_paths[snarl_ref_idx]
            # snarl_ref_path_len = len(snarl_ref_path)
            # 计算snarl中路径节点数量的最大值
            max_len = max(len(lst) for lst in snarl_paths)
            # 如果snarl涉及节点较多，则不进行联合分型
            if max_len > 40:
                # 如果列表中还没有内容，说明这个复杂snarl是此次联合的第一个snarl，需要加到列表中以进行之后的统一分型
                if len(snarl_paths_list) == 0:
                    snarl_paths_list.append(snarl_paths)
                    snarl_paths_dirs_list.append(snarl_paths_dirs)
                    snarl_refs_list.append(snarl_refs)
                    paths_support_reads_list.append(paths_support_reads)
                    snarl_all_reads_list.append(snarl_all_reads)
                    break
                # 如果列表中已经有内容，说明这个复杂snarl不是第一个，此时不用加入列表，skip_num需要回退1
                else:
                    skip_nums -= 1
                    break
            # *********************************************************************
            # 如果前一个snarl的最后一个node和当前snarl的第一个node不同，则不进行联合
            # 如果前一个snarl的最后一个node和当前snarl的第一个node相差超过1，则不进行联合
            now_snarl_first_node = snarl_ref_path[0]
            if first_snarl_flag == 0:  # 不是第一个snarl，不用加入列表，skip_num需要回退1
                # if now_snarl_first_node != last_snarl_last_node:
                if abs(int(now_snarl_first_node) - int(last_snarl_last_node)) > 1:
                    skip_nums -= 1
                    break
            last_snarl_last_node = snarl_ref_path[-1]
            # *********************************************************************

            # 将对应信息存入新列表中，便于扩展snarl
            snarl_paths_list.append(snarl_paths)
            snarl_paths_dirs_list.append(snarl_paths_dirs)
            snarl_refs_list.append(snarl_refs)
            paths_support_reads_list.append(paths_support_reads)
            snarl_all_reads_list.append(snarl_all_reads)

            # 进行right_node_len的更新
            right_node = snarl_paths[0][-1]
            right_node_len = len(node_df[node_df['node'] == int(right_node)].iloc[0, 2])

            # idx进行+1操作
            idx += 1

            # 更新是否是第一个snarl的标识
            first_snarl_flag = 0

        # ---------------------------------------------------------------------------------------------------------
        # 如果list中只包含1个snarl的信息，则进行正常的基因分型过程
        if len(snarl_refs_list) == 1:
            # snarl_ref_idx = snarl_refs_list[0].index('0')
            # snarl_ref_path = snarl_paths_list[0][snarl_ref_idx]
            # snarl_ref_path_len = len(snarl_ref_path)

            max_len = max(len(lst) for lst in snarl_paths_list[0])
            # 如果snarl涉及的节点数量较少，则正常分型
            if max_len <= 40:
                temp_vcf_records, temp_bed_records = normal_genotyping(snarl_paths_list[0], snarl_paths_dirs_list[0],
                                                                       snarl_refs_list[0],
                                                                       paths_support_reads_list[0],
                                                                       snarl_all_reads_list[0])
                for temp_record in temp_vcf_records:
                    vcf_records.append(temp_record)
                for temp_bed_record in temp_bed_records:
                    bed_records.append(temp_bed_record)
            # 如果涉及的节点数量较多，则进行complex_region_genotyping
            else:
                temp_vcf_records, temp_bed_records = complex_region_genotyping(snarl_paths_list[0],
                                                                               snarl_paths_dirs_list[0],
                                                                               snarl_refs_list[0],
                                                                               paths_support_reads_list[0])
                for temp_record in temp_vcf_records:
                    vcf_records.append(temp_record)
                for temp_bed_record in temp_bed_records:
                    bed_records.append(temp_bed_record)
        # 如果list中包含多个snarl的信息，则进行联合分型过程；此时联合分型中已经没有复杂snarl，只可能有一些1/2类型的变异需要单独处理一下
        elif len(snarl_refs_list) > 1:
            temp_vcf_records, temp_bed_records = joint_genotyping(snarl_paths_list, snarl_paths_dirs_list,
                                                                  snarl_refs_list,
                                                                  paths_support_reads_list,
                                                                  snarl_all_reads_list)
            for temp_record in temp_vcf_records:
                vcf_records.append(temp_record)
            for temp_bed_record in temp_bed_records:
                bed_records.append(temp_bed_record)
    return [vcf_records, bed_records]


def file_write_callback(vcf_bed_records):
    vcf_records = vcf_bed_records[0]
    bed_records = vcf_bed_records[1]

    f_vcf = open(vcf_file, 'a')
    for record in vcf_records:
        f_vcf.write(record)
    f_vcf.close()

    # 打开bed文件，将卫星序列的变异位置写入bed文件
    f_bed_all = open(bed_all_file, 'a')
    for record in bed_records:
        f_bed_all.write(record)
    # 关闭bed文件
    f_bed_all.close()

def error_callback(error):
    print(f"Error info: {error}")

# -----------------------------------------------文件处理-----------------------------------------------------------------
parser = argparse.ArgumentParser(description="ref_node_pos")
parser.add_argument('-c', '--chromosome', type=str, default="chr")
parser.add_argument('-l', '--chromosome_len', type=str, default="1")
args = parser.parse_args()
chromosome = args.chromosome
chromosome_len = args.chromosome_len

vcf_file = "./result/sample_all.txt"
log_file = "./result/log.txt"
bed_all_file = "./result/bed_all.txt"
gfa_file = "./job_store/" + chromosome + ".gfa"
ref_nodes_pos_file = "./job_store/ref_nodes_pos.csv"
traversals_with_reads_agg_file = "./job_store/traversals_with_reads_agg_sort_by_pos.csv"
ref_path_file = "./job_store/REF_GRCh38_" + chromosome + ".txt"
edge_pack_file = "./job_store/aln_edge_pack.txt"
chr_fa_file = "./GCA_" + chromosome + ".fa"

# 参考路径ref_path
ref_path_df = pd.read_csv(ref_path_file, sep='\t', header=None, names=['0', '1', '2', '3'])
ref_path_str = ref_path_df.iloc[0, 2]
ref_path = re.findall(r"\d+", ref_path_str)
# print("参考路径：\n", ref_path)

# 从GFA文件中提取node相关信息
node_df = gfa_transform(gfa_file)
# print("node序列：\n", node_df)

# 从ref_nodes_pos.csv中加载每个ref_node的POS位置
ref_nodes_pos_df = pd.read_csv(ref_nodes_pos_file, sep='\t', header=None, names=['node', 'pos'])
ref_nodes_pos_df['pos'] = ref_nodes_pos_df['pos'] + 1  # 准换为1-base的坐标
ref_nodes_pos_df['node'] = ref_nodes_pos_df['node'].astype(str)  # 将node列的数据类型转换为str

# 读取edge_pack文件信息，得到每条边的覆盖度；之后可在涉及节点较多的snarl处进行分解判断
df_edge_pack = pd.read_csv(edge_pack_file, sep='\t')

# 读取chr21的fasta文件，用于判断变异是否是卫星序列
with open(chr_fa_file, 'r') as f:
    chr21_fa_lines = f.readlines()[1:]
    chr21_fa_lines = [line.strip() for line in chr21_fa_lines]


# 记录处理的snarl数量，存入log文件中
logging.basicConfig(filename=log_file, level=logging.INFO)
print('q---------------------------')

if __name__ == '__main__':
    # 向vcf文件写入header信息
    f_vcf = open(vcf_file, 'w')
    header_str = vcf_header(chromosome, chromosome_len)
    f_vcf.write(header_str)
    f_vcf.close()
    
    print(time.localtime())
    start_time = time.time()
    # 将dataframe进行分块
    chunks = pd.read_csv(traversals_with_reads_agg_file, sep='\t', chunksize=5000)

    # 每个chunk一个进程
    pool = mp.Pool(80)
    for chunk in chunks:
        pool.apply_async(chunk_work, (chunk,), callback=file_write_callback, error_callback=error_callback)
    pool.close()
    pool.join()

    end_time = time.time()
    print(time.localtime())
    print(end_time - start_time)
