import math

from Bio import pairwise2

# 控制得分；匹配得1分，不匹配得0分；其中A-T和G-C的匹配稍微得负分，因为它们彼此之间的转换比较困难
match_dict = {
    ('A', 'A'): 1,
    ('A', 'T'): 0,
    ('A', 'G'): 0,
    ('A', 'C'): 0,

    ('T', 'T'): 1,
    ('T', 'G'): 0,
    ('T', 'C'): 0,

    ('G', 'G'): 1,
    ('G', 'C'): 0,

    ('C', 'C'): 1,

    ('N', 'A'): 0,
    ('N', 'T'): 0,
    ('N', 'G'): 0,
    ('N', 'C'): 0,

    ('N', 'N'): 0
}


# 注：这里当变异发生在开头且为indel时(变异形式为'ATT'--'')，pos位置为0；开头是SNP时，pos为1
# 在被调用的时候，如果开头是indel，需要将对应的pos加上1
def variant_call(ref, alt):
    # 如果序列太长，pairwise2比对不了
    ref_len = len(ref)
    alt_len = len(alt)
    if ref_len > 10000 and alt_len > 10000:
        return [], [], []
    
    # 进行比对
    alignments = pairwise2.align.globalds(ref, alt, match_dict, -2, 0)

    left_variant1 = []
    left_variant2 = []
    left_pos = []
    pos_min = math.inf

    # 可以遍历每一个比对，然后将POS位置相加最小的比对进行输出
    for aln in alignments:
        # 初步格式化之后的序列
        ref_seq = aln.seqA
        alt_seq = aln.seqB

        # 从初步格式化的序列中，检测变异并保存到列表中，便于直接写入VCF文件
        variant1 = []
        variant2 = []
        pos = []  # 变异位点的位置
        gap_count = 0  # ref某位置之前的gap数

        # 记录上一个处理的字符
        last_a = ""
        last_b = ""
        for n, (a, b) in enumerate(zip(ref_seq, alt_seq)):
            # match
            if a == b:
                pass
            # ins
            elif a == '-':
                if last_a == '-':
                    variant2[-1] = variant2[-1] + b
                else:
                    variant1.append(last_a)
                    variant2.append(last_a + b)
                    pos.append(n - gap_count)
                gap_count += 1
            # del
            elif b == '-':
                if last_b == '-':
                    variant1[-1] = variant1[-1] + a
                else:
                    variant1.append(last_a + a)
                    variant2.append(last_a)
                    pos.append(n - gap_count)
            # SNP/MNP
            else:
                if last_a == last_b or last_a == '-' or last_b == '-':
                    variant1.append(a)
                    variant2.append(b)
                    pos.append(n + 1 - gap_count)
                else:
                    variant1[-1] = variant1[-1] + a
                    variant2[-1] = variant2[-1] + b
            last_a = a
            last_b = b

        pos_add = 0
        for temp_pos in pos:
            pos_add += temp_pos
        if pos_add < pos_min:
            pos_min = pos_add
            left_variant1 = variant1
            left_variant2 = variant2
            left_pos = pos

    return left_variant1, left_variant2, left_pos
