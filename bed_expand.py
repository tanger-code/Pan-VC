import pandas as pd

bed_file = "./result/bed_all.bed"
bed_expand_file = "./result/bed_all_expand.bed"

# 读取输入bed文件
input_bed = pd.read_csv(bed_file, sep='\t', header=None)

# 将第二列的整数全部减去10
input_bed[1] -= 20

# 将第三列的数值加上10
input_bed[2] += 20

# 将结果写入输出bed文件
input_bed.to_csv(bed_expand_file, sep='\t', header=None, index=None)
