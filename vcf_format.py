# 将vcf文件的header信息存储到一个string中
def vcf_header(chromosome, chromosome_len):
    # 定义一个空字符串，用于存储VCF格式的数据
    header = ""
    # 添加VCF文件的头部信息，包括版本号，日期，字段说明等
    header += "##fileformat=VCFv4.2\n"
    header += "##contig=<ID=" + chromosome + ",length=" + chromosome_len + ">\n"
    header += "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">\n"
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
    header += "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
    header += "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood\">\n"
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    return header


# 定义一个类来存储vcf的一行，即一个记录record，并定义to_string方法将信息以string格式返回
class VcfRecord:
    def __init__(self, CHROM="", POS=0, ID="", REF="", ALT="", QUAL=".", FILTER=".", INFO="", FORMAT="GT:DP:AD:GL",
                 SAMPLE=""):
        self.CHROM = CHROM
        self.POS = POS
        self.ID = ID
        self.REF = REF
        self.ALT = ALT
        self.QUAL = QUAL
        self.FILTER = FILTER
        self.INFO = INFO
        self.FORMAT = FORMAT
        self.SAMPLE = SAMPLE

    def to_string(self):
        variant_line = f"{self.CHROM}\t{self.POS}\t{self.ID}\t{self.REF}\t{self.ALT}\t{self.QUAL}\t{self.FILTER}\t{self.INFO}\t{self.FORMAT}\t{self.SAMPLE}\n"
        return variant_line


# 如果将record数据存储到了列表中，则可以用此方法将列表数据转换为vcf的字符串（很大，包含全部信息），之后直接将字符串写入文件即可
def list_to_vcf(my_list):
    # 定义一个空字符串，用于存储VCF格式的数据
    vcf_string = ""
    # 添加VCF文件的头部信息，包括版本号，日期，字段说明等
    vcf_string += "##fileformat=VCFv4.3\n"
    vcf_string += "##fileDate=20210723\n"
    vcf_string += "##source=MyProgram\n"
    vcf_string += "##reference=GRCh38\n"
    vcf_string += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
    vcf_string += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    # 遍历列表中的每一行数据，将其转换为VCF格式的一行
    for row in my_list:
        # 获取染色体号，位置，参考碱基和变异碱基
        chrom = row[0]
        pos = row[1]
        ref = row[2]
        alt = row[3]
        # 生成一个随机的ID，可以根据需要修改
        id = "rs" + str(pos)
        # 假设质量和过滤条件都是PASS，可以根据需要修改
        qual = "."
        filter = "PASS"
        # 假设深度都是100，可以根据需要修改
        info = "DP=100"
        # 拼接VCF格式的一行，并添加到字符串中
        vcf_string += f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\n"
    # 返回VCF格式的字符串
    return vcf_string


if __name__ == "__main__":
    # record = VcfRecord('x', 9, '>1>6', 'GC', 'GT,GA', 'AT=>1>3>5>6,>1>3>4>6;DP=2', 'GT', '1/1')
    record = VcfRecord()
    record.CHROM = 'x'
    record.POS = 9
    record.ID = '>1>6'

    header_str = vcf_header()
    print(header_str)
    print(record.to_string())
