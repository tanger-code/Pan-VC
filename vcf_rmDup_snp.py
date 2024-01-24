import vcf

vcf_reader = vcf.Reader(filename="z_total_snps.vcf")
vcf_writer = vcf.Writer(open("sample_snps.vcf", "w"), vcf_reader)
header = vcf_reader.metadata

hash_dict = {}
for record in vcf_reader:
    pos = record.POS
    ref = record.REF
    alt = str(record.ALT[0])
    s = str(pos) + ref + alt
    h = hash(s)
    if h in hash_dict:
        pass
    else:
        hash_dict[h] = True
        vcf_writer.write_record(record)
vcf_writer.close()
