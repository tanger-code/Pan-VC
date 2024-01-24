import vcf

vcf_file = "sample_all.vcf"
vcf_rmDup_file = "sample_all_rmDup.vcf"

vcf_reader = vcf.Reader(filename=vcf_file)
vcf_writer = vcf.Writer(open(vcf_rmDup_file, "w"), vcf_reader)
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
