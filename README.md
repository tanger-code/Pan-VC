# Pan-VC 

一个基于三代测序数据和泛基因组图进行小变异的检测和分型的工具。



# 输入文件

在本工具所在的代码目录下，需要提供四个输入文件（以21号染色体为例）：

（1）染色体的泛基因组图文件：`chr21.vg`

（2）样本的测序数据比对到泛基因组图的比对文件：`chr21_sort.gam`

（3）比对文件的索引：`chr21_sort.gam.gai`

（4）染色体的线性参考基因组文件：`GCA_chr15.fa`

注意，以上文件名除了染色体号可以改变之外，其他的要保持一样（比如chr1.vg, GCA_chr1.fa等）。



# 运行工具

本工具的主要流程已经集成到z_pipeline.sh中，进入到工具所在目录即可运行此脚本：

```
sh z_pipeline.sh
```

注意：脚本中用”*“标记的区域是所用的染色体及几种工具的地址，**需要手动修改**。



# 输出文件

程序运行结束之后，会产生三个目录：chunks, job_store, result。其中chunks是对比对文件进行分区之后的一系列结果，job_store是程序运行过程中产生的临时

文件，这两个目录都可以删除。result目录中，log.txt用于跟踪程序的执行进度，bed_all.txt是存在假重复序列变异的区域，这两个文件都可以删除；

sample_all.txt是最终的vcf结果。下面是一些注意事项：

- sample_all.txt中会存在很多重复的记录，需要进行过滤，方法是：将sample_all.txt改名为sample_all.vcf，然后将vcf_rmDup_1.py脚本复制到result目录

  中，运行vcf_rmDup_1.py脚本，之后就可生成sample_all_rmDup.vcf文件。

  ```
  cd result
  mv sample_all.txt sample_all.vcf
  cp ../vcf_rmDup_1.py ./
  python vcf_rmDup_1.py
  ```

- 如果需要将结果分为规范化之后的snp结果和indel结果，可以将res_process.sh复制到result目录中，其中”*“标记的区域是几种工具的地址和线性参考文件的地址，**需要手动修改**。之后运行脚本文件，就可得到分开的vcf文件。

  ```
  cd result
  cp ../res_process.sh ./
  sh res_process.sh
  ```

  