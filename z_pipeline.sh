#!/bin/bash



# **************************��Ҫ������Ⱦɫ��**************************
chromosome="chr15"
# *****************************���ߵ�ַ*******************************
bcftools=/usr/local/bin/bcftools
bedtools=/data/software/bedtools2/bin/bedtools
python=/data/home/tangen/.conda/envs/tee01/bin/python
vg=/data/home/tangen/.conda/envs/vg/bin/vg
# ********************************************************************

date +"%Y-%m-%d %H:%M:%S"

# ���������ļ��д�����ɵ��ļ�
mkdir job_store
mkdir result
mkdir chunks

# ͼת��
vg_graph=./${chromosome}.vg
$vg convert $vg_graph -f -t 10 > ./job_store/${chromosome}.gfa &
$vg convert $vg_graph -f -W -t 10 > ./job_store/${chromosome}_with_P.gfa &
$vg index $vg_graph -x ${chromosome}.xg -t 30 &
wait

gfa_graph=./job_store/${chromosome}.gfa
gfa_with_P_graph=./job_store/${chromosome}_with_P.gfa
aln_file=./${chromosome}_sort.gam

# pack��Ϣ
$vg pack -x $gfa_graph -g $aln_file -t 30 -o ./job_store/graph.pack &
$vg pack -x $gfa_graph -g $aln_file -D -t 30 > ./job_store/aln_edge_pack.txt &
wait

# ��ͼ����ȡsnarl�����Ϣ
$vg snarls $gfa_graph -t 30 > ./job_store/${chromosome}.snarls &
$vg call -t 30 -k ./job_store/graph.pack -T $gfa_graph -S GRCh38 > ./job_store/${chromosome}.vgCallTraversals.gaf &
wait

# ��ȡ�ο�·�������Ϣ
less $gfa_with_P_graph | grep "GRCh38#0#${chromosome}" > ./job_store/REF_GRCh38_${chromosome}.txt
$python ref_node_pos.py -c ${chromosome}
$vg find -x ${chromosome}.xg -N ./job_store/ref_nodes.txt -P GRCh38#0#${chromosome} > ./job_store/ref_nodes_pos.csv

# ��vgCallTraversals�ļ�����ȡsnarl���ơ�·��������֧�ֶȡ��Ƿ��ǲο�·�������ݲο�·������ȡ��snarl���н���
$python vgCallTraversals_to_newTraversals_with_ref.py -c ${chromosome}

# ��gam�ļ�ת��Ϊjson��ʽ�ļ�
$vg view --threads 40 -a $aln_file > ./job_store/graphaligner_aln_sort.json

# -------------------------------------------------------chunk-------------------------------------------------------
# ��gam������Ȼ��snarl��Ӧ��������
total_len_str=$($vg paths -x ${chromosome}.xg -E | grep "GRCh38#0#${chromosome}" | grep -o -E '[0-9]+' | tail -1)
total_len_pre=$((total_len_str))
total_len=$((total_len_pre-1))
step=100000

$vg chunk -x ${chromosome}.xg -a $aln_file -s $step -S ./job_store/${chromosome}.snarls -b ./chunks/chunk -t 40 -p GRCh38#0#${chromosome}

for file in chunks/*
do
    while [ $(jobs | wc -l) -ge 60 ]
    do
        sleep 1
    done
    $vg view -a $file > ${file%.*}.json &
done
wait
# -------------------------------------------------------------------------------------------------------------------

# ��reads��Ӧ��snarl·����
$python reads_to_traversals_linux.py -c ${chromosome} -l $((total_len+1))
$python traversals_agg_sort_by_pos.py
$python Genotyping_ALL_v4.py -c ${chromosome} -l $((total_len+1))

date +"%Y-%m-%d %H:%M:%S"