### 1. create_cisTarget_databases
git clone https://ghproxy.com/https://github.com/aertslab/create_cisTarget_databases
cd create_cisTarget_databases


### 2. virtual environment
conda create -n create_cistarget_databases \
    'python=3.10' \
    'numpy=1.21' \
    'pandas>=1.4.1' \
    'pyarrow>=7.0.0' \
    'numba>=0.55.1' \
    'python-flatbuffers'
conda activate create_cistarget_databases


### 3. install cbust
cd "${CONDA_PREFIX}/bin"
# 下载
git clone -b change_f4_output https://ghproxy.com/https://github.com/ghuls/cluster-buster/
cd cluster-buster
make cbust


### 4. install liftOver and bigWigAverageOverBed
cd "${CONDA_PREFIX}/bin"
# liftOver, 不同物种参考基因组坐标转换或同一物种不同版本参考基因组
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
# bigWigAverageOverBed，计算调控区域信号强度
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
# 修改权限
chmod a+x liftOver bigWigAverageOverBed


### 5. reference genome
genomefa='/data/hanjiangang/sc/pycistarget/database_sheep/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa'
consensdir='/data/hanjiangang/sc/pycistarget/database_sheep' # consensus_regions.bed


### 6. chromosome format
# >1 to > chr1

# 显示染色体名称
grep "^>" Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa | head
grep "^>" Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa | tail

# 添加 chr (1 to chr1)，注意如果修改了 .fa 文件需要删除 .fai 文件，否则会报错
sed -i 's/^>/>chr/g' Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa
grep "^>" Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa | head
grep "^>" Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa | tail

# 仅保留chr1，去掉之后的字符
# >chr 1 ....................  to >chr1
sed -i 's/>\([^ ]*\).*/>\1/' Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa
grep "^>" Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa | head
grep "^>" Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa | tail


### 7. extract fasta file based on consensus_region
bedtools getfasta -fi $genomefa -bed $consensdir/consensus_regions.bed -fo $outdir/consensus_regions.fa


### 1. 转录因子数据集
# Cluster-Buster motifs
# cluster-buster (.cb)格式 motifs 数据存储路径
wget https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
unzip v10nr_clust_public.zip
# 共10249个
# v10nr_clust_public/singletons/ 中的转录因子就是各种转录因子数据库的大杂烩，各种各种，人的、鼠的、果蝇的、以及只关注某些转录因子的数据库
# 同时数据库之间存在重叠，甚至有些数据库比如cispb、jaspa、matacluster直接会包含其他小的转录因子数据库，直接根本无法挑选
# 因此后边直接选了所有的转录因子进行后续分析
cbdir='/data/hanjiangang/sc/pycistarget/database_sheep/v10nr_clust_public/singletons/'

# -m motif_list
# 选择 cbdir 路径中要进行分析的 motifs，名称一定要完全对应，不然会报错
# 例如cbdir路径中为：dbtfbs__ZNF26_HEK293_ENCSR028EGI_merged_N1.cb，all_motifs.tsv 文件中为 dbtfbs__ZNF26_HEK293_ENCSR028EGI_merged_N1
cd v10nr_clust_public/singletons/ # 选择了全部
ls > all_motifs.tsv # 删除 .cb 后缀
motif_list='/data/hanjiangang/sc/pycistarget/database_sheep/v10nr_clust_public/motifs_subset.tsv' # 运行时间太长，仅选择个别数据库

'''
transfac_public
tfdimers
taipale
swissregulon
hocomoco
metacluster
'''



### 2. 其他数据和软件路径
# consensus_regions fasta 文件
consensdir='/data/hanjiangang/sc/pycistarget/database_sheep/consensus_regions.fa'

# CLUSTER_BUSTER_PATH
cbpath='/home/hanjiangang/anaconda3/envs/create_cistarget_databases/bin/cbust'

# 输出文件路径
outdir='/data/hanjiangang/sc/pycistarget/database_sheep'

# 输出文件前缀
tag='cluster_V10_DPCL_sheep'


### 3. create_cistarget_motif_databases
# 指定虚拟环境调用的 python 版本
/home/hanjiangang/anaconda3/envs/create_cistarget_databases/bin/python3.10 \
/data/hanjiangang/sc/pycistarget/create_cisTarget_databases/create_cistarget_motif_databases.py \
-f $consensdir \
-M $cbdir \
-m $motif_list \
-o $outdir/$tag \
-t 50 \
-c $cbpath \
-l
'''
-M $cbdir \                            #Path to directory with Cluster-Buster motifs.
-m $motif_list \                       #motif IDs to be scored from directory specified by "--motifs_dir".
-o $outdir/$tag \                      #输出文件前缀
-t 35 \                                #cpu
-c $cbpath\                            #CLUSTER_BUSTER_PATH
-l                                     #mask，fasta文件中小写字母转换为N
-p ${current_part} 1 \                 #内存不够可以分多（p）批处理
'''

'''命令日志，550980 regions x 10249 motifs
Initialize dataframe (550980 regions x 10249 motifs) for storing CRM scores for each regions per motif.
Adding Cluster-Buster CRM scores (1 of 358) for motif "cisbp__M00011" took 1.578365 seconds. 
Adding Cluster-Buster CRM scores (2 of 358) for motif "c2h2_zfs__M0373" took 0.600930 seconds. 
Adding Cluster-Buster CRM scores (3 of 358) for motif "bergman__slbo" took 0.426880 seconds. 
Adding Cluster-Buster CRM scores (4 of 358) for motif "cisbp__M00001" took 0.394991 seconds. 
Adding Cluster-Buster CRM scores (5 of 358) for motif "bergman__shn-ZFP2" took 0.389542 seconds. 
Adding Cluster-Buster CRM scores (6 of 358) for motif "c2h2_zfs__M0462" took 0.528468 seconds. 
Adding Cluster-Buster CRM scores (7 of 358) for motif "cisbp__M00006" took 0.588859 seconds. 
Adding Cluster-Buster CRM scores (8 of 358) for motif "bergman__tll" took 0.488174 seconds. 

Scoring 358 motifs with Cluster-Buster took: 2397.142409 seconds 

Writing cisTarget regions vs motifs scores db: "/data/hanjiangang/sc/pycistarget/database_sheep/cluster_V10_DPCL_sheep.motifs_v"r
Writing cisTarget regions vs motifs scores db took: 3.325903 seconds 

Writing cisTarget motifs vs regions scores db: "/data/hanjiangang/sc/pycistarget/database_sheep/cluster_V10_DPCL_sheep.regions_"r
Writing cisTarget motifs vs regions scores db took: 552.440004 seconds 

Create rankings from "/data/hanjiangang/sc/pycistarget/database_sheep/cluster_V10_DPCL_sheep.motifs_vs_regions.scores.feather" .2
Creating cisTarget rankings db from cisTarget scores db took: 57.358011 seconds 

Writing cisTarget motifs vs regions rankings db: "/data/hanjiangang/sc/pycistarget/database_sheep/cluster_V10_DPCL_sheep.region"r
Writing cisTarget motifs vs regions rankings db took: 115.978658 seconds 
'''

# 生成三个文件
cluster_V10_DPCL_sheep.motifs_vs_regions.scores.feather
cluster_V10_DPCL_sheep.regions_vs_motifs.rankings.feather
cluster_V10_DPCL_sheep.regions_vs_motifs.scores.feather


### 4. convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs
# 没做这一步，应该是不需要做，因为已经生成了 ranking.features
/home/hanjiangang/anaconda3/envs/create_cistarget_databases/bin/python3.10 \
/data/hanjiangang/sc/pycistarget/create_cisTarget_databases/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py \
-i $outdir/$tag.motifs_vs_regions.scores.feather -s 555
'''
Reading cisTarget motifs vs regions scores db: "/data/hanjiangang/sc/pycistarget/database_sheep/cluster_V10_DPCL_sheep.motifs_v"rehtaef.serocs.snoiger_s 

Reading cisTarget motifs vs regions scores db took: 437.246259 seconds 

Create rankings from "/data/hanjiangang/sc/pycistarget/database_sheep/cluster_V10_DPCL_sheep.motifs_vs_regions.scores.feather" with random seed set to 555. 
Creating cisTarget rankings db from cisTarget scores db took: 316.676783 seconds 

Convert motifs vs regions cisTarget rankings db to regions vs motifs cisTarget rankings db. 
Writing cisTarget motifs vs regions rankings db: "/data/hanjiangang/sc/pycistarget/database_sheep/cluster_V10_DPCL_sheep.regions_vs_motifs.rankings.feather" 
Writing cisTarget motifs vs regions rankings db took: 358.638110 seconds 
'''









