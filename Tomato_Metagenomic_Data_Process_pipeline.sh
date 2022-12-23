[TOC]
#Citation: Zhou et al., Cross-kingdom synthetic microbiota supports tomato suppression of Fusarium wilt disease. Nat Commun 13, 7890 (2022). https://doi.org/10.1038/s41467-022-35452-6
#Author Zhouxin
#Email: zhouxin5518@163.com
#print working directory
pwd 

export PYTHONPATH=/data/serverSoft/serverSoft/Anaconda2/lib/python2.7/site-packages
export PYTHONUSERBASE=/data/serverSoft/serverSoft/Anaconda2
#kneaddata location：/data/serverSoft/serverSoft/Anaconda2/bin/kneaddata
#humann2：/data/serverSoft/serverSoft/Anaconda2/bin/humann2

export PATH=/data/serverSoft/serverSoft/Anaconda2/bin/:$PATH

kneaddata -h

# 可选miniconda https://conda.io/miniconda.html
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda3-5.3.0-Linux-x86_64.sh

wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda2-5.3.0-Linux-x86_64.sh

conda config --set ssl_verify no
bash Anaconda3-5.3.0-Linux-x86_64.sh


unset PYTHONPATH

echo 'PATH=/opt/sysoft/git-2.17.0/bin/:$PATH' >> ~/.bashrc
source ~/.bashrc
#######################################################3
(metawrap) [Cailei@node22 ~]$ cat ~/.condarc 
channels:
  - ursky
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
  - bioconda
  - r
  - conda-forge
show_channel_urls: true
ssl_verify: true
###########################################################################

export PATH=/public/home/Cailei/anaconda3/bin:$PATH

export PYTHONPATH= /public/home/Cailei/anaconda3/lib/python3.7/site-packages

export PYTHONPATH=/data/serverSoft/serverSoft/Anaconda2.5.3/anaconda2/lib/python2.7/site-packages
export PYTHONUSERBASE=/data/serverSoft/serverSoft/Anaconda2.5.3/anaconda2


wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh

db=~/db

# fastqc
conda install fastqc
# multiqc
conda install multiqc

# kneaddata
conda install kneaddata
# 
kneaddata_database

##MetaPhlAn2、HUMAnN2

# install MetaPhlAn2、HUMAnN2
conda install humann2
# 
humann2_test

humann2_databases 
wd=~/db/humann2
mkdir -p $wd 
# 5.37G
humann2_databases --download chocophlan full $wd 

humann2_databases --download uniref uniref90_diamond $wd

humann2_config --print
#
humann2_config --update run_modes threads 8
humann2_config --update database_folders protein $wd/uniref
humann2_config --update database_folders nucleotide $wd/chocophlan

##kraken2: http://ccb.jhu.edu/software/kraken/MANUAL.html#custom-databases

# https://ccb.jhu.edu/software/kraken/
conda install kraken2

DBNAME=/db/kraken2
#kraken2-build --standard --threads 24 --db $DBNAME

kraken2-build --download-library bacteria --threads 24 --db db/kraken2/

kraken-build --db $DBNAME --clean

# megahit 
conda install megahit
# QUEST 
conda install quest 

# prokka 
conda install prokka
# cd-hit 
conda install cd-hit
# emboss transeq
conda install emboss
# salmon
conda install salmon


### COG/eggNOG http://eggnogdb.embl.de
#eggnog
conda install eggnog-mapper
#
mkdir -p ~/db/eggnog && cd ~/db/eggnog
download_eggnog_data.py --data_dir ./ -y -f euk bact arch viruses
#eggNOG
cp eggnog.db /dev/shm
# 
wget -c wget http://bailab.genetics.ac.cn/share/COG.anno
# 
wget -c wget http://bailab.genetics.ac.cn/share/KO.anno

###dbCAN2 http://cys.bios.niu.edu/dbCAN2/
mkdir -p ~/db/dbCAN2 && cd ~/db/dbCAN2
wget -c http://cys.bios.niu.edu/dbCAN2/download/CAZyDB.07312018.fa
wget -c http://cys.bios.niu.edu/dbCAN2/download/CAZyDB.07312018.fam-activities.txt
diamond makedb --in CAZyDB.07312018.fa --db CAZyDB.07312018
# 
grep -v '#' CAZyDB.07312018.fam-activities.txt|sed 's/  //'| \
  sed '1 i ID\tDescription' > fam_description.txt

### Resfam http://dantaslab.wustl.edu/resfams
mkdir -p resfam && cd ~/db/resfam
#
wget http://bailab.genetics.ac.cn/share/Resfams-proteins.dmnd
wget http://bailab.genetics.ac.cn/share/Resfams-proteins_class.tsv


# https://github.com/bxlab/metaWRAP
conda create -n metawrap python=2.7
source activate metawrap
conda config --add channels ursky
conda install -c ursky metawrap-mg 
 

mkdir -p ~/db

## CheckM for Bin
cd ~/db
mkdir checkm && cd checkm

wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz
rm *.gz

checkm data setRoot

conda install -c conda-forge numpy

cd ~/db
mkdir kraken
kraken-build --standard --threads 24 --db kraken
kraken-build --db kraken --clean

cd ~/db
mkdir NCBI_nt && cd NCBI_nt
wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
for a in nt.*.tar.gz; do tar xzf $a; done

cd ~/db
mkdir NCBI_tax
cd NCBI_tax
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz

##
which config-metawrap


# samtools
conda install samtools


cd ~/bin
wget http://bailab.genetics.ac.cn/share/mat_gene2ko.R
chmod +x mat_gene2ko.R

# VizBin
sudo apt-get install libatlas3-base libopenblas-base default-jre
curl -L https://github.com/claczny/VizBin/blob/master/VizBin-dist.jar?raw=true > VizBin-dist.jar
mv VizBin-dist.jar /usr/local/bin

###bowtie2-build
#Solanum lycopersicum
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.3_SL2.50/GCF_000188115.3_SL2.50_genomic.fna.gz

gunzip GCF_000188115.3_SL2.50_genomic.fna.gz

#bowtie2
bowtie2-build GCF_000188115.3_SL2.50_genomic.fna bt2tomato
bowtie2-build GCF_000188115.3_SL2.50_genomic.fna bt2tomato

mkdir -p doc 
mkdir -p seq

#ln -s /db/3metagenome/seq/*.fq seq/
#scp -r Cailei@124.16.144.62:/public/home/Cailei/Xinzhou/Metagenome_tomato_2017/Metagenomic_rawdata .

tree 
# design file
head -n3 doc/design.txt

head -n4 seq/p136C_1.fq
mkdir -p temp 

# Qaulity control and remove host contamination

# kneaddata_database
# kneaddata_database --download human_genome bowtie2 ./ # 

#reference_db = /public/home/projectuser/Xinzhou/DATABASE/bowtie2/tomato/bt2tomato

##ERROR: Unable to list files in trimmomatic directory: /data/serverSoft/serverSoft/Anaconda2/bin/trimmomati
  
# parallel --citation 
time parallel -j 3 --xapply \
  'kneaddata -i {1} -i {2} \
  -o temp/qc -v -t 20 --remove-intermediate-output \
  --trimmomatic /data/serverSoft/serverSoft/Anaconda2/share/trimmomatic-0.36-3 --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
  --bowtie2-options "--very-sensitive --dovetail" -db /public/home/projectuser/Xinzhou/DATABASE/bowtie2/tomato/bt2tomato' \
 ::: seq/*_R1_001.fastq ::: seq/*_R2_001.fastq

nohup gzip seq/*.fastq

kneaddata_read_count_table --input temp/qc --output temp/kneaddata_read_counts.txt
cat temp/kneaddata_read_counts.txt

mkdir -p temp/concat
for i in `tail -n+2 doc/design.txt | cut -f 1`;do \
  cat temp/qc/${i}_R1_001_kneaddata_paired* > temp/concat/${i}.fq; done
ls -lsh temp/concat/*.fq

humann2_config 
# metaphlan2 db_v20 and databases 



time humann2 --input temp/concat/Fan6h_S29_L007.fq  \
  --output temp/ \
  --nucleotide-database /data/serverSoft/serverSoft/Anaconda2/bin/db_v20/
  --threads 40 &

time parallel -j 3 \
  'humann2 --input {}  \
  --output temp/ ' \
  ::: temp/concat/*.fq > log 

mkdir -p metaphlan2

merge_metaphlan_tables.py temp/*_humann2_temp/*_metaphlan_bugs_list.tsv | \
  sed 's/_metaphlan_bugs_list//g' > metaphlan2/taxonomy.tsv

cp /db/bak/taxonomy.tsv metaphlan2/

metaphlan_to_stamp.pl metaphlan2/taxonomy.tsv > metaphlan2/taxonomy.spf


metaphlan_hclust_heatmap.py --in metaphlan2/taxonomy.tsv \
  --out metaphlan2/heatmap.pdf \
  -c bbcry --top 25 --minv 0.1 -s log 

# metaphlan2 to graphlan
export2graphlan.py --skip_rows 1,2 -i metaphlan2/taxonomy.tsv\
  --tree temp/merged_abundance.tree.txt \
  --annotation temp/merged_abundance.annot.txt \
  --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 \
  --annotations 5,6 --external_annotations 7 --min_clade_size 1
# graphlan annotation
graphlan_annotate.py --annot temp/merged_abundance.annot.txt temp/merged_abundance.tree.txt \
  temp/merged_abundance.xml
# output PDF figure, annoat and legend
graphlan.py temp/merged_abundance.xml metaphlan2/graphlan.pdf --external_legends --dpi 300 


sed '1 s/p[0-9]*//g' metaphlan2/taxonomy.tsv | grep -v '#' > metaphlan2/lefse.txt
# 
lefse-format_input.py metaphlan2/lefse.txt temp/input.in -c 1 -o 1000000
# 
run_lefse.py temp/input.in temp/input.res
# 
lefse-plot_cladogram.py temp/input.res metaphlan2/lefse_cladogram.pdf --format pdf --dpi 600 
# 
lefse-plot_res.py temp/input.res metaphlan2/lefse_res.pdf --format pdf --dpi 600
# 
# sort -k3,3n temp/input.res |less -S
lefse-plot_features.py -f one --feature_name "k__Bacteria.p__Firmicutes" --format pdf \
  temp/input.in temp/input.res metaphlan2/Firmicutes.pdf 
# 
lefse-plot_features.py -f diff --archive none --format pdf \
  temp/input.in temp/input.res metaphlan2/features

mkdir -p humann2


humann2_join_tables --input temp/ --file_name pathabundance --output humann2/pathabundance.tsv
sed -i 's/_Abundance//g' humann2/pathabundance.tsv

cp /db/bak/pathabundance.tsv humann2/


humann2_renorm_table --input humann2/pathabundance.tsv --units relab \
  --output humann2/pathabundance_relab.tsv

humann2_split_stratified_table --input humann2/pathabundance_relab.tsv \
  --output humann2/

head -n 1 humann2/pathabundance_relab_stratified.tsv | \
  sed 's/# Pathway/MetaCyc_pathway/' \
  > humann2/pathabundance_relab_stratified_LACTOSECAT-PWY.spf
grep "LACTOSECAT-PWY" humann2/pathabundance_relab_unstratified.tsv \
  >> humann2/pathabundance_relab_stratified_LACTOSECAT-PWY.spf


#MEGAHIT De novo Metagenomic Pipeline

## kraken2 based on NCBI
mkdir -p temp/kraken2
time parallel -j 3 \
  'kraken2 --db /db/kraken2 --paired temp/qc/{1}_1_kneaddata_paired*.fastq \
  --threads 3 --use-names --use-mpa-style --report-zero-counts \
  --report temp/kraken2/{1}_report \
  --output temp/kraken2/{1}_output' \
  ::: `tail -n+2 result/design.txt | cut -f 1`


mkdir -p result/kraken2
parallel -j 6 \
  'cut -f 2 temp/kraken2/{1}_report | sed "1 s/^/{1}\n/" > temp/kraken2/{1}_count ' \
  ::: `tail -n+2 result/design.txt | cut -f 1`
header=`tail -n 1 result/design.txt | cut -f 1`
cut -f 1 temp/kraken2/${header}_report | sed "1 s/^/Taxonomy\n/" > temp/kraken2/0header_count
paste temp/kraken2/*count > result/kraken2/taxonomy_count.txt


## Megahit Assembly
 
rm -rf temp/megahit
time megahit -t 9 \
  -1 `ls temp/qc/*_1_kneaddata_paired_1.fastq|tr '\n' ','|sed 's/,$//'` \
  -2 `ls temp/qc/*_1_kneaddata_paired_2.fastq|tr '\n' ','|sed 's/,$//'` \
  -o temp/megahit # --k-min 27 --k-max 191 --k-step 20

head temp/megahit/final.contigs.fa

# quast
mkdir -p result/megahit/
ln temp/megahit/final.contigs.fa result/megahit/
quast.py result/megahit/final.contigs.fa -o result/megahit/


## Gene annotation & quantitfy


ll temp/megahit/final.contigs.fa # 47Mb
time prokka temp/megahit/final.contigs.fa --outdir temp/prokka \
  --prefix mg --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses \
  --force --cpus 9


### cd-hit
mkdir -p temp/NR

time cd-hit-est -i temp/prokka/mg.ffn -o temp/NR/mg.ffn.nr \
  -aS 0.9 -c 0.95 -G 0 -M 0 -T 9 -n 5 -d 0 -g 1

mkdir -p result/NR
ln temp/NR/mg.ffn.nr result/NR/nucleotide.fa

transeq -sequence result/NR/nucleotide.fa -outseq result/NR/protein.fa

sed -i 's/_1 / /' result/NR/protein.fa

###salmon
mkdir -p temp/salmon
#
salmon index -t result/NR/nucleotide.fa -p 9 --type quasi -k 31 \
  -i temp/salmon/index 
# 
time parallel -j 3 \
  'salmon quant -i temp/salmon/index -l A -p 3 --meta \
  -1 temp/qc/{1}_1_kneaddata_paired_1.fastq -2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
  -o temp/salmon/{1}.quant' \
  ::: `tail -n+2 result/design.txt | cut -f 1`
mkdir -p result/salmon
salmon quantmerge --quants temp/salmon/*.quant -o result/salmon/gene.TPM
salmon quantmerge --quants temp/salmon/*.quant --column NumReads -o result/salmon/gene.count
sed -i '1 s/.quant//g' result/salmon/gene.*


##COG/eggNOG/KEGG

mkdir -p temp/eggnog
time emapper.py -m diamond --no_annot --no_file_comments --data_dir /db/eggnog \
  --cpu 9 -i result/NR/protein.fa -o temp/eggnog/protein.fa --override
# 比对结果功能注释, real 6s, user 44s 
time emapper.py --annotate_hits_table temp/eggnog/protein.fa.emapper.seed_orthologs --no_file_comments \
		-o temp/eggnog/output --cpu 9 --data_dir /db/eggnog --override

mkdir -p result/eggnog
sed '1 i Name\teggNOG\tEvalue\tScore\tGeneName\tGO\tKO\tBiGG\tTax\tOG\tBestOG\tCOG\tAnnotation' \
  temp/eggnog/output.emapper.annotations > temp/eggnog/output

cut -f 1,12 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$'>temp/eggnog/1cog.list

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/eggnog/1cog.list result/salmon/gene.count | \
	sed '/\t$/d' | sed '1 s/COG/KO/' > temp/eggnog/gene_cog.count
# COG,count,RPM
mat_gene2ko.R -i temp/eggnog/gene_cog.count -o result/eggnog/cogtab -n 1000000

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
  /db/eggnog/COG.anno result/eggnog/cogtab.count > result/eggnog/cogtab.count.spf

cut -f 1,7 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$' > temp/eggnog/2ko.list
# KO
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/eggnog/2ko.list \
  result/salmon/gene.count | sed '/\t$/d' > temp/eggnog/gene_ko.count

mat_gene2ko.R -i temp/eggnog/gene_ko.count -o result/eggnog/kotab -n 1000000
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' /db/eggnog/KO.anno \
  result/eggnog/kotab.count | sed 's/^\t/Undescription\t/' > result/eggnog/kotab.count.spf

###dbCAN2, CAZy
mkdir -p temp/dbcan2
time diamond blastp --db /db/dbcan2/CAZyDB.07312018 --query result/NR/protein.fa \
	--outfmt 6 --threads 9 --max-target-seqs 1 --quiet \
	--out temp/dbcan2/gene_diamond.f6

mkdir -p result/dbcan2

cut -f 1,2 temp/dbcan2/gene_diamond.f6 | uniq | sed 's/|/\t/g' | cut -f 1,3 | \
	cut -f 1,2 -d '_' |sed '1 i Name\tKO' > temp/dbcan2/gene_fam.list

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/dbcan2/gene_fam.list \
  result/salmon/gene.count | sed '/\t$/d' > temp/dbcan2/gene_fam.count

mat_gene2ko.R -i temp/dbcan2/gene_fam.count -o result/dbcan2/cazytab

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$$2} NR>FNR{print a[$1],$$0}' /db/dbcan2/fam_description.txt \
	result/dbcan2/cazytab.count > result/dbcan2/cazytab.count.spf


###ResFam
mkdir -p temp/resfam

time diamond blastp --db /db/resfam/Resfams-proteins.dmnd --query result/NR/protein.fa \
	--outfmt 6 --threads 9 --max-target-seqs 1 --quiet \
	--out temp/resfam/gene_diamond.f6

mkdir -p result/resfam
cut -f 1,2 temp/resfam/gene_diamond.f6 | uniq | \
  sed '1 i Name\tResGeneID' > temp/resfam/gene_fam.list

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
  temp/resfam/gene_fam.list result/salmon/gene.count | \
	sed '/^\t/d' > result/resfam/resfam.count

wc -l result/23salmon_gene/gene.count
wc -l result/resfam/resfam.count # 172/7734=2.2%
# stamp
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4"\t"$3"\t"$2} NR>FNR{print a[$1],$0}' \
  /db/resfam/Resfams-proteins_class.tsv  result/resfam/resfam.count \
  > result/resfam/resfam.count.spf

source deactivate
conda info -e
conda env list

conda env remove -n humann

db=/data/db

soft=/data/zhouxin/miniconda2

wd=/data/zhouxin/Tomato_project/Tomato_Metagenome
Rscript='/home/user/anaconda2/bin/Rscript --vanilla'
#
mkdir -p $wd && cd $wd
soft=~/miniconda2

conda activate humann
conda activate metawrap
conda activate metagenome_env

humann3_databases 

cd ${db}
mkdir -p ${db}/humann3


humann2_databases --download chocophlan full ${db}/humann2 --update-config yes

#To upgrade your protein database:
humann2_databases --download uniref uniref90_diamond ${db}/humann2 --update-config yes

humann2_databases --download uniref uniref90_ec_filtered_diamond ${db}/humann2/uniref --update-config yes

humann_databases --download uniref90_ec_filtered_diamond ${db}/humann3/uniref --update-config yes

##To upgrade your annotations database:
humann2_databases --download utility_mapping full ${db}/humann3 --update-config yes
humann_config --print
humann_config --update database_folders nucleotide ${db}/humann3/chocophlan/
humann_config --update database_folders protein ${db}/humann3/new_uniref/
humann_config --update database_folders utility_mapping ${db}/humann3/utility_mapping

humann_config --update

humann2_config --update database_folders nucleotide ${db}/humann3/chocophlan
humann2_config --update database_folders protein ${db}/humann3/new_uniref
humann2_config --update database_folders utility_mapping ${db}/humann3/utility_mapping
humann --threads 20 --input ./p143C_R1.fastq  --output humann143_out/

metaphlan ./p143C_R1.fastq --input_type fastq > SRS014476-Supragingival_plaque_profile.txt

metaphlan --install --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db /data/zhouxin/miniconda2/envs/humann3/lib/python3.7/site-packages/metaphlan/metaphlan_databases

############################

conda create -n humann -c bioconda metaphlan=3.0=pyh5ca1d4c_4
##https://github.com/biobakery/MetaPhlAn/issues/103

mkdir -p ${db}/metaphlan2 && cd ${db}/metaphlan2
wget -c http://210.75.224.110/share/meta/metaphlan2/mpa_v20_m200.tar
tar xvf mpa_v20_m200.tar
bzip2 -d mpa_v20_m200.fna.bz2
bowtie2-build mpa_v30_CHOCOPhlAn_201901.fna mpa_v20_m200

db=~/db
soft=~/miniconda2
mkdir -p ${soft}/bin/db_v20
ln -s ${db}/metaphlan2/* ${soft}/bin/db_v20/
mkdir -p /data/zhouxin/miniconda2/envs/biobakery3/bin/db_v20
ln -s ${db}/metaphlan2/* /data/zhouxin/miniconda2/envs/biobakery3/bin/db_v20/
mkdir -p /data/zhouxin/miniconda2/envs/humann/bin/metaphlan_databases
ln -s /data/zhouxin/db/metaphlan_databases/* /data/zhouxin/miniconda2/envs/humann3/lib/python3.7/site-packages/metaphlan/metaphlan_databases/


mkdir -p temp result
fastqc seq/*.gz -t 40
multiqc -d seq/ -o result/qc

# Solanum lycopersicum
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.3_SL2.50/GCF_000188115.3_SL2.50_genomic.fna.gz

gunzip GCA_002954035.1_ASM295403v1_genomic.fna.gz

bowtie2-build GCA_002954035.1_ASM295403v1_genomic.fna bt2tomato -t 40

     time kneaddata -i seq/HI4526_5_R1.fastq.gz -i seq/HI4526_5_R2.fastq.gz \
      -o temp/qc -v -t 40 --remove-intermediate-output \
       --trimmomatic /data/zhouxin/miniconda2/share/trimmomatic --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:50' \
       --bowtie2-options '--very-sensitive --dovetail' -db ${db}/Tomato_genome/bt2tomato

    nohup parallel -j 6 --xapply \
      "kneaddata -i {1} -i {2} \
      -o temp/qc -v -t 10 --remove-intermediate-output \
      --trimmomatic ${soft}/share/trimmomatic/ --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:50' \
      --bowtie2-options '--very-sensitive --dovetail' -db ${db}/Tomato_genome/bt2tomato" \
      ::: seq/*_R1.fastq.gz ::: seq/*_R2.fastq.gz &

    kneaddata_read_count_table --input temp/qc --output result/01kneaddata_sum.txt
    cat result/01kneaddata_sum.txt

    fastqc temp/qc/*_R1_kneaddata_paired_* -t 40
    multiqc -d temp/qc/ -o result/qc/

# HUMAnN2 Reference-based Metagenomic Pipeline
# https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(Humann2)
    cd ${wd}

    mkdir -p temp/concat
    # for
    for i in `tail -n+2 result/metadata1.tsv | cut -f 1`;do \
      cat temp/qc/${i}*_R1_kneaddata_paired_?.fastq > temp/concat/${i}.fq; done
    # 
    ll temp/concat/*.fq 

mkdir -p temp/humann2
nohup parallel -j 14 \
    'humann --input {}  \
      --output temp/humann3/ --threads 8' \
      ::: temp/concat/*.fq > temp/log &
    cat temp/log
ls -shl temp/humann3/
  
mkdir -p result/metaphlan2

merge_metaphlan_tables.py temp/humann3/*_humann_temp/*_metaphlan_bugs_list.tsv | \
   sed 's/_metaphlan_bugs_list//g' > result/metaphlan2/taxonomy.tsv

${db}/script/metaphlan_to_stamp.pl result/metaphlan2/taxonomy.tsv > result/metaphlan2/taxonomy.spf

metaphlan_hclust_heatmap.py --in result/metaphlan2/taxonomy.tsv \
        --out result/metaphlan2/heatmap.pdf \
        -c jet --top 40 --minv 0.1 -s log


mkdir -p result/humann3
humann2_join_tables --input temp/humann3 --file_name pathabundance \
    --output result/humann3/pathabundance.tsv
sed -i 's/_Abundance//g' result/humann3/pathabundance.tsv
    
humann2_renorm_table --input result/humann3/pathabundance.tsv --units relab \
      --output result/humann3/pathabundance_relab.tsv

humann2_split_stratified_table --input result/humann2/pathabundance_relab.tsv \
      --output result/humann3/


# metaphlan2 to graphlan
export2graphlan.py --skip_rows 1,2 -i result/metaphlan2/taxonomy.tsv \
      --tree temp/merged_abundance.tree.txt \
      --annotation temp/merged_abundance.annot.txt \
      --most_abundant 1000 --abundance_threshold 20 --least_biomarkers 10 \
      --annotations 3,4 --external_annotations 7
# graphlan annotation
graphlan_annotate.py --annot temp/merged_abundance.annot.txt \
      temp/merged_abundance.tree.txt  temp/merged_abundance.xml
# output PDF figure, annoat and legend
graphlan.py temp/merged_abundance.xml result/metaphlan2/graphlan.pdf \
    --external_legends 

sed '1 s/*[0-9]//g' result/metaphlan2/taxonomy.tsv | grep -v '#' > result/metaphlan2/lefse.txt

lefse-format_input.py result/metaphlan2/lefseDay1v15.txt temp/input.in -c 1 -o 1000000
#lefse-format_input.py result/metaphlan2/BAC/lefse.txt temp/input.in -c 1 -o 1000000
run_lefse.py temp/input.in temp/input.res

lefse-plot_cladogram.py temp/input.res result/metaphlan2/lefse_cladogram.pdf --format pdf --dpi 600
#lefse-plot_cladogram.py temp/input.res result/metaphlan2/BAC/lefse_cladogram.pdf --format pdf --dpi 600

lefse-plot_res.py temp/input.res result/metaphlan2/lefse_res.pdf --format pdf --dpi 600
#lefse-plot_res.py temp/input.res result/metaphlan2/BAC/lefse_res.pdf --format pdf --dpi 600

grep -v '-' temp/input.res | sort -k3,3n
lefse-plot_features.py -f one --feature_name "k__Eukaryota.p__Ascomycota" --format pdf \
    temp/input.in temp/input.res result/metaphlan2/Ascomycota.pdf

lefse-plot_features.py -f diff --archive none --format pdf \
temp/input.in temp/input.res result/metaphlan2/BAC/features


## 2.7 kraken2
mkdir -p temp/kraken2
parallel -j 3 \
    "kraken2 --db ${db}/kraken2 --paired temp/qc/{1}_R1_kneaddata_paired*.fastq \
    --threads 40 --use-names --use-mpa-style --report-zero-counts \
    --report temp/kraken2/{1}_report \
    --output temp/kraken2/{1}_output" \
    ::: `tail -n+2 result/metadata.tsv | cut -f 1`


mkdir -p result/kraken2
parallel -j 10 \
    'sort temp/kraken2/{1}_report | cut -f 2 | sed "1 s/^/{1}\n/" > temp/kraken2/{1}_count ' \
      ::: `tail -n+2 result/metadata.tsv | cut -f 1`
header=`tail -n 1 result/metadata.tsv | cut -f 1`
sort temp/kraken2/${header}_report | cut -f 1 | sed "1 s/^/Taxonomy\n/" > temp/kraken2/0header_count
paste temp/kraken2/*count > result/kraken2/taxonomy_count.txt
    
Rscript ${db}/script/kraken2alpha.R -i result/kraken2/taxonomy_count.txt

Rscript ${db}/script/alpha_boxplot.R -i result/kraken2/taxonomy_count.alpha.txt -t shannon \
      -d result/metadata.tsv -n group -w 8 -e 5

for i in richness chao1 ACE shannon simpson invsimpson;do
Rscript ${db}/script/alpha_boxplot.R -i result/kraken2/taxonomy_count.alpha.txt -t ${i} \
    -d result/metadata.tsv -n group1 -w 8 -e 5
done

cd ${wd}

#Megahit Assembly
rm -rf temp/megahit

nohup megahit -t 120 \
   -1 `ls temp/qc/*_1_kneaddata_paired_1.fastq|tr '\n' ','|sed 's/,$//'` \
   -2 `ls temp/qc/*_1_kneaddata_paired_2.fastq|tr '\n' ','|sed 's/,$//'` \
   -o temp/megahit --k-min 27 --k-max 191 --k-step 20 &

ll temp/megahit/final.contigs.fa

head -n2 temp/megahit/final.contigs.fa | cut -c1-60

mkdir -p temp/metaspades
     time metaspades.py -t 140 -m 500 \
      `ls temp/qc/*_1_kneaddata_paired_1.fastq|sed 's/^/-1 /'| tr '\n' ' '` \
       `ls temp/qc/*_1_kneaddata_paired_2.fastq|sed 's/^/-2 /'| tr '\n' ' '` \
       -k 21,33,55,77,99 \
	   -o temp/metaspades 
	   
ll metaspades/contigs.fasta  

mkdir -p result/megahit/
ln temp/megahit/final.contigs.fa result/megahit/
quast.py result/megahit/final.contigs.fa -o result/megahit/


##Gene annotation & quantitfy
# https://2017-cicese-metagenomics.readthedocs.io/en/latest/prokka_tutorial.html
#Prokka: rapid prokaryotic genome annotation https://www.ncbi.nlm.nih.gov/pubmed/24642063

#Prokka http://www.vicbioinformatics.com/software.prokka.shtml
#Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics. 2014 Jul 15;30(14):2068-9. PMID:24642063
ll temp/megahit/final.contigs.fa # 31Mb
conda activate metawrap
#export PERL5LIB=$PERL5LIB:${soft}/envs/metawrap/lib/site_perl/5.26.2:/disk2/home/liuyongxin/miniconda2/lib/site_perl/5.26.2
prokka temp/megahit/final.contigs.fa --outdir temp/prokka \
  --prefix mg --metagenome --kingdom Bacteria,Viruses \
  --force --cpus 30

conda deactivate

### 3.2.2 cd-hit 
mkdir -p temp/NR

grep -c '>' temp/prokka/mg.ffn 

nohup cd-hit-est -i temp/prokka/mg.ffn -o temp/NR/mg.ffn.nr \
        -aS 0.9 -c 0.95 -G 0 -M 0 -T 80 -g 1 &
   
grep -c '>' temp/NR/mg.ffn.nr

mkdir -p result/NR
ln temp/NR/mg.ffn.nr result/NR/nucleotide.fa

transeq -sequence result/NR/nucleotide.fa -outseq result/NR/protein.fa -trim Y

sed -i 's/_1 / /' result/NR/protein.fa

### salmon
mkdir -p temp/salmon

${soft}/envs/metagenome_env/share/salmon/bin/salmon index -t result/NR/nucleotide.fa -p 9 -k 31 \
      -i temp/salmon/index # --type quasi 

#salmon  v0.9.1
parallel -j 6 \
        "${soft}/envs/metawrap/bin/salmon quant -i temp/salmon/index -l A -p 20 --meta \
        -1 temp/qc/{1}_R1_kneaddata_paired_1.fastq \
        -2 temp/qc/{1}_R1_kneaddata_paired_2.fastq \
        -o temp/salmon/{1}.quant" \
        ::: `tail -n+2 result/metadata.tsv | cut -f 1`

mkdir -p result/salmon
${soft}/envs/metawrap/bin/salmon quantmerge \
        --quants temp/salmon/*.quant \
        -o result/salmon/gene.TPM
${soft}/envs/metawrap/bin/salmon quantmerge \
        --quants temp/salmon/*.quant \
        --column NumReads -o result/salmon/gene.count
sed -i '1 s/.quant//g' result/salmon/gene.*

head -n3 result/salmon/gene.*


### eggNOG/COG/KEGG/GO

# diamond eggNOG, emapper-0.12.7
mkdir -p temp/eggnog
nohup emapper.py -m diamond --no_annot --no_file_comments \
    --data_dir ${db}/eggnog --cpu 40 -i result/NR/protein.fa \
    -o temp/eggnog/protein --override &

time emapper.py --annotate_hits_table \
     temp/eggnog/protein.emapper.seed_orthologs --no_file_comments \
    -o temp/eggnog/output --cpu 40 --data_dir ${db}/eggnog --override

mkdir -p result/eggnog
sed '1 i Name\teggNOG\tEvalue\tScore\tGeneName\tGO\tKO\tBiGG\tTax\tOG\tBestOG\tCOG\tAnnotation' \
    temp/eggnog/output.emapper.annotations > temp/eggnog/output

cut -f 1,12 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$' \
    >temp/eggnog/1cog.list

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' \
temp/eggnog/1cog.list result/salmon/gene.count | \
sed '/\t$/d' | sed '1 s/COG/KO/' > temp/eggnog/gene_cog.count

Rscript ${db}/script/mat_gene2ko.R -i temp/eggnog/gene_cog.count \
   -o result/eggnog/cogtab -n 1000000

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
    ${db}/eggnog/COG.anno result/eggnog/cogtab.count > \
    result/eggnog/cogtab.count.spf

cut -f 1,7 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$' \
    > temp/eggnog/2ko.list

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' \
    temp/eggnog/2ko.list \
    result/salmon/gene.count | sed '/\t$/d' > temp/eggnog/gene_ko.count

Rscript ${db}/script/mat_gene2ko.R -i temp/eggnog/gene_ko.count \
    -o result/eggnog/kotab -n 1000000

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
    ${db}/eggnog/KO.anno result/eggnog/kotab.count | \
sed 's/^\t/Undescription\t/' > result/eggnog/kotab.count.spf
      
head result/eggnog/kotab.count.spf

mkdir -p temp/dbcan2
nohup diamond blastp --db ${db}/dbCAN2/CAZyDB.07312018 --query result/NR/protein.fa \
    	--outfmt 6 --threads 40 --max-target-seqs 1 --quiet -e 1e-5 --sensitive \
    	--out temp/dbcan2/gene_diamond.f6 &

mkdir -p result/dbcan2
cut -f 1,2 temp/dbcan2/gene_diamond.f6 | uniq | sed 's/|/\t/g' | cut -f 1,3 | \
    	cut -f 1,2 -d '_' |sed '1 i Name\tKO' > temp/dbcan2/gene_fam.list

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/dbcan2/gene_fam.list \
      result/salmon/gene.count | sed '/\t$/d' > temp/dbcan2/gene_fam.count

# /anaconda2/bin/Rscript --vanilla ${db}/script/mat_gene2ko.R -i temp/dbcan2/gene_fam.count -o result/dbcan2/cazytab
Rscript ${db}/script/mat_gene2ko.R -i temp/dbcan2/gene_fam.count -o result/dbcan2/cazytab

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' ${db}/dbCAN2/fam_description.txt \
    	result/dbcan2/cazytab.count > result/dbcan2/cazytab.count.spf


###ResFam
mkdir -p temp/resfam
nohup diamond blastp --db ${db}/resfam/Resfams-proteins --query result/NR/protein.fa \
    	--outfmt 6 --threads 20 --max-target-seqs 1 --quiet -e 1e-5 --sensitive \
    	--out temp/resfam/gene_diamond.f6 &
    
mkdir -p result/resfam

cut -f 1,2 temp/resfam/gene_diamond.f6 | uniq | \
sed '1 i Name\tResGeneID' > temp/resfam/gene_fam.list
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
      temp/resfam/gene_fam.list result/salmon/gene.count | \
    	sed '/^\t/d' > result/resfam/resfam.count

wc -l result/salmon/gene.count
wc -l result/resfam/resfam.count
    

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4"\t"$3"\t"$2} NR>FNR{print a[$1],$0}' \
      ${db}/resfam/Resfams-proteins_class.tsv  result/resfam/resfam.count \
      > result/resfam/resfam.count.spf


### CARD Binning MetaWRAP
cd ${wd}
mkdir -p binning && cd binning
mkdir -p temp && cd temp

mkdir -p qc
cd qc

ln -s /data/zhouxin/Metagenomic_rawdata/TomatoProject/temp/qc/*paired_1.fastq ./
	ln -s /data/zhouxin/Metagenomic_rawdata/TomatoProject/temp/qc/*paired_1.fastq ./
cd ..

mkdir -p megahit
cd megahit
# gunzip *.gz
ln -s /data/zhouxin/Metagenomic_rawdata/TomatoProject/temp/megahit/final.contigs.fa ./

wd=/data/zhouxin/Metagenomic_rawdata/TomatoProject/
cd ${wd}/binning
conda activate metawrap


## bin
nohup metawrap binning -o temp/binning -t 80 -a temp/megahit/final.contigs.fa \
      --metabat2 --maxbin2 --concoct temp/qc/*R1_kneaddata_paired*.fastq &
cd ${wd}/binning
# /bin/rm -rf temp/bin_refinement
##prodig
nohup metawrap bin_refinement -o temp/bin_refinement -t 80 \
      -A temp/binning/maxbin2_bins/ \
      -c 70 -x 5 &
# cat temp/bin_refinement/metawrap_bins.stats | awk '$2>70 && $3<5' | wc -l
wc -l temp/bin_refinement/binsA.stats
    
head temp/bin_refinement/binsA.stats

##Bin salmon
nohup metawrap quant_bins -b temp/bin_refinement/binsA -t 80 \
      -o temp/bin_quant -a /data/zhouxin/Tomato_project/Tomato_Metagenome/temp/megahit/final.contigs.fa temp/qc/*.fastq  &

ls -l temp/bin_quant/bin_abundance_heatmap.png

## Bin Taxator-tk
nohup metawrap classify_bins -b temp/bin_refinement/binsA \
      -o temp/bin_classify -t 80 &

metaWRAP annotate_bins -o temp/bin_annotate -t 40 \
-b temp/bin_refinement/metawrap_50_10_bins


nohup metawrap blobology -a temp/megahit/final.contigs.fa -t 30 \
      -o temp/bloblogy --bins temp/bin_refinement/metawrap_50_10_bins \
      temp/qc/ERR*.fastq &
 
cat /data/zhouxin/Tomato_project/Tomato_Metagenome/temp/qc/*_1.fastq > temp/qc/all_1.fq
cat /data/zhouxin/Tomato_project/Tomato_Metagenome/temp/qc/*_2.fastq > temp/qc/all_2.fq
nohup metawrap reassemble_bins -o temp/bin_reassemble \
      -1 temp/qc/all_1.fq -2 temp/qc/all_2.fq -t 40 -m 500 \
      -c 50 -x 10 -b temp/bin_refinement/binsA &
 conda deactivate
