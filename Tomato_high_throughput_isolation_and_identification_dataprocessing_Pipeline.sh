[TOC]
#Citation: Zhou et al., Cross-kingdom synthetic microbiota supports tomato suppression of Fusarium wilt disease. Nat Commun 13, 7890 (2022). https://doi.org/10.1038/s41467-022-35452-6
# zhouxin5518@163.com
#print working directory
pwd 

# cd dir
cd /E/Tomato_culturable_data/Bateria_barcode/B1/samples
pwd

###########################
## 2.3 extract sequence Index (lane2library)

# check relationship between sequence and Index, Index usually reverse complementarity with sequence
zcat lane1_1.fq.gz | head
zcat lane1_2.fq.gz | head

# Filter dual-ended files of each library by index list file in parallel
# zcat、grep -A 3 search and print 4 lines
mkdir -p lib
for index in `cut -f 1 index.txt`; do \
	echo ${index}
	zcat lane1_1.fq.gz|grep -A 3 "#${index}"|grep -v -P '^--$' |gzip > lib/${index}_1.fq.gz
	zcat lane1_2.fq.gz|grep -A 3 "#${index}"|grep -v -P '^--$' |gzip > lib/${index}_2.fq.gz
done

# rename samples
awk 'BEGIN{OFS=FS="\t"}{system("mv samples/"$2".fq samples/"$1".fq");}' <(tail -n+2 ${l}.txt)


 cd /Processed_data/Tomato_culturable_data/Bateria_barcode/


gunzip  *.gz
l=library1
#a= LM181115-6001_L1_A001

# merge pair-end sequences
vsearch --fastq_mergepairs   LM181115-6004_L1_A004.R1.clean.fastq --reverse    LM181115-6004_L1_A004.R2.clean.fastq --fastqout B4_merged.fq

####################################################

# extract samples according to library.txt
mkdir -p samples
for index in `tail -n+2 ${l}.txt|cut -f 2`; do
	echo ${index}
	grep -A 2 -B 1 -P "^${index}"  B4_merged.fq|grep -v -P '^--$' > samples/${index}.fq
done 

# rename samples
awk 'BEGIN{OFS=FS="\t"}{system("mv samples/"$2".fq samples/"$1".fq");}' <(tail -n+2 ${l}.txt)


l=library42
mkdir -p isolated_strains
for index in `tail -n+2 ${l}.txt|cut -f 2`; do
	echo ${index}
	grep -A 2 -B 1 -P "${index}$"  samples/B42.fq|grep -v -P '^--$' > isolated_strains/${index}.fq
done

# rename each cell of 96-well plantes
awk 'BEGIN{OFS=FS="\t"}{system("mv isolated_strains/"$2".fq isolated_strains/"$1".fq");}' <(tail -n+2 ${l}.txt)


###########################################################################################
# Cut primers and quality filter
cd 
# cut 27bp barcode and reverse 28bp barcode
mkdir -p temp

for i in `tail -n+2 strain_name.txt | cut -f 1`;do
vsearch --fastx_filter isolated_strains/${i}.fq \
  --fastq_stripleft 27 --fastq_stripright 28 \
  --fastq_maxee_rate 0.01 \
  --fastaout temp/${i}.trimed.fa
done &
# 1m1s


#Dereplicate the sequences contained in queries.fas, take into account the abundance information already
#present, write unwrapped fasta sequences to queries_unique.fas with the new abundance information, discard all sequences with an abundance of 1:
vsearch --derep_fulllength queries.fas --sizein --fasta_width 0 --sizeout --output
queries_unique.fas --minuniquesize 8

mkdir -p derep

###derep sequences
for i in `tail -n+2 strain_name.txt | cut -f 1`;do
vsearch --derep_fulllength temp/${i}.trimed.fa --sizein --fasta_width 0 --sizeout --output \
derep/${i}.derep.fa --minuniquesize 100  --relabel ${i}. 
done &

cd undereped
find . -type f -name "*.fa" > 1.txt
## 
sed 's/\.\///g' 1.txt > 2.txt

sed "s/.derep.fa//g" 2.txt > 3.txt

mv 3.txt ../unclassified2.txt

mkdir -p derep_2

for i in `tail -n+2 unclassified2.txt | cut -f 1`;do
vsearch --derep_fulllength temp/${i}.trimed.fa --sizein --fasta_width 0 --sizeout --output \
derep_2/${i}.derep.fa --minuniquesize 8 --relabel ${i}. 
done &

cd /e/Tomato_culturable_data/Bateria_barcode/Isolated_Bacteria

mkdir -p derep/strain_otus
for i in `tail -n+2 strain_name.txt | cut -f 1`;do
 usearch11 -unoise3 derep/${i}.derep.fa \
   -zotus derep/strain_otus/${i}.otus.fa 
done &

mkdir result

cat derep/strain_otus/*.fa > result/all_strain.fa

grep ">" result/all_strain.fa | wc -l


mkdir result/zclusters
vsearch --cluster_size result/all_strain.fa \
 --centroids result/zotus.fa --clusters result/zclusters/zclusted.fas --id 0.99   
   
####################################
##extract all the names from cluster
cd result/
mkdir cluster_name
for i in `tail -n+2 cluster_name.txt | cut -f 1` ; do
grep ">" zclusters/${i} | sed 's/>//g'| cut -d ';' -f 1 > cluster_name/${i}.txt
done &


cut -d '.' -f 1 cluster1.txt > cluster2.txt

find . -type f -name "*.fa" > 1.txt
## find and delete ./
sed 's/\.\///g' 1.txt > 2.txt

sed "s/.derep.fa//g" 2.txt > 3.txt

mv 3.txt ../unclassified2.txt

for i in `tail -n+2 ./cluster_name.txt | cut -f 1` ; do
  grep ">" ./${i} | cat > cluster_name/${i}.txt

done &
######NCBI taxonomy annotation
blastn -max_target_seqs 3 -db /public/PubDatabase/ncbi/nt/nt -out zouts.blast -query zotus.fas \
-num_threads 8 -outfmt "6 qseqid qlen qstart qend stitle sseqid slen sstart send qcovs bitscore evalue pident"
##
usearch11 -sintax zotus.fas -db db/silva_16s_v138.fa -tabbedout zouts.tax --sintax_cutoff 0.9

# 24s, 103Mb
head -n2 zouts.tax

#rank: p c o f g，phylum, class, order, family, genus

mkdir -p tax
for i in p c o f g s; do
  usearch10 -sintax_summary zouts.tax \
  -rank ${i} \
  -output sum_${i}.txt
done &

# Taxonomy 
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' tax/sum_*.txt

cut -f 1,4 zouts.tax \
  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
  > taxonomy2.txt
head -n2 taxonomy2.txt

awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
  taxonomy2.txt > otus.tax
sed 's/;/\t/g;s/.__//g;' otus.tax|cut -f 1-8 | \
  sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
  > taxonomy.txt
head -n3 taxonomy.txt
###
#make OTU ML tree
mafft --maxiterate 1000 --globalpair otus.fa >otus_Aligned.fa

python fasta2phylip.py

raxmlHPC-PTHREADS-SSE3 -f a -s alpha_pep_align_gb.phy -n tre -m GTRGAMMA -x 1234 -# 1000 -T 8 &

# Find unique read sequences and abundances, miniuniqusize 

vsearch --derep_fulllength temp/all_nonchimera.fa \
  --output temp/uniques.fa --relabel Uni --minuniquesize 8 --sizeout 

mkdir -p clean
awk 'BEGIN{OFS=FS="\t"}{system("usearch10 -fastx_truncate samples/"$1".fq \
  -stripleft "$6" -stripright "$7" -fastqout clean/"$1".fq -relabel "$1".");}' ${l}.txt

head clean/S01.fq

rm -r samples library1.fq
gzip *.fq

#################
#obtain  reverse complementary sequence
grep -v ">" reverse_complement_barcode.fas |cat > complement_barcode.fas
