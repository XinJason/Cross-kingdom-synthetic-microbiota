[TOC]
#!/bin/bash
# Xin Zhou

db=../../../../Public_Scripts

    mkdir -p seq
    mkdir -p result
    mkdir -p temp

    ls seq/
  
for i in `tail -n+2 result/metadata.tsv | cut -f 1`;do    
usearch11 -fastq_mergepairs seq/${i}.R1.fq -reverse seq/${i}.R2.fq \
-fastq_trunctail 30 -fastq_minlen 160 \
-fastq_maxdiffs 15 -fastq_pctid 80 \
-fastq_minmergelen 300 -fastq_maxmergelen 500 \
-fastqout temp1/${i}.merged.fq -relabel ${i}.
 done &
 

   nohup cat temp/*.merged.fq > temp/all.fq &

    ls -lsh temp/all.fq

    head -n 6 temp/all.fq | cut -c1-60

usearch11 -fastx_truncate temp/all.fq \
  -stripleft 10 -stripright 10 \
  -fastqout temp/stripped.fq
  
# fastq filter, keep reads error rates less than 1%
usearch11 -fastq_filter temp/all.fq \
  -fastq_maxee_rate 0.01 \
  -fastaout temp/filtered.fa 
##OTU/ASV Dereplicate and cluster/denoise
###Dereplication
# miniuniqusize 10

usearch11 -fastx_uniques  temp/filtered.fa \
-fastaout temp/uniques.fa -minuniquesize 10 -sizeout -relabel Uni   

usearch11 -unoise3 temp/uniques.fa \
  -zotus temp/zotus.fa

head -n 2 temp/otus.fa

db=../../../Public_Scripts
mkdir -p result/raw
nohup usearch11 -uchime2_ref temp/otus.fa -db ${db}/usearch/silva_16s_v123.fa \
     -chimeras temp/otus_chimeras.fa -strand plus -mode high_confidence -threads 40 &

cat temp/otus.fa temp/otus_chimeras.fa| grep '>' | sort | uniq -u \
  | sed 's/>//' > temp/non_chimeras.id

usearch11 -fastx_getseqs temp/otus.fa -labels temp/non_chimeras.id \
    -fastaout result/raw/otus.fa

sed -i 's/\r//g' result/raw/otus.fa

 # Creat Feature table
nohup usearch11 -otutab temp/filtered.fa -otus result/raw/otus.fa \
   -otutabout result/raw/otutab.txt -threads 40 &

usearch11 -sintax result/raw/otus.fa -db ./All_isolated_bacteria_database.fas \
-strand both -tabbedout result/raw/Owndatabase/otus.sintax -sintax_cutoff 0.9 

usearch11 -sintax result/otus.fa -db ${db}/usearch/Newest_UNITE/UNITE_modified_new2021.fasta \
-strand both -tabbedout result/otus1.sintax -sintax_cutoff 0.9

wc -l result/raw/otutab.txt

db=../../../../../Public_Scripts
Rscript ${db}/script/otutab_filter_nonBac.R -h
Rscript ${db}/script/otutab_filter_nonBac.R \
  --input result/raw/otutab.txt \
  --taxonomy result/raw/otus.sintax \
  --output result/otutab.txt\
  --stat result/raw/otutab_nonBac.stat \
  --discard result/raw/otus.sintax.discard

wc -l result/otutab.txt

usearch11 -otutab_trim result/otutab.txt \
    -min_otu_size 10 \
    -output result/otutab_trim.txt

cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
usearch11 -fastx_getseqs result/raw/otus.fa \
    -labels result/otutab.id -fastaout result/otus.fa
#
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
    result/raw/otus.sintax result/otutab.id \
    > result/otus.sintax

sed -i 's/\t$/\td:Unassigned/' result/otus.sintax

## Taxonomy summary

cut -f 1,4 result/otus.sintax \
  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
   > result/taxonomy2.txt
head -n3 result/taxonomy2.txt

awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
result/taxonomy2.txt > temp/otus.tax
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
result/taxonomy2.txt > temp/otus.tax
sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
> result/taxonomy.txt
head -n3 result/taxonomy.txt

cut -d ' ' -f 1-2 BCSC1/Network_Netshift/network_results/igraph_col.txt  > BCSC1/Network_Netshift/BCSC1_case.txt

less BCSC1/Network_Netshift/BCSC1_case.txt
#rank p c o f g，phylum, class, order, family, genus
mkdir -p result/tax
for i in p c o f g;do
  usearch11 -sintax_summary result/otus.sintax \
  -otutabin result/otutab_rare.txt -rank ${i} \
  -output result/tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
# 列出所有文件
ls -sh result/tax/sum_*.txt
head -n3 result/tax/sum_g.txt



usearch11 -otutab_stats result/otutab.txt \
  -output result/otutab.stat
cat result/otutab.stat


mkdir -p result/alpha
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
  --depth 50000 --seed 1 \
  --normalize result/otutab_rare.txt \
  --output result/alpha/vegan.txt
  
usearch11 -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat

###Calculate alpha diversity index

usearch11 -alpha_div result/otutab_rare.txt \
  -output result/alpha/alpha.txt


#Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
usearch10 -alpha_div_rare result/otutab_rare.txt \
  -output result/alpha/alpha_rare.txt -method without_replacement

##core OTUs
usearch11 -calc_distmx result/otus.fa -tabbedout result/distmx.txt 
  
usearch11 -otutab_core result/otutab_rare.txt -distmxin result/distmx.txt \
  -sintaxin result/otus.sintax -tabbedout result/core.txt
 
Rscript ${db}/script/otu_mean.R --input result/otutab_rare.txt \
   --design result/metadata.txt \
   --group group --thre 0 \
   --output result/otutab_mean.txt
head -n3 result/otutab_mean.txt


awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
    else {for(i=2;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
    result/otutab_mean.txt > result/alpha/otu_group_exist.txt
head result/alpha/otu_group_exist.txt

##########################################################

mkdir -p Sub_OTU/Shandong/result
grep "Shangdong" result/metadata.txt | cut -f 1 > Sub_OTU/Shandong/SampleID.txt

usearch11 -otutab_sample_subset result/otutab_rare.txt \
-labels Sub_OTU/Shandong/SampleID.txt -output Sub_OTU/Shandong/otutab.txt

usearch11 -otutab_trim Sub_OTU/Shandong/otutab.txt -min_otu_size 2 -output Sub_OTU/Shandong/result/otutab_rare.txt

usearch11 -otutab_stats Sub_OTU/Shandong/result/otutab_rare.txt \
   -output Sub_OTU/Shandong/result/otutab.stat
cat Sub_OTU/Shandong/result/otutab.stat


cut -f 1 Sub_OTU/Shandong/result/otutab_rare.txt | sed '1 s/#OTU ID/OTUID/' > Sub_OTU/Shandong/result/otu.id

usearch11 -fastx_getseqs result/otus.fa -labels Sub_OTU/Shandong/result/otu.id \
-fastaout Sub_OTU/Shandong/result/otus.fa


awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/otus.sintax \
Sub_OTU/Shandong/result/otu.id > Sub_OTU/Shandong/result/otus.sintax

awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/taxonomy.txt \
Sub_OTU/Shandong/result/otu.id > Sub_OTU/Shandong/result/taxonomy.txt


usearch11 -otutab_stats Sub_OTU/Shandong/result/otutab_rare.txt \
   -output Sub_OTU/Shandong/result/otutab.stat
cat Sub_OTU/Shandong/result/otutab.stat

#awk '{print $NF}' result/tax/test.txt
awk 'BEGIN{IFS='\t'}{$NF="";print $0}'  result/tax/sum_g.txt > result/tax/abundant_genus.txt

sed -i '2d' result/tax/abundant_genus.txt
cat result/tax/abundant_genus.txt | tr -s ' ' '\t' > result/tax/abundant_genus1.txt
cat result/tax/abundant_genus1.txt

mv result/tax/abundant_genus1.txt result/tax/abundant_genus.txt

less result/tax/abundant_genus.txt

cut -d ' ' -f 1-2 Network_Netshift/network_results/igraph_col.txt  > ../GH_case.txt

cut -d ' ' -f 1-2 Network_Netshift/network_results/igraph_col.txt  > ../NF_control.txt


cat ../GH_case.txt | tr -s ' ' '\t' > ../GH_case1.txt
cat ../NF_control.txt | tr -s ' ' '\t' > ../NF_control1.txt
cut -f 1-2,5 result/metadata.tsv > temp/group.txt
mkdir -p result_for_Netshift
cut -f 1-2 script/HLJNF_genus_edge.txt  > result_for_Netshift/HLJNF_Case.txt
sed -i '1d' result_for_Netshift/HLJNF_Case.txt
less result_for_Netshift/HLJNF_Case.txt
q

cd /mnt/sdd/ZhouXin/Tomato_project/DataProcess4Publication/Field_Bacteria_16S/Sub_OTU/Nature_field
mkdir -p result/OTU_tree
# Muscle
nohup muscle -in result/otus.fa -out result/OTU_tree/otus_aligned.fas &


trimal -in result/OTU_tree/otus_aligned.fas -out result/OTU_tree/otus_aligned_trimed.fa -gt 0.95

### IQ-TREE ML
nohup iqtree -s result/OTU_tree/otus_aligned_trimed.fa \
        -bb 1000 -redo -alrt 1000 -nt AUTO -T 40  \
        -pre result/OTU_tree/otus &

#Rscript ${db}/script/betaNTI_linux.R \
 --input result/otutab_rare.txt --tree result/OTU_tree/otus.contree \
  --output result/GMbeta_NTI
  
#
Rscript ${db}/script/RaupCrick_linux.R \
 --input result/otutab_rare.txt --output result/GMbeta_NTI
#
Rscript ${db}/script/summary_NTI_RC_linux.R \
 --bNTI result/bNTI.csv --RaupCrick result/GMbeta_NTI/table.csv \
 --output result/GMbeta_NTI


Rscript ${db}/script/alpha_boxplot.R -h

Rscript ${db}/script/alpha_boxplot.R --alpha_index richness \
  --input result/alpha/vegan.txt --design result/metadata.tsv \
  --group group --output result/alpha/ \
  --width 180 --height 120

for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
  Rscript ${db}/script/alpha_boxplot.R --alpha_index ${i} \
     --input result/alpha/vegan.txt --design result/metadata.tsv \
     --group group --output result/alpha/ \
    --width 180 --height 120
done
alpha_boxplot(alpha_div, index = "Shannon", metadata, groupID = group)
### 
Rscript ${db}/script/alpha_rare_curve.R \
  --input result/alpha/alpha_rare.txt --design result/metadata.tsv \
  --group group --output result/alpha/ \
  --width 180 --height 120



bash ${db}/script/sp_pheatmap.sh \
  -f result/beta/bray_curtis.txt \
  -H 'TRUE' -u 5 -v 5

cut -f 1-2,5 result/metadata.tsv > temp/group.txt

bash ${db}/script/sp_pheatmap.sh \
  -f result/beta/bray_curtis.txt \
  -H 'TRUE' -u 16 -v 12 \
  -P temp/group.txt -Q temp/group.txt

Rscript ${db}/script/beta_pcoa.R \
  --input result/beta/bray_curtis.txt --design result/metadata.csv \
  --group group --output result/beta/bray_curtis.txt.pcoa.pdf \
  --width 180 --height 120



## Taxonomy Stackplot

Rscript ${db}/script/tax_stackplot.R \
  --input result/tax/sum_p.txt --design result/metadata.tsv \
  --group Group --output result/tax/sum_p.stackplot \
  --legend 7 --width 180 --height 120

for i in p c o f g; do
Rscript ${db}/script/tax_stackplot.R \
  --input result/tax/sum_${i}.txt --design result/metadata.tsv \
  --group Group --output result/tax/sum_${i}.stackplot \
  --legend 13 --width 200 --height 150; done

i=p
Rscript ${db}/script/tax_circlize.R \
  --input result/tax/sum_${i}.txt --design result/metadata.tsv \
  --group Group1 --legend 5

mv circlize.pdf result/tax/sum_${i}.circlize.pdf
mv circlize_legend.pdf result/tax/sum_${i}.circlize_legend.pdf


Rscript ${db}/script/tax_maptree.R \
  --input result/otutab.txt --taxonomy result/taxonomy.txt \
  --output result/tax/tax_maptree.pdf \
  --topN 100 --width 183 --height 118

# Difference comparison
mkdir -p result/compare/

db=../../../../Public_Scripts
compare="HLJNF-HLJGH"
Rscript ${db}/script/compare.R \
  --input result/otutab.txt --design result/metadata.tsv \
  --group group --compare ${compare} --threshold 0.04 \
  --method edgeR --pvalue 0.05 --fdr 0.2 \
      --output result/compare/


Rscript ${db}/script/compare_volcano.R \
  --input result/compare/${compare}.txt \
  --output result/compare/${compare}.txt.volcano.pdf
  --width 89 --height 59

    bash ${db}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
       -d result/metadata.tsv -A group \
       -t result/taxonomy.txt \
       -w 8 -h 5 -s 7 \
       -o result/compare/${compare}.txt

    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_p.txt \
       -w 183 -v 59 -s 7 -l 10 \
       -o result/compare/${compare}

    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 183 -v 59 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.txt

cut -f 1 result/compare/SDNF-SDGH.txt > result/compare/SDID.txt 

awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' result/taxonomy.txt \
    result/compare/SDID.txt > result/compare/SD_genus.tax
cut -f 1,7 result/compare/SD_genus.tax > result/compare/SD_genus1.tax   

    Rscript ${db}/script/format2stamp.R -h
    mkdir -p result/stamp
    Rscript ${db}/script/format2stamp.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --threshold 0.1 \
      --output result/stamp/tax


Rscript ${db}/script/format2lefse.R -h
Rscript ${db}/script/format2lefse.R --input result/otutab.txt \
  --taxonomy result/taxonomy.txt --design result/metadata.tsv \
  --group Group1 --threshold 0.1 \
  --output result/LEfSe

########################################################

setwd("D:/DATA_PROCESS/Tomato_Project/Metagenomic_data/result")

db=../../../../../Public_Scripts


Rscript ${db}/script/format2lefse.R -h
Rscript ${db}/script/format2lefse.R --input resfam/CK_otutab.txt \
--taxonomy resfam/resfam_taxonomy.txt --design ./CKmetadata.txt \
--group group --threshold 0.05 \
--output resfam/CKLEfSe
##
Rscript ${db}/script/format2lefse.R --input dbcan2/dbcan2_otutab.txt \
--taxonomy dbcan2/dbcan2_taxonomy.txt --design ./CKmetadata.txt \
--group group --threshold 0.05 \
--output dbcan2/CKLEfSe


    bugbase=${db}/script/BugBase
    rm -rf result/bugbase/
    Rscript ${bugbase}/bin/run.bugbase.r -L ${bugbase} \
      -i result/gg/otutab.txt -m result/metadata.tsv -c Group1 -o result/bugbase/


mkdir -p result/tree
cd result/tree

tail -n+2 ../otutab_rare.txt | wc -l

usearch11 -otutab_trim ../otutab_rare.txt \
    -min_otu_freq 0.001 \
    -output otutab.txt

tail -n+2 otutab.txt | wc -l

cut -f 1 otutab.txt | sed '1 s/#OTU ID/OTUID/' > otutab_high.id


usearch11 -fastx_getseqs ../otus.fa -labels otutab_high.id \
    -fastaout otus.fa
head -n 2 otus.fa


awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../taxonomy.txt \
    otutab_high.id > otutab_high.tax

awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean.txt otutab_high.id \
    | sed 's/#OTU ID/OTUID/' > otutab_high.mean
head -n3 otutab_high.mean


cut -f 2- otutab_high.mean > temp
paste otutab_high.tax temp > annotation.txt
head -n 3 annotation.txt


time muscle -in otus.fa -out otus_aligned.fas

#trimAL
trimal -in otus_aligned.fas -out otus_aligned_trimed.fa -gt 0.95

mkdir -p iqtree
time ../../${db}/win/iqtree -s  otus_aligned_trimed.fa \
        -bb 1000 -redo -alrt 1000 -nt AUTO \
        -pre iqtree/otus

# cd ${wd}/result/tree
Rscript ../../${db}/script/table2itol.R -a -c double -D plan1 -i OTUID -l Genus -t %s -w 0.5 annotation.txt


## 
Rscript ../../${db}/script/table2itol.R -a -d -c none -D plan2 -b Phylum -i OTUID -l Genus -t %s -w 0.5 annotation.txt

## 
Rscript ../../${db}/script/table2itol.R -c keep -D plan3 -i OTUID -t %s otutab.txt

##
Rscript ${db}/script/table2itol.R -a -c factor -D plan4 -i OTUID -l Genus -t %s -w 0 annotation.txt

cd ${wd}

#MachineLearning
tsv-utils  transpose OTU_table.txt            \
| tsv-utils annotation                        \
       <(tsv-utils  add_headline              \
           "#OTU\tCatalog" metadata.txt)  -   \
 | sed 's/#OTU ID/Samples/'                   \
 > feature_table.txt

 usearch11  -forest_train     \
   feature_table.txt   -trees 100 -forestout forest.txt

grep -w "^varimp" forest.txt                  \
|cut -f3,4                                    \
|sort -rgk2                                   \
| tsv-utils annotation sintax.txt -           \
| tsv-utils add_headline                      \
     "#OTU\tImportant\tLineage" -             \
>rank.txt 

    mkdir -p ~/amplicon/24Compare/LEfSe
    cd ~/amplicon/24Compare/LEfSe

    lefse-format_input.py LEfSe.txt input.in -c 1 -o 1000000
    #lefse
    run_lefse.py input.in input.res
 
    lefse-plot_cladogram.py input.res cladogram.pdf --format pdf

    lefse-plot_res.py input.res res.pdf --format pdf

    head input.res 
    lefse-plot_features.py -f one --feature_name "Bacteria.Firmicutes.Bacilli.Bacillales.Planococcaceae.Paenisporosarcina" \
       --format pdf input.in input.res Bacilli.pdf

    mkdir -p features
    lefse-plot_features.py -f diff --archive none --format pdf \
      input.in input.res features/



cut -d ' ' -f 1-2 Network_Netshift/network_results/igraph_col.txt  > ../HLJGH_case.txt
#
cat ../HLJGH_case.txt | tr -s ' ' '\t' > ../HLJGH_case1.txt

cut -d ' ' -f 1-2 Network_Netshift/network_results/igraph_col.txt  > ../HLJNF_control.txt
cat ../HLJNF_control.txt | tr -s ' ' '\t' > ../HLJNF_control1.txt

cut -d ' ' -f 1-2 Network_Netshift/network_results/igraph_col.txt  > ../SDGH_case.txt
cat ../SDGH_case.txt | tr -s ' ' '\t' > ../SDGH_case1.txt

cut -d ' ' -f 1-2 Network_Netshift/network_results/igraph_col.txt  > ../SDNF_control.txt
cat ../SDNF_control.txt | tr -s ' ' '\t' > ../SDNF_control1.txt


cd ~/amplicon/result

biom convert -i gg/otutab.txt -o otutab_gg.biom --table-type="OTU table" --to-json
sed '1 s/^/#/' metadata.tsv > MappingFile.txt

