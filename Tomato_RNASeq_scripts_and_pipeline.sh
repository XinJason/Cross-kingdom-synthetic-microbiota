
##Author：#Citation: Zhou et al., Cross-kingdom synthetic microbiota supports tomato suppression of Fusarium wilt disease. Nat Commun 13, 7890 (2022). https://doi.org/10.1038/s41467-022-35452-6

# Download genome of Solanum lycopersicum in NCBI
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.3_SL2.50/GCF_000188115.3_SL2.50_genomic.fna.gz

soft=/data/zhouxin/Tomato_project/RNASeqCleanData/Tomato_transcriptome/soft

conda activate metagenome_env

cd ~/transcriptome/data

ls -sh *

#cd ~/transcriptome/data
#fastq-dump -v --split-3 --gzip SRR1039508
#rename "SRR1039508"  "untrt_N61311"  SRR1039508*

#/bin/rm ~/ncbi/public/sra/*.sra

cd /data/zhouxin/Tomato_project/RNASeqCleanData/Tomato_transcriptome/data

zcat AllT1_R1.fq.gz | head -n 8

echo "`zcat AllT1_R1.fq.gz | wc -l` / (4*1000000)" | bc -l

zcat AllT1_R1.fq.gz | awk '{if(FNR%4==0) base+=length}END{print base/10^9,"G";}'

fastqc AllT1_R1.fq.gz

fastqc *.fq.gz -t 20

multiqc -d . -o multiqc

cd data
trimmomatic PE -phred33 \
    trt_N061011_1.fq.gz trt_N061011_2.fq.gz \
    trt_N061011_1.trimmomatic.fq.gz trt_N061011_1.trimmomatic.unpaired.fq.gz \
    trt_N061011_2.trimmomatic.fq.gz trt_N061011_2.trimmomatic.unpaired.fq.gz \
    LEADING:20 TRAILING:20 MINLEN:36

trimmomatic.sh PE -phred33 \
    input_forward.fq input_reverse.fq \
    output_forward_paired.fq output_forward_unpaired.fq \
    output_reverse_paired.fq output_reverse_unpaired.fq \
    ILLUMINACLIP:adaptor-PE.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:36

#conda install -c bioconda ucsc-bedgraphtobigwig 
#conda install gtfToGenePred
#conda install -c bioconda ucsc-bedsort

cd ~/transcriptome/data/genome

gffread genome/Slycopersicum_514_ITAG3.2.gene.gff3 -T -o genome/Slycopersicum.gtf

#gtf2gff
#gffread merged.gtf -o- > merged.gff3

../../Tomato_transcriptome/soft/gtf2bed12.sh -f genome/Slycopersicum.gtf

gtfToGenePred -ignoreGroupsWithoutExons tomato.gtf \
    GRCh38.gtf.50505050.pred
genePredToBed GRCh38.gtf.50505050.pred GRCh38.gtf.bed12
rm -f GRCh38.gtf.50505050.pred
awk '$3-$2>1000 && $3-$2<2000' GRCh38.gtf.bed12 >GRCh38.model.gtf.bed12
head GRCh38.gtf.bed12

conda install -c anaconda mysql
# # cp ../bak/GRCh38.chromsize .

awk 'BEGIN{OFS="\t"}{if($0~/>/) {if(size>0) print chrname, size; size=0; chrname=$0; sub(">","", chrname);} else size+=length;}END{print chrname,size}' genome/Slycopersicum_514SL.fa
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg38.chromInfo" | tail -n +2 >GRCh38.chromsize
head GRCh38.chromsize

conda install -c bioconda hisat2 

hisat2_extract_splice_sites.py genome/Slycopersicum.gtf >genome/Slycopersicum.splice_site.bed
hisat2_extract_exons.py genome/Slycopersicum.gtf >genome/Slycopersicum.exon.bed

hisat2-build -p 10 --ss genome/Slycopersicum.splice_site.bed \
        --exon genome/Slycopersicum.exon.bed --seed 1 genome/Slycopersicum_514SL.fa Slycopersicum
ls Slycopersicum.*


#########################
# mkdir
mkdir -p Tomato_Hisat2

for i in `tail -n+2 ./metadata.tsv | cut -f 1`;do \
hisat2 -x genome/Slycopersicum -1 seq/${i}_R1.fq.gz -2 seq/${i}_R2.fq.gz \
    --known-splicesite-infile genome/Slycopersicum.splice_site.bed \
    --novel-splicesite-outfile Tomato_Hisat2/${i}.splicesiteNew.bed \
    -p 10 --mm -S Tomato_Hisat2/${i}.sam; done
	
nohup hisat2 -x genome/Slycopersicum -1 seq/FunT3_R1.fq.gz -2 seq/FunT3_R2.fq.gz \
    --known-splicesite-infile genome/Slycopersicum.splice_site.bed \
    --novel-splicesite-outfile Tomato_Hisat2/FunT3.splicesiteNew.bed \
    -p 6 --mm -S Tomato_Hisat2/FunT3.sam &


for i in `tail -n+2 ./metadata.tsv | cut -f 1`;do \
samtools view --threads 80 -F4 Tomato_Hisat2/${i}.sam -b \
        >Tomato_Hisat2/${i}.bam; done
		
		
samtools view --threads 20 -F4 Tomato_Hisat2/FunT3.sam -b \
        >Tomato_Hisat2/FunT3.bam 

for i in `tail -n+2 ./metadata.tsv | cut -f 1`;do \
samtools sort --threads 80 Tomato_Hisat2/${i}.bam \
        -o Tomato_Hisat2/${i}.sortP.bam; done

nohup samtools sort --threads 10 Tomato_Hisat2/FunT3.bam \
        -o Tomato_Hisat2/FunT3.sortP.bam & 
		
#######################
# built index
for i in `tail -n+2 ./metadata.tsv | cut -f 1`;do \
samtools index Tomato_Hisat2/${i}.sortP.bam; done

nohup bam2wig.py -i Tomato_Hisat2/AllT3.sortP.bam \
    -s genome/Slycopersicum.chromsize -o Tomato_Hisat2/AllT3 -t 1000000000 -q 0 &

ls -ltr Tomato_Hisat2/



############################################
# (salmon)

gffread tomato.gtf -g GCF_000188115.3_SL2.50_genomic.fna -w tomato.transcript.fa.tmp

cut -f 1 -d ' ' tomato.transcript.fa.tmp >tomato.transcript.fa
head tomato.transcript.fa

## 

salmon index -t genome/Slycopersicum_514_ITAG3.2.transcript.fa -i genome/Slycopersicum_Transcript_salmon

ls genome/Slycopersicum_Transcript_salmon

# duplicate_clusters.tsv  hash.bin  header.json  indexing.log  quasi_index.log  
# refInfo.json  rsd.bin  sa.bin  txpInfo.bin  versionInfo.json

head tomato.transcript.salmon/duplicate_clusters.tsv 

# RetainedTxp	DuplicateTxp
# ENST00000486475	ENST00000629881

cat <<END | sed 's/|/\t/' >sampleFile
Samp|conditions
BacT1|BAC
BacT2|BAC
BacT3|BAC
FunT1|FUN
FunT2|FUN
FunT3|FUN
AllT1|All
AllT2|All
AllT3|All
BlankT1|CK
BlankT2|CK
BlankT3|CK
END

cat -A sampleFile

 cd ~/transcriptome/data

# -g genome/GRCh38.gtf
for i in `tail -n+2 ./metadata.tsv | cut -f 1`
do
  salmon quant --gcBias -l A -1 seq/${i}_R1.fq.gz -2 seq/${i}_R2.fq.gz  -i genome/Slycopersicum_Transcript_salmon \
    -g genome/Slycopersicum.gtf -o ${i}/${i}.salmon.count -p 80 >${i}.salmon.log 2>&1
  # 
  cut -f 1,5 ${i}/${i}.salmon.count/quant.genes.sf | sed "s/NumReads/${i}/" >${i}/${i}.salmon.gene.count.tab
done
#############################################################

head BlankT3/BlankT3.salmon.gene.count.tab

find . -name quant.sf

tail -n+2 ./metadata.tsv | cut -f 1 | xargs -i echo -e "{}\t{}/{}.salmon.count/quant.sf" >salmon.output
head salmon.output
# Samp    Samp/Samp.salmon.count/quant.sf
# untrt_N61311    untrt_N61311/untrt_N61311.salmon.count/quant.sf
# untrt_N052611   untrt_N052611/untrt_N052611.salmon.count/quant.sf
######################################################################

sed 's/"/\t/g' genome/Slycopersicum.gtf | \
  awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) print "TXname\tGene"; if($3=="transcript") print $14, $10}' >genome/tomato.tx2gene
head genome/tomato.tx2gene
# TXname  Gene
# ENST00000608838 ENSG00000178591
# library("tximport")
# library("readr")
# salmon_file <- read.table("salmon.output", header=T,  row.names=1, sep="\t")
# tx2gene <- read.table("genome/GRCh38.tx2gene", header=T, row.names=NULL, sep="\t")
# txi <- tximport(salmon_file, type = "salmon",  tx2gene = tx2gene)
# dds <- DESeqDataSetFromTximport(txi,  sample,  ~conditions)


#cd ~/transcriptome/data
paste `find . -name *.salmon.gene.count.tab` | \
  awk 'BEGIN{OFS=FS="\t" }{line=$1; \
    for(i=2;i<=NF;i++) if(i%2==0) {if(FNR==1) count=$i; else count=int($i+0.5); line=line"\t"count;} print line;}' \
  >tomato_trans.Count_matrix.xls
head tomato_trans.Count_matrix.xls

#############################################################################
#STAR
gtfToGenePred -ignoreGroupsWithoutExons genome/Slycopersicum.gtf genome/Slycopersicum.gtf.pred

genePredToBed genome/Slycopersicum.gtf.pred genome/Solanum_lycopersicum.gtf.bed12

rm -f genome/Slycopersicum.gtf.pred

##../soft/gtf2bed12.sh -f genome/Solanum_lycopersicum.SL3.0.47.gtf
#
awk '$3-$2>1000 && $3-$2<2000' genome/Solanum_lycopersicum.gtf.bed12 >genome/Solanum_lycopersicum.model.gtf.bed12
head genome/Solanum_lycopersicum.model.gtf.bed12

# # cp ../bak/GRCh38.chromsize .

# awk 'BEGIN{OFS="\t"}{if($0~/>/) {if(size>0) print chrname, size; size=0; chrname=$0; sub(">","", chrname);} else size+=length;}END{print chrname,size}' GRCh38.fa
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg38.chromInfo" | tail -n +2 >GRCh38.chromsize
head genome/tomato_genome/GRCh38.chromsize


#conda install -c bioconda star 
mkdir -p genome/indexCn
nohup STAR --runMode genomeGenerate \
  --runThreadN 40 \
  --genomeDir genome/indexCn \
  --genomeFastaFiles genome/Slycopersicum_514SL.fa \
  --sjdbGTFfile genome/Slycopersicum.gtf \
  --sjdbOverhang 134 \
  --genomeSAindexNbases 11 &
mv Log.out STAR_index.log

# STAR gene
wc -l genome/indexCn/geneInfo.tab


ls -sh star_GRCh38

# 4.0K chrLength.txt      368K exonInfo.tab          1.5G SAindex
# 4.0K chrNameLength.txt   24K geneInfo.tab          204K sjdbInfo.txt
# 4.0K chrName.txt         64M Genome                204K sjdbList.fromGTF.out.tab
# 4.0K chrStart.txt       4.0K genomeParameters.txt  204K sjdbList.out.tab
# 732K exonGeTrInfo.tab   516M SA                    224K transcriptInfo.tab

wc -l star_GRCh38/geneInfo.tab

grep -cP '\tgene\t' GRCh38.gtf


mkdir -p trt_N061011

max_intron_size=5000

star_p=" --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
       --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
       --alignIntronMin 20 --alignIntronMax ${max_intron_size} \
                --alignMatesGapMax ${max_intron_size} \
                --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
                --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
                --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
                --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
                --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts"
          
nohup STAR --runMode alignReads --runThreadN 10 \
        --readFilesIn seq/FunT3_R1.fq.gz seq/FunT3_R2.fq.gz \
        --readFilesCommand zcat --genomeDir genome/indexCn \
        --outFileNamePrefix FunT3/FunT3. ${star_p} &


# Aug 03 11:44:27 ..... started STAR run
# Aug 03 11:44:27 ..... loading genome
# Aug 03 11:44:30 ..... started mapping
# Aug 03 11:44:48 ..... finished successfully

ls -sh trt_N061011/*

# using the --outFilterType Normal vs BySJout options will create slightly different SJ.out.tab counts for the following reason.
# Imagine that you have a read two junctions, with only 1st junction passing the filter.
# Then with the Normal option, the 1st junction will be counted in the SJ.out.tab.
# If the BySJout option is used, the entire read alignment may be prohibited because of the 2nd junction, and then the 1st junction will be not counted in the SJ.out.tab.
mkdir -p tmp
nohup samtools sort -@ 12 -T tmp/BlankT3\
  -o BlankT3/BlankT3.Aligned.sortedByCoord.out.bam \
  BlankT3/BlankT3.Aligned.out.bam &
  
for i in `tail -n+2 ./metadata.tsv | cut -f 1`
do
 samtools index ${i}/${i}.Aligned.sortedByCoord.out.bam &
done



for i in `tail -n+2 ./metadata.tsv | cut -f 1`
do
STAR --runMode inputAlignmentsFromBAM \
        --inputBAMfile ${i}/${i}.Aligned.sortedByCoord.out.bam \
        --outWigType bedGraph --outFileNamePrefix ${i}/${i}. \
        --outWigNorm RPM --outWigStrand Unstranded
bedSort ${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
        ${i}/${i}.Signal.UniqueMultiple.str1.out.bg
bedGraphToBigWig ${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
        genome/indexCn/chrNameLength.txt \
        ${i}/${i}.Signal.UniqueMultiple.str1.out.bw
done


multiqc -f -d . -o multiqc

geneBody_coverage2.py -i \
  trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bw \
  -r genome/GRCh38.model.gtf.bed12 -o trt_N061011/trt_N061011.geneBody_coverage

read_distribution.py -i trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
  -r genome/GRCh38.gtf.bed12 >trt_N061011/trt_N061011.read_distrib.xls

cat trt_N061011/trt_N061011.read_distrib.xls

RPKM_saturation.py -i \
  trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
  -r genome/GRCh38.gtf.bed12 -s 10 -q 0 -o trt_N061011/trt_N061011.RPKM_saturation

ls -ltr trt_N061011

max_intron_size=5000

star_p=" --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
       --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
       --alignIntronMin 20 --alignIntronMax ${max_intron_size} \
                --alignMatesGapMax ${max_intron_size} \
                --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
                --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
                --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
                --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
                --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts"

for i in `tail -n+2 ./metadata.tsv | cut -f 1`; do 
	mkdir -p ${i}
	mkdir -p tmp
	STAR --runMode alignReads --runThreadN 80 \
        --readFilesIn seq/${i}_R1.fq.gz seq/${i}_R2.fq.gz \
        --readFilesCommand zcat --genomeDir genome/indexCn \
        --outFileNamePrefix ${i}/${i}. ${star_p}
    samtools sort -@ 80 -T ${i}.tmp \
                -o ${i}/${i}.Aligned.sortedByCoord.out.bam \
                ${i}/${i}.Aligned.out.bam
	samtools index ${i}/${i}.Aligned.sortedByCoord.out.bam
	STAR --runMode inputAlignmentsFromBAM \
        --inputBAMfile ${i}/${i}.Aligned.sortedByCoord.out.bam \
        --outWigType bedGraph --outFileNamePrefix ${i}/${i}. \
        --outWigNorm RPM --outWigStrand Unstranded
    bedSort ${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
        ${i}/${i}.Signal.UniqueMultiple.str1.out.bg
    bedGraphToBigWig ${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
        genome/indexCn/chrNameLength.txt \
        ${i}/${i}.Signal.UniqueMultiple.str1.out.bw 
done &

ls -sh AllT1/*

#awk '$3-$2>1000 && $3-$2<2000' genome/Solanum_lycopersicum.gtf.bed12 >genome/Solanum_lycopersicum.model.gtf.bed12

nohup for i in `tail -n+2 ./metadata.tsv | cut -f 1`; do
  geneBody_coverage2.py -i \
    ${i}/${i}.Signal.UniqueMultiple.str1.out.bw \
    -r genome/Solanum_lycopersicum.gtf.bed12 -o ${i}/${i}.geneBody_coverage
	
  read_distribution.py -i ${i}/${i}.Aligned.sortedByCoord.out.bam \
    -r genome/Solanum_lycopersicum.gtf.bed12 >${i}/${i}.read_distrib.xls
	
  RPKM_saturation.py -i \
    ${i}/${i}.Aligned.sortedByCoord.out.bam \
    -r genome/Solanum_lycopersicum.gtf.bed12 -s 10 -q 0 -o ${i}/${i}.RPKM_saturation 
done &

   
# sed '5 i\Gene\ttrt_N061011\ttrt_N061011\ttrt_N061011' trt_N061011/trt_N061011.ReadsPerGene.out.tab trt_N061011/trt_N061011.ReadsPerGene.out.tab.ehbio

for i in `tail -n+2 ./metadata.tsv | cut -f 1`;  do 
  sed "5 i\Gene\t${i}\t${i}\t${i}" ${i}/${i}.ReadsPerGene.out.tab >${i}/${i}.ReadsPerGene.tomato
done

paste `find . -name *.ReadsPerGene.tomato` | tail -n +5 | \
  awk 'BEGIN{OFS=FS="\t" }{line=$1; \
    for(i=2;i<=NF;i++) if(i%2==0 && i%4!=0) line=line"\t"$i; print line;}' \
  >tomato_trans.Count_matrix.xls
head tomato_trans.Count_matrix.xls

for i in `tail -n +2 ./metadata.tsv | cut -f 1`; do echo -e "$i | $i\tSTAR_${i}"; done >star_map
for i in `tail -n +2 ./metadata.tsv | cut -f 1`; do echo -e "$i | $i.salmon.count | aux_info | ${i}.salmon.count\tSalmon_${i}"; done >salmon_map


nohup for i in `tail -n +2 metadata.tsv  | cut -f 1`; do 
	htseq-count -f bam -r pos -a 10 -t exon -s no -i gene_id -m union ${i}/${i}.Aligned.sortedByCoord.out.bam genome/Slycopersicum.gtf >${i}/${i}.readsCount 
	done &
	
	grep -v '^__' ${i}/${i}.readsCount | sed "1 iGene\t${i}" >${i}/${i}.readsCount2
done

grep -v '^__' ${i}/${i}.readsCount | sed "1 iGene\t${i}" >${i}/${i}.readsCount2

paste `find . -name *.readsCount2` | awk 'BEGIN{OFS=FS="\t" }{line=$1; for(i=2;i<=NF;i++) if(i%2==0) line=line"\t"$i; print line;}' >tomato_trans.Count_matrix.xls



DESeq2.sh -f tomato_trans.Count_matrix.xls -s sampleFile -P 0.01

R
source("https://bioconductor.org/biocLite.R")
BiocManager::install("AnnotationHub")
BiocManager::install("biomaRt")

library(AnnotationHub)
library(biomaRt)

hub <- AnnotationHub::AnnotationHub()

query(hub, "Solanum")  
Solanum.OrgDb <- hub[["AH59087"]]

saveDb(Solanum.OrgDb, file = "Solanum.orgdb")

loadDb(file = "Solanum.orgdb")

columns(Solanum.OrgDb)

ah <- AnnotationHub()
org <- ah[ah$rdataclass == "OrgDb",]
hm <- query(org, "Solanum")

hm
hm_org <- ah[["AH59087"]]
hm_org <- hm[[1]]
hm_org

DEG.gene_symbol = as.character(output.gene_id$gene_id) #获得基因 symbol ID
DEG.entrez_id = mapIds(x = Solanum.OrgDb,
                       keys = DEG.gene_symbol,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

head(keys(hm_org, keytype = "REFSEQ"))

keys(hm_org, keytype = "ENTREZID", column = "REFSEQ", pattern = "NM_001246833.2")

head(select(hm_org, keys = keys(hm_org), columns = c("ENTREZID", "GO"), keytype = "ENTREZID"))

gene.df <- bitr(gene, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = Solanum.OrgDb)

tomato <- bitr(x, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="Solanum.OrgDb")

source("https://bioconductor.org/biocLite.R")

BiocManager::install("clusterProfiler") 
BiocManager::install("topGO")  #GO
BiocManager::install("Rgraphviz")
BiocManager::install("pathview") #KEGG pathway
BiocManager::install("org.Hs.eg.db") 

library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(Solanum.OrgDb)


DEG.gene_symbol = as.character(output.gene_id$gene_id) 
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = DEG.gene_symbol,
                       keytype = "SYMBOL",
                       column = "ENTREZID")


sp_enrichmentPlot.sh -f ehbio.DESeq2.all.DE.entrez.untrt._higherThan_.trt.KEGG.xls \
	-o GeneRatio -T numeric -v Description -c qvalue -s Count -l qvalue 
	
sp_enrichmentPlot.sh -f ehbio.DESeq2.all.DE.entrez.all.KEGG.xls \
	-o GeneRatio -T numeric -v Description -c qvalue -s Count -l qvalue -S Group

awk 'BEGIN{OFS=FS="\t"}ARGIND==1{entrez[$1]=$3;}ARGIND==2{if(entrez[$1]!="") print entrez[$1],$2;}' genome/tomato_genome/tomato.idmap tomato_trans.Count_matrix.xls.DESeq2.all.DE >tomato.DESeq2.all.DE.entrez
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{symbol[$1]=$2;}ARGIND==2{if(symbol[$1]!="") print symbol[$1],$2;}' genome/tomato_genome/tomato.idmap tomato_trans.Count_matrix.xls.DESeq2.all.DE >tomato.DESeq2.all.DE.symbol


clusterProfileGO.sh -f ehbio.DESeq2.all.DE.entrez -r org.Hs.eg.db -s hsa
clusterProfileKEGG.sh -f ehbio.DESeq2.all.DE.entrez -s hsa -r org.Hs.eg.db

head -n 11 tomato_trans.Count_matrix.xls.DESeq2.All._higherThan_.BAC.xls >tomato_trans.Count_matrix.xls.DESeq2.All._higherThan_.BAC.top10.xls 
sp_enrichmentPlot.sh -f tomato_trans.Count_matrix.xls.DESeq2.All._higherThan_.BAC.top10.xls  -o GeneRatio -T numeric -v Description -c qvalue -s Count -l qvalue 

sp_enrichmentPlot.sh -f ehbio.DESeq2.all.DE.entrez.untrt._higherThan_.trt.KEGG.xls \
	-o GeneRatio -T numeric -v Description -c qvalue -s Count -l qvalue 

sp_enrichmentPlot.sh -f ehbio.DESeq2.all.DE.entrez.all.KEGG.xls \
	-o GeneRatio -T numeric -v Description -c qvalue -s Count -l qvalue -S Group


cut -f 1,5 ehbio_trans.Count_matrix.xls.DESeq2.untrt._vs_.trt.results.xls | tail -n +2 >ehbio.DESeq2.untrt._vs_.trt.gsea
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{entrez[$1]=$3;}ARGIND==2{if(entrez[$1]!="") print entrez[$1],$2;}' genome/GRCh38.idmap ehbio.DESeq2.untrt._vs_.trt.gsea >ehbio.DESeq2.untrt._vs_.trt.gsea.entrez

clusterProfileGOGSEA.sh -f ehbio.DESeq2.untrt._vs_.trt.gsea.entrez -r org.Hs.eg.db 



WGCNA.sh -f WGCNA/LiverFemaleClean.txt -t WGCNA/TraitsClean.txt


for i in `tail -n +2 sampleFile | cut -f 1`; do 
	htseq-count -f bam -r pos -a 10 -t exon -s no -i gene_id -m union ${i}/${i}.Aligned.sortedByCoord.out.bam ehbio_trans.gtf >${i}/${i}.readsCount
	grep -v '^__' ${i}/${i}.readsCount | sed " 1 iGene\t${i}"  >${i}/${i}.readsCount2
done

cd ~/data
gffcompare -R -r genome/GRCh38.gtf -o assembeCompare2Ref ehbio_trans.gtf


cut -f 3 assembeCompare2Ref.ehbio_trans.gtf.tmap | tail -n +2 | sort | uniq -c
