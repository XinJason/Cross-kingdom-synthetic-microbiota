###Analysis commonds for the fungal community based on PacBio sequencing
#set working directory
#setwd("D:/DATA_PROCESS/Tomato_Project/DataProcess4Publication/Field_Bacteria_16S")

setwd("D:/DATA_PROCESS/Tomato_Project/Field_Amplicon/16SV3V4/result/gg")
##Install and loading of R packages
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list <- c("reshape2","ggplot2","devtools","bindrcpp",
                  "ggthemes","agricolae","dplyr","stringr","ggsci","vegan","randomForest")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 安装bioconductor常用包
package_list <- c("digest","ggrepel")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


if (!requireNamespace("BiocManager", quietly = TRUE)) # TRUE or FALSE 均可。
  install.packages("BiocManager")

BiocManager::install(version = "3.11") # 指定3.11版本。
library("BiocManager")
install("edgeR") #此时使用的是BiocManager的命令，因此不能用install.packages(“edgeR”)。
# 安装Github常用包
# 参数解析、数据变换、绘图和开发包安装
package_list <- c("phyloseq")
for(p in package_list){
  q=unlist(strsplit(p,split = "/"))[2]
  if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install_github(p)
    suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

###############################
#alpha 多样性分析
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析参数-h显示帮助信息dominance, richness, simpson, shannon, ACE
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="alpha/vegan.txt",
                help="Input table file to read; Alpha [default %default]"),
    make_option(c("-t", "--type"), type="character", default="ACE",
                help="type of alpha index; [default %default]"),
    make_option(c("-d", "--design"), type="character", default="./metadata.csv",
                help="design file; [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=8,
                help="Width of figure; 宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=5,
                help="Height of figure; 高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 有统计表txt和矢量图pdf [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  # 如果没有设置输出，根据默认文件夹和类型参数设置输出位置
  if (opts$output==""){
    opts$output=paste("alpha/",opts$type, sep = "")}
  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("Type of alpha index is ", opts$type,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}
# 读取 alpha文件
alpha <- read.table("alpha/vegan.txt", header=T, row.names=1, sep="\t", comment.char="") 
# 读取实验设计
design <- read.table("./metadata.txt", header=T, row.names=1, sep="\t", comment.char="") 
# 提取样品组信息,默认为group可指定
sampFile = data.frame(group=design[,opts$group],
                      sample=row.names(design), 
                      row.names = row.names(design))
# 数据筛选，筛选两文件中共有
idx = rownames(sampFile) %in% rownames(alpha) # match design with alpha
sampFile = sampFile[idx,]
alpha = alpha[rownames(sampFile),] 
## richness index
# add design to alpha
#index = cbind(alpha[rownames(design),]$richness, sampFile) 
index = cbind(alpha[rownames(sampFile),][[opts$type]], sampFile) 
colnames(index) = c(opts$type,"group","sample") # add richness colname is value
# 统计各组间差异
#model = aov(richness ~ group, data=index)
model = aov(index[[opts$type]] ~ group, data=index)
# 计算Tukey显著性差异检验
Tukey_HSD <- TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
# 提取比较结果
Tukey_HSD_table <- as.data.frame(Tukey_HSD$group) 
Tukey_HSD_table
# LSD检验，添加差异组字母
out <- LSD.test(model,"group", p.adj="none")
stat = out$groups
# 分组结果添入Index
index$stat=stat[as.character(index$group),]$groups
# 设置分组位置为各组y最大值+高的3%
max=max(index[,c(opts$type)])
min=min(index[,opts$type])
x = index[,c("group",opts$type)]
# 下方的richness如何替换为变量
# y = x%>% group_by(group)%>% summarise(Max=max(richness))
y = x %>% group_by(group) %>% summarise_(Max=paste('max(',opts$type,')',sep=""))
y=as.data.frame(y)
rownames(y)=y$group
index$y=y[as.character(index$group),]$Max + (max-min)*0.05
##自己加的，调整分组顺序命令
index$group <- factor(index$group, levels=c("HLJNF", "HLJGH", "SDNF", "SDGH",ordered=TRUE))
#                      
p = ggplot(index, aes(x=group, y=index[[opts$type]], color=group)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste(opts$type, "index")) + theme_classic() +
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)
p= p + scale_color_npg()
p
# 5. 保存图表
# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, ".pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, ".png", sep=""), p, width = opts$width, height = opts$height)
# 保存一个制表符，解决存在行名时，列名无法对齐的问题
write.table("\t", file=paste(opts$output,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
write.table(Tukey_HSD_table, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)

######################################
# Bacterial Beta Diversity
########PCoA########
# 主坐标轴分析，可选距离矩阵bray_curtis、unifrac、unifrac_binary、jaccard、manhatten、euclidean
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析命令行
if (TRUE){
  option_list = list(
    make_option(c("-t", "--type"), type="character", default="manhatten",
                help="Distance type; 距离类型, 可选bray_curtis, bray_curtis_binary, euclidean, jaccard, jaccard_binary, manhatten, unifrac, unifrac_binary [default %default]"),   
    make_option(c("-i", "--input"), type="character", default="",
                help="Input beta distance; 距离矩阵,默认beta目录下与t同名，可指定 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="./metadata.tsv",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Province",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=6,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=4,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 有txt和矢量图pdf [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$input==""){opts$input=paste("beta/",opts$type, ".txt", sep = "")}
  if (opts$output==""){opts$output=paste("beta/pcoa_",opts$type, sep = "")}
  # 显示输入输出确认是否正确
  print(paste("The distance matrix file is ", opts$input,  sep = ""))
  print(paste("Type of distance type is ", opts$type,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}
# 3. 读取输入文件
# 读取距离矩阵文件
dis = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="") 
# 读取实验设计
#design <- read.csv(opts$design,header=T, row.names=1) 
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char="") 
# 提取样品组信息,默认为group可指定
# sampFile = as.data.frame(design[,opts$group],row.names = row.names(design))
# colnames(sampFile)[1] = "group" # 单列数据框筛选会变为list，改为双列
sampFile = data.frame(group=design[,opts$group],
                      sample=row.names(design), 
                      row.names = row.names(design))
# 数据筛选，筛选两文件中共有
idx = rownames(sampFile) %in% colnames(dis) # match design with alpha
sampFile = sampFile[idx,]
dis = dis[rownames(sampFile),rownames(sampFile)] 
# 4. 统计与绘图
# vegan:cmdscale计算矩阵矩阵中主坐标轴坐标，取前3维
pcoa = cmdscale(dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
##summary(points)
eig = pcoa$eig
points = cbind(points, sampFile[rownames(points),])
colnames(points) = c("x", "y", "z","group") 
#各 PCoA 轴的解释量
pcoa_exp <- pcoa$eig/sum(pcoa$eig)
##查看PCOA各轴解释率
head(pcoa_exp)
#或
site <- scores(pcoa)
##自己加的，调整分组顺序命令
design$group <- factor(design$group, levels=c("HLJNF", "HLJGH", "SDNF", "SDGH",ordered=TRUE))
# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y,  color=design$group)) + geom_point(alpha=.7, size=4) +
    ## 在aes中添加颜色命令aes(x=x, y=y, color=design$group, shape =design$Province) 
  #aes(x=x, y=y, color=design$group, shape =design$Year)+
  aes(x=x, y=y, color=design$group)+
  # # geom用于设置填充形状，alpha设置透明度。不设置则为实心填充，遮盖椭圆中的点, levels设置confidence ellipses的置信区间, 在0-1范围内。levels越小椭圆面积越小，涵盖的点越集中。
  # 不需要填充的时候去掉fill及对fill的补充参数geom,alpha等即可
stat_ellipse(aes(x = x, y = y, fill = design$group),type = "norm", geom = "polygon",color=NA, alpha = 0.2, levles = 0.9) + 
  # stat_ellipse(level = 0.9, aes(fill = group), type = "norm", geom = "polygon",alpha= 0.2,color=NA) + #添加椭圆形置信区间
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=paste(opts$type," PCoA",sep=""))  + 
  theme_classic() + scale_color_npg()
# # 颜色可以自己设置，或者直接用scale_color_brewer()
# scale_color_brewer(palette="Set3")
#scale_shape_manual(values = c(design$Group1))
p
# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, ".pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, ".png", sep=""), p, width = opts$width, height = opts$height)
print(paste(opts$output, ".pdf finished.", sep = ""))
# 添加样品标签
p=p+geom_text_repel(label=paste(rownames(points)),colour="black",size=3)
p
# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, "_label.pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, "_label.png", sep=""), p, width = opts$width, height = opts$height)
print(paste(opts$output, "_label.pdf finished.", sep = ""))
# Compare each group beta by vegan adonis in bray_curtis
da_adonis = function(sampleV){
  sampleA = as.matrix(sampleV$sampA)
  sampleB = as.matrix(sampleV$sampB)
  design2 = subset(sampFile, group %in% c(sampleA,sampleB))
  if (length(unique(design2$group))>1) {
    sub_dis_table = dis_table[rownames(design2),rownames(design2)]
    sub_dis_table = as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
    adonis_table = adonis(sub_dis_table~group, data=design2, permutations = 10000) 
    adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
    print(paste("In ",opts$type," pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
    adonis_pvalue = paste(opts$type, sampleA, sampleB, adonis_pvalue, sep="\t")
    write.table(adonis_pvalue, file=paste(opts$output, ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
  }
}

# loop for each group pair
dis_table = as.matrix(dis)
if (TRUE) {
  compare_data = as.vector(unique(design[[opts$group]]))
  len_compare_data = length(compare_data)
  for(i in 1:(len_compare_data-1)) {
    for(j in (i+1):len_compare_data) {
      tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
      print(tmp_compare)
      da_adonis(tmp_compare)
    }
  }
}else {
  compare_data = read.table("doc/compare.txt", sep="\t", check.names=F, quote='', comment.char="")
  colnames(compare_data) = c("sampA", "sampB")
  for(i in 1:dim(compare_data)[1]){da_adonis(compare_data[i,])}
}	 
# 5. 保存图表
# 提示工作完成
print(paste("Adnois statistics result in",opts$output, ".txt is finished.", sep = ""))

##计算群落之间的相似性
#若是已经提供好了距离矩阵，则直接使用现有的距离矩阵进行分析即可
anosim_result_dis <- anosim(dis, design$Province, permutations = 999) 
###查看结果
summary(anosim_result_dis)
anosim_result_dis
#anosim_result_otu <- anosim(otu, group$site, permutations = 999, distance = 'bray')        

pdf(file = "anosim.locations.pdf",width = 10)
plot(anosim_result_dis, col = c( 'red', 'green', 'blue', 'orange', 'purple'))
dev.off()

###############################
##网络性质柱形图绘制
colours <- c("#D55E00","#426600","#66CC99", "#009E73", "#0072B2","#990000","#FF0010","#2BCE48");



data<-read.table("BacteriaFungi_Network_character.txt",header=T,sep="\t", comment.char="")
data1 = data[ ,1:15]

data1$sample <- factor(data1$sample, levels=c("FunGH","FunNF",ordered=TRUE))
data1<-na.omit(data1)
p <- ggplot(data = data1, mapping = aes(x = sample, y =Centralization_closeness,fill=sample)) + geom_bar(stat = 'identity',position="dodge")+
  xlab('')+ylab('Number of centralization closeness')+theme_classic()+
scale_fill_manual(values=colours[1:6])
p= p+ scale_color_npg()
# 保存pdf和png格式方便查看和编辑
ggsave("Fun Centralization_closeness.pdf", p, width = 4, height = 2.5)



#随机森林结果比较
Genus <- read.table("./Genus_compare.txt", header=T, row.names=1, sep="\t", comment.char="") 
# 读取实验设计
design <- read.table("./design.txt", header=T, row.names=1, sep="\t", comment.char="") 

imp = read.table("family_imp_phylum.txt", header=T, row.names= 1, sep="\t") 
imp = tail(imp, n = optimal)
imp$Family = factor(rownames(imp), levels = rownames(imp))



DA_list = read.table("NF-GH.txt", header=T, row.names=1, sep="\t") # , comment.char=""
DA_list = DA_list[rownames(imp) ,1:5]
DA_list$level = ifelse(DA_list$logFC>0, "Enriched","Depleted")
otu_error_bar = merge(otu_error_bar, DA_list, by.x="Family", by.y = "row.names", all.x = T)

### b1. Plot japonica (TEJ) enriched family

TEJ = rownames(DA_list[DA_list$level == "Enriched",])
p = ggplot(otu_error_bar[otu_error_bar$Family %in% TEJ, ], aes(x=Family, y=value, fill=subspecies)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.5,                    # Width of the error bars
                position=position_dodge(.9)) + main_theme
p=p+theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave(paste("Family_barplot_TEJ",".pdf", sep=""), p, width=89 * 0.6, height=50 * 1.5, unit='mm')
p
##########################################################################################
##进行BetaNTI以及RCBray的可视化
###展示betaNTI值###
library(reshape2)
library(ggplot2)
library(stringr)
meta.data <- read.csv("metadata_samples.csv", row.names = 1)
bnti <- read.csv("beta_NTI/clbi_weighted_bNTI.csv", row.names = 1)
bnti.1 <- data.frame(sample1 = rep(colnames(bnti)[1],dim(bnti)[1]-1), 
                     sample2 = rownames(bnti)[2:dim(bnti)[1]],
                     bnti = bnti[2:dim(bnti)[1],1])
for(i in 2:(dim(bnti)[1]-1)){
  bnti.2 <- data.frame(sample1 = rep(colnames(bnti)[i],dim(bnti)[1]-i), 
                       sample2 = rownames(bnti)[(i+1):dim(bnti)[1]],
                       bnti = bnti[(i+1):dim(bnti)[1],i])
  bnti.1 <- rbind(bnti.1, bnti.2)
}
bnti.1$loc1 <- str_sub(bnti.1$sample1,1,nchar(bnti.1$sample1)-3)
bnti.1$loc2 <- str_sub(bnti.1$sample2,1,nchar(bnti.1$sample2)-3)
bnti.1$Group <- "All"
locs <- levels(factor(meta.data$group))
bnti.0 <- subset(bnti.1, loc1 == loc2)
for(i in locs){
  bnti.2 <- subset(bnti.0, loc1 == i)
  bnti.2$Group <- i
  print(c(i, dim(bnti.2)))
  bnti.1 <- rbind(bnti.1, bnti.2)
}
bnti.locs <- subset(bnti.1, Group != "All")
P2<- ggplot(bnti.locs, aes(x=Group, y=bnti,fill=Group)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  #  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  scale_y_continuous(limits = c(min(bnti.1$bnti)-0.2,max(bnti.1$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  #  scale_y_continuous(breaks = c(-4,-2,0,2,4)) +
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",size=12), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 14,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", colour="black",size=12),  #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black",size=14), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("betNTI")+xlab("Mangrove") #设置x轴和y轴的标题
P2
ggsave("beta_NTI/bnti_locations.violin.pdf", height = 5.5, width = 7)

#betaNTI与环境参数差异
env.dif <- as.matrix(t(sapply(1:nrow(bnti.1), function(x) meta.data[bnti.1$sample1[x],10:21]-meta.data[bnti.1$sample2[x],10:21])))
rownames(env.dif) <- rownames(bnti.1)
write.csv(env.dif,"beta_NTI/environment.differenc.csv", row.names = T)
env.dif <- read.csv("beta_NTI/environment.differenc.csv", row.names = 1)
bnti.env <- cbind(bnti.1, env.dif)
colnames(env.dif)

(sum.lm.MAT <- summary(lm(bnti.env$bnti~bnti.env$MAT)))
(cortest.res.MAT <- cor.test(bnti.env$bnti,bnti.env$MAT, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.MAT <- ggplot(data = bnti.env, aes(MAT, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "MAT", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$MAT)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.mat$estimate, 3), ", P= ",sum.lm.mat$coefficients[2,4], sep =""), size = 3)

(sum.lm.MAP <- summary(lm(bnti.env$bnti~bnti.env$MAP)))
(cortest.res.MAP <- cor.test(bnti.env$bnti,bnti.env$MAP, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.MAP <- ggplot(data = bnti.env, aes(MAP, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "MAP", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$MAP)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.MAP$estimate, 3), ", P= ",sum.lm.MAP$coefficients[2,4], sep =""), size = 3)

(sum.lm.gravel_proportion <- summary(lm(bnti.env$bnti~bnti.env$gravel_proportion)))
(cortest.res.gravel_proportion <- cor.test(bnti.env$bnti,bnti.env$gravel_proportion, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.gravel_proportion <- ggplot(data = bnti.env, aes(gravel_proportion, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "gravel_proportion", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$gravel_proportion)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.gravel_proportion$estimate, 3), ", P= ",sum.lm.gravel_proportion$coefficients[2,4], sep =""), size = 3)

(sum.lm.Salinity <- summary(lm(bnti.env$bnti~bnti.env$Salinity)))
(cortest.res.Salinity <- cor.test(bnti.env$bnti,bnti.env$Salinity, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.Salinity <- ggplot(data = bnti.env, aes(Salinity, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "Salinity", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$Salinity)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.Salinity$estimate, 3), ", P= ",sum.lm.Salinity$coefficients[2,4], sep =""), size = 3)

(sum.lm.pH <- summary(lm(bnti.env$bnti~bnti.env$pH)))
(cortest.res.pH <- cor.test(bnti.env$bnti,bnti.env$pH, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.pH <- ggplot(data = bnti.env, aes(pH, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "pH", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$pH)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.pH$estimate, 3), ", P= ",sum.lm.pH$coefficients[2,4], sep =""), size = 3)

(sum.lm.TC <- summary(lm(bnti.env$bnti~bnti.env$TC)))
(cortest.res.TC <- cor.test(bnti.env$bnti,bnti.env$TC, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.TC <- ggplot(data = bnti.env, aes(TC, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "TC", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$TC)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.TC$estimate, 3), ", P= ",sum.lm.TC$coefficients[2,4], sep =""), size = 3)

(sum.lm.TOC <- summary(lm(bnti.env$bnti~bnti.env$TOC)))
(cortest.res.TOC <- cor.test(bnti.env$bnti,bnti.env$TOC, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.TOC <- ggplot(data = bnti.env, aes(TOC, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "TOC", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$TOC)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.TOC$estimate, 3), ", P= ",sum.lm.TOC$coefficients[2,4], sep =""), size = 3)

(sum.lm.TN <- summary(lm(bnti.env$bnti~bnti.env$TN)))
(cortest.res.TN <- cor.test(bnti.env$bnti,bnti.env$TN, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.TN <- ggplot(data = bnti.env, aes(TN, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "TN", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$TN)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.TN$estimate, 3), ", P= ",sum.lm.TN$coefficients[2,4], sep =""), size = 3)

(sum.lm.N.NH4. <- summary(lm(bnti.env$bnti~bnti.env$N.NH4.)))
(cortest.res.N.NH4. <- cor.test(bnti.env$bnti,bnti.env$N.NH4., alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.N.NH4. <- ggplot(data = bnti.env, aes(N.NH4., bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "N.NH4.", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$N.NH4.)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.N.NH4.$estimate, 3), ", P= ",sum.lm.N.NH4.$coefficients[2,4], sep =""), size = 3)

(sum.lm.N.NO3. <- summary(lm(bnti.env$bnti~bnti.env$N.NO3.)))
(cortest.res.N.NO3. <- cor.test(bnti.env$bnti,bnti.env$N.NO3., alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.N.NO3. <- ggplot(data = bnti.env, aes(N.NO3., bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "N.NO3.", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$N.NO3.)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.N.NO3.$estimate, 3), ", P= ",sum.lm.N.NO3.$coefficients[2,4], sep =""), size = 3)

(sum.lm.TP <- summary(lm(bnti.env$bnti~bnti.env$TP)))
(cortest.res.TP <- cor.test(bnti.env$bnti,bnti.env$TP, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.TP <- ggplot(data = bnti.env, aes(TP, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "TP", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$TP)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.TP$estimate, 3), ", P= ",sum.lm.TP$coefficients[2,4], sep =""), size = 3)

(sum.lm.TS <- summary(lm(bnti.env$bnti~bnti.env$TS)))
(cortest.res.TS <- cor.test(bnti.env$bnti,bnti.env$TS, alternative = "greater", method = "pearson", conf.level = 0.95)) 
p.TS <- ggplot(data = bnti.env, aes(TS, bnti)) + 
  geom_hline(yintercept = c(-2,0,2), color = "gray", linetype = "longdash")+
  geom_point(alpha = 1/5) + geom_smooth(size = 0.5, method = "lm", se = T)+
  labs(title = NULL,x = "TS", y = "BNTI") + 
  scale_y_continuous(limits = c(min(bnti.env$bnti)-0.2,max(bnti.env$bnti)+0.3),breaks = c(-4,-2,0,2,4)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = max(bnti.env$TS)*0.75, y = max(bnti.env$bnti)*0.9, 
           label = paste( "R= ", round(cortest.res.TS$estimate, 3), ", P= ",sum.lm.TS$coefficients[2,4], sep =""), size = 3)

pdf("beta_NTI/bnti_env.dif.pdf", width = 16,height = 12)
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 4)))
vplayout= function(x, y)viewport(layout.pos.row = x, layout.pos.col = y)
print(p.MAT, vp=vplayout(1,1))
print(p.MAP, vp=vplayout(1,2))
print(p.gravel_proportion, vp=vplayout(1,3))
print(p.Salinity, vp=vplayout(1,4))
print(p.pH, vp=vplayout(2,1))
print(p.TC, vp=vplayout(2,2))
print(p.TOC, vp=vplayout(2,3))
print(p.TN, vp=vplayout(2,4))
print(p.N.NH4., vp=vplayout(3,1))
print(p.N.NO3., vp=vplayout(3,2))
print(p.TP, vp=vplayout(3,3))
print(p.TS, vp=vplayout(3,4))
dev.off()

###展示RCbray值###
library(reshape2)
library(ggplot2)
library(stringr)
meta.data <- read.csv("metadata_samples.csv", row.names = 1)
bnti <- read.csv("beta_NTI/RC_bray.csv", row.names = 1)
bnti.1 <- data.frame(sample1 = rep(colnames(bnti)[1],dim(bnti)[1]-1), 
                     sample2 = rownames(bnti)[2:dim(bnti)[1]],
                     bnti = bnti[2:dim(bnti)[1],1])
for(i in 2:(dim(bnti)[1]-1)){
  bnti.2 <- data.frame(sample1 = rep(colnames(bnti)[i],dim(bnti)[1]-i), 
                       sample2 = rownames(bnti)[(i+1):dim(bnti)[1]],
                       bnti = bnti[(i+1):dim(bnti)[1],i])
  bnti.1 <- rbind(bnti.1, bnti.2)
}
bnti.1$loc1 <- str_sub(bnti.1$sample1,1,nchar(bnti.1$sample1)-3)
bnti.1$loc2 <- str_sub(bnti.1$sample2,1,nchar(bnti.1$sample2)-3)
bnti.1$Group <- "All"
locs <- levels(factor(meta.data$group))
bnti.0 <- subset(bnti.1, loc1 == loc2)
for(i in locs){
  bnti.2 <- subset(bnti.0, loc1 == i)
  bnti.2$Group <- i
  print(c(i, dim(bnti.2)))
  bnti.1 <- rbind(bnti.1, bnti.2)
}
bnti.locs <- subset(bnti.1, Group != "All")
bnti.all <- subset(bnti.1, Group == "All")
P2<- ggplot(bnti.1, aes(x=Group, y=bnti,fill=Group)) + 
  #  geom_violin(trim=FALSE,color="white", position = "dodge") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.8,position=position_dodge(0.9))+ #绘制箱线图
  #  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  scale_y_continuous(limits = c(min(bnti.1$bnti)-0.1,max(bnti.1$bnti)+0.1),breaks = c(-1,-0.5,0,0.5,1)) +
  #  scale_y_continuous(breaks = c(-4,-2,0,2,4)) +
  geom_hline(yintercept = c(-0.95,0,0.95), color = "gray", linetype = "longdash")+
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",size=12), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=12,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 14,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", colour="black",size=12),  #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black",size=14), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("RCbray")+xlab("Mangrove") #设置x轴和y轴的标题
P2
ggsave("beta_NTI/RCbray_locations.violin.pdf", height = 5.5, width = 10)







#Fungal analysis
######################################################################################
setwd("D:/DATA_PROCESS/Tomato_Project/DataProcess4Publication/TomatoITS_Final")

#群落数据co-occurrence分析
###安装需要的包，默认不安装，没安装过的请取消如下注释
#BiocManager::install(c("igraph", "psych", "AnnotationDbi", "impute", "preprocessCore", "GO.db", "WGCNA", "multtest"))
#BiocManager::install("Hmisc")
###加载包
library(igraph)
library(WGCNA)
library(multtest)
library(psych)

###封装到自定义函数中
matrix2igraph<-function(matr,r.threshold,p.threshold){
  #数据量小时可以用psych包corr.test求相关性矩阵，数据量大时，可应用WGCNA中corAndPvalue, 但p值需要借助其他函数矫正
  #occor<-corAndPvalue(matr,method = c( "spearman"))
  # multiple test the p values
  #mtadj<-mt.rawp2adjp(unlist(occor$p),proc="BH")
  #adpcor<-mtadj$adjp[order(mtadj$index),2]
  #occor.p<-matrix(adpcor,dim(matr)[2])

  #用psych包corr.test求相关性矩阵
  occor = corr.test(matr,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
  occor.r = occor$r # 取相关性矩阵R值
  occor.p = occor$p # 取相关性矩阵p值
  # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold] = 0   
  # 构建igraph对象
  igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
  igraph
  # NOTE:可以设置weighted=NULL,但是此时要注意此函数只能识别相互作用矩阵内正整数，所以应用前请确保矩阵正确。
  # 可以按下面命令转换数据
  # occor.r[occor.r!=0] <- 1
  # igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)
  # 是否去掉孤立顶点，根据自己实验而定
  # remove isolated nodes，即去掉和所有otu均无相关性的otu 可省略，前期矩阵已处理过
  bad.vs <- V(igraph)[degree(igraph) == 0]
  igraph <- delete.vertices(igraph, bad.vs)
  igraph
}

###文件读取
otu <- read.csv("otumap.110samples.fungi.csv",row.names=1)
otu_tax <- read.csv("otumap.110samples.taxa.csv",row.names=1)[,1:7]
##otu筛选,可选步骤，如果otu数量太多，建议筛选一下
otu<-otu[rowSums(otu)/sum(otu)>=(0.01/100),]
##文件处理,由于otu经过标准化后减少，但tax未减少，需重新与匹配二者数据
otu_tax<-otu_tax[rownames(otu),]
#otu_tax<-strsplit(as.vector(otu_tax[,1]),split=";")
otu_phylum <- as.character(otu_tax$Phylum)
otu_abundance<-rowSums(otu)
otu_pro<-cbind(otu_phylum,otu_abundance)

#从OTU table数据获得network关系
igraph<-matrix2igraph(t(otu),0.6,0.05)
# 将igraph weight属性赋值到igraph.weight,用于后边做图
(igraph.weight <- E(igraph)$weight)
# 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
E(igraph)$weight <- NA
#igraph<-remove.edge.attribute(igraph,"weight")#把边值删除

#按相关类型设置边颜色
# 如果构建网络时，weighted=NULL,此步骤不能统计
sum(igraph.weight>0)# number of postive correlation
sum(igraph.weight<0)# number of negative correlation

# set edge color，postive correlation 设定为red, negative correlation设定为blue
E.color <- igraph.weight
E.color <- ifelse(E.color>0, "red",ifelse(E.color<0, "blue","grey"))
E(igraph)$color <- as.character(E.color)
E(igraph)$width = abs(igraph.weight)*4

# 添加OTU注释信息，如分类单元和丰度
# 另外可以设置vertices size, vertices color来表征更多维度的数据
#根据otu的分类地位、多度等性质设置顶点颜色、大小、线粗细等
# set vertices size
igraph.otu <- as.data.frame(otu_pro[V(igraph)$name,]) # 筛选对应OTU属性
igraph.size <- log(as.numeric(igraph.otu$otu_abundance),10)*2
V(igraph)$size <- igraph.size
# set vertices color
igraph.col<-igraph.otu$otu_phylum
colour.s <- c("#2BCE48","#990000","#740AFF","#FFFF00","#FF0010","#FFA405", "#0075DC","#9DCC00","#5EF1F2","#808080","#8F7C00")
colour.ss <- colour.s[1:length(levels(factor(igraph.col)))]
for(i in 1:length(levels(factor(igraph.col)))){
  igraph.col[igraph.col == levels(factor(igraph.col))[i]] <- colour.ss[i]
}
V(igraph)$color <- igraph.col

# 创建存放结果文件夹
#dir.create("network")
#pdf(file = "network/co-occurrence_network_phylum.pdf",width = 13)

set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,0.6, legend=levels(as.factor(igraph.otu$otu_phylum)),col=colour.ss, pch=16,cex=1,bty="n")
dev.off()

#更改layout:  component_wise(), layout_as_bipartite(), layout_as_star(), layout_as_tree(), layout_in_circle(), 
#layout_nicely(), layout_on_grid(), layout_on_sphere(), layout_randomly(), layout_with_dh(), layout_with_fr(), 
#layout_with_gem(), layout_with_graphopt(), layout_with_kk(), layout_with_lgl(), layout_with_mds(), 
#layout_with_sugiyama(), merge_coords(), norm_coords(), normalize()
#set.seed(123)
#plot(igraph,main="Co-occurrence network",layout=layout_nicely,vertex.frame.color=NA,vertex.label=NA,
#     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
###按模块着色
# 模块性 modularity
fc = cluster_fast_greedy(igraph,weights =NULL)# cluster_walktrap cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
modularity = modularity(igraph,membership(fc))
# 按照模块为节点配色
comps = membership(fc)
colbar = rainbow(max(comps))
V(igraph)$color = colbar[comps]
#pdf(file = "network/co-occurrence_network_modular.pdf",width = 13)
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
legend(1.2,0.6, legend=paste("module", levels(as.factor(comps)), sep=" "),col=colbar, pch=16,cex=1,bty="n")
#dev.off()

###对不同group绘制network###
###文件读取
otu <- read.csv("otumap.110samples.fungi.csv",row.names=1)
otu_tax <- read.csv("otumap.110samples.taxa.csv",row.names=1)[,1:7]
meta.data <- read.csv("metadata_samples.csv", row.names = 1)
otu <- otu[,rownames(meta.data)]
otu.t1 <- otu[,meta.data$Depth == "20-30cm"]
otu.t1 <- otu.t1[rowSums(otu.t1) > 0,]

##otu筛选,可选步骤，如果otu数量太多，建议筛选一下
otu.t1<-otu.t1[rowSums(otu.t1)/sum(otu.t1)>=(0.01/100),]
##文件处理,由于otu经过标准化后减少，但tax未减少，需重新与匹配二者数据
otu_tax<-otu_tax[rownames(otu.t1),]
#otu_tax<-strsplit(as.vector(otu_tax[,1]),split=";")
otu_phylum <- as.character(otu_tax$Phylum)
otu_abundance<-rowSums(otu.t1)
otu_pro<-cbind(otu_phylum,otu_abundance)

#输出为graphml文件，使用Gephi进行可视化
igraph<-matrix2igraph(t(otu.t1),0.6,0.05)
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
(igraph.weight <- E(igraph)$weight)
E(igraph)$correlation <- igraph.weight
E(igraph)$weight <- abs(igraph.weight)
E.cor <- igraph.weight
E.cor <- ifelse(E.cor>0, 1,ifelse(E.cor<0, -1, 0))
E(igraph)$cor <- as.character(E.cor)
#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
tax <- otu_tax[as.character(V(igraph)$name), ]
V(igraph)$kingdom <- tax$Kingdom
V(igraph)$phylum <- tax$Phylum
V(igraph)$class <- tax$Class
V(igraph)$order <- tax$Order
V(igraph)$family <- tax$Family
V(igraph)$genus <- tax$Genus
V(igraph)$species <- tax$Species
#根据丰度设置节点size
igraph.otu <- as.data.frame(otu_pro[V(igraph)$name,]) # 筛选对应OTU属性
igraph.size <- log(as.numeric(igraph.otu$otu_abundance),10)*2
V(igraph)$size <- igraph.size
igraph
plot(igraph)
##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
adj_matrix <- as.matrix(get.adjacency(igraph, attr = 'correlation'))
write.csv(data.frame(adj_matrix, check.names = FALSE), 'network/network.adj_matrix.notXMD.csv')

#边列表
edge <- data.frame(as_edgelist(igraph))    #igraph 的邻接列表转为边列表
edge_list <- data.frame(source = edge[[1]], target = edge[[2]], weight = E(igraph)$weight,
  correlation = E(igraph)$correlation, cor =  E(igraph)$cor)
head(edge_list)
write.csv(edge_list, 'network/network.edge_list.notXMD.csv', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(label = names(V(igraph)), kingdom = V(igraph)$kingdom, phylum = V(igraph)$phylum,
  class = V(igraph)$class, order = V(igraph)$order, family = V(igraph)$family,
  genus = V(igraph)$genus, species = V(igraph)$species, size = V(igraph)$size)
head(node_list)
write.csv(node_list, 'network/network.node_list.notXMD.csv', row.names = FALSE, quote = FALSE)
#write.graph(igraph, 'network/network.graphml', format = 'graphml')

#节点属性 node property
igraph.node.pro <- cbind(igraph.degree = igraph::degree(igraph),   #节点度
                         igraph.closeness = centralization.degree(igraph)$res,  ##节点度中心性
                         igraph.betweenness = centralization.betweenness(igraph)$res, ##节点介数中心性
                         igraph.cen.degree = centralization.closeness(igraph)$res) ##节点中心性
##写出节点性质到文件
write.csv(igraph.node.pro,"network/network.node_list.infor.csv")

nodeinfo <- read.csv("network/network.all.node_infor.csv",row.names = 1)
taxa.tab <- read.csv("otumap.110samples.taxa.csv", row.names = 1)
for(i in 1:dim(nodeinfo)[1]){
  nodeinfo$phylum[i] <- taxa.tab$Phylum[rownames(taxa.tab) == rownames(nodeinfo)[i]]
}
write.csv(nodeinfo,"network/network.all.node_infor.csv")


