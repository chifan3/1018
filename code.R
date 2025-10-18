rm(list = ls())
# setwd("")
##################
library(msigdbr)
all_gene_sets<-msigdbr(species = "Homo sapiens")
all_gene_sets<-as.data.frame(all_gene_sets)
Hallmarker_gene_set<-all_gene_sets[which(all_gene_sets$gs_cat==""),]
H <- c()

Hallmarker_gene_set<-Hallmarker_gene_set[which(Hallmarker_gene_set$gs_name %in% H),]
length(table(Hallmarker_gene_set$gs_name))

Hallmarker_gene_set

geneset_description<-unique(Hallmarker_gene_set$gs_name)#25
geneset_description
write.csv(geneset_description,file = "")

diff <- read_csv("")
diff <- diff[,c(1,3)]
exp <- read_csv("")
names(exp)[1] <- "Transcript_id"
exp <- merge(diff,exp,by = "Transcript_id")
exp <- exp[,-1]

library(dplyr)
library(tibble)
#去重
exp <- exp %>% 
  mutate(rowMean =rowMeans(.[grep("s", names(.))])) %>% #求出平均数  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(Gene_Symbol,.keep_all = T) %>% # gene_name留下第一个
  select(-rowMean)
  
gene_list<-exp$Gene_Symbol
gene_list<-gene_list[gene_list%in%Hallmarker_gene_set$gene_symbol]## 
genelist <- unique(gene_list)
sample_hallmarker_expre<-exp[exp$Gene_Symbol %in% gene_list,]
sample_hallmarker_expre <- na.omit(sample_hallmarker_expre)
rownames(sample_hallmarker_expre) <- sample_hallmarker_expre$Gene_Symbol
sample_hallmarker_expre <- sample_hallmarker_expre[,-1]
# sample_hallmarker_expre <- mutate_all(sample_hallmarker_expre, as.numeric)
sample_hallmarker_expre <- log2(sample_hallmarker_expre+1)
sample_hallmarker_expre <- sample_hallmarker_expre[which(rowSums(sample_hallmarker_expre)>0),]

Hallmarker_gene_set<-Hallmarker_gene_set[which(Hallmarker_gene_set$gene_symbol %in% rownames(sample_hallmarker_expre)),]
# 提取出gene_symbol, gs_description
Hallmarker_gene_set<-Hallmarker_gene_set[,c(3,4)]
Hallmarker_gene_set<-as.data.frame(Hallmarker_gene_set)

mat <- as.matrix(sample_hallmarker_expre)

gen<-Hallmarker_gene_set
head(gen)

module=levels(as.factor(gen$gs_name))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,1]==module[y])[,2]#提取术语相同的基因名
}

library("GSVA")
es <- gsva(mat, gs, method="ssgsea", verbose=FALSE, parallel.sz=2)

######## a, Overall distribution #####################
#standarization by z-score
library(matrixStats)
GSEA_sd<-rowSds(es)
es<-as.data.frame(es)
GSEA_mean<-rowMeans(es)
es$GSEA_mean<-GSEA_mean
es$GSEA_sd<-GSEA_sd
#head(es,2)
es[,-c(7,8)]<-(es[,-c(7,8)]-es$GSEA_mean)/es$GSEA_sd

#加载自己整理过后的的术语

final_enriched_scores<-as.data.frame(matrix(NA,nrow=25))
cluster <- read_excel("picture/cluster.xlsx")
cluster$clustering <- ifelse(cluster$cluster=="SAAS",1,2)

for (i in 1:2) {
  tmp_sample<-cluster[which(cluster$clustering==i),]$sample
  Mean_cluster_score<-rowMeans(as.matrix((es[,tmp_sample])))
  final_enriched_scores[,i]<-Mean_cluster_score
}
rownames(final_enriched_scores)<-rownames(es)
colnames(final_enriched_scores)<-c("SAAS","SAAD")

#final_enriched_scores<-final_enriched_scores[,c(1:2)]
bk<-c(seq(min(final_enriched_scores),0,length.out = 25),seq(0.00001,max(final_enriched_scores),length.out = 25))
###先z-score 后 scale
# pdf(file="./Hallmark_GSEA.pdf",width = 8,height = 10)
pheatmap::pheatmap(as.matrix(final_enriched_scores),cluster_cols = F,cluster_rows = T,color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk,main="",scale ="row")
# dev.off()
write.csv(es,file = "picture/result.csv")
pheatmap::pheatmap(es[,-c(7,8)],
         cluster_cols = F,
         cluster_rows = T,
         color=colorRampPalette(c("#2A37EE","#FDF5E6","#EE5731"))(100),
         show_colnames = T,#
         border_color = NA,scale = "row",
         show_rownames =T,
         height = 1,          # 
         width = 3,            # 
         aspect_ratio = 0.6,
         fontsize = 8,
         fontsize_row=8,
         fontsize_col=8
)
###############################################################
library(dplyr)
library(tibble)
########### 提取生存信息
dput(colnames(TCGA_clinical))

coxdata <- TCGA_clinical[,c("TCGA_id","OS", "OS.Status")]

data <- as.data.frame(t(expr_TCGA_TPM_tumor))
data <- data %>% 
  rownames_to_column("TCGA_id")

coxdata <- coxdata %>% 
  inner_join(data, by = "TCGA_id") %>% 
  column_to_rownames("TCGA_id")

########## colnames(coxdata) <- gsub("-","_",colnames(coxdata)) ##处理一下，不然会报错
genes <- colnames(coxdata)[-c(1:2)]

library(survival)
res <- data.frame()
for (i in 1:length(genes)) {
  #i=1
  print(i)
  surv =as.formula(paste('Surv(OS, OS.Status)~', genes[i]))
  x = coxph(surv, data = coxdata)
  x = summary(x)
  p.value=signif(x$wald["pvalue"], digits=2)
  HR =signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
  CI <- paste0("(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res[i,1] = genes[i]
  res[i,2] = HR
  res[i,3] = CI
  res[i,4] = p.value
}
names(res) <- c("ID","HR","95% CI","p.value")
save(res,file = "./output1/01exercise.univariate_Cox.Rdata")

res <- res %>% 
  filter(p.value < 0.05)
###################################################################
####
NMF_mat<-tumor_tpm[inter,]
NMF_mat<-NMF_mat[which(rowSums(NMF_mat)>0),]#数据的过滤

#res.multirun_extraH <- nmf(NMF_mat, rank=2:6, nrun=5)

mads<-apply(NMF_mat, 1, mad)
NMF_mat<-NMF_mat[rev(order(mads)),]
NMF_mat= sweep(NMF_mat,1, apply(NMF_mat,1,median,na.rm=T))
library(ConsensusClusterPlus)
############################### ###############################
results = ConsensusClusterPlus(as.matrix(NMF_mat),maxK=5,reps=500,pItem=0.8,pFeature=1, 
                               title="consensus_cluster",clusterAlg="km",distance="euclidean",plot="pdf")
icl = calcICL(results,title="consensus_cluster",plot="pdf")
icl[["clusterConsensus"]]
##################------------------------获取最佳k值-----------------------------
maxK =5
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK

for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i

# The optimal K
optK = Kvec[which.min(PAC)]
optK## [1] 2
################-------------------------------------------------------------------
k= 2
consensus_clustering<-icl[["itemConsensus"]]
consensus_clustering<-consensus_clustering[which(consensus_clustering$k==2),]
consensus_clustering<-consensus_clustering[which(consensus_clustering$itemConsensus>0.5),]

cluster_information<-data.frame(sample_name=consensus_clustering$item,clustering=consensus_clustering$cluster)
library(dplyr)
cluster_information <- cluster_information %>% 
   distinct(sample_name,.keep_all = T)  #保留唯一样本
  
Merged_data<-merge(TCGA_clinical_tumor,cluster_information,by.x = "TCGA_id",by.y="sample_name")
Merged_data<-Merged_data[which(!is.na(Merged_data$clustering)),]#选出非缺失值clustering的行
##################################分型热图-------------------------------------------------------------
Merged_data=Merged_data[order(Merged_data$clustering),]#按照聚类排序
Merged_data <- Merged_data[!duplicated(Merged_data[,1]), ]
exp<-tumor_tpm[,Merged_data$TCGA_id]

rownames(Merged_data) <- Merged_data[,1]
data <- as.data.frame(Merged_data[,18,drop = FALSE])#clustering
#drop = FALSE可以保留原数据结构
head(data)


######################################################

library(pheatmap)
annotation_col <- da[,c("clustering","TCGA_id","age_at_index","gender")]
                                
# colnames(Merged_data)
# annotation_col <- Merged_data[,c("clustering","TCGA_id","age_at_index","gender",
#                           "ajcc_pathologic_stage","ajcc_pathologic_t",
#                           "ajcc_pathologic_n","ajcc_pathologic_m")]
#                                  
colnames(annotation_col) <- c("Cluster","TCGA_id","Age","Gender")
                              # ,"Stage","T","N","M")#改列名
row.names(annotation_col) <- annotation_col$TCGA_id

annotation_col <- annotation_col[,c(1,3:4),drop=F]
# ############
#########
annotation_col$Age<- ifelse(annotation_col$Age > 65, "> 65", "<= 65")
# annotation_col$Stage <- gsub("Stage ", "", annotation_col$Stage)#使内容简单化
annotation_col$Cluster<- ifelse(annotation_col$Cluster == 1, "C1", "C2")
table(annotation_col$Cluster)

tRNAClustercolor <- c("#B93234","#399FA3") 
names(tRNAClustercolor) <- c("C1","C2") #类型颜色
table(annotation_col$Cluster)

Agecolor <- c("#2D29EF","#B5ABAB")
names(Agecolor) <- c("> 65","<= 65")

Gendercolor <- c("#FBEDB8","#B007FD")
names(Gendercolor) <- c("female","male") #类型颜色
ann_colors <- list("Cluster"=tRNAClustercolor,"Age"=Agecolor,"Gender"=Gendercolor)
library(tibble)
annotation_col <- annotation_col[order(annotation_col$Cluster), ]  data2 <- da[,c(1,1：30)]
data2 <-column_to_rownames(data2,var = "TCGA_id")
data3 <- t(data2)
data3 <- data3[,rownames(annotation_col)]#样本顺序
library(pheatmap)
pheatmap(
  data3,
  cluster_cols = F,
  cluster_rows = F,
  color=colorRampPalette(c("#111DA7","blue","white","#D72463","#760707"))(30),
  show_colnames = F,#显示样本名
  border_color = NA,scale = "row",
  show_rownames =T,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  height = 0.2,          # 增加热图的高度，例如设为12英寸
  width = 6,            # 可以调整宽度以保持合适的长宽比
  aspect_ratio = 0.6,
  fontsize = 8,
  fontsize_row=7.5,
  fontsize_col=6)


###################################---亚型的生存分析----------------------------
colnames(mydata)[3] <- "Cluster"
mydata$Cluster<- ifelse(mydata$Cluster == 1, "C1","C2")

library(survival)
fit <- survfit(Surv(OS,OS.Status) ~ Cluster,data = mydata)
#设置主题
library(survminer)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
library(ggsci)
#绘制KM曲线
ggsurvplot(fit,
           palette= c(pal_nejm()(2)),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = F,
           legend.labs=c("C1","C2"), 
           legend.title="Cluster",
           xlab="Time (years)",
           ylab='Survival probability',
           risk.table=TRUE,
           break.time.by = 2,
           risk.table.title="Number at risk",
           risk.table.height=.4,
           risk.table.y.text = FALSE,
           surv.median.line = "hv",
           ggtheme = mytheme)
##########################################################
library(tidyestimate)

scores <-filter_common_genes(df =tpm_exp, id = "hgnc_symbol", tidy = FALSE, tell_missing = T, find_alias = F)
#tidy：逻辑值，指示基因ID是作为rownames还是数据框的第一列。
#若基因ID为rownames，则设为FALSE；若在第一列，则设为TRUE
scores<-estimate_score(df = scores,is_affymetrix = F)#RNA-seq,illumina
Merged_data<-merge(TCGA_clinical_tumor,cluster_information,by.x = "TCGA_id",by.y="sample_name")
Merged_data<-Merged_data[which(!is.na(Merged_data$clustering)),]#选出非缺失值clustering的行

Merged_data=Merged_data[order(Merged_data$clustering),]#按照聚类排序
names(Merged_data)[1] <- "SampleID"
Merged_data<-merge(Merged_data,scores,all.x = T,all.y = T,by.x = "SampleID",by.y="sample")
Merged_data<-Merged_data[which(!is.na(Merged_data$clustering)),]
Merged_data$clustering<- ifelse(Merged_data$clustering == 1, "C1", "C2")
Merged_data=Merged_data[order(Merged_data$clustering),]
library(ggpubr)

# pdf(file = paste0("./immune cell/stromal score",".pdf"),width = 10,height=8)
A=ggboxplot(Merged_data, x = "clustering",y = "stromal",color = "clustering",
            ylab = "Stromal score", xlab = "Molecular subtype")+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16,size = 1, position=position_jitter(0.3),alpha=0.3)+
  stat_compare_means(label.y = 2000)+
  ggtitle(label="ESTIMATE:stromal score vs molecular subtype")
print(A)
#03.stromal.c1.c2.7.6
################
ggboxplot(Merged_data, 
             x = "clustering", y = "immune", 
             color = "clustering", 
             ylab = "Immune score", xlab = "Molecular subtype")+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16,size = 1, position=position_jitter(0.3),alpha=0.3)+
  stat_compare_means(label.y = 5000)+
  ggtitle(label="ESTIMATE:immune score vs molecular subtype")
print(A)
##################

#############------------------------------------------------------library(CIBERSORT)
#如果是RNA-Seq数据，使用FPKM和TPM都很合适
#读取包自带的LM22文件（免疫细胞特征基因文件）
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
data(LM22)
LM22[1:4,1:4]
#样本表达矩阵文件
load(file="D:/works/R/Subtype prediction/0407chol/output/002tpm_exp.not.log.Rdata")
exp[1:4,1:4]#[1] 19493    151
load(file = "D:/works/R/Subtype prediction/0407chol/output/002TCGA_clinical_tumor.Rdata")
tumor_tpm <- exp[,TCGA_clinical_tumor$entity_submitter_id]

tumor_matrix <- as.matrix(tumor_tpm)


results <- cibersort(sig_matrix = LM22, mixture_file = tumor_matrix,perm = 1000,QN = T)
save(results,file = "output1/03CIBERSORT.Rdata")
#------------------------------------------------------------------
cluster_information <- read.csv("output1/02cluster_information_trna.csv")
cluster_information$clustering<- ifelse(cluster_information$clustering == 1, "C1", "C2")

re <- as.data.frame(results[,-(23:25)])
re$SampleNames = rownames(re)

merge.re<-merge(re,cluster_information,by.x = "SampleNames",by.y="sample_name")
merge.re <- merge.re[,-24]

library(reshape2)
dat <- melt(merge.re, id = c('SampleNames','clustering'))

colnames(dat)[3:4] <- c("Cell_type","Proportion")
#----------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
theme <- theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=14), 
        axis.text.y = element_text(size = 14), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14), 
        axis.line = element_line(size = 1),
        plot.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent'), 
        legend.position = "top",
        legend.title = element_blank(),
        panel.border = element_blank(),    
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.background = element_blank())
p <- ggplot(dat,aes(x = Cell_type,y = Proportion,fill = factor(clustering))) +
  # geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape=NA, 
               width = 0.6, alpha = 0.8, 
               # show.legend = FALSE
  ) +
  scale_fill_manual(values = alpha(brewer.pal(8,"Set1")[1:2],0.6)) +
  scale_y_continuous(expand = c(0,0),limits = c(-0.05,1), breaks = seq(0,1,0.2),
                     labels = seq(0,1,0.2)) +
  theme +
  labs(y = 'Immune Cell Relative Proportion',fill = NULL)

###添加显著性检验 
Data_summary <- as.data.frame(compare_means(Proportion~clustering, dat, method = "wilcox.test", 
                                            paired = FALSE,group.by = "Cell_type"))    

stat.test <- Data_summary[,-c(2,9)]
stat.test$xmin=c(1:22)-0.2  
stat.test$xmax=c(1:22)+0.2 
stat.test$p.signif[stat.test$p.signif == "ns"] <- NA
p2 <- p + geom_signif(xmin = stat.test$xmin,xmax = stat.test$xmax,
                      annotations = stat.test$p.signif,margin_top = 0.00,
                      y_position = 0.9,size = 0.5,textsize = 7.5,
                      tip_length = 0)#横线两侧折线长度
p2
###################################################################

sample_hallmarker_expre<-tumor_tpm[gene_list,]

Hallmarker_gene_set<-Hallmarker_gene_set[which(Hallmarker_gene_set$gene_symbol %in% row.names(sample_hallmarker_expre)),]
# 提取出gene_symbol, gs_description
Hallmarker_gene_set<-Hallmarker_gene_set[,c(4,15)]
Hallmarker_gene_set<-as.data.frame(Hallmarker_gene_set)

mat<-sample_hallmarker_expre
gen<-Hallmarker_gene_set
head(gen)

module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]#提取术语相同的基因名
}
mat=as.matrix(mat)

library("GSVA")
es <- gsva(mat, gs, method="ssgsea", verbose=FALSE, parallel.sz=2)

######## a, Overall distribution #####################
#standarization by z-score
library(matrixStats)
GSEA_sd<-rowSds(es)
es<-as.data.frame(es)
GSEA_mean<-rowMeans(es)
es$GSEA_mean<-GSEA_mean
es$GSEA_sd<-GSEA_sd
#head(es,2)
es[,-c(131,132)]<-(es[,-c(131,132)]-es$GSEA_mean)/es$GSEA_sd

#加载自己整理过后的的术语
geneset_description <- read.csv("D:/works/R/Subtype prediction/0407chol/output/009geneset_description.csv",row.names = 2)
geneset_description<-geneset_description[rownames(es),]
merge.es <- merge(geneset_description,es,by="row.names")

rownames(merge.es)<-merge.es$function.
es <- merge.es[,-c(1:3)]
final_enriched_scores<-as.data.frame(matrix(NA,nrow=50))
cluster_information <- read.csv("output1/02cluster_information_trna.csv")
for (i in 1:2) {
  tmp_sample<-cluster_information[which(cluster_information$clustering==i),]$sample_name
  Mean_cluster_score<-rowMeans(as.matrix((es[,tmp_sample])))
  final_enriched_scores[,i]<-Mean_cluster_score
}
rownames(final_enriched_scores)<-rownames(es)
colnames(final_enriched_scores)<-c("C1","C2")

#final_enriched_scores<-final_enriched_scores[,c(1:2)]
bk<-c(seq(min(final_enriched_scores),0,length.out = 50),seq(0.00001,max(final_enriched_scores),length.out = 50))
###先z-score 后 scale
# pdf(file="./Hallmark_GSEA.pdf",width = 8,height = 10)
pheatmap::pheatmap(as.matrix(final_enriched_scores),cluster_cols = F,cluster_rows = T,color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk,main="pheatmap of hallmarker enrichment",scale ="row")
# dev.off()
###########################################################

##########---------------------------enrich-------------------------------------
library(clusterProfiler)
library(tidyverse)

library(org.Hs.eg.db)
entrezid = mapIds(x = org.Hs.eg.db,
                  keys = gene,
                  keytype = "SYMBOL",
                  column = "ENTREZID")
length(entrezid)#471
entrezid = na.omit(entrezid)
length(entrezid)#458
entre = data.frame(entrezid)
head(entre)


save(entre,file = "")
####################################################
######GO富集分析####################################
GO_enrich = enrichGO(gene = entre[,1],#待富集的基因列表
                         
                         OrgDb = org.Hs.eg.db,#指定物种的基因数据库
                         keyType = 'ENTREZID',#输入数据的类型
                         ont = 'ALL',#可以指定BP\MF\CC\ALL
                         pAdjustMethod = 'fdr',#指定p值校正方法
                         pvalueCutoff = 0.05,#指定p值阈值（指定1以输出全部）
                         qvalueCutoff = 0.05,#同上
                         readable = FALSE)

GO_en = data.frame(GO_enrich)
save(GO_en,file = "")


#####################GO/KEGG富集结果可视化###############


highcol="#0000FF"     
lowcol="red3" 
pvalueFilter=0.05      
qvalueFilter=0.05
showNum=6
GO=GO_en[(GO_en$pvalue<pvalueFilter),]

GO_BP_count = sum(GO$ONTOLOGY == "BP")
GO_CC_count = sum(GO$ONTOLOGY == "CC")
GO_MF_count = sum(GO$ONTOLOGY == "MF")
# 选择每个类别的结果
GO = GO[c(
  which(GO$ONTOLOGY == "BP")[1:min(shownum, GO_BP_count)],
  which(GO$ONTOLOGY == "CC")[1:min(shownum, GO_CC_count)],
  which(GO$ONTOLOGY == "MF")[1:min(shownum, GO_MF_count)]
), ]
# GO=GO[c(which(GO$ONTOLOGY=="BP")[1:shownum],which(GO$ONTOLOGY=="CC")[1:shownum],which(GO$ONTOLOGY=="MF")[1:shownum]),]

# rt=GO[,c(1,2,3,10,4,6,8)] 
 rt=GO[,c(2,3,4,11,5,7,9)]       

names(rt)=c("Ontology","ID","Term","Count","Ratio","pvalue","qvalue")

# rt$ID=gsub(":","",rt$ID)
# for (i in 1:nrow(rt)) {
#   rt[i,3]=paste0(rt[i,2],":",rt[i,3])
# }

split_b<-str_split(rt$Ratio,"/")
b<-sapply(split_b,"[",1)
c<-sapply(split_b,"[",2)
rt$Ratio=as.numeric(rt$Count)/as.numeric(c[1])

##################################################################
#尝试让term拐弯的操作
library(ggplot2)
library(ggtext)

rt$Term_wrapped <- stringr::str_wrap(rt$Term, width = 60) # 调整宽度到适合的字符数

labels=rt[order(rt$Ratio,decreasing =F),"Term_wrapped"]
rt$Term_wrapped = factor(rt$Term_wrapped,levels=labels)
# 对于Term标签的换行，仍需预先用str_wrap处理

# 更新aes映射
p <- ggplot(rt, aes(Ratio, Term_wrapped)) + 
  geom_point(aes(size = Count, color = pvalue))

###############################################################
# p = 
#   ggplot(rt,aes(Ratio, Term)) +
#   geom_point(aes(size=Count, color=pvalue))+
p1 = p +
  scale_colour_gradient(high=highcol, low = lowcol) + 
  # labs(color="pvalue",size="Count",x="Gene ratio",y="Term")+
  labs(color="pvalue",size="Count",x="Gene ratio",y="")+
  theme_bw()+
  theme(axis.text.x=element_text(color="black", size=16),axis.text.y=element_text(color="black", size=16))+
  scale_size_continuous(range=c(6,12))
  # facet_grid(Ontology~., scale = 'free_y', space = 'free_y')

  p1
# facet_grid( Ontology~. ,scales="free")
ggsave(".pdf", width=10, height=8) #按照ratio进行的排序   
################################柱状图###########
rt=GO[,c(1,2,3,10,4,6,8)]        
names(rt)=c("Ontology","ID","Term","Count","Ratio","pvalue","qvalue")
# 
# rt$ID=gsub(":","",rt$ID)
# for (i in 1:nrow(rt)) {
#   rt[i,3]=paste0(rt[i,2],":",rt[i,3])
# }
rt$Term_wrapped <- stringr::str_wrap(rt$Term, width = 60) # 调整宽度到适合的字符数

labels=rt[order(rt$Count,decreasing =F),"Term_wrapped"]
rt$Term_wrapped = factor(rt$Term_wrapped,levels=labels)
# 对于Term标签的换行，仍需预先用str_wrap处理


p=ggplot(data=rt)+geom_bar(aes(x=Term_wrapped, y=Count, fill=pvalue), stat='identity')+
  coord_flip() + scale_fill_gradient(high=highcol, low = lowcol) +    
  xlab("") + ylab("Gene count") +  
  theme_bw()+
  theme(axis.text.x=element_text(color="black", size=16),axis.text.y=element_text(color="black", size=16)) + 
  scale_y_continuous(expand=c(0.02, 0.02)) + scale_x_discrete(expand=c(0.02,0.02))
  # facet_grid(Ontology~., scale = 'free_y', space = 'free_y')
# facet_grid( Ontology~. ,scales="free")
print(p)
ggsave("", width=10, height=8)       
#按照count进行的排序

##-----------------首先进行单因素cox，再训练集和验证集的划分--------------------
rm(list = ls())
setwd("D:/works/R/Subtype prediction/0603sulfur")
options(stringsAsFactors = F)
# train <- read.csv("output/013train.csv")#105
# test <- read.csv("output/013test.csv")
########## 加载R包
library(glmnet)
library(survival)
library(dplyr)
library(tidyr)
library(tibble)
library(survminer)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(timeROC)
library(survivalROC)
library(caret)
library(rms)
library(foreign)
library(My.stepwise)
########## 读取数据,
load(file = "")
gene_use <- diff$gene_id
#gene_use(".")中基因与tumor——tpm(".","-")中的行名连接符号有些不一致
gene_use <- gsub("\\.","_",gene_use)
########## 提取gene_use的表达矩阵
rownames(tumor_tpm) <- gsub("\\.","_",rownames(tumor_tpm))
rownames(tumor_tpm) <- gsub("-","_",rownames(tumor_tpm))

tumor_tpm <- tumor_tpm[gene_use,]
########## 提取生存信息
dput(colnames(TCGA_clinical_tumor))
names(TCGA_clinical_tumor)[1] <- "TCGA_id"
coxdata <- TCGA_clinical_tumor[,c("TCGA_id","OS", "OS.Status")]

data <- as.data.frame(t(tumor_tpm))
data <- data %>% 
  rownames_to_column("TCGA_id")

coxdata <- coxdata %>% 
  inner_join(data, by = "TCGA_id") %>% 
  column_to_rownames("TCGA_id")#130,473

########## 下面开始选择基因和建模
res <- data.frame()
for (i in 1:length(genes)) {
  #i=1
  print(i)
  surv =as.formula(paste('Surv(OS, OS.Status)~', genes[i]))
  x = coxph(surv, data = coxdata)
  x = summary(x)
  p.value=signif(x$wald["pvalue"], digits=2)
  HR =signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
  CI <- paste0("(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res[i,1] = genes[i]
  res[i,2] = HR
  res[i,3] = CI
  res[i,4] = p.value
}
names(res) <- c("ID","HR","95% CI","p.value")
res <- res %>% 
  filter(p.value < 0.05)
######################################################################################
#------------------------------------------------------------------------------
data <- t(tumor_tpm[res$ID,])
set.seed()
index <-  sort(sample(nrow(data), nrow(data)*.8))
train <- data[index,]#
test <-  data[-index,]#
write.csv(train,file = "output1/08train.csv")
write.csv(test,file = "output1/08test.csv")

## 二、lasso回归
# ------------------------------------coxdata----------------------------
dput(colnames(TCGA_clinical_tumor))
names(TCGA_clinical_tumor)[1] <- "TCGA_id"
coxdata <- TCGA_clinical_tumor[,c("TCGA_id","OS", "OS.Status")]
library(dplyr)
library(tibble)
data <- as.data.frame(train)
data <- data %>% 
  rownames_to_column("TCGA_id")

coxdata <- coxdata %>% 
  inner_join(data, by = "TCGA_id") %>% 
  column_to_rownames("TCGA_id")#104，283

## 获取表达数据--------------------------------------------------------------------
x <- as.matrix(coxdata[,3:length(colnames(coxdata))])
## 获取生存数据
y <- as.matrix(Surv(coxdata[,1],coxdata[,2]))
#使用了Surv函数来构建一个生存对象，该对象接受两个参数，
#第一个参数是表示生存时间的向量，第二个参数是表示生存状态的向量
## Part1.初步探索
library(glmnet)
set.seed()
fit <- glmnet(x,y, family = "cox",alpha = 1,nlambda = 1000)

plot(fit)
plot(fit,xvar = "lambda",label = F)
plot(fit,xvar = "dev",label = T)

set.seed() 
cvfit <- cv.glmnet(x, 
                   y,
                   family = "cox",
                   nfolds = 10, #多少折检验，常用10折交叉检验
                   nlambda = 100,#λ数目，会影响图片竖杠数目
                   alpha=1)
plot(cvfit)
library(export)
graph2pdf(file="./output1/picture/08lasso.train2.6.6.pdf",width = 6,height = 6)

## Part3.根据拟合λ值，求系数
fitMin <- glmnet(x,y,family = "cox",lambda = cvfit$lambda.min,alpha = 1)
fitLse <- glmnet(x,y,family = "cox",lambda = cvfit$lambda.1se,alpha = 1)
cbind(fitLse$beta,fitMin$beta)#里面的beta等同于coef()的结果。

## Part4.提取预测变量，预测lasso回归的风险评分
lassoGene <- rownames(fitMin$beta)[as.numeric(fitMin$beta)!=0] #13,提取不为0的系数，对应的基因
# lassoGene <- rownames(fitLse$beta)[as.numeric(fitLse$beta)!=0]
# lasso.coef <- fitMin$beta@x  ##提取不为0的系数
# lasso_result <- data.frame(Gene = lassoGene,coef=lasso.coef) #lambda.min的结果
# TCGA_lasso.Riskscore <- predict(fit,newx = x,s = cvfit$lambda.min,type = "link") #计算lasso回归的风险评分
# colnames(TCGA_lasso.Riskscore) <- "lasso.Riskscore"
# TCGA_lasso.Riskscore <- as.data.frame(TCGA_lasso.Riskscore)


###-------- 
coxdata <- cbind(coxdata[,1:2],coxdata[,lassoGene])#104hang15lie
str(coxdata)
#--------------------------------------------------
cox_data <- as.formula(paste0("Surv(OS, OS.Status)~",paste0(lassoGene,sep = "",collapse = "+")))
cox_more <- coxph(cox_data,data=coxdata)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]
cox_formula <- as.formula(paste("Surv(OS, OS.Status)~",
                                paste(rownames(cox_table)[cox_table[,3]>0.05],
                                      collapse = "+")))
#筛选的基因拿来构建多因素模型
cox_more_2 <- coxph(cox_formula,data=coxdata)
#library(survminer) #载入所需R包
ggforest(cox_more_2, #直接用前面多因素cox回归分析的结果
         main = "Hazard ratio",
         cpositions = c(0.02,-0.15, 0.25), #前三列的位置，第二列是样品数，设了个负值，相当于隐藏了
         fontsize = 1, #字体大小
         refLabel = "reference", 
         noDigits = 2) 

summary_cox <- summary(cox_more_2)

multiCox_HR <- cbind(coef = summary_cox$coefficients[,"coef"],
                     HR=signif(summary_cox$conf.int[,"exp(coef)"],digits = 2), #保留2位置小数
                     HR.confint.lower = signif(summary_cox$conf.int[,"lower .95"],digits = 2),
                     HR.confint.upper = signif(summary_cox$conf.int[,"upper .95"],digits = 2),
                     pvalue = signif(summary_cox$coefficients[,"Pr(>|z|)"],digits = 2))

multiCox_HR <- as.data.frame(multiCox_HR)
multiCox_HR$CI <- paste0("(", multiCox_HR$HR.confint.lower, "-", multiCox_HR$HR.confint.upper, ")")
colnames(multiCox_HR)[6] <- "95% CI"
multiCox_HR <- multiCox_HR %>% 
  rownames_to_column("Gene")
#----------------------------------------
inter <- c("OS","OS.Status",multiCox_HR$Gene)
dt <- coxdata[,inter]#去掉没有经过ph检验的基因
#提取回归系数用于不同队列风险评分计算：
coef <- coef(cox_more_2)
#--------------------------------------------------------------------------------------
#####
cutoff <- median(risk_data$score)#[1] -0.01551213

risk_data$risk_group <- ifelse(risk_data$score > cutoff, "High risk","Low risk")
library(export)
graph2pdf(file="./output1/picture/08train.risk_group.OS6.6.pdf",width = 6,height = 6)
save(risk_data,file = "./output1/08train.risk_data.group.Rdata")

load("./output1/08train.risk_data.group.Rdata")

##########################################################################################
##--------------------------桑基图----------------------------------------------------

df <- to_lodes_form(data,
                    key = "x", value = "stratum", id = "alluvium",
                    axes = 1:3)
ggplot(df, aes(x = x,  fill=stratum,label=stratum,
               stratum = stratum, alluvium  = alluvium))+#数据
  geom_flow(width = 0.3,#连线宽度
            curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
            alpha = 0.7,#透明度
            color = 'white')+#间隔颜色
  geom_stratum(width = 0.28, alpha=0.7)+#图中方块的宽度
  geom_text(stat = 'stratum', size = 4, color = 'black')+
           
  scale_fill_brewer(palette = "Set3")+
  #自定义颜色
  theme_void()+
  theme(legend.position = 'none')+
  scale_x_discrete(
    limits =  c("tRNA metabolic cluster", "Gene cluster", "Risk", "Survival state")
  )

#######---------------------风险评分图-----------------------------------------------
library(ggplot2)
## 样品按score值排序
risk_df <- risk_df[order(risk_df$score),]
p1 <- ggplot(data = risk_df,aes(x=seq(1:104), y=score)) +
  geom_point(aes(fill=risk_group),pch = 21,color = "white", size = 3, stroke = 0.1)+
  scale_fill_manual(values = c("red","blue")) + # 自定义颜色映射
  ggtitle("Training set") +
  theme(plot.title = element_text(hjust = 0.5)) + # title居中
  xlab("") + 
  ylab("Risk Score") +
  theme_bw()+
  geom_vline(aes(xintercept=52.5),colour = "#BB0000",linetype = "dashed")
p1#08.risk.score.train1.6.3
### 3.生存状态散点图
p2 <- ggplot(data = risk_df) +
  geom_point(aes(x=seq(1:104), y=OS, color=OS.Status)) + 
  scale_color_manual(values = c("blue","red")) +
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5)) + # title居中
  xlab("Patients(increasing risk score)") + 
  ylab("Survival time (years)")+
  theme_bw()+
  guides(color = guide_legend(reverse = TRUE))#转换图例顺序
p2#08.risk.score.train2.6.3
### 4.组成图
library(cowplot)
plot_grid(p1,p2,ncol = 1, align = "h",labels = c("A","B"))
####------------------------gene风险热图------------------------------------------------------
rm(list = ls())
load(file = ".Rdata")
rownames(data) <- data$sample_name
heat=data[,6:14]%>%  #t(scale(t(exp[gene,])))
  scale()%>%  
  as.data.frame() %>%                 
  rownames_to_column("sample_name") %>% #行名 提取出来
  mutate(sample_name=factor(sample_name,levels = sample_name))%>%
  gather(key = gene,value = expression,-sample_name)%>%
  inner_join(.,data[,c("sample_name","risk_group")],by=c("sample_name"))%>%
  mutate(risk_group=factor(risk_group,levels = c("High risk","Low risk")))

library(ggh4x)
p_heat=ggplot(heat, aes(sample_name,gene, fill = expression)) +
  geom_tile() +
  scale_fill_gradientn(colors=c('#00FFD4','#015749','#160A0A','#670303','#F90101'))+
  theme( panel.grid = element_blank(),    
         panel.background = element_blank(),
         axis.text.x.bottom = element_blank(),
         axis.ticks.x = element_blank(),
         legend.position = "bottom")+
  xlab("")+ylab("")+
  theme(#panel.grid = element_blank(),
    axis.title = element_text(size=15,color='black'),
    axis.text = element_text(size=15,color='black'))+
  facet_nested(.~risk_group,drop=T,scale="free",space="free",
               strip =strip_nested(background_x =elem_list_rect(fill=c("red","blue")),
                                   
                                   text_x=element_text(size = 14),
                                   by_layer_x = F))


#绘制ROC曲线的方式很多种，这里使用timeROC绘制 1年，3年和5年的ROC曲线
library(timeROC)
library(survival)
with(data,
     ROC_riskscore <<- timeROC(T = OS,
                               delta = OS.Status,
                               marker = score,
                               cause = 1,
                               weighting = "cox",
                               times = c(1,3,5),
                               # times=quantile(data$OS,probs=seq(0.2,0.8,0.1)),
                               ROC = TRUE)
)
# quantile(data$OS,probs=seq(0.2,0.8,0.1))
plot(ROC_riskscore, time = 1, col = "red", add = F,title = "",lwd = 2)
plot(ROC_riskscore, time = 3, col = "blue", add = T,lwd = 2)
plot(ROC_riskscore, time = 5, col = "purple", add = T,lwd = 2)
legend(x=0.45,y=0.35,c(paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)),
                      paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)),
                      paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3))),
       col=c("red","blue","purple"),lty=1,lwd=2,cex = 1,bty = "n")
# text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
# text(0.5,0.15,paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
# text(0.5,0.1,paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3)))
##########################################
library(ggsignif)
library(ggplot2)
my_comparisons <- list(c("<= 65", "> 65"))

p <- ggplot(data1, aes(x=Age, y=score, fill=Age)) +#color=age,只有边框颜色，
  #改为fill=age,可以填充颜色
  geom_boxplot() + 
  # geom_jitter(shape = 21, width = 0.20,
  # size = 5, stroke = 2)+
  theme_minimal() +
  labs(title=paste(), x="Age", y="Risk score")+
  theme(axis.text = element_text (size = 18))+#调整坐标轴字体大小
  theme(axis.title.x=element_text(vjust=0, size=18,face = "plain"))+#调整xlab字体大小
  theme(axis.title.y=element_text(vjust=0, size=18,face = "plain"))+#调整ylab字体大小
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(axis.line = element_line(color = "black",linewidth = 0.8))+#调整坐标轴样式
  theme(legend.title = element_text(size=20),  # 调整图例标题的大小
        legend.text = element_text(size=20),
        legend.position = "top")+   # 调整图例标签的大小
  scale_fill_manual(values=c("#89ABE3FF","#F55665FF" ))
p
p <- p + coord_cartesian(ylim = c(NA, max(data1$score) + 1.5))# 使用coord_cartesian调整y轴范围
p + geom_signif(
  comparisons = my_comparisons,#指定比较对象
  test = "wilcox.test",#指定使用的检验方法
  textsize = 4,#指定标记中文字的大小
  y_position = 4,#指定标记在y轴上的坐标，按照前面指定比较对象的顺序
  color = "black" 
)
##########################------------data2,gender,risk.score---------------------------
my_comparisons <- list(c("female", "male"))

p <- ggplot(data2, aes(x=Gender, y=score, fill=Gender)) +
  geom_boxplot() + 
  # geom_jitter(shape = 21, width = 0.20,
  # size = 5, stroke = 2)+
  theme_minimal() +
  labs(title=paste(), x="Gender", y="Risk score")+
  theme(axis.text = element_text (size = 18))+#调整坐标轴字体大小
  theme(axis.title.x=element_text(vjust=0, size=18,face = "plain"))+#调整xlab字体大小
  theme(axis.title.y=element_text(vjust=0, size=18,face = "plain"))+#调整ylab字体大小
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(axis.line = element_line(color = "black",linewidth = 0.8))+#调整坐标轴样式
  theme(legend.title = element_text(size=20),  # 调整图例标题的大小
        legend.text = element_text(size=20),
        legend.position = "top")+   # 调整图例标签的大小
  scale_fill_manual(values=c("#89ABE3FF","#F55665FF" ))
p
p <- p + coord_cartesian(ylim = c(NA, max(data1$score) + 1.5))# 使用coord_cartesian调整y轴范围
p + geom_signif(
  comparisons = my_comparisons,#指定比较对象
  test = "wilcox.test",#指定使用的检验方法
  textsize = 4,#指定标记中文字的大小
  y_position = 4,#指定标记在y轴上的坐标，按照前面指定比较对象的顺序
  color = "black" 
)
#######################################################
##定义模型参数
# coxFile="output/.uniCox.txt"
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #可视化
  pdf(file=forestFile, width=6, height=4)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  #绘制左边的图
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=3)
  abline(v=1, col="black", lty=2, lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
  axis(1)
  dev.off()
}

###森林图函数
#定义独立预后分析函数
indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
  risk=read.table(riskFile,header=T, sep="\t", check.names=F, row.names=1)    #风险文件
  cli=read.table(cliFile,header=T, sep="\t", check.names=F, row.names=1)      #临床文件
  rt <- merge(risk,cli,by="id")
  
  rt <- rt[,c(1,2,3,13,15,16)]
  rt <- column_to_rownames(rt,var = "id")
  #单因素分析
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(OS, OS.Status) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  
  write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")#可视化
  
  #多因素分析
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("OS", "OS.Status", as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(OS,OS.Status) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
}

#调用函数进行独立预后分析
indep(riskFile=".txt",
      cliFile=".txt",
      uniOutFile=".txt",
      multiOutFile=".txt",
      uniForest=".pdf",
      multiForest=".pdf")
#################################################################
#####单基因相关性分析---目的基因与免疫细胞相关性
rm(list = ls()) 
##定义目的基因
mygene <- c("")  
#免疫细胞丰度数据
inputFile=".txt"
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#基因表达数据
load(file = "")
tpm_exp<-tumor_tpm[,rownames(risk_data)]
#合并
nc = t(rbind(data,tpm_exp[mygene,]))  ;#将你的目的基因匹配到表达矩阵---行名匹配--注意大小写
library(Hmisc)
m = rcorr(nc)$r[1:nrow(data),(ncol(nc)-length(mygene)+1):ncol(nc)]

##计算p值
p = rcorr(nc)$P[1:nrow(data),(ncol(nc)-length(mygene)+1):ncol(nc)]
head(p)

library(dplyr)
tmp <- matrix(case_when(as.vector(p) < 0.001 ~ "***",
                        as.vector(p) < 0.01 ~ "**",
                        as.vector(p) < 0.05 ~ "*",
                        TRUE ~ ""), nrow = nrow(p))

##绘制热图
library(pheatmap)
p1 <- pheatmap(m,
               display_numbers =tmp,
               angle_col =45,
               
               color = colorRampPalette(c("#008B8B", "#FFF2F2", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 30, 
               cellheight = 15,
               width = 20, 
               height=16,
               treeheight_col = 0,
               treeheight_row = 0)
p1#12.gene.cor.immunecell10.10

######################################################################
geneset <- read_excel("output1/16.cells.13.function.immune.xlsx", sheet = 3)
head(geneset)

genesets <- split(geneset$gene,geneset$`function`)
head(genesets)

ssgseaScore=gsva(data, genesets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE,
                 verbose   = F#是否提供每个步骤详细信息
                 # ,min.sz = 10
)

normalize = function(x){
  return((x - min(x)) / (2 * (max(x) - min(x))))
}#将数值缩放到0~0.5
ssgseaScore=normalize(ssgseaScore)

ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="output1/13ssGSEA.result13function.txt",sep="\t",quote=F,col.names=F)

##########################################################################
##############16cells.13function.riskgroup.boxplot########################
library(stringr)

inputFile=".txt"

r=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#风险评分的分组
load("")
p <- risk_data

annotation <- data.frame(p[,27,drop = FALSE])
table(annotation)
head(annotation)

texp <- t(r)
data1 <- merge(annotation,texp,by = "row.names")

rownames(data1) <- data1[,1]
data1 <- data1[,-1]
table(data1$group)
####
dat=data1[order(data1$group),]

ann_colors=list()

crgCluCol <- c("#FF4500","#009ACD")
names(crgCluCol)=levels(factor(dat$group))
ann_colors[["group"]]=crgCluCol

head(dat, 10)

library(reshape2)
da=melt(dat, id.vars=c("group"))
colnames(da)=c("cluster", "cells", "value")
#
library(ggpubr)
p=ggboxplot(da, x="cells", y="value", color = "cluster", 
            ylab="Score",
            xlab="",
            legend.title="Risk",
            palette = crgCluCol,
            width=0.6, add = "none")
p=p+rotate_x_text(60)
#标显著星号
  p+stat_compare_means(aes(group=cluster),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                         symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

##############################################################
#############################################################
####################肿瘤突变负荷tmb与风险评分的相关性分析###################
####maf矩阵提取

#这个函数有两个参数，
#@metadata是从TCGA数据下载的sample sheet
#@path是保存maf文件的路径

merge_maf <- function(metadata, path){
  #通过合并path,还有sample sheet前两列得到每一个文件的完整路径
  filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                         fsep = .Platform$file.sep)
  
  message ('############### Merging maf data ################\n',
           '### This step may take a few minutes ###\n')
  #通过lapply循环去读每一个样本的maf，然后通过rbind合并成矩阵，按行来合并
  #colClasses指定所有列为字符串
  mafMatrix <- do.call("rbind", lapply(filenames, function(fl) 
    read.table(gzfile(fl),header=T,sep="\t",quote="",fill=T,colClasses="character")))
  return (mafMatrix)
}

#定义去除重复样本的函数FilterDuplicate
FilterDuplicate <- function(metadata) {
  filter <- which(duplicated(metadata[,'sample']))
  if (length(filter) != 0) {
    metadata <- metadata[-filter,]
  }
  message (paste('Removed', length(filter), 'samples', sep=' '))
  return (metadata)
}

#读入maf的sample sheet文件
metaMatrix.maf=read.table("data/maf_sample_sheet.2024-07-06.tsv",sep="\t",header=T)
#替换.为下划线，转换成小写，sample_id替换成sample
names(metaMatrix.maf)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.maf))))
#删掉最后一列sample_type中的空格
metaMatrix.maf$sample_type=gsub(" ","",metaMatrix.maf$sample_type)

#删掉重复的样本
metaMatrix.maf <- FilterDuplicate(metaMatrix.maf)#Removed 4 samples,149

#调用merge_maf函数合并maf的矩阵
maf_value=merge_maf(metadata=metaMatrix.maf, 
                    path="data/maf_20240706_data"
)
#查看前三行前十列
maf_value[1:3,1:10]


#保存合并后的maf文件
write.table(file="output1/14combined_maf_value.txt",maf_value,row.names=F,quote=F,sep="\t")
######################################################################
###############################计算TMB值##############################
#BiocManager::install("maftools")
rm(list = ls())
library(maftools)

laml <- read.maf(maf = "output1/14combined_maf_value.txt")

#????tmbֵ
tmb_table_wt_log = tmb(maf = laml)
write.table(tmb_table_wt_log,file="output1/14TMB_log.txt",sep="\t",row.names=F)
############################突变负荷分析###############################
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

tmbFile="output1/14TMB_log.txt"         #????ͻ???????ļ?

#样本的风险评分数据
load("~/0603/data/09whole.risk_data.Rdata")
data <- risk_data[,26,drop=F]
data=avereps(data)
#保留样本前12位id
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\", rownames(data))

#tmb数据
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F)
tmb$Tumor_Sample_Barcode=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",tmb$Tumor_Sample_Barcode)
library(dplyr)
library(tibble)
tmb <- tmb %>% 
  arrange(desc(total_perMB)) %>% #把表达量的平均值按从大到小排序
  distinct(Tumor_Sample_Barcode,.keep_all = T) %>% # 样本留下第一个
  column_to_rownames("Tumor_Sample_Barcode")
tmb=avereps(tmb)
#相同的样本
sameSample=intersect(row.names(data), row.names(tmb))
data=data[sameSample,,drop=F]
tmb=tmb[sameSample,,drop=F]
rt=cbind(data, tmb)

#?????Է???
x=as.numeric(rt[,"score"])
y=log2(as.numeric(rt[,"total_perMB"])+1)
df1=as.data.frame(cbind(x,y))
corT=cor.test(x, y, method="spearman")
p1=ggplot(df1, aes(x, y)) + 
  xlab(paste0("Risk score"))+ylab("Tumor mutation burden")+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
p2
#??????????ͼ??
pdf(file="picture/14cor.tmb.riskscore5.5.pdf",width=5,height=5)
print(p2)
dev.off()
###############################################################################
#########################tmb.riskgroup.boxplot#################################
rm(list = ls())
# setwd("D:/works/R/Subtype prediction/0603sulfur")

load("~/0603/data/09whole.risk_data.Rdata")
data <- risk_data[,27,drop=F]#risk.group
#保留样本前12位id
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\", rownames(data))

##tmb数据
tmbFile="output1/14TMB_log.txt" 
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F)
tmb$Tumor_Sample_Barcode=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",tmb$Tumor_Sample_Barcode)
library(dplyr)
library(tibble)
tmb <- tmb %>% 
  arrange(desc(total_perMB)) %>% #把表达量的平均值按从大到小排序
  distinct(Tumor_Sample_Barcode,.keep_all = T) %>% # 样本留下第一个
  column_to_rownames("Tumor_Sample_Barcode")
tmb=avereps(tmb)
#相同的样本
sameSample=intersect(row.names(data), row.names(tmb))#86
data=data[sameSample,,drop=F]
tmb=tmb[sameSample,,drop=F]
rt=cbind(data, tmb)
rt$tmb <- log2(as.numeric(rt[,"total_perMB"])+1)
###################
data1 <- rt[,c(1,5)]#score,tmb
names(data1) <- c("Cluster","tmb")
data1$Cluster <-factor(data1$Cluster,levels = unique(data1$Cluster))
library(ggsignif)
library(ggplot2)
my_comparisons <- list(c("High risk", "Low risk"))
# data1 <- data1[data1$tmb != 0, ]#因为绘图时的警告，去掉tmb为0的行,还是警告，重复值太多。

p <- ggplot(data1, aes(x=Cluster, y=tmb, fill=Cluster)) +#color=age,只有边框颜色，
  #改为fill=,可以填充颜色
  geom_boxplot() + 
  # geom_jitter(shape = 21, width = 0.20,
  # size = 5, stroke = 2)+
  theme_minimal() +
  labs(title=paste(), x="Risk", y="Tumor mutational burden")+
  theme(axis.text = element_text (size = 18))+#调整坐标轴字体大小
  theme(axis.title.x=element_text(vjust=0, size=18,face = "plain"))+#调整xlab字体大小
  theme(axis.title.y=element_text(vjust=0, size=18,face = "plain"))+#调整ylab字体大小
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(axis.line = element_line(color = "black",linewidth = 0.8))+#调整坐标轴样式
  theme(legend.title = element_text(size=20),  # 调整图例标题的大小
        legend.text = element_text(size=20),
        legend.position = "top")+   # 调整图例标签的大小
  scale_fill_manual(values=c("red","blue" ))
p
p <- p + coord_cartesian(ylim = c(NA, max(data1$tmb) + 1.5))# 使用coord_cartesian调整y轴范围
p + geom_signif(
  comparisons = my_comparisons,#指定比较对象
  test = "wilcox.test",#指定使用的检验方法
  textsize = 4,#指定标记中文字的大小
  y_position = 4,#指定标记在y轴上的坐标，按照前面指定比较对象的顺序
  color = "black" 
)#14.riskgroup.tmb.boxplot6.6.pdf

##############################################################################
#####################肿瘤突变负荷与高低风险分组的生存分析#####################
library(limma)
library(survival)
##tmb数据
tmbFile="output1/14TMB_log.txt" 
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F)
tmb$Tumor_Sample_Barcode=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",tmb$Tumor_Sample_Barcode)
library(dplyr)
library(tibble)
tmb <- tmb %>% 
  arrange(desc(total_perMB)) %>% #把表达量的平均值按从大到小排序
  distinct(Tumor_Sample_Barcode,.keep_all = T) %>% # 样本留下第一个
  column_to_rownames("Tumor_Sample_Barcode")
tmb=avereps(tmb)
#相同的样本
sameSample=intersect(row.names(data), row.names(tmb))data=data[sameSample,,drop=F]
tmb=tmb[sameSample,,drop=F]
rt=cbind(data, tmb)
rt$TMB <- log2(as.numeric(rt[,"total_perMB"])+1)
############################
rt <- rt[,c(1,2,3,7)] 
cutoff <- median(rt$TMB)#[1] 0.1375035

rt$Cluster2 <- ifelse(rt$TMB > cutoff, "H-TMB","L-TMB")
table(rt$Cluster2)
# H-TMB L-TMB
# 42    44

fit <- survfit(Surv(OS, OS.Status) ~ Cluster2, data = rt)
library(survminer)
p1  <-  ggsurvplot(fit,
                   legend.title = "",#定义图例的名称
                   legend.labs = c("H-TMB","L-TMB"),
                   #legend = "top",#图例位置
                   pval = T,
                   pval.method = F,#添加p值的检验方法
                   #conf.int = TRUE,#添加置信区间
                   risk.table = TRUE, #在图下方添加风险表
                   risk.table.col = "strata", #根据数据分组为风险表添加颜色
                   risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
                   #linetype = "strata", #改变不同组别的生存曲线的线型
                   surv.median.line = "hv", #标注出中位生存时间
                   xlab = "Time in years", #x轴标题
                   xlim = c(0,max(risk_data$OS)), #展示x轴的范围
                   ylab = "Overall survival rate",
                   break.time.by = 1, #x轴间隔
                   size = 1, #线条大小
                   #ggtheme = theme_bw(), #为图形添加网格
                   palette = c("red","blue")#图形颜色风格
)
#######################细化分组进行生存分析-riskgroup.tmbgroup#######################
library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)
#
rt$Cluster <- paste(rt$Cluster2, "+", rt$group)

# 生存分析
fitd <- survdiff(Surv(OS, OS.Status) ~ Cluster,
                 data      = rt,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS, OS.Status)~ Cluster,
               data      = rt,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

# 配对生存分析
ps <- pairwise_survdiff(Surv(OS, OS.Status)~ Cluster,
                        data            = rt,
                        p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
mycol <- brewer.pal(n = 10, "Paired")[c(2,4,6,8)]
# 绘制基础图形
## 隐藏类标记
fit$strata
names(fit$strata) <- gsub("Cluster=", "", names(fit$strata))
## 生存曲线图
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE, # 不绘制置信区间
                risk.table        = TRUE, # 生存风险表
                risk.table.col    = "strata",
                risk.table.y.text = FALSE,
                palette           = mycol, # KM曲线颜色
                data              = rt,
                xlim = c(0,max(rt$OS)), #展示x轴的范围
                size              = 1,
                legend.title      = "",
                break.time.by = 1, #x轴间隔
                xlab              = "Time (years)",
                ylab              = "Overall survival rate",
                tables.height     = 0.3) # 风险表的高度
p
## 添加overall pvalue
p.lab <- paste0("P",
                ifelse(p.val < 0.001, " < 0.001", # 若P值<0.001则标记为“<0.001”
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 0, y = 0.55, # 在y=0.55处打印overall p值
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
## 添加配对表格
addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                         round(ps$p.value, 3))))
addTab[is.na(addTab)] <- "-"
df <- tibble(x = 7, y = 1, tb = list(addTab))
p$plot <- p$plot + 
  geom_table(data = df, 
             aes(x = x, y = y, label = tb), 
             table.rownames = TRUE)

save(rt,file = "output1/14.tmb.risk.group.survival.Rdata")
##########################################################################################高低风险分组的体细胞突变特征的瀑布图#####################
library(maftools)           
maf = read.maf(maf = 'output1/14combined_maf_value.txt')

# oncoplot(maf = maf,top = 30, fontSize = 0.8 ,showTumorSampleBarcodes = F )

#整理成数据框
a=maf@data
#提取需要的列
a=a %>% .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>% 
  as.data.frame() %>% 
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

#去重复
gene =as.character(unique(a$Hugo_Symbol))#2414
sample = as.character(unique(a$Tumor_Sample_Barcode))#122个样本

#整理成数据框表头
mat=as.data.frame(matrix("",length(gene),length(sample),
                         dimnames = list(gene,sample)))
mat_0_1=as.data.frame(matrix(0,length(gene),length(sample),
                             dimnames = list(gene,sample)))

#把数据填入表格
for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] = as.character(a[i,2])
} 
#填入变异类型
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] = 1}
#进行变异基因计数并整理成计数表格
library(tibble)
gene_count = data.frame(gene=rownames(mat_0_1),
                        count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))


gene_top = gene_count$gene[1:20] # 修改数字，代表TOP多少
save(mat,mat_0_1,gene_count,file = "output1/14.TMB-laml.rda")
#####################################################################
################分组maf文件的构建
##########把突变样本ID转化为分组样本的ID
load(file = "output1/14.TMB-laml.rda")
gene_top = gene_count$gene[1:20]
a=maf@data
#载入分组的信息
#将两个表的样本ID统一
a$Tumor_Sample_Barcode = substring(a$Tumor_Sample_Barcode,1,12)
load("~/0603/output1/14.tmb.risk.group.survival.Rdata")#rt

group <- rownames_to_column(rt,var="sample")#将rt行名转为"sample"列,整个数据框赋给group
colnames(group)[1]="Tumor_Sample_Barcode"

#把分组文件各自拆分
group_high=group[group$group=="High risk",]#46个
group_low=group[group$group=="Low risk",]#40个

#提取cluster的maf文件
maf_1=a[a$Tumor_Sample_Barcode%in%group_high$Tumor_Sample_Barcode,]
maf_2=a[a$Tumor_Sample_Barcode%in%group_low$Tumor_Sample_Barcode,]

#构建两组的maf文件
maf.coad_high=read.maf(maf=maf_1)
maf.coad_low=read.maf(maf=maf_2)

# 3.2作图
#绘制High risk瀑布图
oncoplot(maf = maf.coad_high,
         genes = gene_top,   #显示前20个的突变基因信息
         fontSize = 0.6,   #设置字体大小
         showTumorSampleBarcodes = F)   
#14.risk.high.pubu.8.8
#绘制Low risk瀑布图
oncoplot(maf = maf.coad_low,
         genes = gene_top,   #显示前20个的突变基因信息
         fontSize = 0.6,   #设置字体大小
         showTumorSampleBarcodes = F)
#14.risk.low.pubu.8.8

########################################################
#######################################
rm(list = ls())
setwd("")
##
library(readxl)
rt1 <- read_excel(path = "DTP_NCI60_ZSCORE.xlsx", skip = 7)

colnames(rt1) <- rt1[1,]
rt1 <- rt1[-1,-c(67,68)]

table(rt1$`FDA status`)

rt1 <- rt1[rt1$`FDA status` %in% c("FDA approved", "Clinical trial"),]
rt1 <- rt1[,-c(1, 3:6)]
write.table(rt1, file = "16.drug.txt",sep = "\t",row.names = F,quote = F)

rt2 <- read_excel(path = "data/nci60_RNA__RNA_seq_composite_expression/output/RNA__RNA_seq_composite_expression.xls", skip = 9)
colnames(rt2) <- rt2[1,]
rt2 <- rt2[-1,-c(2:6)]
write.table(rt2, file = "16.geneExp.txt",sep = "\t",row.names = F,quote = F)

###############################################
library(impute)
library(limma)
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
drug <- rt[,2:ncol(rt)]
dimnames <- list(rownames(drug),colnames(drug))
data <- matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)

mat <- impute.knn(data)
drug <- mat$data
library("tibble")
drug <- avereps(drug) %>% t() %>% as.data.frame()
str(drug)

exp <- read.table("geneExp.txt", sep="\t", header=T, row.names = 1, check.names=F)
dim(exp)
exp[1:4, 1:4]

genelist <- c()
genelist <- intersect(genelist,row.names(exp))
genelist
exp <- exp[genelist,] %>% t() %>% as.data.frame()

identical(rownames(exp),rownames(drug))

##======药物敏感性计算
corTab <-cor(drug,exp,method="pearson")
library(WGCNA)
corPval <- corPvalueStudent(corTab,nSamples = nrow(drug))

##=========筛选有显著性的
##为例
fitercor <- lapply(genelist, function(g){
  index <- abs(corTab[,g])> 0.5 & corPval[,g] < 0.01
  df <- cbind(corTab[index,g],corPval[index,g])
  colnames(df) <- c("pearson","Pvalue")
  write.csv(df,file = paste0("output1/",g,"-cor.csv"))
  df
})
length(fitercor) == length(genelist)#[1] TRUE
names(fitercor) <- genelist

save(fitercor,drug,exp,file = "fitercor.Rdata")

###======可视化========
library(ggplot2)
library(ggpubr)

ifelse(dir.exists("output1/16.opFig"),FALSE,dir.create("output1/16.opFig"))
g <- names(fitercor)[1]
#fitercor$SLC38A10[1]
# [1] 0.3359995
lapply(names(fitercor), function(g){
  data <- fitercor[[g]]
  data <- na.omit(data)
  if(!is.null(data)){
    #dr <- rownames(data)[1]
    for(dr in rownames(data)){
      df <- data.frame(exp = exp[,g],dr = drug[,dr])
      tit <- paste0("R:", round(data[dr, 1], 3))
      if (data[dr, 2] < 0.001) {
        tit <- paste0(tit, ", p value < 0.001")
      } else {
        tit <- paste0(tit, ", p value = ", round(data[dr, 2], 3))
      }
      
      p <- ggplot(data = df, aes(x = exp, y = dr)) + #数据映射
        geom_point(alpha = 0.6,shape = 19,size=3,color="#DC143C") +#散点图，alpha就是点的透明度
        #geom_abline()+
        labs(title = tit)+
        geom_smooth(method = lm, formula = y ~ x,aes(colour = "lm"), size = 1.2,se = T)+
        scale_color_manual(values = c("#808080")) + #手动调颜色c("#DC143C","#00008B", "#808080")
        theme_bw() +#设定主题
        theme(axis.title=element_text(size=15,face="plain",color="black"),
              axis.text = element_text(size=12,face="plain",color="black"),
              legend.position =  "none",
              panel.background = element_rect(fill = "transparent",colour = "black"),
              plot.background = element_blank(),
              plot.title = element_text(size=15, lineheight=.8,hjust=0.5, face="plain"),
              legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
        ylab(paste0("Activity z scores of ",dr)) + #expression的作用就是让log10的10下标
        xlab(paste0("The expression of ",g))
      ggsave(filename = paste0("output1/16.opFig/",g,"-",dr,"-cor.pdf"),plot = p,width = 5,height = 5)
    }
  }
})

#####################################pRRophetic####################
library(pRRophetic)
trace(calcPhenotype, edit = T)
trace(summarizeGenesByMean, edit = T)#包和软件版本不兼容，使用这个对源代码进行修改

rm(list = ls())
library(pRRophetic)
library(ggplot2)
library(cowplot)

load("") #读入表达矩阵dat
load("")#读入分组信息,高低风险分组

dat <- tumor_tpm[,rownames(risk_data)]
table(risk_data$group)

data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
allDrugs=unique(drugData2016$Drug.name)
drugs <- c(
  )
allDrugs <- intersect(allDrugs,drugs)#61

jco <- c("#BC3C29", "#0072B5")
info <- IC50 <- expr <- cvOut <- predictedPtype <- predictedBoxdat <- list()

testMatrix = as.matrix(dat[,rownames(risk_data)])
for (drug in ) {
  set.seed() 
  cat(drug," starts!\n")
 
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix=testMatrix,
                                              drug = "CCT018159",
                                              tissueType = "blood",
                                              selection = 1) 
  
  if(!all(names(predictedPtype[[drug]])==rownames(risk_data))) {stop("Name mismatched!\n")} 
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "group"=risk_data$group, 
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$group <- factor(predictedBoxdat[[drug]]$group,levels = c("High risk","Low risk"))
  
}
#计算药物敏感性，写入csv文件，画出高低风险分组之间的箱线图。
p <- vector()
for (drug in allDrugs[1:59]) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$group %in% "High risk"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$group %in% "Low risk"),"est.ic50"]))$p.value
  p <- append(p,tmp)
}
names(p) <- allDrugs[1:59]
sig.drugs <- names(p[p < 0.05])
p.drugs<- data.frame(drug = names(p), p_value = as.numeric(p))
sig.p.drugs <- p.drugs[p.drugs$drug %in% sig.drugs,]
save(predictedBoxdat,file = "Rdata")
###################################################################sig.p.drugs <- read.csv("output1/16.2pvalue.csv",row.names = 1)
my_comparisons <- list(c("High risk","Low risk"))

plotp <- list()
library(ggpubr)
# for (drug in sig.p.drugs$drug) {
for (drug in allDrugs[1:59]) {
p <- ggplot(data = predictedBoxdat[[drug]], aes(x=group, y=est.ic50))
p <- p + geom_boxplot(aes(fill = group)) + 
  scale_fill_manual(values=c("#DD4CA4","#27ABC6")) + 
  theme_minimal() +
  theme(axis.text = element_text (size = 10))+#调整坐标轴字体大小
  theme(axis.title.x=element_text(vjust=0, size=10,face = "plain"))+#调整xlab字体大小
  theme(axis.title.y=element_text(vjust=0, size=10,face = "plain"))+#调整ylab字体大小
  theme(plot.title = element_text(size = 11, face = "bold"))+
  theme(axis.line = element_line(color = "black",linewidth = 0.8))+#调整坐标轴样式
  theme(legend.position = "none")+   # 调整图例标签的大小
  xlab("") + ylab(paste0(drug," Estimated IC50"))+
  ggtitle(drug)
# +theme_bw()
p <- p + stat_compare_means(comparisons = my_comparisons,
                       method="wilcox.test",
                       # label = "p.format",
                       map_signif_level=function(p) sprintf("%.2g",p),
                       label.x = 1.5)
p
plotp[[drug]] <- p
cat(drug," has been finished!\n")
# write.csv(predictedBoxdat[[drug]], file = paste0(" "_est.ic50.csv"))
# ggsave(paste0("", drug, "_boxplot of predicted IC50.pdf"), width = 4, height = 4)
}



























