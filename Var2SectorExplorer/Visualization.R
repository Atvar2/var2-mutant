suppressPackageStartupMessages(library(MySeuratWrappers))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpmisc))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library("ggpubr"))

#========================= Figure 1 ============================================
#celltype: colors
colors<-c("#E4F4D9","#BAF2AF","#8AF270", "#0BFC11", "#18BF06","#089B43","#13783C","#034F17" ,# 叶肉 0,1,2,3,4,5,6,13
          "#F2E4AA","#DBCC3E", #表皮 7, 12
          "#ABD5EF","#4A91BC", #保卫  10, 14
          "#F2DAE4","#E5389B","#F4ABCD","#F769B7","#D3097D" # 8, 9, 11,15,16
)
colors=c("#BAF2AF","#DBCC3E","#4A91BC","#F2DAE4","#E5389B","#F4ABCD")
# update colors
colors <-c("#E5D4DE","#51AA5D","#F1BC74","#F3B3A2","#D7E8A5",
           "#50AED6","#446C88","#EC5C57","#E59DC5","#1E4129",
           "#AC2E83","#BD966A","#8D529B","#A0A3A9","#E2D5CD",
           "#5D3967","#C6E0BD","#E4C753")
# 处理和对照color
groupcolors=c("#6089A0","#EB5C57")   # WT, color
# 热图颜色
# colorRampPalette(c("steelblue", "white", "darkorange2"))(41)


data<-readRDS("/szrmyy/wangjgLab/scRNA/chenjh/P0006_Arabidopsis/05.Var2geneAnalysis/objAnnoMerge.RDS")
data$Groups<-factor(data$Groups,levels=c("WT","mutant"))

# // UMAP and marker gene
data$seurat_clusters=factor(data$seurat_clusters, levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
data$seurat_clusters=factor(data$seurat_clusters, levels=c(0,1,2,3,4,5,6,13, 7,12, 10,14, 8,9,11,15))
data$Celltype=factor(data$Celltype,levels=c("Mesophyll cell","Epidermal cell", "Guard cell", "Phloem parenchyma","Xylem cell" ,"Compansion cell"))
pdf("1-DimPlotCelltype.pdf")
DimPlot(data,label=F,group.by="Celltype",reduction="umap",pt.size=0.1,cols =c("#C9D69B","#92C5AF","#7F8BAB","#9D629E","#4694A1","#5183B2")) + scale_color_manual(values=colors[1:6])
dev.off()
pdf("1-DimPlotCelltypeGroups.pdf")
DimPlot(data,group.by="Groups",label=F,reduction="umap",pt.size=0.1) + scale_color_manual(values=groupcolors)
dev.off()
pdf("1-DimPlotCelltypeSplit.pdf")
DimPlot(data,group.by="Celltype",split.by="Groups",label=F,reduction="umap",pt.size=0.1) + scale_color_manual(values=colors)
dev.off()

# marker基因 展现
gene=c("AT5G38430","AT2G39470","AT1G70760",# Mesophyll
       "AT2G26250", "AT2G39400",           # 表皮
       "AT3G24140","AT2G46070","AT3G26744", #保卫细胞
       "AT3G48740","AT5G23660", #Phloemparenchyma
       "AT3G25710", 'AT2G29440', #Xylem cell
       "AT1G22710", "AT1G79430" #Companion cell
)

DefaultAssay(data)<-"RNA"
pdf("1-MarkergeneVlnPlotCelltype.pdf")
VlnPlot(data, features = gene,group.by="Celltype",stacked=T,pt.size=0)+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+labs(x="",y="") + scale_color_manual(values=colors)+scale_fill_manual(values=colors[1:6])
dev.off()

marker<-read.csv(file="/szrmyy/wangjgLab/scRNA/chenjh/P0006_Arabidopsis/05.Var2geneAnalysis/1-Celltype_markers.csv",header=T)
dat<-AverageExpression(data)
dat<-dat$SCT
dat<-dat[,c("Mesophyll cell","Phloem parenchyma","Epidermal cell","Xylem cell","Guard cell","Compansion cell")]
pdf("1-Celltype-Top10markers.pdf",height=16,width=8)
pheatmap(dat,scale="row",cluster_rows = FALSE,cluster_cols = FALSE,show_rownames =T,color=colorRampPalette(c(pal_startrek()(2)[2],"white",pal_startrek()(2)[1]))(100),border_color=NA,angle_col =c("45"))
dev.off()

# 比例变化
count<-table(data$Groups,data$Celltype)
percentage<-apply(count,2,function(x){x/rowSums(count)})
p<-percentage[c("mutant","WT"),]
logfc<-percentage["mutant",]/percentage["WT",]
mege<-rbind(p,abs(log2(logfc)))
rownames(mege)<- c("mutant","WT","logfc(abs)")
write.csv(file="CelltypeLogfc.csv",mege)

pdf("1-CelltypeCompositionBarplot.pdf")
dtype<-as.data.frame(table(data@meta.data[,c('Celltype','Groups')]))
ggplot(dtype, aes(fill=Celltype, y=Freq, x=Groups)) +scale_fill_manual(values=colors[1:6])+
  geom_bar(position="fill", stat="identity")+theme_classic()
ggplot(dtype, aes(fill=Celltype, y=Freq, x=Groups)) +
  geom_bar(position="stack", stat="identity")+ scale_fill_manual(values=colors[1:6])+theme_classic()
dev.off()

# // 亚群特异表达基因
markers<-read.csv(file="/data02/backup/02.result/chenjh/P0006_Arabidopsis/05.Var2geneAnalysis/1-Celltype_markers.csv",header=T)
markers$gene<-as.character(markers$gene)
Idents(data)<-"Celltype"
dat<-AverageExpression(data)
dat<-dat$RNA
dat<-dat[,c("Mesophyll cell","Phloem parenchyma","Epidermal cell","Xylem cell","Guard cell","Compansion cell")]
celltype= unique(markers$cluster)
annotation_col = data.frame(
  celltype = factor(celltype,levels = c("Mesophyll cell","Phloem parenchyma","Epidermal cell","Xylem cell","Guard cell","Compansion cell"))
)
rownames(annotation_col)=as.character(colnames(dat))
ann_colors=list(celltype=c("Mesophyll cell"="#E5D4DE", "Phloem parenchyma"="#51AA5D","Epidermal cell"= "#F1BC74", 
                           "Xylem cell"="#F3B3A2","Guard cell"= "#D7E8A5", "Compansion cell"="#50AED6"))
top10 <- markers  %>%  group_by(cluster)  %>%  top_n(n = 15, wt = avg_log2FC)
pdf("1-Celltype-Top15markers.pdf",height=16,width=8)
pheatmap(dat[top10$gene,],scale="row",cluster_rows = FALSE,cluster_cols = FALSE,gaps_col=c(1,2,3,4,5),annotation_col=annotation_col,annotation_colors=ann_colors,show_rownames =F,color=colorRampPalette(c(pal_startrek()(2)[2],"white",pal_startrek()(2)[1]))(100),border_color=NA,angle_col =c("45"))
dev.off()


# // var2 小提琴图 
gene<-"AT2G30950"
pdf("1-Var2geneExpressionBygroups.pdf")
VlnPlot(data, features = gene,group.by="Groups",stacked=T,adjust = 2, pt.size=0)+labs(x="",y="") + scale_color_manual(values=groupcolors)+scale_fill_manual(values=groupcolors)
dev.off()

# // 补光色素柱状图, WT和mutant
data<-read.table(file="calculatedWTAndmutantRatio.xls",sep="\t",header=T)
data<-data.frame(Bin=rownames(data),score=rep(data$AverageScore.1000.,2),type=c(rep("WT",40),rep("mutant",40)),value=c(data$wildRatio,data$mutantRatio))
data$value <-data$value*100
plotdata<-dcast(data,type~score)
plotdata<-plotdata[,-1]
plotdata<-as.matrix(plotdata)
pdf("1-Mesophyllcellsharvestlightw1000.pdf")
ze_barplot <- barplot(plotdata , beside=T , legend.text=T,col=c("#6089A0", "#EB5C57") , ylim=c(0,100) , ylab="Percentage (%)")
legend("topright",
       legend = c("WT", "mutant"),
       fill = c("#6089A0", "#EB5C57"))
dev.off()
p <- data %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="")
colors<-c("#69b3a2", "#404080")
data$Bin<-factor(data$Bin,levels=c(paste("Bin",seq(0,39,1),sep="")))
pdf("1-Mesophyllcellsharvestlightw1000Stacked.pdf",height=8,width=10)
ggplot(data, aes( x = Bin,y=100 * value,fill = type))+
  geom_col(position = 'stack', width = 0.6)+
  #geom_bar(position = "stack", stat = "identity", width = 0.6) 
  theme_bw()+   #设置固定主题为传统的白色背景和深灰色的网格线
  scale_fill_manual(values=colors)+  #自定义的修改填充颜色
  scale_y_continuous(expand = c(0,0))+# 调整y轴属性，使柱子与X轴坐标接触
  theme(
    text=element_text(size=12),
    plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    axis.text.y=element_text(size=12,color = "black"),
    axis.text.x=element_text(size=12,  color = "black",angle = 45, hjust = 0.5,vjust = 0.5),
    legend.title=element_text(size=12), 
    #图例字体设置
    legend.text=element_text(size=12)
    #legend.position = ' none' #删除图例
    #legend.position="bottom" ,#图例位置放中间，可选参数为“left”,“top”, “right”, “bottom”.
    #legend.background = element_rect(fill="lightblue",size=0.5, linetype="solid", colour ="darkblue")
  )+
  theme(panel.grid = element_blank(), #从图中删除非数据元素
        #修改面板背景
        panel.background = element_rect(color = 'black', fill = 'transparent'),
  )+ 
  guides(fill=guide_legend(keywidth = 1, keyheight = 1)) #修改图例的框大小
dev.off()
# // 铁相关基因的表达，待筛选
# var2, 铁相关基因 小提琴图
colors=c("#bababa","#94c4e8")
gene=c("AT1G08830","AT2G28190","AT1G12520","AT1G77490","AT4G08390","AT5G23310")
gene="AT2G30950"
data$Groups<-factor(data$Groups,levels=c("WT","mutant"))
pdf("4-var2Expression.pdf")
VlnPlot(data, features = gene,group.by="Groups",pt.size=0,adjust = .2)+labs(x="",y="") + scale_color_manual(values=colors)+scale_fill_manual(values=colors)
dev.off()


# 补光基因
gene<-c("AT1G29920","AT1G29910","AT1G29930","AT2G34430","AT2G34420","AT2G05100","AT2G05070","AT3G27690","AT5G54270","AT5G01530","AT3G08940","AT2G40100","AT4G10340","AT1G15820","AT3G54890","AT3G61470","AT1G61520","AT3G47470","AT1G45474","AT1G19150")
# 光合作图用
gene<-c("AT1G29920","AT1G29910","AT1G29930","AT2G05100","AT5G54270","AT3G54890","AT3G61470","AT1G61520")
# 光核心
gene<-c("ATCG00020","ATCG00270","ATCG00680","ATCG00280","ATCG00350","ATCG00340")
# ROS清除
gene<-c("AT1G08830","AT2G28190","AT1G07890","AT4G09010")
# 胁迫相关pathway
gene<-c("AT2G31570","AT1G02930","AT5G16980","AT5G01600","AT3G56090")
# // plot gene 
gene<-c("AT1G29910", "AT2G05100", "AT5G54270",  # // 光合作图用
        "ATCG00680", "ATCG00280", "ATCG00350", #// 光核心
        "AT2G31570","AT1G02930","AT3G56090",
        "AT1G08830","AT2G28190"   #// ROS清除
)
DefaultAssay(data)<-"RNA"
median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}

pgene<-gene[1:3]
pdf("2-Keygenesvlnplotmerge.pdf",width=8,height=9)
p<-VlnPlot(data, features = pgene,ncol=3,
           group.by="Groups",pt.size=0,combine = T)+
  stat_summary(fun.y = median.stat, geom='line', size = 10, colour = "darkred") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ scale_color_manual(values=groupcolors)+
  scale_fill_manual(values=groupcolors)
print(p)
dev.off()

# method 2
plotlist=list()
for(i in 1:length(gene)){
  p<-VlnPlot(data, features = gene[i],ncol=3,
             group.by="Groups",pt.size=0,adjust=1.2,combine = T)+
    stat_summary(fun.y = median.stat, geom='line', size = 1, colour = "darkred") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ scale_color_manual(values=groupcolors)+
    scale_fill_manual(values=groupcolors)
  plotlist[[i]]<-p
}
pdf("2-Keygenesvlnplotmerge.pdf",width=8,height=9)
grid.arrange(grobs=plotlist,ncol = 4, padding = unit(1, "line"))
dev.off()

#========================= Figure 2 ============================================
# // Var2 和 Ftsh8 折线图  ---> barplot
Idents(data)<-"Celltype"
WT<-subset(data,Groups=="WT")
mutant<-subset(data,Groups=="mutant")
dat<-AverageExpression(WT)
dat<-dat$RNA
gene<-c("AT2G30950","AT1G06430")
Var2_mutant<-as.data.frame(t(dat[gene,]))
dat<-AverageExpression(mutant)
dat<-dat$RNA
mutantFt<-as.data.frame(dat["AT1G06430",])
colnames(mutantFt)<-"mutantFtsh8"
Averagevar2_Ftsh8<-melt(cbind(Var2_mutant,mutantFt))
Averagevar2_Ftsh8$celltype=rep(c(rownames(mutantFt)),3)
pdf("1-var2_ftsh8Avaerage.pdf",height=6,width=10)
ggplot(data = Averagevar2_Ftsh8,aes(x=celltype,y=value,group = variable,color= variable,shape=variable))+
  geom_point()+
  geom_line()+ scale_color_manual(values=c("#8D93E2","#B398C6","#AA5B9B"))+
  xlab("Celltype")+#横坐标名称
  ylab("Average expression level")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
        legend.position = c(.1,.915),#更改图例的位置，放至图内部的左上角
        legend.box.background = element_rect(color="black"))#为图例田间边框线
dev.off()

# barplot
pdf("1-var2_ftsh8AvaerageBarplot.pdf",height=6,width=10)
ggplot(data = Averagevar2_Ftsh8,aes(x=celltype,y=value,group = variable,fill= variable,shape=variable))+
  geom_bar(stat = 'identity', position = 'dodge')+ scale_color_manual(values=c("#8D93E2","#B398C6","#AA5B9B"))+
  scale_fill_manual(values=c("#8D93E2","#B398C6","#AA5B9B"))+
  xlab("Celltype")+#横坐标名称
  ylab("Average expression level")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
        legend.position = c(.1,.915),#更改图例的位置，放至图内部的左上角
        legend.box.background = element_rect(color="black"))#为图例田间边框线
dev.off()

# // 表达趋势图
wtcells<-sample(colnames(WT), 35000)
mutantcells<-sample(colnames(mutant), 35000)
WT<-subset(WT,cells=wtcells)
mutant<-subset(mutant,cells=mutantcells)

gene20<-list(c("AT1G29920","AT1G29910","AT1G29930","AT2G34430","AT2G34420","AT2G05100","AT2G05070","AT3G27690","AT5G54270","AT5G01530","AT3G08940","AT2G40100","AT4G10340","AT1G15820","AT3G54890","AT3G61470","AT1G61520","AT3G47470","AT1G45474","AT1G19150"))
gene<-list(c("AT1G29920","AT1G29910","AT2G05100","AT2G05070","AT5G54270","AT5G01530","AT3G54890","AT1G19150"))
# var2
gene<-list(c("AT2G30950","AT1G06430"))
DefaultAssay(data) <- "RNA"
WT <- AddModuleScore(
  object = WT,
  features = gene20,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'Lightharvest'
)
WT <- AddModuleScore(
  object = WT,
  features = gene,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'var2'
)
mutant <- AddModuleScore(
  object = mutant,
  features = gene20,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'Lightharvest'
)
mutant <- AddModuleScore(
  object = mutant,
  features = gene,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'var2'
)
WTmet<-WT@meta.data
mutantmet<-mutant@meta.data
WTmet  %>%  arrange(Lightharvest1) ->WTmet
mutantmet %>% arrange(Lightharvest1)->mutantmet

wtkeygene<-WT@assays$RNA@data[c("AT2G30950","AT1G06430"),]             # var2: AT2G30950 FtsH8: AT1G06430; var1:
wtkeygene<-as.data.frame(t(as.data.frame(wtkeygene)))
WTmet$var2epr <-wtkeygene$AT2G30950
WTmet$FtsH8  <-wtkeygene$AT1G06430
WTmet$merge<- WTmet$var2epr + WTmet$FtsH8
mutantkeygene<-mutant@assays$RNA@data[c("AT2G30950","AT1G06430"),]
mutantkeygene<-as.data.frame(t(as.data.frame(mutantkeygene)))
mutantmet$var2epr <-mutantkeygene$AT2G30950
mutantmet$FtsH8  <-mutantkeygene$AT1G06430
write.table(file="1-wtmetsorted.xls",WTmet,sep="\t")
write.table(file="1-mutantmetsorted.xls",mutantmet,sep="\t")

#  补充补光色素基因
WTmet<-WT@meta.data
mutantmet<-mutant@meta.data
WTmet  %>%  arrange(Lightharvest1) ->WTmet
mutantmet %>% arrange(Lightharvest1)->mutantmet

wtkeygene<-WT@assays$RNA@data[c("AT3G56940","AT4G25080","AT4G13250","ATCG00020","AT1G20620","AT1G07890","AT1G08830"),]             # var2: AT2G30950 FtsH8: AT1G06430; var1:
wtkeygene<-as.data.frame(t(as.data.frame(wtkeygene)))
colnames(wtkeygene)<-c("CRD1","CHLM","NYC1","PSBA","CAT3","APX1","CSD1")
WTmet<-cbind(WTmet,wtkeygene)

mutantkeygene<-mutant@assays$RNA@data[c("AT3G56940","AT4G25080","AT4G13250","ATCG00020","AT1G20620","AT1G07890","AT1G08830"),]
mutantkeygene<-as.data.frame(t(as.data.frame(mutantkeygene)))
colnames(mutantkeygene)<-c("CRD1","CHLM","NYC1","PSBA","CAT3","APX1","CSD1")
mutantmet<-cbind(mutantmet,mutantkeygene)
write.table(file="1-wtmetsorted.xls",WTmet,sep="\t")
write.table(file="1-mutantmetsorted.xls",mutantmet,sep="\t")

# ========  merge  展现 ======================================

# // WT Var2 expression
wtbinExr<-read.table("wtexprvar2_Ftsh8.txt",header=T,sep="\t")
wtvioPlot<-wtbinExr[which(wtbinExr$var2epr>0),]
wtline<-read.table("1-var2Bin.xls",col.names=c("Bin_1000","Ratio"),header=F,sep="\t")
wtvioPlot$Bin_1000<-factor(wtvioPlot$Bin_1000,levels=unique(wtvioPlot$Bin_1000))
wtline$Bin_1000 <-factor(wtline$Bin_1000,level=unique(wtline$Bin_1000))
p<-ggplot(data =wtvioPlot ,aes(x = Bin_1000,y = var2epr)) + geom_violin(trim = FALSE,fill="#bababa") +  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                                                                                                                     geom="pointrange", color = "red") +
  geom_line(data = wtline,aes(x = Bin_1000,y = Ratio,group = 1),size=1,color="#bababa") +
  geom_point(data = wtline,aes(x = Bin_1000,y = Ratio),shape=21,fill="white",size=2)+
  scale_y_continuous(breaks=pretty_breaks(10),sec.axis = sec_axis( ~./1,breaks=seq(0,80,5),labels=paste(seq(0,80,5),"%",sep=""),name = "Percentage of Var2 expression"))+
  #scale_color_manual(label = c("Categroy1", "Categroy2"),values = c("#ee8f71","#C10534")) +
  labs(
    title="WT Var2 expression levels",
    subtitle="Percantage of WT Var2 expression",
  )+
  theme_minimal(base_size=16) %+replace%
  theme(
    axis.text.x=element_text(angle=30, hjust=1),
    plot.caption = element_text(hjust=0),
    plot.margin = unit(c(1,0.5,1,0.5), "lines")
  )+theme_classic()
pdf("1-var2ExpressionDistribution.pdf",height=6,width=12)
gg.gap(plot=p,
       segments=c(4,60),
       tick_width = c(1,2),
       ylim=c(0,75))
dev.off()

# // WT  Ftsh8 expression
wtbinExr<-read.table("wtexprvar2_Ftsh8.txt",header=T,sep="\t")
wtvioPlot<-wtbinExr[which(wtbinExr$FtsH8>0),]
wtline<-read.table("1-Ftsh8wtBin.xls",col.names=c("Bin_1000","Ratio"),header=F,sep="\t")
wtvioPlot$Bin_1000<-factor(wtvioPlot$Bin_1000,levels=unique(wtvioPlot$Bin_1000))
wtline$Bin_1000 <-factor(wtline$Bin_1000,level=unique(wtline$Bin_1000))
p<-ggplot(data =wtvioPlot ,aes(x = Bin_1000,y = FtsH8)) + geom_violin(trim = FALSE,fill="#bababa") +  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                                                                                                                   geom="pointrange", color = "red") +
  geom_line(data = wtline,aes(x = Bin_1000,y = Ratio,group = 1),size=1,color="#bababa") +
  geom_point(data = wtline,aes(x = Bin_1000,y = Ratio),shape=21,fill="white",size=2)+
  scale_y_continuous(breaks=pretty_breaks(10),sec.axis = sec_axis( ~./1,breaks=seq(0,25,5),labels=paste(seq(0,25,5),"%",sep=""),name = "Percentage of Ftsh8 expression"))+
  #scale_color_manual(label = c("Categroy1", "Categroy2"),values = c("#ee8f71","#C10534")) +
  labs(
    title="WT Ftsh8 expression levels",
    subtitle="Percantage of WT Ftsh8 expression",
  )+
  theme_minimal(base_size=16) %+replace%
  theme(
    axis.text.x  = element_text(angle=30, vjust=0.5),
    plot.caption = element_text(hjust=0),
    plot.margin = unit(c(1,0.5,1,0.5), "lines")
  ) + theme_classic()
pdf("1-WtFtsh8ExpressionDistribution.pdf",height=6,width=12)
gg.gap(plot=p,
       segments=c(3,12),
       tick_width = c(1,2),
       ylim=c(0,20))
dev.off()
# // mutant  Ftsh8 expression
wtbinExr<-read.table("mutantexprvar2_Ftsh8.txt",header=T,sep="\t")
wtvioPlot<-wtbinExr[which(wtbinExr$FtsH8>0),]
wtline<-read.table("1-Ftsh8mutantBin.xls",col.names=c("Bin_1000","Ratio"),header=F,sep="\t")
wtvioPlot$Bin_1000<-factor(wtvioPlot$Bin_1000,levels=unique(wtvioPlot$Bin_1000))
wtline$Bin_1000 <-factor(wtline$Bin_1000,level=unique(wtline$Bin_1000))
p<-ggplot(data =wtvioPlot ,aes(x = Bin_1000,y = FtsH8)) + geom_violin(trim = FALSE,fill="#94c4e8") +  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                                                                                                                   geom="pointrange", color = "red") +
  geom_line(data = wtline,aes(x = Bin_1000,y = Ratio,group = 1),size=1,color="#94c4e8") +
  geom_point(data = wtline,aes(x = Bin_1000,y = Ratio),shape=21,fill="white",size=2)+
  scale_y_continuous(breaks=pretty_breaks(10),sec.axis = sec_axis( ~./1,breaks=seq(0,25,5),labels=paste(seq(0,25,5),"%",sep=""),name = "Percentage of Ftsh8 expression"))+
  #scale_color_manual(label = c("Categroy1", "Categroy2"),values = c("#ee8f71","#C10534")) +
  labs(
    title="mutant Ftsh8 expression levels",
    subtitle="Percantage of WT Ftsh8 expression",
  )+
  theme_minimal(base_size=16) %+replace%
  theme(
    plot.caption = element_text(hjust=0),
    plot.margin = unit(c(1,0.5,1,0.5), "lines")
  )+theme_classic()
pdf("1-mutantFtsh8ExpressionDistribution.pdf",height=6,width=12)
gg.gap(plot=p,
       segments=c(3,10),
       tick_width = c(1,2),
       ylim=c(0,20))
dev.off()


# // mutant  Ftsh8 expression
wtbinExr<-read.table("wtexprvar2_Ftsh8.txt",header=T,sep="\t")
wtvioPlot<-wtbinExr[which(wtbinExr$var2epr>0),]
wtline<-read.table("1-var2Bin.xls",col.names=c("Bin_1000","Ratio"),header=F,sep="\t")
wtvioPlot$Bin_1000<-factor(wtvioPlot$Bin_1000,levels=unique(wtvioPlot$Bin_1000))
wtline$Bin_1000 <-factor(wtline$Bin_1000,level=unique(wtline$Bin_1000))

pdf("1-mutantFtsh8Line.pdf",height=5,width=14)
ggplot()+ geom_line(data = wtline,aes(x = Bin_1000,y = Ratio,group = 1),size=2,color="#4148A3") +
  geom_point(data = wtline,aes(x = Bin_1000,y = Ratio),shape=23,fill="white",size=2)+
  scale_y_continuous(n.breaks=7,limits=c(0,70),position ="right",name = "Percentage of Ftsh8 expression")+
  labs(
    title="mutant Ftsh8  expression levels",
    subtitle="Percantage of mutant Ftsh8 expression",
  )+
  theme_minimal(base_size=16) %+replace%
  theme(
    plot.caption = element_text(hjust=0),
    plot.margin = unit(c(1,0.5,1,0.5), "lines")
  )+theme_classic()+theme(axis.text.x  = element_text(angle=30, vjust=0.5))
dev.off()
pdf("1-mutantFtsh8Vioplot.pdf",height=5,width=14)
p<-ggplot(data =wtvioPlot ,aes(x = Bin_1000,y = FtsH8)) + geom_violin(trim = F,color="#C8B1DD") +  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),geom="pointrange", color = "red")+ scale_y_continuous(n.breaks=4,limits=c(0,4),position ="left",name = "Percentage of var2 expression")+
  theme_classic()+
  theme(axis.text.x  = element_text(angle=30, vjust=0.5))
print(p)
dev.off()


#//  拟合图  data 所有的细胞
ggplot(data = wtmet, aes(x = Lightharvest1, y = merge)) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(grp.label), "*\": \"*",
                                 after_stat(eq.label), "*\", \"*",
                                 after_stat(rr.label), sep = "")))
formula <- y ~ poly(x, 2, raw = TRUE)
#formula<- y ~ log(x)
my.formula=y ~ x
# 拟合图形
setwd("/szrmyy/wangjgLab/scRNA/chenjh/P0006_Arabidopsis/09.Figures//Figures/Figure2/")
wtdf<-data.frame(group="var2+Ftsh8",LightharvestScore=WTmet$Lightharvest1, expressionLevel=WTmet$merge)
mutantdf<-data.frame(group="Ftsh8",LightharvestScore=mutantmet$Lightharvest1,expressionLevel=mutantmet$FtsH8)
test<-rbind(wtdf,mutantdf)
pdf("2-var2_Ftsh8_geom_smoothmerge.pdf")
ggplot(test, aes(LightharvestScore, expressionLevel, colour=group)) +
  geom_smooth(method = "auto", level=0.90)+
  scale_colour_manual(values=groupcolors)+
  #stat_fit_deviations(formula = formula, colour = "red")
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) + theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),legend.justification=c(0.85,0.1),legend.position=c(0.85,0.1))
dev.off()

#========================= Figure 3 ============================================
# // Me 聚类补光色素相关分析
gene20<-list(c("AT1G29920","AT1G29910","AT1G29930","AT2G34430","AT2G34420","AT2G05100","AT2G05070","AT3G27690","AT5G54270","AT5G01530","AT3G08940","AT2G40100","AT4G10340","AT1G15820","AT3G54890","AT3G61470","AT1G61520","AT3G47470","AT1G45474","AT1G19150"))
gene<-c("AT1G29920","AT1G29910","AT2G05100","AT2G05070","AT5G54270","AT5G01530","AT3G54890","AT1G19150")
DefaultAssay(data) <- "RNA"
pbmc <- AddModuleScore(
  object = data,
  features = gene20,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'Lightharvest'
)
# seurat clusters
pdf("1-LightharvestFeaturePlot.pdf")
p2<-FeaturePlot(pbmc, features = "Lightharvest1",pt.size=0.1,reduction="umap")
print(p2)
dev.off()

Lh<-aggregate(metadata$Lightharvest1, by=list(type=metadata$seurat_clusters),mean)
Lh$type=factor(Lh$type,levels=Lh$type)
pbmc$seurat_clusters<-factor(pbmc$seurat_clusters,levels=Lh$type)
Idents(pbmc)="seurat_clusters"
pdf("4-seuratclustersorderByLightharvest1.pdf")
VlnPlot(object = pbmc, features = 'Lightharvest1',pt.size=0)+stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="black")
dev.off()

pbmc$cluster<-factor(pbmc$cluster,levels=c("Mc1", "Mc6","Mc5","Mc2","Mc4", "Mc3"))
Idents(pbmc)<-"cluster"
pdf("1-clustersorderByLightharvest1.pdf")
VlnPlot(object = pbmc, features = 'Lightharvest1',pt.size=0)+
  scale_color_manual(values=c(Mc1="#EC5C57",Mc2="#E59DC5",Mc3="#1E4129",Mc4="#AC2E83",Mc5="#BD966A",Mc6="#8D529B"))+
  scale_fill_manual(values=c(Mc1="#EC5C57",Mc2="#E59DC5",Mc3="#1E4129",Mc4="#AC2E83",Mc5="#BD966A",Mc6="#8D529B")) +
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="black")
dev.off()

#// Me merge cluster聚类分析
# old version
data$cluster<-"Mc1"
data$cluster[data$seurat_clusters %in% c(4,8)]<-"Mc2"
data$cluster[data$seurat_clusters %in% c(0,2,6,12)]<-"Mc3"
data$cluster[data$seurat_clusters %in% c(3,11)]<-"Mc4"
data$cluster[data$seurat_clusters %in% c(9)]<-"Mc5"
data$cluster[data$seurat_clusters %in% c(15)]<-"Mc6"
Mc2<-subset(data,cluster=="Mc2")
umap_integrated<-as.data.frame(Mc2@reductions$umap_integrated@cell.embeddings)
cells<-rownames(umap_integrated[umap_integrated$umapintegrated_1> 0.1,])
data@meta.data[cells,]$cluster="Mc3"
colors=pal_jco()(6)

colors<-colors[8:13]
data<-readRDS(file="/szrmyy/wangjgLab/scRNA/chenjh/P0006_Arabidopsis/07.subType/Mesophyll_cell/02.updateanno/leafAnnocluster.RDS")
data$cluster=factor(data$cluster,levels=c("Mc1","Mc2", "Mc3", "Mc4", "Mc5", "Mc6"))
pdf("1-DimplotAnnocluster.pdf")
Seurat::DimPlot(data,group.by="cluster",reduction="umap",cols=colors,label=T,label.box =T)  #+scale_color_manual(values=colors)
dev.off()
data$Groups<-factor(data$Groups,levels=c("WT","mutant"))
pdf("1-DimplotAnnoclustersplitBygroups.pdf")
Seurat::DimPlot(data,group.by="Groups",split.by="Groups",cols=groupcolors,reduction="umap") #+scale_color_manual(values=colors)
dev.off()
pdf("1-DimplotAnnoclusterGroups.pdf")
Seurat::DimPlot(data,group.by="Groups",reduction="umap",label=T,label.box =F)+scale_color_manual(values=groupcolors)
dev.off()

#// 每个亚群的功能  markers 直接读入文件 1-cluster_markers.csv
markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <-subset(markers, p_val_adj<0.05)
markers<-read.csv("/szrmyy/wangjgLab/scRNA/chenjh/P0006_Arabidopsis/07.subType/Mesophyll_cell/02.updateanno/1-cluster_markers.csv",header=T)
dat<-AverageExpression(data)
dat<-dat$RNA
anno_col<-data.frame(orig.ident=colnames(dat))
rownames(anno_col)<-factor(colnames(dat),levels=c("Mc1","Mc2", "Mc3", "Mc4", "Mc5", "Mc6"))
ann_colors = list(
  orig.ident = c("Mc1"="#EC5C57","Mc2"="#E59DC5", "Mc3"="#1E4129", "Mc4"="#AC2E83", "Mc5"="#BD966A", "Mc6"="#8D529B")
)
dat<-dat[,c("Mc1","Mc2","Mc3","Mc4","Mc5","Mc6")]

pdf("1-leafcluster-markersgene.pdf",height=20,width=10)
pheatmap(dat[markers$gene,],scale="row",show_rownames=F,gaps_col=c(1,2,3,4,5),annotation_colors = ann_colors,annotation_col = anno_col,cluster_rows = FALSE,cluster_cols = FALSE,color=colorRampPalette(c("#594867","#219186","#e5ea36"))(100),border_color=NA,angle_col=45)  
dev.off()

# //细胞比例图
count<-table(data$Groups,data$cluster)
percentage<-apply(count,2,function(x){x/rowSums(count)})
p<-percentage[c("mutant","WT"),]
logfc<-percentage["mutant",]/percentage["WT",]
mege<-rbind(p,abs(log2(logfc)))
rownames(mege)<- c("mutant","WT","logfc(abs)")
write.csv(file="CelltypeLogfc.csv",mege)

mutant<-subset(data,Groups=="mutant")
WT<-subset(data,Groups=="WT")
mutantpie<-data.frame(celltype=rownames(table(mutant$cluster)),(round(table(mutant$cluster)/sum(table(mutant$cluster))*100,2)))
WTpie<-data.frame(celltype=rownames(table(WT$cluster)),(round(table(WT$cluster)/sum(table(WT$cluster))*100,2)))
mutantpie$celltype<-factor(mutantpie$celltype,levels=levels(Idents(data)))
WTpie$celltype<-factor(WTpie$celltype,levels=levels(Idents(data)))

pdf("1-AnnoclusterPiechatstat.pdf")
p3 <- ggpie(mutantpie, 'Freq',  #绘图，只用写频数就行，切记不用再写分组
            fill = 'celltype', palette = colors,legend = 'right') #按照Cylinders填充，颜色板为jco.  ggsci
print(p3)
p4 <- ggpie(WTpie, 'Freq',  #绘图，只用写频数就行，切记不用再写分组
            fill = 'celltype', palette = colors,legend = 'right') #按照Cylinders填充，颜色板为jco.  ggsci
print(p4)
dev.off()

pdf("1-MesophyllCompositionBarplot.pdf")
dtype<-as.data.frame(table(data@meta.data[,c('cluster','Groups')]))
ggplot(dtype, aes(fill=cluster, y=Freq, x=Groups)) +scale_fill_manual(values=colors[8:13])+
  geom_bar(position="fill", stat="identity")+theme_classic()
ggplot(dtype, aes(fill=cluster, y=Freq, x=Groups)) +
  geom_bar(position="stack", stat="identity")+ scale_fill_manual(values=colors[8:13])+theme_classic()
dev.off()

#  monocle
#// 拟时各个亚类的分布
load("4-ClusterMarkersGenes.rData")
plotdf2=as.data.frame(t(cds@reducedDimS))
colnames(plotdf2)=c("component1","component2")
plotdf2$Pseudotime=cds$Pseudotime
#plotdf2$geneset<-pbmc$mesophyll1[rownames(plotdf2)]
plotdf2$cluster<-pbmc$cluster[rownames(plotdf2)]
plotdf2$AT2G25080<-AT2G25080Expr[rownames(plotdf2),]

pdf("2-PseudotimeDestributionByCluster.pdf")
ggplot(plotdf2, aes(x=Pseudotime,y=cluster,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
dev.off()


# 展示
pdf("4-Pseudotime-monocle.pdf")
plot_cell_trajectory(cds,show_branch_points=T,show_tree=F)
plot_cell_trajectory(cds,color_by="Pseudotime",show_branch_points=F,show_tree=F)
plot_cell_trajectory(cds,color_by="Groups",show_branch_points=F,show_tree=F)+scale_color_startrek(alpha=0.6)
plot_cell_trajectory(cds,color_by="cluster",show_branch_points=F,show_tree=F)+scale_color_manual( values=c(pal_d3()(6)[4], pal_d3()(6)[5], pal_d3()(6)[6]))
dev.off()

BEAM_res <- BEAM(cds, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("fd", "pval", "qval")]

b<-plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                    qval < 1e-4)),],
                               branch_point = 1,
                               num_clusters = 3,
                               cores = 1,
                               use_gene_short_name = T,
                               show_rownames = F,
                               return_heatmap=T)
ic1<-subset(b$annotation_row,Cluster==1)
c2<-subset(b$annotation_row,Cluster==2)
c3<-subset(b$annotation_row,Cluster==3)
c1<-as.character(row.names(c1))
genelist<-unique(c1)
c1_go<-enrichGO(genelist,OrgDb = "org.At.tair.db", keyType="TAIR",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
write.table(c1_go,"c1_go.xls",quote=F,sep="\t")
c2<-as.character(row.names(c2))
genelist<-unique(c2)
c2_go<-enrichGO(genelist,OrgDb = "org.At.tair.db", keyType="TAIR",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
write.table(c2_go,"c2_go.xls",quote=F,sep="\t")
c3<-as.character(row.names(c3))
genelist<-unique(c3)
c3_go<-enrichGO(genelist,OrgDb = "org.At.tair.db", keyType="TAIR", ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
write.table(c3_go,"c3_go.xls",quote=F,sep="\t")

# // 补光映射
# Ros, lightharvest
DefaultAssay(data) <- "RNA"
pbmc <- AddModuleScore(
  object = data,
  features = gene20,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'Lightharvest'
)
ligthHarvestScore <- pbmc$Lightharvest1[colnames(cds)]
pData(cds)$ligthHarvestScore<-ligthHarvestScore
pdf("4-monocleAnnoClusterlightharvest.pdf")
plot_cell_trajectory(cds, color_by = "ligthHarvestScore") +
  scale_color_gradient2(low = muted("#e8efe1"), mid = "#c4dda4", high = muted("#2f6818"), midpoint = 0, guide = "colourbar") # the code does not really use my color bar
dev.off()

pbmc <- AddModuleScore(
  object = data,
  features = gene,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'ros'
)
rosScore <- pbmc$ros1[colnames(cds)]
pData(cds)$rosScore<-rosScore
pdf("4-monocleAnnoClusterRosscore.pdf")
plot_cell_trajectory(cds, color_by = "ligthHarvestScore") +
  scale_color_gradient2(low = muted("#e8efe1"), mid = "#c4dda4", high = muted("#2f6818"), midpoint = 0, guide = "colourbar") # the code does not really use my color bar
dev.off()

# // strss 基因上调
# /szrmyy/wangjgLab/scRNA/chenjh/P0006_Arabidopsis/09.Figures/01.Lightharvest/monocleHeatmap.R  热图


# // 轨迹终点的确定
gene<- c("AT3G48430","AT3G44680","AT1G79000","AT4G34060","AT5G13170","AT1G66580","AT3G10985","AT4G02380","AT2G29350","AT5G14930","AT1G71190","AT1G20620","AT4G23810","AT1G62300","AT2G40750","AT4G01250","AT5G13080","AT3G56400","AT1G34180","AT1G69490","AT3G29035","AT5G39610","AT2G43000","AT1G32640","AT5G46760","AT4G17880","AT4G16430","AT1G01260","AT4G00870","AT2G46510","AT3G18520","AT5G03740","AT4G13250","AT5G04900","AT4G22920","AT5G13800","AT2G23140","AT2G03670")

load("/szrmyy/wangjgLab/scRNA/chenjh/P0006_Arabidopsis/07.subType/Mesophyll_cell/02.updateanno/4-ClusterMarkersGenesAddaging.rData")
load("/szrmyy/wangjgLab/scRNA/chenjh/P0006_Arabidopsis/07.subType/Mesophyll_cell/02.updateanno/4-ClusterMarkersGenes.rData")

pData(Agingcds)$Pseudotime <-pData(cds)$Pseudotime
pData(Agingcds)$State  <-pData(cds)$State
pData(Agingcds)$Size_Factor  <-pData(cds)$Size_Factor
Agingcds@reducedDimS <- cds@reducedDimS
gene<-gene[gene  %in% fData(Agingcds)$gene_short_name]
annocol<-as.data.frame(pData(cds)$cluster)
colnames(annocol)<-"cluster"
rownames(annocol)<-rownames(pData(cds))
#// 热图注释
Binner <- function(cds_object){
  df <- data.frame(pData(cds_object))
  df <- df[,c("Pseudotime", "cluster")]
  df <- df[order(df$Pseudotime, decreasing = F),]
  len <- length(df$Pseudotime)
  bin<-round(len/100)
  State <- c()
  value <- c()
  for(i in 0:99){
    if(i < 99){
      start <- 1+(bin*i)
      stop <- bin+(bin*i)
      index<-median(c(start:stop))
      #value <- median(as.numeric(as.vector(df$cluster[c(start:stop)])))
      value <- as.vector(df$cluster[index])
      State <- c(State, value)
    }
    else{
      State <- c(State, value)
    }
  }
  return(as.data.frame(cluster=State))
}

bin <- Binner(Agingcds)
pdf("2-ExpressionExtrendHeatmap_2.pdf")
plot_pseudotime_heatmap(Agingcds[gene,],
                        cluster_rows = T,
                        cores = 1, add_annotation_col = bin,
                        show_rownames = F, return_heatmap = T, num_clusters = 5)
dev.off()
pdf("2-ExpressionExtrendHeatmap.pdf")
plot_pseudotime_heatmap(Agingcds[gene,],
                        num_clusters = 4,
                        cores = 1, add_annotation_col=annocol,
                        show_rownames = T)
dev.off()
gene<-c("AT3G44680","AT3G10985","AT4G02380","AT1G62300","AT4G01250","AT4G13250","AT4G22920","AT5G13800")
pdf("2-plot_genes_in_pseudotime_eightGenes.pdf",height=6,width=8)
plot_genes_in_pseudotime(Agingcds[gene,], ncol =4,color_by = "cluster")+scale_colour_manual(values=c("#EC5C57","#E59DC5", "#1E4129"))
dev.off()

pdf("Psedo-time_ByGroups.pdf",width=12,height=8)
pData(cds)$Groups<-factor(pData(cds)$Groups,levels=c("WT","mutant"))
plot_cell_trajectory(cds,color_by="Groups")+facet_wrap (~Groups)+scale_color_manual(values=c("#395659","#958A67"))
dev.off()

