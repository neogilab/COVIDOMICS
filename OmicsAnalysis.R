############################################# NEW PLAN
DAT=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Input.txt",row.names = 1,check.names = FALSE)
L=log2(DAT)
write.table(L,file = "/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Log2Data.txt",sep="\t",col.names = NA,quote = FALSE)

library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Metadata.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(HCvsMild=Mild - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_HCvsMild.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsSevere=Severe - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_HCvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsConv=Conv - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsConv",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_HCvsConv.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(MildvsSevere=Severe - Mild,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="MildvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_MildvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(ConvvsMild=Mild - Conv,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ConvvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_ConvvsMild.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(ConvvsSevere=Severe - Conv,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ConvvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_ConvvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

##########################################

library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Metadata.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(HCvsMild=Mild - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_HCvsMild.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsSevere=Severe - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_HCvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(MildvsSevere=Severe - Mild,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="MildvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_MildvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)


######################## MUVR
library(doParallel)
library(MUVR)

nCore=10
nRep=25
nOuter=8
varRatio=0.8
method='RF'

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/X.txt",check.names = FALSE)
Y=c(rep("Mild",29),rep("Severe",12))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar

vip=getVIP(classModel, model='min')

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/Model_Validation.pdf")
plotVAL(classModel)
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/SwimLanePlot.pdf")
plotMV(classModel, model='min')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/StabilityPlot.pdf")
plotStability(classModel, model='mid')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/VIP.pdf")
plotVIP(classModel, model='mim')
dev.off()

vip
write.table(vip,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/Variables.txt",sep="\t",col.names = NA,quote = FALSE)

######################## HC-Mild

nCore=10
nRep=25
nOuter=8
varRatio=0.8
method='RF'

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/XX.txt",check.names = FALSE)
Y=c(rep("Mild",29),rep("HC",31))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar

vip=getVIP(classModel, model='min')

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/Model_Validation.pdf")
plotVAL(classModel)
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/SwimLanePlot.pdf")
plotMV(classModel, model='min')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/StabilityPlot.pdf")
plotStability(classModel, model='mid')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/VIP.pdf")
plotVIP(classModel, model='mim')
dev.off()

vip
write.table(vip,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/Variables.txt",sep="\t",col.names = NA,quote = FALSE)
##################  HC-Severe

nCore=10
nRep=25
nOuter=8
varRatio=0.8
method='RF'

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/XX.txt",check.names = FALSE)
Y=c(rep("Severe",12),rep("HC",31))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar

vip=getVIP(classModel, model='min')

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/Model_Validation.pdf")
plotVAL(classModel)
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/SwimLanePlot.pdf")
plotMV(classModel, model='min')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/StabilityPlot.pdf")
plotStability(classModel, model='mid')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/VIP.pdf")
plotVIP(classModel, model='mim')
dev.off()

vip
write.table(vip,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/Variables.txt",sep="\t",col.names = NA,quote = FALSE)

######################## Heatmap


############################### Volcano

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_MildvsSevere.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Volcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = BM),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,box.padding=0,nudge_x = 0.55,nudge_y = 0.25)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_HCvsMild.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Volcano_HCvsMild.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.55,nudge_y = 1.2)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_HCvsSevere.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Volcano_HCvsSevere.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 1,nudge_y = 0.7)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_MildvsSevere.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Volcano_MildvsSevere.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0.2,nudge_x = -0.2,nudge_y = 0.24)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/M_S.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/M_S_Volcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()
#########################


library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/MetsL.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Heatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Temp/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Zscore.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Metadata.txt",row.names = 1)

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("HC"="#48690E","Conv"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))

ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("Healthy Control (HC)"="#48690E","HC (CoV-2 Ab+)"="#C6E2FF","Hospitalized-mild"="#ffd700","Hospitalized-severe"="#ffa500")))
colnames(sampleinfo)
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_split = 2,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 14),height  = unit(41, "cm"),width  = unit(22, "cm"),
           column_split =c(rep("a_HC",21),rep("b_Antibody_positive",10),rep("c_Covid_mild",29),rep("d_Covid_severe",12)))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Heatmap.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "bottom", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

C1=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/C1.txt")
C=melt(C1)
head(C)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Cluster1.pdf")
ggplot(C, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c("HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-3,3)+
  theme(plot.margin = margin(3.5,3.5,3.5,3.5, "cm"),legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_text(size=15,color="black"),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 15,colour = "black"))+
  labs(x="Zscore")
dev.off()

C2=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/C2.txt")
C=melt(C2)
head(C)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Cluster2.pdf")
ggplot(C, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c("HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-3,3)+
  theme(plot.margin = margin(3.5,3.5,3.5,3.5, "cm"),legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_text(size=15,color="black"),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 15,colour = "black"))+
  labs(x="Zscore")
dev.off()


###################### Correlation

library(psych)
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Metabolom.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Protein.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Results.txt",sep="\t",col.names = NA,quote = FALSE)





##################### Bio Marker UMAP

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/BioM.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+stat_ellipse(aes(colour=Group),level=0.90,linetype = 2)+
  scale_color_manual(labels = c("Hospitalized-mild","Hospitalized-severe"),values=c(Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Hospitalized-mild","Hospitalized-severe"),values=c(Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),axis.text=element_text(size = 12),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()

############# PCA Biomarker
library(PCAtools)
count=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/BioM.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/Metdata.txt",row.names = 1,check.names = FALSE)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/PCA.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_SigniPCA.txt",row.names = 1,check.names = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/PCA.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+stat_ellipse(aes(colour=Group),level=0.95)+
  scale_color_manual(labels = c("Hospitalized-mild","Hospitalized-severe"),values=c(Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Hospitalized-mild","Hospitalized-severe"),values=c(Mild="#ffd700",Severe="#ffa500"))+
  labs(x="PC1, 61.06%",y="PC2, 12.27%")+theme(axis.title = element_text(size=13),axis.text=element_text(size = 12),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()
?stat_ellipse
##################

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/Sev.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()



data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/24_IP.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/24_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/24_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/24_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()



data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/58_IP.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/58_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/58_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/58_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()


################## Biomarker Proteins MIld-Severe

library(doParallel)
library(MUVR)

nCore=10
nRep=25
nOuter=8
varRatio=0.8
method='RF'

XX=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/XX.txt",check.names = FALSE)
Y=c(rep("Mild",29),rep("Severe",12))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar

vip=getVIP(classModel, model='min')

pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Model_Validation.pdf")
plotVAL(classModel)
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/SwimLanePlot.pdf")
plotMV(classModel, model='min')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/StabilityPlot.pdf")
plotStability(classModel, model='mid')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/VIP.pdf")
plotVIP(classModel, model='mim')
dev.off()

vip
write.table(vip,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Variables.txt",sep="\t",col.names = NA,quote = FALSE)

###############


library(psych)
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Metabolom.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Protein.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Results.txt",sep="\t",col.names = NA,quote = FALSE)

#################### IPA figures
detach("package:randomForest", unload=TRUE)
data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_pwy.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_pwy.pdf")
ggplot(data, aes(x=Zscore, y=pval)) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),aes(x=Zscore, y=pval,size=pval),fill="#99aab5",color="#7a8890",pch=21)+
  geom_point(data=subset(data, Zscore <= -2),aes(x=Zscore, y=pval,size=pval),fill="#667BC4",color="#51629c",pch=21)+
  geom_point(data=subset(data, Zscore >=  2),aes(x=Zscore, y=pval,size=pval),fill="#A83131",color="#862727",pch=21)+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,box.padding=0,nudge_x = 0.55,nudge_y = 0.25)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_disease.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_disease.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.3,nudge_x = 0.5,nudge_y = 0.75)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_Toxicity.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_Toxicity.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.3,nudge_x = 0.5,nudge_y = 0.1)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_disease.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_disease.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.5,nudge_x = 0.5,nudge_y = 1)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_Pwy.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_Pwy.pdf")
ggplot(data, aes(x=Zscore, y=pval)) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),aes(x=Zscore, y=pval,size=pval),fill="#99aab5",color="#7a8890",pch=21)+
  geom_point(data=subset(data, Zscore <= -2),aes(x=Zscore, y=pval,size=pval),fill="#667BC4",color="#51629c",pch=21)+
  geom_point(data=subset(data, Zscore >=  2),aes(x=Zscore, y=pval,size=pval),fill="#A83131",color="#862727",pch=21)+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,box.padding=0,nudge_x = 0.55,nudge_y = 0.25)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_Toxicity.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_Toxicity.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.5,nudge_x = 0.5,nudge_y = 0.2)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/MIld_Severe_disease.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/MIld_Severe_disease.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.5,nudge_x = 0.5,nudge_y = 0.1)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


#########################

library(psych)
library(reshape2)

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Metabolome.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Protein.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)

corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Results.txt",sep="\t",col.names = NA,quote = FALSE)
DAT=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Melted.txt")
head(DAT)
Re <- acast(DAT,Var1~Var2)
head(Re)
write.table(Re,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/For_HeatMap.txt",sep="\t",col.names = NA,quote = FALSE)


library(ComplexHeatmap)
library(circlize)
col_corr1= colorRamp2(c(-1, -0.75,-0.5,-0.3, 0.3,0.5,0.75,1), c("#1919ff","#4c4cff","#6f6fff","#9393ff" ,"#ff9999","#ff3232" ,"#ff0000","#cc0000"))
corr=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/For_HeatMap.txt",check.names = FALSE,row.names = 1)
head(corr)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/For_HeatMap.pdf",width = 7,height = 10)
Heatmap(as.matrix(corr),cluster_rows=FALSE,cluster_columns = FALSE,name="Spearman",width = unit(8, "cm"),show_row_names = TRUE,border = FALSE,
        row_names_gp =gpar(fontsize = 8),height  = unit(10, "cm"),column_names_gp =gpar(fontsize = 8),na_col = "#e6e6e6")
dev.off()

#####################

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Metabolome.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Protein.txt",check.names = FALSE)
mbl=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/MBL.txt",check.names = FALSE)

Res1=corr.test(as.matrix(Met),as.matrix(mbl),use = "pairwise",method="spearman",adjust="none")
Res2=corr.test(as.matrix(Prot),as.matrix(mbl),use = "pairwise",method="spearman",adjust="none")

corr=melt(Res1$r)
pval=melt(Res1$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/MBL_Met.txt",sep="\t",col.names = NA,quote = FALSE)

##########################   Correlation Heatmap
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Covid_Metabolome.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Covid_Protein.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")
head(Res$r)
corrplot(Res$r, p.mat = Res$p, sig.level = 0.05,insig = "blank",order="AOE",is.corr = FALSE)
library(corrplot)
warnings()
library(ggcorrplot)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Covid.pdf",width = 10,height = 10)
ggcorrplot(Res$r,p.mat = Res$p,sig.level = 0.1,hc.order = FALSE,insig = "blank",lab_col="black",colors=c("#0000cc","white","#cc0000"),
              tl.cex=5,tl.srt=90,legend.title="Correlation",outline.color="#e6e6e6")
dev.off()
head(DD1)
write.table(DD1$data,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Correlation.txt",sep="\t",col.names = NA,quote = FALSE)

DAT=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Correlation.txt")
head(DAT)
library(reshape2)
Re <- acast(DAT,factor(Var1,levels=unique(Var1))~factor(Var2,levels=unique(Var2)))

head(Re)
library(ComplexHeatmap)
library(circlize)
col_corr= colorRamp2(c(-1, -0.75,-0.5,-0.3,0, 0.3,0.5,0.75,1), c("#1919ff","#4c4cff","#6f6fff","#9393ff" ,"white","#ff9999","#ff3232" ,"#ff0000","#cc0000"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Correlation.pdf",width = 7,height = 10)
Heatmap(as.matrix(Re),col=col_corr,cluster_rows=FALSE,cluster_columns = FALSE,name="Spearman",width  = unit(8, "cm"),show_row_names = TRUE,border = FALSE,
        row_names_gp =gpar(fontsize = 8),height  = unit(23, "cm"),column_names_gp =gpar(fontsize = 8),na_col = "#e6e6e6")
dev.off()

############ For soham

Ol=read.delim("/home/anoop/Desktop/Proteo-Transcriptomics/New_Analysis/PLS/PLS_NEW/FinalFigures/Gender_Deseq2/For_Soham/RNA.txt",check.names=FALSE)
Ref=read.delim("/home/anoop/Desktop/Proteo-Transcriptomics/New_Analysis/PLS/PLS_NEW/FinalFigures/Gender_Deseq2/For_Soham/parp.txt",check.names = FALSE)
Res=corr.test(as.matrix(Ol),as.matrix(Ref),use = "pairwise",method="spearman",adjust="none")
head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/Proteo-Transcriptomics/New_Analysis/PLS/PLS_NEW/FinalFigures/Gender_Deseq2/For_Soham/res_RNA.txt",sep="\t",col.names = NA,quote = FALSE)

################  Sankey Plot Mild-Severe

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/pwy.txt")
library(ggplot2)
library(ggalluvial)
head(XX)
geom_stratum(fill = c("#975d72","#2675a0","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2"),width = 0.2) +
  
col=c(rep("black",24),rep("#92bacf",28))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Pathway.pdf",height = 10,width = 7.5)
ggplot(XX,
       aes(axis1 = Pathway, axis2 = Metabolites)) +
  geom_alluvium(aes(),fill="#ec8e56")+#scale_fill_manual(values = c("#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56"))+
  geom_stratum(fill = c("#98b3c2","#98b3c2","#98b3c2","#98b3c2","#98b3c2",
                        "#98b3c2","#98b3c2","#98b3c2","#98b3c2","#98b3c2",
                        "#98b3c2","#98b3c2","#98b3c2","#98b3c2","#A83131",
                        "#667BC4","#667BC4","#667BC4","#A83131","#A83131",
                        "#A83131","#667BC4","#A83131"),width = 0.2) +
  geom_label(stat = "stratum", infer.label = TRUE,size=4) +
  scale_x_discrete(limits = c("Pathway", "Metabolites"), expand = c(.3, .3)) +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.grid = element_blank())+theme(legend.position = "none",plot.margin = margin(3,0.5,3,0.5, "cm"),
                                                                          axis.text.x = element_text(size=15,color="black"))
dev.off()

############# Violin
Bean=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/BioMarker.txt",check.names = FALSE)
head(Bean)
P1=ggplot(Bean,aes(x=Mets,y=`eicosanedioate (C20-DC)`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "eicosanedioate (C20-DC)")+
  theme_bw()+scale_y_continuous(breaks = seq(-2, 2, by = 1))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P2=ggplot(Bean,aes(x=Mets,y=mannose,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "mannose")+
  theme_bw()+scale_y_continuous(breaks = seq(-2, 2.5, by = 1))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))

P3=ggplot(Bean,aes(x=Mets,y=`eicosenedioate (C20:1-DC)*`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "eicosenedioate (C20:1-DC)*")+
  theme_bw()+scale_y_continuous(breaks = seq(-3, 3, by = 1.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P4=ggplot(Bean,aes(x=Mets,y=`hydantoin-5-propionate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "hydantoin-5-propionate")+
  theme_bw()+scale_y_continuous(breaks = seq(-3.5, 3, by = 2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P5=ggplot(Bean,aes(x=Mets,y=`6-bromotryptophan`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "6-bromotryptophan")+
  theme_bw()+scale_y_continuous(breaks = seq(-3, 2, by = 1.2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

#P6=ggplot(Bean,aes(x=Mets,y=`1-palmitoyl-2-oleoyl-GPC (16:0||18:1)`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "1-palmitoyl-2-oleoyl-GPC")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 1.5, by = 0.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P7=ggplot(Bean,aes(x=Mets,y=`4-hydroxyphenylacetate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "4-hydroxyphenylacetate")+
  theme_bw()+scale_y_continuous(breaks = seq(-2.5, 4, by = 1.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P8=ggplot(Bean,aes(x=Mets,y=`6-oxopiperidine-2-carboxylate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "6-oxopiperidine-2-carboxylate")+
  theme_bw()+scale_y_continuous(breaks = seq(-2.5, 3, by = 1.2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position =c(1.7,0.5),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

#P6=ggplot(Bean,aes(x=Mets,y=`1-palmitoyl-2-oleoyl-GPC (16:0||18:1)`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "1-palmitoyl-2-oleoyl-GPC")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 1.5, by = 0.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(1.4,1.1),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/BioMarker.pdf",height = 10,width =8)
ggarrange(P2,P1,P5,P3,P7,P4,P8,nrow = 4,ncol=2)+theme(plot.margin = margin(0.5,7,0.5,0.5, "cm"))
dev.off()

############### Violin pathway

Bean=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/pathway.txt",check.names = FALSE)
head(Bean)
P1=ggplot(Bean,aes(x=Mets,y=`alpha-ketobutyrate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "alpha-ketobutyrate")+
  theme_bw()+scale_y_continuous(breaks = seq(-3.5, 2.5, by = 1.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=20),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))




P3=ggplot(Bean,aes(x=Mets,y=`glucose`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "D-glucose")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 2, by = 0.6))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P4=ggplot(Bean,aes(x=Mets,y=`formiminoglutamate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "formiminoglutamic acid")+
  theme_bw()+scale_y_continuous(breaks = seq(-2.5, 2.5, by = 1.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))


P5=ggplot(Bean,aes(x=Mets,y=`hydantoin-5-propionate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "hydantoin-5-propionate")+
  theme_bw()+scale_y_continuous(breaks = seq(-3.5, 3, by = 2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P6=ggplot(Bean,aes(x=Mets,y=`indoleacetate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "indoleacetic acid")+
  theme_bw()+scale_y_continuous(breaks = seq(-3, 5, by = 2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P7=ggplot(Bean,aes(x=Mets,y=`alanine`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "L-alanine")+
  theme_bw()+scale_y_continuous(breaks = seq(-1.5, 1, by = 0.6))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P8=ggplot(Bean,aes(x=Mets,y=`tryptophan`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "L-tryptophan")+
  theme_bw()+scale_y_continuous(breaks = seq(-1.5, 1.5, by = 0.6))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P9=ggplot(Bean,aes(x=Mets,y=`quinolinic`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "quinolinic acid")+
  theme_bw()+scale_y_continuous(breaks = seq(-2, 2.5, by = 1))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P2=ggplot(Bean,aes(x=Mets,y=`citrate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "citric acid")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 1, by = 0.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position ="none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))
library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/Pathway.pdf",height = 10,width =30)
ggarrange(P1,P4,P5,P3,P9,P2,P6,P7,P8,nrow = 1)+theme(plot.margin = margin(5,0,5,0, "cm"))
dev.off()


#################################


Bean=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Data.txt",check.names = FALSE)


P1=ggplot(Bean,aes(x=Group,y=HGF,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+
  labs(y="NPX",title = "HGF")+
  theme_bw()+scale_y_continuous(breaks = seq(4, 14, by = 2))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=9),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P2=ggplot(Bean,aes(x=Group,y=PTN,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "PTN")+
  theme_bw()+scale_y_continuous(breaks = seq(0, 8, by = 4))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P3=ggplot(Bean,aes(x=Group,y=CXCL13,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "CXCL13")+
  theme_bw()+scale_y_continuous(breaks = seq(6, 14, by = 2))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P4=ggplot(Bean,aes(x=Group,y=`MCP-3`,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "MCP-3")+
  theme_bw()+scale_y_continuous(breaks = seq(0, 10, by = 3))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P5=ggplot(Bean,aes(x=Group,y=IL12,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "IL12")+
  theme_bw()+scale_y_continuous(breaks = seq(4, 12, by = 1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")




P9=ggplot(Bean,aes(x=Group,y=CXCL12,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "CXCL12")+
  theme_bw()+scale_y_continuous(breaks = seq(0, 2, by = 0.4))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),
        legend.position = c(1.6,1.3))




library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Violin.pdf",height = 5,width =9)
ggarrange(P1,P2,P3,P4,P5,P9,nrow = 2,ncol = 3)+theme(plot.margin = margin(2,5,2,2, "cm"))
dev.off()

##############################################

ip=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/Rank.txt")
head(ip)
ip$Mets <- factor(ip$Mets, levels = ip$Mets)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/Rank.pdf")
ggplot(ip, aes(y=Mets)) + 
  geom_point(data=ip,aes(x=1,y=Mets,size=rank,color=rank))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#191978",high="#7f7fb4",breaks=c(1,10,50,200,400))+
  scale_y_discrete(position = "right")+theme_bw()+coord_fixed(ratio = 0.95)+scale_size(range = c(5,8),breaks=c(1,10,50,200,400))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(5,2,5.4,2.2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(3, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Rank",nrow = 1),color=guide_legend(title = "Rank",nrow=1))
dev.off()



ip=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Rank.txt")
head(ip)
ip$Mets <- factor(ip$Mets, levels = ip$Mets)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Rank.pdf")
ggplot(ip, aes(y=Mets)) + 
  geom_point(data=ip,aes(x=1,y=Mets,size=rank,color=rank))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#65880f",high="#a2b76f",breaks=c(1,5,10,20,50))+
  scale_y_discrete(position = "right")+theme_bw()+coord_fixed(ratio = 0.99)+scale_size(range = c(5,8),breaks=c(1,5,10,20,50))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(5,3.5,7,2.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(1.8, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Rank",nrow = 1),color=guide_legend(title = "Rank",nrow=1))
dev.off()

#################### Donut

data <- data.frame(
  category=c("A", "B", "C"),
  count=c(10, 60, 30)
)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Donut.txt")
head(data)
data$Level <- factor(data$Level, levels = data$Level)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Donut.pdf",width = 10,height = 10)
ggdonutchart(data, size=0.7,"Percent", label = "Label",fill="Category",color = "white")+theme(plot.margin = margin(4,10,4,0, "cm"),axis.text = element_text(size=22),
                                                                                              legend.position = c(1.35,0.5),legend.text = element_text(size = 15),
                                                                                              legend.title = element_blank())
dev.off()


#####
Bean=read.delim("/home/anoop/Desktop/COVID_Omics/Maike/Input.txt")
head(Bean)
P1=ggplot(Bean,aes(x=Group,y=MBL,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "MBL",y="ng/mL")+
  theme_bw()+
  theme(axis.title.x = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = c(1.8, 0.5))+scale_y_continuous(breaks = seq(0.5, 12, by =3))

P2=ggplot(Bean,aes(x=Group,y=Neopterin,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "Neopterin",y="nmol/L")+
  theme_bw()+
  theme(axis.title.x = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = c(1.4, 0.5))+scale_y_continuous(breaks = seq(4, 55, by =15))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Maike/Violin.pdf",height = 5,width =10)
ggarrange(P1,P2,nrow = 1)+theme(plot.margin = margin(2,6,2,0.5, "cm"))
dev.off()

Bean=read.delim("/home/anoop/Desktop/COVID_Omics/Maike/MBL.txt")
pdf("/home/anoop/Desktop/COVID_Omics/Maike/MBL_newActualValue.pdf",height = 5,width =10)
ggplot(Bean,aes(x=Group,y=MBL,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "MBL",y = expression(mu*"g/mL"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),legend.text = element_text(size = 9),
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),plot.margin = margin(2,13,2,0.5, "cm"),
        panel.border = element_blank(),legend.position = c(1.25, 0.5))+scale_y_continuous(breaks = seq(0, 5, by =1))
dev.off()


####################### NEW Limma with Co-factors

library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/MetaData.txt")

design <- model.matrix(~0+Age+BMI+Gender+Group,data = sampleinfo)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
head(design)
cont.matrix <- makeContrasts(HCvsConv=GroupConv - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsConv",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_HCvsConv.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsMild=GroupMild - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_HCvsMild.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsSevere=GroupSevere - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_HCvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(ConvvsMild=GroupMild - GroupConv,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ConvvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_CONVvsMild.txt",sep="\t",quote=FALSE,col.names = NA)


cont.matrix <- makeContrasts(ConvvsSevere=GroupSevere - GroupConv,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ConvvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_CONVvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)


cont.matrix <- makeContrasts(MildvsSevere=GroupSevere - GroupMild,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="MildvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_MildvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

################################
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/MetaDataHC_Covid.txt")

design <- model.matrix(~0+Age+BMI+Gender+Group,data = sampleinfo)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
head(design)
cont.matrix <- makeContrasts(HCvsCovid=GroupCovid - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsCovid",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_HCvsCovid.txt",sep="\t",quote=FALSE,col.names = NA)

########################  3 Group with co-variates

sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/MetaData.txt")

design <- model.matrix(~0+Age+BMI+Gender+Group,data = sampleinfo)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
head(design)
cont.matrix <- makeContrasts(HCvsMild=GroupMild - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_HCvsMild.txt",sep="\t",quote=FALSE,col.names = NA)


cont.matrix <- makeContrasts(HCvsSevere=GroupSevere - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_HCvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(MildvsSevere=GroupSevere - GroupMild,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="MildvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_MildvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)


############# Volcano plots NEW

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_HCvsMild2.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Volcano_HCvsMild.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0.2,nudge_x = -0.2,nudge_y = 0.24)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_HCvsSevere2.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Volcano_HCvsSevere.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0.2,nudge_x = 0.1,nudge_y = 0.15)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_MildvsSevere2.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Volcano_MildvsSevere.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = BM),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0.2,nudge_x = 0.2,nudge_y = 0.1)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Adjusted/MetaDataHC_Covid.txt")

design <- model.matrix(~0+Age+BMI+Gender+Group,data = sampleinfo)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Adjusted/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
head(design)
cont.matrix <- makeContrasts(HCvsCovid=GroupCovid - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsCovid",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Adjusted/Limma_HCvsCovid.txt",sep="\t",quote=FALSE,col.names = NA)


####################################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/MetsL.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Heatmap1.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Zscore.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Design.txt",row.names = 1)

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("HC"="#48690E","Conv"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))

head(Zscore)
ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("Healthy Control (HC)"="#48690E","HC (CoV-2 Ab+)"="#8BB443","Hospitalized-mild"="#ffd700","Hospitalized-severe"="#ffa500")))

MET=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Mets.txt",row.names = 1)
ha2 = rowAnnotation(df = MET,show_annotation_name = FALSE,annotation_legend_param = list(Super_Pathway = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Super_Pathway=c("Amino Acid"="#ed8d18","Carbohydrate"="#017399","Cofactors and Vitamins"="#f2af5d","Energy"="#02ace5",
                                                  "Lipid"="#8e540e","Nucleotide"="#e5dfa4","Peptide"="#7f7c5b")))


ha3 = rowAnnotation(foo = anno_mark(at = c(38,55,100,119,123,150,32,56,58,61,47,111,37,101,122,152,50),
                                    labels_gp = gpar(fontsize=20),lines_gp = gpar(col="black"),link_height = unit(35, "mm"),link_width=unit(12, "mm"),
                                    labels = c("aspartate","glutamate","phenylalanine","fructose","mannose","alpha-ketoglutarate","alanine","glutamine",
                                               "glycine","histidine","cysteine","tryptophan","arginine","proline","glycerate","citrate","cystine")))

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_split = 2,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           right_annotation = ha3,left_annotation= ha2,top_annotation = ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 4),height  = unit(41, "cm"),width  = unit(22, "cm"),
           column_split =c(rep("a_HC",21),rep("b_Antibody_positive",10),rep("c_Covid_mild",27),rep("d_Covid_severe",12)))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Heatmap.pdf",height = 20,width =25)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

packageVersion("GGally")
Cls=row_order(H1)
row.names(Zscore)

??annotation_legend_param
clu = t(t(row.names(Zscore[row_order(H1)[[2]],])))
write.table(as.data.frame(clu),file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/CL2.txt",sep="\t",quote = FALSE,col.names = NA)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Bargraph/superCnt.txt")

library(ggplot2)
data$Super <- factor(data$Super, levels = data$Super)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Bargraph/SuperPathway.pdf")
ggplot(data, aes(x = Type,y = Perc, fill = Super)) +theme_bw()+labs(y="Percetage")+guides(fill=guide_legend(title="Super pathway"))+
  geom_bar(stat="identity")+scale_fill_manual(values = c("Lipid"="#8e540e","Amino Acid"="#ed8d18","Cofactors and Vitamins"="#f2af5d","Peptide"="#7f7c5b",
                                                         "Nucleotide"="#e5dfa4","Carbohydrate"="#017399","Energy"="#02ace5","Partially Characterized Molecules"="#35ccff"))+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),plot.margin = margin(2.5,3,2.5,3, "cm"))


dev.off()

lipid=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Bargraph/Lipid.txt")
head(lipid)
lipid$Type <- factor(lipid$Type, levels = lipid$Type)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Bargraph/Lipid.pdf")
ggplot(lipid, aes(x = Group,y = Perc, fill = Type)) +theme_bw()+labs(y="Percetage")+guides(fill=guide_legend(title="Sub pathway"))+
  geom_bar(stat="identity")+scale_fill_manual(values = c("Fatty Acid, Dicarboxylate"="#994500","Lysophospholipid"="#e56800","Androgenic Steroids"="#ff8f32",
                                                         "Sphingomyelins"="#ffb97f","Monoacylglycerol"="#665600","Diacylglycerol"="#998100","Plasmalogen"="#ccac00",
                                                         "Phosphatidylcholine (PC)"="#fdd017",
                                                         "Long Chain Polyunsaturated Fatty Acid (n3 and n6)"="#746776","Phosphatidylethanolamine (PE)"="#baa5bc",
                                                         "Hexosylceramides (HCER)"="#e9cfec","Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)"="#4e599e",
                                                         "Secondary Bile Acid Metabolism"="#6473cb","Fatty Acid, Monohydroxy"="#8c99e7",
                                                         "Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)"="#0aa091",
                                                         "Fatty Acid Metabolism (Acyl Choline)"="#25cebd",
                                                         "Other"="#9ee9e1"))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),plot.margin = margin(2.5,2,2.5,1.5, "cm"))
dev.off()


############# Sankey  HC-Covid After adjusted

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Sankey/HC_Covid.txt")
library(ggplot2)
library(ggalluvial)
head(XX)

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Sankey/HC_Covid_pwy.pdf",height = 10,width = 10)
ggplot(XX,
       aes(axis1 = Pathway, axis2 = Metabolite)) +
  geom_alluvium(aes(),fill="#ec8e56")+#scale_fill_manual(values = c("#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56"))+
  geom_stratum(width = 0.1) +
  geom_label(stat = "stratum", infer.label = TRUE,size=2) +
  scale_x_discrete(limits = c("Pathway", "Metabolites"), expand = c(.5, .5)) +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.grid = element_blank())+theme(legend.position = "none",plot.margin = margin(1,0.5,1,0.5, "cm"),
                                                                                                       axis.text.x = element_text(size=15,color="black"))
dev.off()



################## NEW Experiment Soham

library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Caco2_Meta.txt",row.names = 1)


library(ComplexHeatmap)
library(circlize)


ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("UnTreated"="#48690E","Treated"="orange")))
colnames(sampleinfo)
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 14),height  = unit(10, "cm"),width  = unit(10, "cm"),
           column_split =c(rep("a_UT",3),rep("b_T",3)))
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

############ ISG Volcano plot
library(ggrepel)
data=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/ISG/Caco2_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/ISG/Caco2_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + xlim(-1.5, 2.5)+
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "#c2f542", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/ISG/Calu3_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/ISG/Calu3_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + xlim(-2, 2)+
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "#c2f542", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/293T_Huh7/ISG/293T_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/293T_Huh7/ISG/293T_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + xlim(-1, 1)+
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "#c2f542", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/293T_Huh7/ISG/Huh7_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/293T_Huh7/ISG/Huh7_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + xlim(-0.5, 0.5)+
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "#c2f542", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()



ip=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Calu3/NonMetabolic.txt")
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Calu3/NonMetabolic.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Ratio,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#191978",high="#7f7fb4",breaks=c(1E-20,1E-10,0.000001,0.0001,0.001,0.05))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=7,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(2,6,4,6.2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Ratio(%)",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()


ip=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Calu3/Metabolic.txt")
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Calu3/Metabolic.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Ratio,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#664c7f",high="#ead6ff",breaks=c(1E-20,1E-10,0.000001,0.0001,0.001,0.05))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=7,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(0,6.2,0,5.7, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(15,0.5),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Ratio(%)",nrow = 2),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()


ip=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Caco2/Input.txt")
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Caco2/Caco2.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Ratio,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#05a167",high="#9bd9c2",breaks=c(0.0001,0.001,0.05,0.1,0.3,0.5))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(3,6,5,4.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(12.5,0.5),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Ratio(%)",nrow = 2),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()

########## pairs plot

library(ggplot2)
library(GGally)
data <- data.frame( var1 = 1:100 + rnorm(100,sd=20), v2 = 1:100 + rnorm(100,sd=27), v3 = rep(1, 100) + rnorm(100, sd = 1)) 
data$v4 = data$var1 ** 2 
data$v5 = -(data$var1 ** 2)
?ggpairs
both=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/ggpair/Data.txt")
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Both.pdf")
ggpairs(both, title="correlogram with ggpairs()") 
dev.off()
head(flea)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/ggpair/Corr.pdf")
ggpairs(both, mapping = aes(fill=Group,color = Group,alpha = 0.5),showStrips=NULL,lower = list(continuous = "smooth",size=5),
        upper = list(continuous = wrap("cor", size = 3)),diag = list(continuous = "barDiag")
        )+theme_bw()+theme(panel.grid = element_blank(),axis.text = element_text(size=9))+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))
dev.off()

####### Correlation Biomarkers 3 groups

library(psych)
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Meta_HC.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Prot_HC.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="none")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_HC_1.txt",sep="\t",col.names = NA,quote = FALSE)

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Meta_Mild.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Prot_Mild.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="none")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_Mild_1.txt",sep="\t",col.names = NA,quote = FALSE)

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Meta_Severe.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Prot_Severe.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="none")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_Severe_1.txt",sep="\t",col.names = NA,quote = FALSE)


Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Severe.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Severe_MBL.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_Severe_MBL.txt",sep="\t",col.names = NA,quote = FALSE)

###

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Covid.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Covid_MBL.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="none")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value
?corr.test
write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_Covid_MBL_1.txt",sep="\t",col.names = NA,quote = FALSE)


################## New violin
library(ggplot2)
Bean=read.delim("/home/anoop/Desktop/COVID_Omics/Maike/NewViolin.txt")
head(Bean)

P1=ggplot(Bean,aes(x=Group,y=Glucose,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glucose")+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-4, 2, by =0.65))

P2=ggplot(Bean,aes(x=Group,y=Mannose,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Mannose")+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 3, by =1))


P3=ggplot(Bean,aes(x=Group,y=MBL,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "MBL",y = expression(mu*"g/mL"))+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(0, 5, by =1.25))

P4=ggplot(Bean,aes(x=Group,y=Alanine,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Alanine")+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.65))


P5=ggplot(Bean,aes(x=Group,y=Tryptophan,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Tryptophan")+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.6))


P6=ggplot(Bean,aes(x=Group,y=Glutamate,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glutamate")+
  theme_bw()+ stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 1, by =1))

my_comparisons = list( c("Mild", "Severe"))
?stat_compare_means
library(ggsignif)
library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Maike/NewViolin.pdf",height = 5,width =9)
ggarrange(P1,P2,P3,P4,P5,P6,nrow = 2,ncol = 3)+theme(plot.margin = margin(1.5,5,1.5,2, "cm"))
dev.off()


#############

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Sankey/Input.txt",header = TRUE)
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Sankey/HC_Covid.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=Term,size=pval,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#fce79a",high="#b3962d")+
  scale_y_discrete(position = "right")+theme_bw()+scale_size(range = c(5,7))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),legend.position = c(6, -0.1),
        legend.box="vertical",legend.text = element_text(size=10,colour = "black"),
        legend.title = element_text(size=10),
        panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(4.5,5,4.5,5.2, "cm"))+
  guides(color=guide_legend(nrow =1 ,title = "-log10(Adj.Pvalue)"), size = guide_legend(nrow = 2,title = "-log10(Adj.Pvalue)"))
dev.off()
################## Supplement Box
Bean=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/Suppl.txt",header = TRUE)  
head(Bean)
P1=ggplot(Bean,aes(x=Group,y=glycine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glycine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 9),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 2, by =0.75))


P2=ggplot(Bean,aes(x=Group,y=proline,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Proline",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.65))

P3=ggplot(Bean,aes(x=Group,y=tryptophan,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Tryptophan",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.55))


P4=ggplot(Bean,aes(x=Group,y=alanine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Alanine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 1, by =0.55))


P5=ggplot(Bean,aes(x=Group,y=histidine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Histidine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.4))


P6=ggplot(Bean,aes(x=Group,y=glutamine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glutamine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.35))


P7=ggplot(Bean,aes(x=Group,y=arginine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Arginine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.75))





P9=ggplot(Bean,aes(x=Group,y=aspartate,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Aspartate",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.75))

P10=ggplot(Bean,aes(x=Group,y=cystine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Cystine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.5))

P11=ggplot(Bean,aes(x=Group,y=phenylalanine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Phenylalanine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.45))


P12=ggplot(Bean,aes(x=Group,y=fructose,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Fructose",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 3, by =1.25))


P13=ggplot(Bean,aes(x=Group,y=mannose,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Mannose",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 3, by =1))


P14=ggplot(Bean,aes(x=Group,y=glucose,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glucose",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 2, by =0.65))

P8=ggplot(Bean,aes(x=Group,y=glutamate,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glutamate",y="Log2(Measurement)")+
  theme_bw()+guides(color=guide_legend(nrow =1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 9),legend.position = c(2.5,-0.2),legend.title=element_blank(),
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 1, by =0.85))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/Suppl.pdf",height = 5,width =12)
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,nrow = 2,ncol = 7)+theme(plot.margin = margin(1,0,1.5,0, "cm"))
dev.off()

############## Facs
library(reshape2)
C1=read.delim("/home/anoop/Desktop/COVID_Omics/ubha/Data.txt",header = TRUE)
head(C1)
C=melt(C1)
library(ggplot2)
library(ggridges)

D1=subset(C,variable=="CD8_GLUT")
P1=ggplot(D1, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(90,102), breaks = seq(90,102, by = 5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "CD8 GLUT")
?stat_density_ridges
D2=subset(C,variable=="IM_GLUT")
P2=ggplot(D2, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(-2,9), breaks = seq(-2,9, by = 5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "IM GLUT")

D3=subset(C,variable=="NCM_GLUT")
P3=ggplot(D3, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(-2,9), breaks = seq(-2,9, by = 5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "NCM GLUT")

D4=subset(C,variable=="NCM_XCT")
P4=ggplot(D4, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(0,120), breaks = seq(0,120, by = 50))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "NCM XCT")

D5=subset(C,variable=="CM_XCT")
P5=ggplot(D5, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(99,100.2), breaks = seq(99,100.2, by = 0.5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "CM XCT")

D6=subset(C,variable=="IM_XCT")
P6=ggplot(D6, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(99,100.2), breaks = seq(99,100.2, by = 0.5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "IM XCT")

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/ubha/Facs.pdf",width = 10,height = 10)
ggarrange(P1,P4,P2,P5,P3,P6,nrow = 3,ncol = 2)+theme(plot.margin = margin(1,1,1,1, "cm"))
dev.off()



C1=read.delim("/home/anoop/Desktop/COVID_Omics/ubha/Data.txt",header = TRUE)
data=read.table("/home/anoop/Desktop/COVID_Omics/ubha/Data.txt",sep="\t",header=TRUE)
rnames <- data[,1] 
mat_data <- data.matrix(data[,2:ncol(data)])
S=scale(mat_data, center = TRUE, scale = TRUE)
rownames(S) <- rnames

C=melt(S)
library(ggplot2)
library(ggridges)
head(C)
D1=subset(C,Var2=="CD8_GLUT")
P1=ggplot(D1, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "bottom",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "CD8 GLUT")

D2=subset(C,Var2=="IM_GLUT")
P2=ggplot(D2, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "IM GLUT")

D3=subset(C,Var2=="NCM_GLUT")
P3=ggplot(D3, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "NCM GLUT")


D4=subset(C,Var2=="NCM_XCT")
P4=ggplot(D4, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "NCM XCT")

D5=subset(C,Var2=="CM_XCT")
P5=ggplot(D5, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "CM XCT")


D6=subset(C,Var2=="IM_XCT")
P6=ggplot(D6, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "IM XCT")


pdf("/home/anoop/Desktop/COVID_Omics/ubha/FacsScaled.pdf",width = )
ggarrange(P1,P4,P2,P5,P3,P6,nrow = 3,ncol = 2)+theme(plot.margin = margin(1,1,1,1, "cm"))
dev.off()


###########################


library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Input.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Heatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Zscore.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Design.txt",row.names = 1)

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)


head(Zscore)
ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("Healthy Control (HC)"="#48690E","HC (CoV-2 Ab+)"="#8BB443","Hospitalized-mild"="#ffd700","Hospitalized-severe"="#ffa500")))

MET=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Pwy.txt",row.names = 1)
ha2 = rowAnnotation(df = MET,show_annotation_name = FALSE,annotation_legend_param = list(Pwy = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                    col = list(Pwy=c("Glycolysis, Gluconeogenesis, and Pyruvate Metabolism"="#0284d0","Fructose, Mannose and Galactose Metabolism"="#f07aaa","TCA Cycle"="#045463")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           left_annotation= ha2,top_annotation = ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 20),height  = unit(18, "cm"),width  = unit(40, "cm"),
           column_split =c(rep("a_HC",21),rep("b_Antibody_positive",10),rep("c_Covid_mild",29),rep("d_Covid_severe",12)))


col_fun_lfc = colorRamp2(c(-1, -0.05,0, 0.5,1), c("#065535","#508871" ,"white", "#ff7f7f","#e50000"))
LFC=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/LFC.txt",row.names = 1)
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(2, "cm"),show_row_names = TRUE,
           heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           row_names_gp =gpar(fontsize = 20),height  = unit(18, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 
t=H1+H2


pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Heatmap.pdf",height = 20,width =25)
draw(t,heatmap_legend_side = "bottom", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

####### Dotplot

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/DotPlot.txt",header = TRUE,check.names = FALSE)
head(data)
P1=ggplot(data, aes(x=Group, y=Fructose,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Fructose")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits=c(-2,4),breaks = seq(-2, 4, by =1.75))


P2=ggplot(data, aes(x=Group, y=Lactate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Lactate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 1.5),breaks = seq(-2, 1, by =1))


P3=ggplot(data, aes(x=Group, y=Pyruvate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Pyruvate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-0.9, 1.5),breaks = seq(-1, 1, by =.65))


P4=ggplot(data, aes(x=Group, y=Citrate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Citrate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-1, 1.5),breaks = seq(-1, 1, by =.65))

P5=ggplot(data, aes(x=Group, y=`Aconitate [cis or trans]`,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Aconitate [cis or trans]")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 1.5),breaks = seq(-2, 1, by =.8))


P6=ggplot(data, aes(x=Group, y=`Alpha-ketoglutarate`,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = expression(alpha*"-ketoglutarate"))+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 2.5),breaks = seq(-2, 3, by =1.3))

P7=ggplot(data, aes(x=Group, y=Fumarate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Fumarate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-5, 1.6),breaks = seq(-5, 1, by =1.7))


P8=ggplot(data, aes(x=Group, y=Malate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Malate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-1, 1.2),breaks = seq(-1, 1, by =.55))


P9=ggplot(data, aes(x=Group, y=Succinate,fill=Group,color=Group,size=6)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Succinate")+
  theme(legend.position = c(1.5,0.5),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-5, 2),breaks = seq(-5, 1.5, by =2))+
  stat_summary(fun = "median", color="blue", fill = "blue", size = 2, geom = "point")

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/DotPlot.pdf",width = 13,height = 10)
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,nrow = 2,ncol = 5)+theme(plot.margin = margin(6,0,6,0, "cm"))
dev.off()

###########


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/DotPlot.txt",header = TRUE,check.names = FALSE)
head(data)
P1=ggplot(data, aes(x=Group, y=Fructose,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Fructose")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits=c(-2,4),breaks = seq(-2, 4, by =1.75))


P2=ggplot(data, aes(x=Group, y=Lactate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Lactate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 1.5),breaks = seq(-2, 1, by =1))


P3=ggplot(data, aes(x=Group, y=Pyruvate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Pyruvate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-0.9, 1.5),breaks = seq(-1, 1, by =.65))


P4=ggplot(data, aes(x=Group, y=Citrate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Citrate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-1, 1.5),breaks = seq(-1, 1, by =.65))

P5=ggplot(data, aes(x=Group, y=`Aconitate [cis or trans]`,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Aconitate [cis or trans]")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 1.5),breaks = seq(-2, 1, by =.8))


P6=ggplot(data, aes(x=Group, y=`Alpha-ketoglutarate`,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = expression(alpha*"-ketoglutarate"))+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 2.5),breaks = seq(-2, 3, by =1.3))

P7=ggplot(data, aes(x=Group, y=Fumarate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Fumarate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-5, 1.6),breaks = seq(-5, 1, by =1.7))


P8=ggplot(data, aes(x=Group, y=Malate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Malate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-1, 1.2),breaks = seq(-1, 1, by =.55))


P9=ggplot(data, aes(x=Group, y=Succinate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Succinate")+
  theme(legend.position = c(1.5,0.5),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-5, 2),breaks = seq(-5, 1.5, by =2))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/BoxPlot.pdf",width = 13,height = 10)
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,nrow = 2,ncol = 5)+theme(plot.margin = margin(6,0,6,0, "cm"))
dev.off()


Data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/Significant_Result.txt",header = TRUE,row.names = 1)
head(Data)
Pval=Data$padj
myPval <- Data["padj"]
myLFC<-Data["log2FoldChange"]

head(myLFC)
BiocManager::install("rsbml")
n
library(libSBML)
library(piano)
install.packages("/home/anoop/Tools/libSBML_5.19.0.tar.gz", repos = NULL, type="source")
library(rsbml)
myGsc <- loadGSC("/home/anoop/Desktop/COVID_Omics/Transcriptomics/Human-GEM.xml")


library(psych)
library(reshape2)
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/ForShuba/Severe.txt",check.names=FALSE)

Res=corr.test(as.matrix(Met),use = "pairwise",method="spearman",adjust="none")
melt(Res$r)
melt(Res$p)


