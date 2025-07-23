library(dplyr)
library(ggplot2)
library(ggsci)
library(corrplot)
library(RColorBrewer)
library(scales)

cols <- pal_npg("nrc")(9)
show_col(cols)



###########################plotPCA################################

rawdata <- read.csv('../data/Phosphosite_data/Phosphosite_Diffanalysis(p.adjust0.05_fc2).csv',
                    check.names = F, stringsAsFactors = F)

length(unique(rawdata$ProteinID))


an_col <- read.csv('../data/Raw_data/an_col.csv')
an_col$Type = factor(an_col$Type, levels = unique((an_col$Type)))

exprs <- dplyr::select(rawdata, an_col$Sample_id)

an_cols <- an_col
plotPCA <- function(data, an_cols, sc=FALSE){
  an_col <- an_cols[which(an_cols$Sample_id %in% colnames(data)),]
  an_col$Type <- factor(an_col$Type)
  exprs<- data[an_col$Sample_id]
  pca <-prcomp(t(exprs), scale. = sc)
  PC = round(pca$sdev **2 / sum(pca$sdev ** 2) *100, 2)
  df <- data.frame(Type = an_col$Type, 
                   PC1 = pca$x[, 1], 
                   PC2 = pca$x[, 2],
                   PC3 = pca$x[, 3], 
                   PC4 = pca$x[, 4])
  
  gpca = ggplot(df,aes(x=PC1, y=PC2, color=Type)) +
    geom_point(size=2) + 
    scale_color_npg() +
    theme_bw()+ theme(panel.grid=element_blank()) +
    theme(plot.title=element_text(hjust=0.5))+ 
    theme(plot.title = element_text(color="black", size=12),
          axis.title.x = element_text(color="black", size=12),
          axis.title.y = element_text(color="black", size=12))+
    theme(axis.text.x = element_text(color="black",size=12),
          axis.text.y = element_text( color="black",size=12))+
    theme(legend.title = element_text(size= 12),
          legend.text = element_text(size=12))+
    # theme(legend.position=c(0.8, 0.8)) + 
    labs(# title = 'PCA', 
      x = paste0('PC1 (', PC[1], '%)'),
      y = paste0('PC2 (', PC[2], '%)')) +
    theme(aspect.ratio = 1)
  
  return(gpca)
}
gpca <- plotPCA(exprs, an_col, sc = T)
gpca

ggsave(gpca, filename = '../data/Phosphosite_data/Phosphosite_PCA.pdf',
       width = 5, height = 4)

######################Intensity###################################

an_col$Raw_Sample_id <- paste0('Raw_', an_col$Sample_id)
rawexprs <- dplyr::select(rawdata, an_col$Raw_Sample_id)

logexprs <- log10(rawexprs)
lexprs <- NULL
for (i in 1:ncol(logexprs)){
  lexprs <- c(lexprs, logexprs[,i])
}

lexprs1 <- data.frame(Value = lexprs)



gi <- ggplot(lexprs1, aes(x=Value)) +
  geom_histogram( binwidth=0.1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw()+ theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_text( color="black",size=12),
        legend.position = "none")+
  labs(x = 'Log10(Intensity)', y = 'Number of phosphosites') +
  theme(aspect.ratio = 0.9)
gi

ggsave(gi, file = '../data/Phosphosite_data/Phosphosite_Intensity.pdf',
       height = 4, width = 5)



###############CV################
exprs <- dplyr::select(rawdata, an_col$Sample_id)

calc_cv <- function(df)(apply(df, 1, sd, na.rm = TRUE) / rowMeans(df, na.rm = TRUE) * 100)

df_CV <- NULL
for (i in unique(an_col$Type)) {
  df_CV <- rbind(df_CV, data.frame(Type = i,
                                   CV = calc_cv(exprs[,an_col[which(an_col$Type == i),
                                                              'Sample_id']])))
}

df_CV$Type <- factor(df_CV$Type, levels = levels(an_col$Type))
gcv <- ggplot(df_CV, aes(x = Type, y=CV, fill=Type)) + 
  geom_boxplot(alpha = 0.7, outlier.size=0.5) + 
  scale_fill_npg() +scale_y_continuous(limits=c(0,50)) +
  theme_bw()+ theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_text( color="black",size=12),
        legend.position = "none")+
  labs(y = 'Coefficient of variation (%)', x = 'Type') +
  theme(aspect.ratio = 1.2)

gcv

ggsave(gcv, filename = '../data/Phosphosite_data/Phosphosite_CV.pdf',
       width = 5, height = 4)




##################Corr_All#########################

exprs <- rawdata[,an_col$Sample_id]
corr <- cor(t(scale(t(exprs))))
pdf(file="../data/Phosphosite_data/Phosphosite_Corr_AllPhos.pdf",
    height =11, width = 11)

hmcol<- rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))

corrplot(corr, tl.col='black',is.corr = FALSE,
         addCoef.col=1,
         rect.col=0,
         method = "color",
         # order = 'hclust',
         # hclust.method = 'ward.D2', 
         # addrect = length(unique(an_col$Type)),
         # hclust.method = "complete",
         col= hmcol,
         tl.pos="lt", number.cex = 1,
         col.lim=c(min(corr),max(corr)),
         # type = "lower",
         addgrid.col = 'gray'
)

dev.off()



##################Corr_Diff#########################

exprs <- rawdata[which(rawdata$Diff == 'YES'),an_col$Sample_id]
corr <- cor(t(scale(t(exprs))))
pdf(file="../data/Phosphosite_data/Phosphosite_Corr_DiffPhos.pdf",
    height =11, width = 11)

hmcol<- rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))

corrplot(corr, tl.col='black',is.corr = FALSE,
         addCoef.col=1,
         rect.col=0,
         method = "color",
         # order = 'hclust',
         # hclust.method = 'ward.D2', 
         # addrect = length(unique(an_col$Type)),
         # hclust.method = "complete",
         col= hmcol,
         tl.pos="lt", number.cex = 1,
         col.lim=c(min(corr),max(corr)),
         # type = "lower",
         addgrid.col = 'gray'
)

dev.off()



#####################pro2phos_count##############################

pros <- data.frame(table(table(rawdata$ProteinID)))

pros1 <- pros[which(pros$Var1 %in% c(1:6)), ]
pros1$Var1 <- as.character(pros1$Var1)
pros1 <-rbind(pros1,
              c('>6', sum(pros[-which(pros$Var1 %in% c(1:6)), 'Freq'])))
pros1$Freq <- as.numeric(pros1$Freq)

pros1$Var1 <- paste0(pros1$Var1, ' (', pros1$Freq, ')')

pros1$Var1 <- factor(pros1$Var1, levels = pros1$Var1)

png = pal_npg()(10)
png = pal_npg()(10)[-8]
png[1] <- png[2]
png[2] <- '#EB6F5D'
pdf('../data/Phosphosite_data/Pro2Phoscount_pie.pdf',
    height = 5, width = 5)
pie(pros1$Freq, labels = pros1$Var1,
    col = png, border = "#e9ecef")
dev.off()


#####################STY precent#######################
num <- table(rawdata$`Amino acid`)  # 统计 S/T/Y 的数量

percent <- round(num / sum(num) * 100, 1)

labels <- paste0(num, "\n(", percent, "%)")

df <- data.frame(
  AminoAcid = paste0('p', names(num)),
  Count = as.numeric(num),
  Label = labels)

p <- ggplot(df, aes(x = AminoAcid, y = Count)) +  # 去掉 fill 分组
  geom_bar(stat = "identity", width = 0.6, fill = "grey30") +  # 固定颜色 grey30
  geom_text(aes(label = Label), vjust = -0.5, size = 4) +
  theme_minimal() +
  labs(x = "Amino Acid", y = "Count") +
  ylim(0, max(df$Count) * 1.2)+
  theme_bw()+ theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12))+
  theme(axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_text( color="black",size=12))


p
ggsave(p, filename ="../data/Phosphosite_data/Phosphosite_STY_barplot.pdf", 
       width = 4, height = 4)







