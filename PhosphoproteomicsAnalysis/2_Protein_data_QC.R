library(dplyr)
library(ggplot2)
library(ggsci)
library(corrplot)
library(RColorBrewer)
library(scales)

# install.packages("corrplot")
# install.packages("RColorBrewer")

# library(stringr)
# library(ggridges)
# library(hrbrthemes)
# library(DescTools)
cols <- pal_npg("nrc")(9)
show_col(cols)

###########################plotPCA################################
rawdata <- read.csv('../data/Protein_data/Protein_Diffanalysis(p.adjust0.05_fc2).csv',
                    check.names = F, stringsAsFactors = F)

an_col <- read.csv('../data/Raw_data/an_col.csv')
an_col$Type <- factor(an_col$Type, levels = unique(an_col$Type))

exprs <- rawdata[, an_col$Sample_id]

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
    theme(plot.title = element_text(color="black", size=11),
          axis.title.x = element_text(color="black", size=11),
          axis.title.y = element_text(color="black", size=11))+
    theme(axis.text.x = element_text(color="black",size=9),
          axis.text.y = element_text( color="black",size=9))+
    # theme(legend.position=c(0.8, 0.8)) + 
    labs(# title = 'PCA', 
      x = paste0('PC1 (', PC[1], '%)'),
      y = paste0('PC2 (', PC[2], '%)')) +
    theme(aspect.ratio = 1)
  
  return(gpca)
}
gpca <- plotPCA(exprs, an_col, sc = F)
gpca

ggsave(gpca, filename = '../data/Protein_data/Protein_PCA.pdf',
       width = 5, height = 5)

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
  geom_histogram( binwidth=0.1, fill='#F39B7FFF', color="#e9ecef", alpha=0.9) +
  theme_bw()+ theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(color="black", size=11),
        axis.title.x = element_text(color="black", size=11),
        axis.title.y = element_text(color="black", size=11))+
  theme(axis.text.x = element_text(color="black",size=9),
        axis.text.y = element_text( color="black",size=9),
        legend.position = "none")+
  labs(x = 'Log10(Intensity)', y = 'Number of Protein') +
  theme(aspect.ratio = 0.9)

gi
ggsave(gi, file = '../data/Protein_data/Protein_Intensity.pdf',
       height = 5, width = 6)


######################Sample Intensity##############################
# pdf('../data/Protein_data/ProteinRawExpression.pdf',
#     height = 5, width = 6)
# boxplot(log2(rawdata[, paste0('Raw_', an_col$Sample_id)] + 1),
#         ylab='Log2(X+1)',
#         xlab='Sample_id',
#         col=an_col$Type)
# dev.off()


pdf('../data/Protein_data/ProteinNormaExpression.pdf',
    height = 5, width = 10)
boxplot(rawdata[,an_col$Sample_id],
        ylab='ProteinNormaExpression',
        xlab='Sample_id',
        col = an_col$Type)
dev.off()


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
  geom_boxplot(alpha = 0.8, outlier.size=0.5) + 
  scale_fill_npg() +scale_y_continuous(limits=c(0,50)) +
  theme_bw()+ theme(panel.grid=element_blank()) +
  theme(plot.title = element_text(color="black", size=11),
        axis.title.x = element_text(color="black", size=11),
        axis.title.y = element_text(color="black", size=11))+
  theme(axis.text.x = element_text(color="black",size=9),
        axis.text.y = element_text( color="black",size=9),
        legend.position = "none")+
  labs(y = 'Coefficient of variation (%)', x = 'Type') +
  theme(aspect.ratio = 1.2)

gcv

ggsave(gcv, filename = '../data/Protein_data/Protein_CV.pdf',
       width = 5, height = 5)


##################Corr#########################
exprs <- rawdata[which(rawdata$Diff == "YES"),an_col$Sample_id]

corr <- cor(t(scale(t(exprs))))
pdf(file="../data/Protein_data/Protein_Corr.pdf",
    height =10, width = 10)
hmcol<- rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))

corrplot(corr, tl.col='black',is.corr = FALSE,
         addCoef.col=1,
         rect.col=0,
         method = "color",
         # order = 'hclust', 
         # hclust.method = 'ward.D2', 
		 addrect = length(unique(an_col$Type)),
         # hclust.method = "complete",
         col= hmcol,
         tl.pos="lt", number.cex = 1,
         col.lim=c(min(corr),max(corr)),
         # type = "lower",
         addgrid.col = 'gray')

dev.off()






##################Corr_All#########################

exprs <- rawdata[,an_col$Sample_id]
corr <- cor(t(scale(t(exprs))))
pdf(file="../data/Protein_data/Protein_Corr_Allpro.pdf",
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
pdf(file="../data/Protein_data/Protein_Corr_Diffpro.pdf",
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







