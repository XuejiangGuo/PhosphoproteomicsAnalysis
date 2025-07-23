library(dplyr)
# library(openxlsx)
# install.packages('DescTools')
library(DescTools)


##################an_col############################
an_col <- read.csv('../data/Raw_data/an_col.csv')
an_col$Type <- factor(an_col$Type, levels = Types)
an_col <- an_col[order(an_col$Type), ]



##############protein####################
rawdata <-  read.csv('../data/Raw_data/TMT18_spermatocyteannotation.txt', sep='\t',
                     check.names = F,stringsAsFactors = F)
rawdata <- dplyr::select(rawdata, -Protein.ID)
names(rawdata)[names(rawdata) == 'Protein.ISOID'] <- 'ProteinID'



data <- rawdata[which(rawdata$`Potential contaminant` != '+' &
                        rawdata$Reverse != '+' & 
                        rawdata$`Unique peptides` >= 1 & 
                        apply(rawdata[, an_col$Channle_id], 1, sum) != 0), ]

IMF_id <- c("ProteinID", "ENSEMBL", "ENTREZID","SYMBOL", "GENENAME",
            'Protein IDs', 'Majority protein IDs',
            'Peptides', 'Razor + unique peptides',
            'Unique peptides', 'Sequence coverage [%]', 
            'Unique + razor sequence coverage [%]',
            'Unique sequence coverage [%]','id')
imf <- data[,IMF_id]



rawexprs <- data[, an_col$Channle_id]
colnames(rawexprs) <- paste0('Raw_', an_col$Sample_id)
# row norma
exprs <- sweep(rawexprs, 1, rowMeans(rawexprs), FUN = '/')
colnames(exprs) <- an_col$Sample_id


##################MeanExpresion&FoldChange#######################

getmean <- function(name)(rowMeans(exprs[, an_col[which(an_col$Type==name), 'Sample_id']]))

result <- NULL
for (i in levels(an_col$Type)){
  result[paste0(i, '_mean')] <- list(getmean(i))
}
result <-data.frame(result, check.names = F)


result$`FoldChange(Max/Min)` = apply(result, 1, max) / apply(result, 1, min)


######anova##########

anova <- function(x){
  df <- data.frame(value = x, group=an_col$Type)
  fit <- aov(value~group,data=df)
  out <- PostHocTest(fit, method="hsd")[["group"]][, 'pval']
  names(out) <- as.character(lapply(names(out), function(x)(sprintf('PostHoc(%s)', x))))
  out['ANOVA_pvalue'] <- summary(fit)[[1]][["Pr(>F)"]][1]
  return(out)
}


result <- cbind(result, data.frame(t(apply(exprs, 1, anova)), check.names = F))
result$`ANOVA_p.adjust` = p.adjust(result$ANOVA_pvalue, method = 'BH')


qval <- 0.05
fc <- 2
result$Diff = 'NO'
result[which(result$`FoldChange(Max/Min)` > fc &
               result$ANOVA_p.adjust < qval), 'Diff'] = 'YES'

table(result$Diff)

imf$`>>ProteinRawExpression` <- '>>ProteinRawExpression'
output<- cbind(imf, rawexprs)
output$`>>ProteinNormaExpression` <- '>>ProteinNormaExpression'
output<- cbind(output, exprs)
output$`>>ProteinDiffResult` <- '>>ProteinDiffResult'
output<- cbind(output, result)


write.csv(output, '../data/Protein_data/Protein_Diffanalysis(p.adjust0.05_fc2).csv',
          na = '',row.names = F)

write.csv(output, '../data/Raw_data/Protein_Diffanalysis(p.adjust0.05_fc2).csv',
          na = '',row.names = F)



