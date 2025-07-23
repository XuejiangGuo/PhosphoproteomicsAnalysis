library(dplyr)
library(DescTools)


###########################Diff################################
rawdata <- read.csv('../data/Raw_data/Phosphosites_map.csv',
                    stringsAsFactors = F,check.names = F)

an_col <- read.csv('../data/Raw_data/an_col.csv')
an_col$Type = factor(an_col$Type, levels = unique((an_col$Type)))


exprs <- dplyr::select(rawdata, an_col$Sample_id)


getmean <- function(name)(rowMeans(exprs[, an_col[which(an_col$Type==name), 'Sample_id']]))

result <- NULL
for (i in levels(an_col$Type)){
  result[paste0(i, '_mean')] <- list(getmean(i))
}
result <-data.frame(result, check.names = F)

levels(an_col$Type)


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
result$Diff <- 'NO'
result[which(result$`FoldChange(Max/Min)` > fc &
               result$ANOVA_p.adjust < qval), 'Diff'] <- 'YES'
table(result$Diff)


rawdata$`>>PhosphositeDiffResult` <- '>>PhosphositeDiffResult'
output<- cbind(rawdata, result)

output$Peptide <- gsub("_", "*", substr(output$`Sequence window`, 9, 23))


write.csv(output, '../data/Phosphosite_data/Phosphosite_Diffanalysis(p.adjust0.05_fc2).csv',
          na = '',row.names = F)
write.csv(output, '../data/Raw_data/Phosphosite_Diffanalysis(p.adjust0.05_fc2).csv',
          na = '',row.names = F)


####################GPS_innput###################
writeLines(paste0(">", output$PhosphositeID, "\n", output$Peptide), 
           "../data/KS_Net/GPS/GPS_input.txt")



######WebLogo######
library(stringr)
first_part <- str_split(output$`Sequence window`, pattern = ";", simplify = TRUE)[, 1]
writeLines(first_part, "../data/WebLogo/STY_WebLogo_input.txt")




