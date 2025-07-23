library(dplyr)


############phos########################
rawdata <- read.csv('../data/Raw_data/phos_TMT18_spermatocyteannotation.txt', sep='\t',
                    check.names = F,stringsAsFactors = F)

rawdata <- dplyr::select(rawdata, -Protein.ID)
names(rawdata)[names(rawdata) == 'Protein.ISOID'] <- 'ProteinID'


an_col <- read.csv('../data/Raw_data/an_col.csv', check.names = F)

IMF_id <- c("ProteinID", "ENSEMBL", "ENTREZID","SYMBOL", "GENENAME",
             'Proteins',	'Positions within proteins',
            'Localization prob', 'Score diff',
            'PEP','Score', 'Delta score',
            'Number of Phospho (STY)',
            'Sequence window', 'Amino acid', 'Position',
            'Positions','id')


RIC_id <- paste0(an_col$Channle_id, '___', 1)
data <- rawdata[which(
    rawdata$`Potential contaminant` != '+' &
    rawdata$Reverse != '+' &
    rawdata$`Localization prob` >= 0.75 &
    apply(rawdata[,RIC_id], 1, sum) != 0), ]


imf <- data[,IMF_id]
imf <- imf %>%
  mutate(Position = ifelse(
    is.na(Position),
    as.integer(str_split(`Positions within proteins`, ";", simplify = TRUE)[, 1]),
    Position
  ))

imf <- imf %>%
  mutate(PhosphositeID = paste(ProteinID, `Amino acid`, Position, sep = "_")) %>%
  select(PhosphositeID, everything())



rawexprs <- data[, RIC_id]
colnames(rawexprs) <- paste0('Raw_', an_col$Sample_id)

# row norma
exprs <- sweep(rawexprs, 1, rowMeans(rawexprs), FUN = '/')
colnames(exprs) <- paste0(an_col$Sample_id, '_Phos')


imf$`>>PhosphositeRawExpression` <- '>>PhosphositeRawExpression'
output <- cbind(imf, rawexprs)
output$`>>PhosphositeNormaExpression` <- '>>PhosphositeNormaExpression'
output<- cbind(output, exprs)
output$`>>ProteinNormaExpression` <- '>>ProteinNormaExpression'


dfpro <- read.csv('../data/Protein_data/Protein_Diffanalysis(p.adjust0.05_fc2).csv',
                  check.names = F)

pro <- dfpro[,an_col$Sample_id]
colnames(pro) <- paste0(an_col$Sample_id, '_Pro')
pro$ProteinID <- dfpro$ProteinID
dfphos <- left_join(output, pro, by='ProteinID')

exppro <- dfphos[, paste0(an_col$Sample_id, '_Pro')]
# 将 NA 填充为 1
exppro[is.na(exppro)] <- 1
# 将 0 替换为 1e-6
# exppro[exppro == 0] <- 1e-6

# 找出非零的最小值
min_nonzero <- min(exppro[exppro != 0], na.rm = TRUE)
# 将0替换为该非零最小值
exppro[exppro == 0] <- min_nonzero


expphos <- dfphos[, paste0(an_col$Sample_id, '_Phos')] / exppro
colnames(expphos) <- an_col$Sample_id


dfphos$`>>PhosphositeAdjustExpression` <- '>>PhosphositeAdjustExpression'
result <- cbind(dfphos, expphos)

write.csv(result, '../data/Raw_data/Phosphosites_map.csv', row.names = F)
