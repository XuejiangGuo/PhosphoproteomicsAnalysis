library(dplyr)

# 读数据
phos <- read.csv("../data/Phosphosite_data/Phosphosite_Diffanalysis(p.adjust0.05_fc2).csv", 
                 check.names=FALSE)

GPS <- read.csv("../data/KS_Net/KSNet_GPS.csv", check.names =FALSE)
GPS <- KS_net[which(KS_net$Peptide %in% phos$Peptide),]
GPS <- unique(KS_net[,c('PhosphositeID', 'Peptide','Kinase')])

# 统计激酶计数

KScount <- function(GPS, cols) {
  
  df_tab <- as.data.frame(table(GPS$Kinase))
  
  df_tab$rest <- length(unique(GPS$PhosphositeID)) - df_tab$Freq
  colnames(df_tab) <-  c('Kinase', cols)
  return(df_tab)
}

fisher_test <- function(in_GPS, out_GPS) {
  # 假设 KScount 是返回一个数据框的函数，行名为 Kinase
  KS1 <- KScount(in_GPS, c('a', 'b'))  # 比如包含列 a 和 b
  KS2 <- KScount(out_GPS, c('c', 'd'))  # 比如包含列 c 和 d
  
  # 合并两个 data.frame，按行名（KinaseGene）对齐
  KS <- left_join(KS1, KS2, by = "Kinase")
  KS[is.na(KS)] <- 0


  # 进行 Fisher 精确检验
  fisher_results <- apply(dplyr::select(KS, -Kinase), 1, function(x) {
    mat <- matrix(as.numeric(x), nrow = 2, byrow = TRUE)
    fisher.test(mat, alternative = "greater")
  })
  
  # 提取 OR 和 P 值
  KS$oddsratio <- sapply(fisher_results, function(x) x$estimate)
  KS$pvalue <- sapply(fisher_results, function(x) x$p.value)
  
  # 排序和调整 P 值
  KS <- KS[order(KS$pvalue), ]
  KS$p.adjust <- p.adjust(KS$pvalue, method = "BH")
  
  # 输出筛选后显著结果的行数
  print(dim(KS[KS$pvalue < 0.05, ]))
  
  return(KS)
}


clus <- phos$Mfuzz_cluster[phos$Mfuzz_cluster != "NO"]
results <- NULL
for (clu in sort(unique(clus))){
  
  in_ps <- phos$PhosphositeID[phos$Mfuzz_cluster == clu]
  out_ps <- phos$PhosphositeID[phos$Mfuzz_cluster != clu]
  
  in_GPS <- GPS[GPS$PhosphositeID %in% in_ps, ]
  out_GPS <- GPS[GPS$PhosphositeID %in% out_ps, ]
  
  result <- fisher_test(in_GPS, out_GPS)
  result$Mfuzz_cluster <- clu
  results <- rbind(results, result)
  
} 


write.csv(results, "../data/Phosphosite_Mfuzz_KS_enrich/Phosphosite_Mfuzz_KS_enrich.csv", 
          row.names=FALSE)
