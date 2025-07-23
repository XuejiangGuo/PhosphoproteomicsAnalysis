library(openxlsx)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
# library(export)
# library(Hmisc)



R.utils::setOption( "clusterProfiler.download.method",'auto')

data <- read.csv('../data/Phosphosite_data/Phosphosite_Diffanalysis(p.adjust0.05_fc2).csv',
                 check.names = F, stringsAsFactors = F)


enrichtest <- function(entrez_id){
  allGO <- enrichGO(na.omit(entrez_id),
                    OrgDb = org.Mm.eg.db,
                    ont='ALL',
                    pAdjustMethod = 'BH',
                    pvalueCutoff = 1, 
                    keyType = 'ENTREZID',
                    minGSSize = 10,
                    maxGSSize = 500,
                    readable = TRUE,
                    pool = FALSE)
  goresult <- allGO@result
  
  kegg <- enrichKEGG(na.omit(entrez_id), 
                     organism = "mmu", 
                     keyType = "ncbi-geneid", 
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH", 
                     minGSSize = 3, 
                     maxGSSize = 3000,
                     use_internal_data = F)
  
  
  keggresult <- setReadable(kegg, org.Mm.eg.db, 
                            keyType = "ENTREZID")@result
  keggresult$ONTOLOGY <- "KEGG"
  
  enrichresult <- bind_rows(goresult, keggresult)
  # c = 3
  # p = 0.05
  # enrichresult <- enrichresult[which(enrichresult$Count >= c &
  #                                      enrichresult$pvalue < p),]
  # q = 0.05
  # enrichresult <- enrichresult[which(enrichresult$Count >= c &
  #                                      enrichresult$p.adjust < q),]
  
  print(table(enrichresult$ONTOLOGY))
  return(enrichresult)
}

clus <- data$Mfuzz_cluster[data$Mfuzz_cluster != "NO"]
results <- NULL

for (clu in sort(unique(clus))){
  entrez_id <- data$ENTREZID[data$Mfuzz_cluster ==clu]
  result <- enrichtest(entrez_id)
  result$Mfuzz_cluster <- clu
  
  results <-  bind_rows(results, result)
  results$qvalue <- NULL
}

write.csv(results, "../data/Phosphosite_Mfuzz_enrich/Phosphosite_Mfuzz_enrich.csv", 
          row.names=FALSE)





