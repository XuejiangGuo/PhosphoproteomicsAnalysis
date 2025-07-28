library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(export)
library(openxlsx)


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
  c = 3
  # p = 0.05
  # enrichresult <- enrichresult[which(enrichresult$Count >= c &
  #                                      enrichresult$pvalue < p),]
  q = 0.05
  enrichresult <- enrichresult[which(enrichresult$Count >= c &
                                       enrichresult$p.adjust < q),]
  
  print(table(enrichresult$ONTOLOGY))
  return(enrichresult)
}
entrez_id <- data$ENTREZID
enrichresult <- enrichtest(data$ENTREZID)
enrichresult$qvalue <- NULL


write.csv(enrichresult, file = '../data/Phosphoproteins_enrich/phosphoproteins_enrich.csv',
          row.names = F)

table2excel(x=enrichresult,
            file='../data/Phosphoproteins_enrich/phosphoproteins_enrich.xlsx',
            sheetName='enrichresult',
            append = FALSE,
            add.rownames=FALSE,
            fontName="Calibri",
            border=c("top", "bottom"))




################plot#######################

enrich <- read.xlsx('../data/Phosphoproteins_enrich/phosphoproteins_enrich_select.xlsx', 1)


enrich$GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,'/'), 
                      function(x)(as.numeric(x[1])/as.numeric(x[2]))))

enrich$Description <- capitalize(enrich$Description)
# enrich$Description <- gsub("MRNA","mRNA",enrich$Description)


unique(enrich$ONTOLOGY)


plotEnrich <- function(enrich, types) {
  df <- enrich[which(enrich$ONTOLOGY %in% types),]
  df$ONTOLOGY <- factor(df$ONTOLOGY, 
                        levels = unique(df$ONTOLOGY),ordered = T)
  df <- df[order(df$ONTOLOGY,df$GeneRatio),]
  df$Description <- factor(df$Description, 
                           levels = unique(df$Description),ordered = T)
  
  p1 <- ggplot(df, aes(GeneRatio, Description)) + 
    geom_point(aes(size=Count, color = -log10(p.adjust))) +
    scale_color_gradient(low = "blue",  high = "red") + 
    # scale_color_gradient(low="#4DBBD5FF", mid='#91D1C2FF', high = "#E64B35FF") +
    theme_bw() + theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(color="black", size=14),
          axis.title.x = element_text(color="black", size=12),
          axis.title.y = element_text(color="black", size=12))+
    theme(axis.text.x = element_text(color="black",size=10),
          axis.text.y = element_text( color="black",size=10)) +
    labs(
      y = NULL)
  
  if (length(types) > 1){
    p2 <- p1 + facet_grid(ONTOLOGY~.,scales= "free") 
  }
  else{
    p2 <- p1 + facet_grid(ONTOLOGY~.)
  }
  return(p2)
}

(gp <- plotEnrich(enrich, c('BP', 'CC', 'MF')))
ggsave(gp, filename = '../data/Phosphoproteins_enrich/Phosphoproteins_GO.pdf',
       height = 6, width = 6.8)


(kp <- plotEnrich(enrich, 'KEGG'))
ggsave(kp, filename = '../data/Phosphoproteins_enrich/Phosphoproteins_KEGG.pdf',
       height = 4, width = 6)





