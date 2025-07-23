
library(dplyr)
library(stringr)


phos <- read.csv('../data/Phosphosite_data/Phosphosite_Diffanalysis(p.adjust0.05_fc2).csv',
                check.names = FALSE)

# writeLines(paste0(">", phos$PhosphositeID, "\n", phos$Peptide), 
#            "../data/KS_Net/GPS/GPS_input.txt")


data <- read.csv('../data/KS_Net/GPS/GPS_input.gps.txt', sep = '\t',
                 check.names = FALSE)

result <- left_join(phos[,c("PhosphositeID", "Peptide")], data, by = "Peptide") %>%
  distinct()


result$KinaseGene <- sapply(strsplit(result$Kinase, "/"), function(x) tail(x, 1))



write.csv(result, '../data/KS_Net/KSNet_GPS.csv', row.names = FALSE)
