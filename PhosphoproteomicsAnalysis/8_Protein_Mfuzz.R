
library(Mfuzz)
library(ComplexHeatmap)
library(ggsci)
library(circlize)
library(RColorBrewer)


npg <- pal_npg()(10)
colmap <- c(npg[1],npg[2],"#808080")


set.seed(1234)

rawdata <- read.csv('../data/Raw_data/Protein_Diffanalysis(p.adjust0.05_fc2).csv',
                    check.names = F, stringsAsFactors = F)

rownames(rawdata) <- rawdata$ProteinID

an_col <- read.csv('../data/Raw_data/an_col.csv')


# 提取差异表达数据
diff_data <- rawdata[which(rawdata$Diff == 'YES'),
                paste0(unique(an_col$Type), '_mean')]

# 创建 ExpressionSet
eset <- new('ExpressionSet',exprs = as.matrix(diff_data))
eset <- standardise(eset)


cluster_num <- 4
mfuzz_cluster <- mfuzz(eset,
                       c = cluster_num, 
                       m = mestimate(eset))
mfuzz.plot(eset, mfuzz_cluster,mfrow=c(2,3),
           new.window=FALSE,min.mem=0.35,
           time.labels = unique(an_col$Type))


pdf('../data/Protein_Mfuzz/Protein_Mfuzz.pdf', 
    height = 8, width = 12)
mfuzz.plot(eset, mfuzz_cluster,mfrow=c(2,3),
           new.window=FALSE,min.mem=0.35,
           time.labels = unique(an_col$Type))
dev.off()


cluster_df <- data.frame(ProteinID = names(mfuzz_cluster$cluster),
                         Mfuzz_cluster = paste0('Cluster', mfuzz_cluster$cluster),
                         stringsAsFactors = FALSE)

data <- rawdata %>%
  left_join(cluster_df, by = "ProteinID") %>%
  mutate(Mfuzz_cluster = ifelse(Diff == 'YES', Mfuzz_cluster, 'NO'))



table(data$Mfuzz_cluster)


write.csv(data,file = '../data/Protein_data/Protein_Diffanalysis(p.adjust0.05_fc2).csv',
          row.names = F)


######################heatmap##################################
# ---------- 1. 数据读取 ----------
rawdata <- read.csv('../data/Protein_data/Protein_Diffanalysis(p.adjust0.05_fc2).csv',
                    check.names = FALSE, stringsAsFactors = FALSE)
rownames(rawdata) <- rawdata$ProteinID

an_col <- read.csv('../data/Raw_data/an_col.csv')
an_col$Type <- factor(an_col$Type)
an_col$Type <- factor(an_col$Type, levels = unique((an_col$Type)))
an_col <- an_col[order(an_col$Type), ]

# -------- 2. 表达矩阵标准化 --------
df <- rawdata[which(rawdata$Diff == 'YES'), ]
df$Cluster <- df$Mfuzz_cluster
df <- df[order(df$Cluster),]
dfexprs <- df[,an_col$Sample_id]
dfscal = t(scale(t(dfexprs)))


# -------- 3. 注释颜色定义 --------
col_colors <- pal_npg("nrc")(length(unique(an_col$Type)))
names(col_colors) <- unique(an_col$Type)

row_colors <- pal_nejm()(length(unique(df$Cluster)))
names(row_colors) <- unique(df$Cluster)

# -------- 4. 构建注释 --------
ha_col <- HeatmapAnnotation(Type = an_col$Type,
                            col = list(Type = col_colors))

ha_row <- rowAnnotation(Cluster = df$Cluster,
                        col = list(Cluster = row_colors))

# -------- 5. 设置热图色板 --------
hm_colors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))

# -------- 6. 绘制热图 --------
hp <- Heatmap(
  as.matrix(dfscal),
  name = "Z-score",
  col = colorRamp2(seq(-2, 2, length.out = 100), hm_colors),
  border = TRUE,
  cluster_rows = FALSE,           # 关闭行聚类
  cluster_columns = FALSE,        # 关闭列聚类
  column_order = an_col$Sample_id, # 手动设置列顺序
  row_split = df$Mfuzz_cluster,   # 按 Mfuzz_cluster 分组，顺序由因子levels控制
  column_split = an_col$Type,     # 按类型分组，顺序由因子levels控制
  show_column_dend = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  top_annotation = ha_col,
  left_annotation = ha_row,
  row_names_gp = gpar(fontsize = 6, col = "black")
)

hp

# -------- 7. 输出 PDF --------
pdf('../data/Protein_Mfuzz/Protein_Mfuzz_heatmap.pdf', width = 8, height = 8)
draw(hp)
dev.off()
