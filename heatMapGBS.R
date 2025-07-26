# 加载必要的包
library(limma)
library(pheatmap)
library(ggplot2)

# 读取数据
data <- read.delim("expr_agg.txt", row.names = 1, check.names = FALSE)

# 查看数据结构
head(data)

# 假设分组：前7列（GSM768539-541）为Group1，后7列（GSM768542-545）为Group2
group <- factor(c(rep("GBS", 7), rep("Normal", 7)))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 使用limma进行差异分析
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(Group2vs1 = GBS - Normal, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 提取差异分析结果
results <- topTable(fit2, number = Inf, coef = "Group2vs1")
results$Gene <- rownames(results)
write.csv(results, "GSE34308results.txt", row.names=FALSE)
# 标记显著差异基因（假设p<0.05且|logFC|>1）
results$Significant <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1, "Yes", "No")
write.csv(results, "GSE34308results.txt", row.names=FALSE)
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_bw() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(adj.P.Val)")

# 选择差异最显著的20个基因
top_genes <- rownames(results)[order(results$adj.P.Val)[1:50]]
# 提取表达矩阵并标准化（Z-score）write.csv(heatmap_data, "GSE34308top_genes.csv", row.names=TRUE)


# 读取样本名称
col_names <- read.csv("GSE34308top_genes.csv", skip = 1, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
sample_names <- as.character(col_names[1, -1])
heatmap_data <- t(scale(t(data[top_genes, ])))

# 分组处理（保持原始分组逻辑不变）
h1n1_cols <- which(sample_groups == "H1N1")
hc_cols <- which(sample_groups == "Normal")  # 这里自动适应新标签

library(pheatmap)

# 读取分组信息
groups <- read.csv("com_samples.csv", nrows = 1, header = FALSE, stringsAsFactors = FALSE)
sample_groups <- as.character(groups[1, -1])

# 将HC标签改为Normal
sample_groups <- ifelse(sample_groups == "HC", "Normal", sample_groups)  # 关键修改

# 读取样本名称
col_names <- read.csv("com_samples.csv", skip = 1, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
sample_names <- as.character(col_names[1, -1])

# 读取基因数据并处理重复
data <- read.csv("com_samples.csv", skip = 2, header = FALSE)
colnames(data) <- c("Symbol", sample_names)

# 移除空Symbol和重复项
data <- data[data$Symbol != "", ]
data <- data[!duplicated(data$Symbol), ]

# 设置行名并转换为数值矩阵
rownames(data) <- data$Symbol
data_matrix <- as.matrix(data[, -1])
mode(data_matrix) <- "numeric"

# 分组处理（保持原始分组逻辑不变）
h1n1_cols <- which(sample_groups == "H1N1")
hc_cols <- which(sample_groups == "Normal")  # 这里自动适应新标签

# 计算差异最大的前50个基因
mean_diff <- rowMeans(data_matrix[, h1n1_cols]) - rowMeans(data_matrix[, hc_cols])
abs_mean_diff <- abs(mean_diff)
top_50_genes <- names(sort(abs_mean_diff, decreasing = TRUE))[1:50]
top_50_data <- data_matrix[top_50_genes, ]

# 创建注释（标签已修改）
annotation_col <- data.frame(Group = factor(sample_groups))
rownames(annotation_col) <- sample_names

# 绘制热图
pheatmap(top_50_data,
         annotation_col = annotation_col,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols ="euclidean",
         cluster_cols = FALSE,  
         fontsize_row = 8)

