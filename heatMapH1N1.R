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
