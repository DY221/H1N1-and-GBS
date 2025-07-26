# 读取数据
data <- read.delim("GSE111368_Non-normalized_data.txt", header = TRUE, sep = "\t", check.names = FALSE)

# 提取所有AVG_Signal列
avg_cols <- grep("^AVG_Signal", colnames(data), value = TRUE)

# Z-score标准化
data[avg_cols] <- scale(data[avg_cols])

# 查看前几行归一化后的数据
head(data[avg_cols])

# 方法二定义Min-Max函数
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# 应用归一化
data[avg_cols] <- lapply(data[avg_cols], normalize)

# 查看前几行归一化后的数据
head(data[avg_cols])

write.table(data, "GSE111368_normalized_data.txt", sep = "\t", row.names = FALSE)
