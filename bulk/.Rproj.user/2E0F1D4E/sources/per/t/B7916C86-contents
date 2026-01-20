
#加载R包
library(ggplot2)
library(BiocManager)
library(GEOquery)
library(tidyverse)
library(pheatmap)
library(limma)
library(genefilter)
library(dplyr)
library(DESeq2)
library(venn)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Hs.eg.db)  # 人类数据库，小鼠请替换为org.Mm.eg.db
library(enrichplot)    # 用于高级可视化
library(VennDiagram)
library(glmnet)
library(e1071)
library(kernlab)
library(caret)



getwd()
# 2.数据下载----

# 输入数据信息（）

GEO_id2 <- "GSE80178"
GPL_num2 <- "GPL16686"
SMname2 <- "GSE80178_series_matrix.txt.gz"
GEO_id3 <- "GSE29221"
GPL_num3 <- "GPL6947"
SMname3 <- "GSE29221_series_matrix.txt.gz"

options(timeout = 600)


chip_project2 <- getGEO(GEO = GEO_id2, destdir = '.', getGPL = F)
ann_info2 <- getGEO(GEO = GPL_num2, destdir = ".")
chip_project3 <- getGEO(GEO = GEO_id3, destdir = '.', getGPL = F)
ann_info3 <- getGEO(GEO = GPL_num3, destdir = ".")




# 3.提取临床信息和表达矩阵----


Input_sm2 <- getGEO(filename = SMname2, getGPL = F)#芯片数据提取信息到Input_sm
pd_info2 <- pData(Input_sm2)      # Input_sm中提取临床信息
expr2 <- exprs(Input_sm2)
Input_sm3 <- getGEO(filename = SMname3, getGPL = F)#芯片数据提取信息到Input_sm
pd_info3 <- pData(Input_sm3)      # Input_sm中提取临床信息
expr3 <- exprs(Input_sm3)


#它的作用是用getGEO函数从本地文件中读取数据集的信息，并保存到一个变量中：Input_sm。
#pd_info列名是样本id，可以找其他感兴趣的临床表达数据
#expr 列名是探针id，不可以直接用于差异分析，需要把探针id转换成基因id
#基因id在注释信息中

#提取平台文件的注释信息

Meta(ann_info2)$title
anno_geoquery2 <- Table(ann_info2)
Meta(ann_info3)$title
anno_geoquery3 <- Table(ann_info3)
#把注释信息对应到探针id


# 数据对应前对应前需要进行数据评估
library(RColorBrewer)                               # R配色的一个包
colors = rep(brewer.pal(6, "Set2")[1:2], each = 4)  # 选择颜色

boxplot(expr2, las = 2, outline = F, col = colors)
boxplot(expr3, las = 2, outline = F, col = colors)
R.version

#它的作用是用boxplot函数绘制表达矩阵的箱线图，用来http://127.0.0.1:46819/graphics/plot_zoom_png?width=1186&height=619显示每个样本的基因表达分布。
#它先导入了一个包：RColorBrewer，它提供了一些美观的配色方案。
#然后，它用brewer.pal函数选择了一组颜色，并用rep函数重复了两次，用来区分不同的样本组。
#最后，它用boxplot函数绘制了箱线图，并设置了一些参数，如las，outline和col，来调整图形的样式。

#是否需要对数据进行标准化


expr3 <- normalizeBetweenArrays(expr3)
boxplot(expr3, las = 2, outline = F, col = colors)
boxplot(expr2, las = 2, outline = F, col = colors)
expr2 <- normalizeBetweenArrays(expr2)





probe2symbol2 <- anno_geoquery2[, c("ID", "GB_ACC")] 
colnames(probe2symbol2) <- c("PROBE_ID", "SYMBOL_ID")
probe2symbol3 <- anno_geoquery3[, c("ID", "GB_ACC")] 
colnames(probe2symbol3) <- c("PROBE_ID", "SYMBOL_ID")

# 定义一个函数用于 ID 转换

p2g <- function(expr, probe2symbol) { 
  library(dplyr)
  library(tibble)
  library(tidyr)
  expr <- as.data.frame(expr)
  p2g_expr <- expr %>% 
    rownames_to_column(var = "PROBE_ID") %>%                    # 合并探针的信息
    inner_join(probe2symbol, by = "PROBE_ID") %>%               # 去掉多余信息
    dplyr::select(-PROBE_ID) %>%                                       # 重新排列
    dplyr::select(SYMBOL_ID, everything()) %>%                  # 求出平均数(这边的点号代表上一步产出的数据)
    mutate(rowMean = rowMeans(.[grep("GSM", names(.))])) %>%    # 去除symbol中的NA
    filter(!is.na(SYMBOL_ID))%>%                               # 把表达量的平均值按从大到小排序
    arrange(desc(rowMean)) %>%                                  # symbol留下第一个
    distinct(SYMBOL_ID, .keep_all = T) %>%                      # 反向选择去除rowMean这一列
    dplyr::select(-rowMean) %>%                                 # 列名变成行名
    column_to_rownames(var = "SYMBOL_ID")
  p2g_expr <- p2g_expr[-1, ]
  save(p2g_expr, file = "p2g_expr.Rdata")
  return(p2g_expr)
}


# 4.芯片注释---- 


expr_anno3 <- p2g(expr3, probe2symbol3)
expr_anno3.1 <- p2g(expr3, probe2symbol3)
expr_anno2 <- p2g(expr2, probe2symbol2)



expr_matrix <- expr_anno3

# 2. 检查数据格式
cat("原始数据类型：", class(expr_matrix), "\n")
cat("重复基因数：", sum(duplicated(rownames(expr_matrix))), "\n")

# 3. 处理负数
# 如果所有值都是负数，需要特殊处理
if(any(expr_matrix > 0, na.rm = TRUE)) {
  min_positive <- min(expr_matrix[expr_matrix > 0], na.rm = TRUE)
  expr_matrix[expr_matrix < 0] <- min_positive / 10
} else {
  # 如果没有正值，可以设置为一个小的正数
  warning("所有表达值都为负数或零，将设置一个默认小值")
  expr_matrix[expr_matrix < 0] <- 1e-10
}

# 4. 去除NA值
expr_matrix[is.na(expr_matrix)] <- 0

# 5. 按基因名合并重复行
# 确保是数据框格式进行dplyr操作
if(is.matrix(expr_matrix)) {
  expr_df <- as.data.frame(expr_matrix)
} else {
  expr_df <- expr_matrix
}

expr_unique <- expr_df %>%
  rownames_to_column("GeneSymbol") %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("GeneSymbol")

# 6. 转换为矩阵
expr_matrix <- as.matrix(expr_unique)

# 7. 验证
cat("最终数据类型：", class(expr_matrix), "\n")
cat("最终行数：", nrow(expr_matrix), "\n")
cat("最终列数：", ncol(expr_matrix), "\n")
cat("重复基因数：", sum(duplicated(rownames(expr_matrix))), "\n")

# 8. 写入文件（适合CIBERSORT的格式）
write.table(expr_matrix, "cibersort_input.txt", sep = "\t", quote = FALSE, col.names = NA)





# 安装并加载 AnnotationDbi 和注释包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))  # 小鼠需改为 org.Mm.eg.db
library(AnnotationDbi)
library(org.Hs.eg.db)

# 1. 处理行名（去除RefSeq的版本号）
ids <- rownames(expr_anno2)
ids_no_version <- gsub("\\..*", "", ids)  # 例如将"NM_001234.5"转换为"NM_001234"

# 2. 转换RefSeq ID为Gene Symbol
symbol <- mapIds(
  org.Hs.eg.db,           # 人类基因注释库
  keys = ids_no_version,  # 输入去除版本号的RefSeq ID
  column = "SYMBOL",      # 目标转换类型：基因符号
  keytype = "REFSEQ",     # 输入ID类型：RefSeq mRNA
  multiVals = "first"     # 多匹配时取第一个结果
)

# 3. 处理未匹配的ID（保留原始RefSeq ID）
symbol[is.na(symbol)] <- ids_no_version[is.na(symbol)]

# 4. 对重复基因符号的表达值取平均值
# 将表达矩阵转换为数据框（方便分组操作）
expr_df <- as.data.frame(expr_anno2)
# 添加基因符号列（用于分组）
expr_df$gene_symbol <- symbol
# 按基因符号分组，对每个样本列计算平均值
expr_mean <- aggregate(
  . ~ gene_symbol,  # 按gene_symbol分组，对所有样本列（.）计算
  data = expr_df,
  FUN = mean        # 聚合函数：平均值
)

# 5. 整理结果（将基因符号设为行名，删除多余列）
rownames(expr_mean) <- expr_mean$gene_symbol  # 基因符号作为行名
expr_mean <- expr_mean[, -1]  # 移除gene_symbol列（已作为行名）

# 查看处理结果（去重后并取平均的表达矩阵）
head(expr_mean)
expr_anno1<- expr_mean










group_info2 <- as.data.frame(cbind(pd_info2$geo_accession, pd_info2$title)) # 获取分组信息
group_info3 <- as.data.frame(cbind(pd_info3$geo_accession, pd_info3$title)) # 获取分组信息


colnames(group_info2) <- c("Id", "group")  # 修改列名
group_info2 <- group_info2[order(group_info2[, 2]), ]
group_info <- rbind(group_info1, group_info2)
colnames(group_info3) <- c("Id", "group")  # 修改列名
group_info3 <- group_info3[order(group_info3[, 2]), ]

#手动整理 Group_info.csv 并得到 sample1 和 sample2 文件
DF<- as.data.frame(read.table(file = "sample1.txt"))
control<- as.data.frame(read.table(file = "sample2.txt"))


#它的作用是从两个文件中读取正常和肿瘤样本的基因表达数据，
#并保存到两个数据框中，命名为normal和tumor。
colnames(DF) <- c("DF")
colnames(control) <- c("control")

# 并转换为数据框类的对象
#然后，它用colnames函数给这两个数据框的列重新命名为******normal和turmor。
DF_num <- length(DF$DF)
control_num <- length(control$control)
#接着，它用length函数计算出正常和肿瘤样本的数量，
#并保存到两个变量中，命名为normal_num和tumor_num。
Total_info <- c(DF$DF, control$control)
#最后，它用c函数把normal和tumor两个数据框中的数据合并为一个向量，
#并保存到一个变量中，命名为Total_info。

dim(expr_anno3)[1]

out_exp_counts <- as.data.frame(c(1:38770))


# for循环对表达矩阵重新排序
for (each in Total_info) {
  out_exp_counts <- cbind(out_exp_counts, subset(expr_anno3, select = c(each)))
}

#接着，它用for循环遍历Total_info向量中的每个元素，它表示样本的编号。
#在每次循环中，它用subset和cbind函数
#从expr_anno数据框中选取对应样本的基因表达数据，并添加到新的数据框中。
out_exp_counts <- out_exp_counts[, -1]
#最后，它用方括号运算符去除新数据框中的第一列，它是基因的编号。

# 对数据框的每个元素处理：≤0设为NaN，再log2转换
out_exp_counts[out_exp_counts <= 0] <- NaN  # 直接用逻辑索引，无需which()
out_exp_counts_log2 <- log2(out_exp_counts)
# 保留无NaN的行（去掉含NaN的行）
out_exp_counts_log2_clean <- out_exp_counts_log2[complete.cases(out_exp_counts_log2), ]
out_exp_counts<-out_exp_counts_log2_clean







sample_info <- data.frame(
  row.names = colnames(expr_matrix),  # 行名必须与表达矩阵的列名完全一致
  condition = c(rep("Disease", 12), rep("Control", 12))
)

# 检查样本顺序是否匹配
all(rownames(sample_info) == colnames(expr_matrix))
# 如果返回 TRUE，则说明匹配成功，可以继续

# 加载绘图包
library(ggplot2)

# 1. 转置矩阵
# prcomp函数要求样本在行，基因在列，所以需要用 t() 转置
counts_for_pca <- t(expr_matrix)

# 2. 运行PCA
# scale. = TRUE 是非常重要的一步，它会对每个基因（现在是列）进行标准化
# 这可以确保表达量高的基因不会主导整个PCA结果
pca_result <- prcomp(counts_for_pca, scale. = T)

# 3. 提取PCA结果用于绘图
# pca_result$x 包含了每个样本在各个主成分上的坐标
pca_data <- as.data.frame(pca_result$x[, 1:2]) # 我们只取前两个主成分
pca_data$condition <- sample_info$condition    # 添加分组信息
pca_data$sample_name <- rownames(pca_data)     # 添加样本名
library(ggpubr)

# 计算每个主成分的解释方差比例
percent_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = sample_name)) +
  geom_point(size = 4) + 
  # 绘制样本点
  
  # geom_text(hjust = 0.5, vjust = -0.8, size = 3) + # 如果需要，可以显示样本名
  labs(
    title = "PCA Analysis of GSE29221",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition"
  ) +
  theme_bw() +  # 使用一个简洁的主题
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )



# 假设 pca_data 和 percent_var 已经存在于你的R环境中

# 加载必要的包
library(ggplot2)
library(ggpubr)

# 绘制PCA图 (PC1 vs PC2)，风格与示例代码一致
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  
  # ① 添加置信椭圆（按分组填充，透明度0.2，去掉边框）
  stat_ellipse(
    aes(fill = condition),
    type = "norm",
    alpha = 0.2,
    geom = "polygon",
    level = 0.99,
    colour = NA  # 关键修改1：设置colour=NA去除椭圆边框
  ) +
  
  # ② 绘制样本点（调小尺寸，从3改为1.5）
  geom_point(size = 1.5, alpha = 0.8) +  # 关键修改2：样本点大小调小（可根据需求改1/2等）
  
  # ③ 添加参考虚线（x=0和y=0）
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  # ④ 标签与标题设置（自动提取方差解释率）
  labs(
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)"),
    color = "Condition",
    shape = "Condition",
    fill = "Condition",
    title = "PCA Analysis of GSE29221"
  ) +
  
  # ⑤ 自定义颜色、形状，并设置图例标签
  scale_color_manual(
    values = c("Disease" = "#e74c3c", "Control" = "#1f77b4"),
    labels = c("Disease", "Control")
  ) +
  scale_shape_manual(
    values = c("Disease" = 17, "Control" = 16), # 17=三角形, 16=实心圆
    labels = c("Disease", "Control")
  ) +
  scale_fill_manual(
    values = c("Disease" = "#e74c3c", "Control" = "#1f77b4"),
    labels = c("Disease", "Control")
  ) +
  
  # ⑥ 主题调整
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


# 6.差异基因表达分析（limma）----
# 输入参数
logFoldChange <- 0.5
adjustP <- 0.05
conNum1 <- control_num
treatNum1 <- DF_num
#定义四个变量：logFoldChange，adjustP，conNum和treatNum。
#这些变量分别表示：
#logFoldChange：用于筛选差异表达基因的对数倍数变化阈值，它表示正常和帕金森病样本之间的基因表达量的比值的对数，为1。
#adjustP：用于筛选差异表达基因的校正后的p值阈值，它表示差异表达基因的显著性水平，经过多重检验校正，为0.05。
#conNum：表示正常样本的数量，为normal_num变量的值。
#treatNum：表示tumor样本的数量，为tumor_num变量的值。
#这些变量可以用于后续的差异表达基因分析和热图绘制。
modType1 <- factor(
  c(rep("DF", treatNum1), 
    rep("control", conNum1)),
  levels = c("DF", "control")  # 显式指定水平顺序
)
design1 <- model.matrix(~0 + modType1)  # 直接用modType1（已显式控制水平）
colnames(design1) <- c("DF", "control")
cont.matrix1 <- makeContrasts(DF-control, levels = design1)

#构建 对比矩阵，进行比较肿瘤在前正常在后
#modType1 <- c(rep("Diabetic_peripheral_neuropathy", treatNum1), rep("diabetes_mellitus", conNum1))
#design1 <- model.matrix(~0 + factor(modType1))
#colnames(design1) <- c("Diabetic_peripheral_neuropathy", "diabetes_mellitus")
#它的作用是构建一个线性模型和一个对比矩阵，用于差异表达基因分析。
#它先用c和rep函数创建一个向量，命名为modType，它包含了每个样本的分组信息，
#分为con（正常）和tumor（肿瘤）。
#然后，它用model.matrix函数根据modType向量创建一个设计矩阵，
#命名为design，它表示每个样本属于哪个分组。
#用colnames函数给设计矩阵的列重新命名为con和treat。
#它用makeContrasts函数根据设计矩阵创建一个对比矩阵，命名为cont.matrix，它表示要比较的两个分组之间的差异，即treat-con。
# 直接对整个矩阵转换



#模型拟合
fit.1 <- lmFit(out_exp_counts, design1)     # 线性模型拟合
#用lmFit函数对out_exp_counts矩阵中的表达量数据
#和design矩阵中的实验设计进行线性模型拟合，然后赋值给fit变量。
#lmFit函数是属于limma包的。
fit2.1 <- contrasts.fit(fit.1, cont.matrix1)
#用contrasts.fit函数对fit对象中的线性模型进行对比分析，
#使用cont.matrix矩阵中的对比矩阵，然后赋值给fit2变量。
#contrasts.fit函数也是属于limma包的，它是用来设置对比系数的。
fit2.1 <- eBayes(fit2.1)         # 贝叶斯模型拟合
#这一行是用eBayes函数对fit2对象中的对比分析结果
#进行经验贝叶斯方法（EBM）的校正，然后重新赋值给fit2变量。
#eBayes函数也是属于limma包的，它是用来计算调整后的p值和log-odds值的。



#【所有的】差异基因(结果输出)
res_DEG <- topTable(fit2.1, adjust = 'fdr', number = Inf)
#提取差异表达基因（DEG）的结果，然后赋值给res_DEG变量。
#topTable函数是属于limma包的，
#。adjust = 'fdr'参数是指定用假阳性发现率（FDR）来校正p值
#以控制多重检验的错误。
#number = Inf参数是指定返回所有的基因，而不是默认的前6个。
write.csv(res_DEG, file = "res_DEG.csv")


#【显著】差异表达的基因(结果输出)
res_DEG_sig <- res_DEG[with(res_DEG, (abs(logFC) > 1 & adj.P.Val < 0.05 )), ]
res_DEG_sig1 <- res_DEG[with(res_DEG, (abs(logFC) > 0.5 & adj.P.Val < 0.05 )), ]
#用with函数从res_DEG数据框中筛选出满足条件的基因，
res_DEG_up_sig1 <- res_DEG[with(res_DEG, (logFC > 0.5 & adj.P.Val < 0.05 )), ]

res_DEG_down_sig1 <- res_DEG[with(res_DEG, (logFC < -0.5 & adj.P.Val < 0.05 )), ]






#然后赋值给res_DEG_sig变量。
#条件是基因的对数倍数变化（logFC）的绝对值大于logFoldChange变量的值，
#且校正后的p值（adj.P.Val）小于adjustP变量的值
#这些变量应该是在之前定义过的。
res_DEG_sig <- rbind(id = colnames(res_DEG_sig), res_DEG_sig)

#用rbind函数把res_DEG_sig数据框的列名作为一行添加到数据框的最前面，
#并用"id"作为第一列的名字，然后重新赋值给res_DEG_sig变量。
res_DEG_sig <- res_DEG_sig[-1, ]
#用负索引去掉第一列，也就是"id"列，然后重新赋值给res_DEG_sig变量。
write.csv(res_DEG_sig, file = "res_DEG_sig.csv")
write.table(res_DEG_sig, file = "res_DEG_sig.txt", sep = "\t", quote = F, col.names = F)
#把res_DEG_sig变量中的数据写入一个名为"res_DEG_sig.txt"的文件中，
#以制表符分隔，并且不加引号和列名。
#txt文件方便分析
#csv方便自己查看

sig_genes <- rownames(res_DEG_sig)

# 2. 筛选expr_anno3中存在的基因（避免NA/空值）
sig_genes_filtered <- sig_genes[sig_genes %in% rownames(expr_anno3)]

# 3. 提取表达量矩阵（核心步骤）
immu <- expr_anno3[sig_genes_filtered, ]
## 方式2：四舍五入后转整数（推荐有小数的标准化表达量）
# 步骤1：先将immu转为矩阵（无论原始是数据框/矩阵，统一格式）
immu_mat <- as.matrix(immu)

# 步骤2：四舍五入后转整数（核心修正：先矩阵操作，再转整数）
immu_int <- round(immu_mat)  # 对矩阵整体四舍五入（保留维度）
mode(immu_int) <- "integer"  # 直接将矩阵类型改为整数型（避免list问题）
write.table(
  immu_int,
  file = "immu_integer.txt",
  sep = "\t",
  row.names = TRUE,    # 保留基因名
  col.names = NA,      # 列名对齐
  quote = FALSE        # 无多余引号
)
#处理整合应激反应基因
ISRgene1 <- read.csv("ISRgene.csv", stringsAsFactors = FALSE)
ISRgene2 <- read.csv("ISRgene2.csv", stringsAsFactors = FALSE)
#合并基因
ISRgene <- c(ISRgene1$gene, ISRgene2$gene)
#去重
ISRgene <- unique(c(ISRgene1$gene, ISRgene2$gene))
ISRgene <- data.frame(
  gene = unique(c(ISRgene1$gene, ISRgene2$gene))  # 若保留重复则去掉unique()
)
# 提取ISR基因（去重）
ISRgene_ven <- unique(ISRgene$gene)  # 替换"gene"为实际列名

# 提取显著差异基因（res_DEG_sig中的基因，去重）
deg_gene <- unique(rownames(res_DEG_sig1))  # 替换"gene"为实际列名



#计算交集
common_genes <- intersect(ISRgene_ven, deg_gene)
cat("交集基因数量：", length(common_genes), "\n")
# 1. 筛选出common_genes中存在于expr_anno3行名的基因（避免提取时报错）
common_genes_valid <- common_genes[common_genes %in% rownames(expr_anno3)]

# 2. 提取核心基因的表达矩阵
expr_core <- expr_anno3[common_genes_valid, ]

# 将核心基因表达矩阵输出为txt文件（保存到当前工作目录）
write.table(
  x = expr_core,                # 要输出的核心基因表达矩阵
  file = "core_genes_expr.txt", # 输出文件名称（可自定义，如改路径："D:/data/core_genes_expr.txt"）
  sep = "\t",                   # 分隔符设为制表符（\t），避免数据错位
  row.names = TRUE,             # 保留行名（核心基因名）
  col.names = NA,               # 让列名（样本名）与数据对齐（避免行名占用第一列导致列名偏移）
  quote = FALSE                 # 不添加引号（避免基因名/数值被引号包裹）
)

# 绘制韦恩图（同上）
library(VennDiagram)
library(grid)

venn_plot <- draw.pairwise.venn(
  area1 = length(ISRgene_ven),
  area2 = length(deg_gene),
  cross.area = length(common_genes),
  category = c("ISR Genes", "DEG Genes"),
  fill = c("#377eb8", "#ff7f00"),
  lty = "blank",
  col = "black",
  cat.col = c("#377eb8", "#ff7f00"),
  cat.cex = 1.2,
  main = "ISR Genes vs DEG Genes"
)

grid.draw(venn_plot)
write.csv(common_genes, file = "common_genes.csv", row.names = TRUE)


ISRgene_ven<- as.character(ISRgene_ven)
deg_gene<- as.character(deg_gene)
length(deg_gene)
# 3. 绘制总差异基因的韦恩图

total_venn <- venn.diagram(
  x = list(
    "ISR related genes" = ISRgene_ven,  # 第一组：GSE95849所有差异基因
    "DEGs" = deg_gene   # 第二组：GSE80178所有差异基因
  ),
  filename = "total_deg_venn.tiff",  # 输出文件
  col = "transparent",              # 圆圈边框透明
  fill = c("#4DBBD5FF", "#E64B35FF"),  # 填充色（蓝色和红色）
  alpha = 0.5,                      # 透明度
  label.col = "black",              # 数量标签颜色
  cex = 1.5,                        # 数量标签大小（放大更清晰）
  cat.col = c("#4DBBD5FF", "#E64B35FF"),  # 数据集名称颜色
  cat.cex = 1.2,                    # 名称字体大小
  cat.pos = c(0, 0),               
  cat.dist = c(0.05, 0.05) ,# 名称与圆圈的距离
  scaled = FALSE  # 固定圆圈大小一致（关键参数）
)




BiocManager::install("ReactomePA", dependencies = TRUE, force = TRUE)
library(clusterProfiler)   # 用于基因ID转换
library(ReactomePA)       # 用于GSEA分析
library(org.Hs.eg.db)     # 人类基因注释包（若研究其他物种，需替换为对应包，如小鼠用org.Mm.eg.db）
library(dplyr)            # 数据处理辅助

GSEAanalyse <- data.frame(
  symbol = rownames(res_DEG),
  logfc = res_DEG$logFC,
  stringsAsFactors = FALSE  # 避免因子类型干扰
)

# 执行ID转换（fromType=“SYMBOL”，toType=“ENTREZID”）
gene_entrez <- bitr(
  GSEAanalyse$symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = FALSE  # 保留未匹配到的基因（后续过滤）
)

# 合并原始logFC与ENTREZ ID
# 先修改GSEAanalyse的列名（symbol -> SYMBOL）
colnames(GSEAanalyse)[colnames(GSEAanalyse) == "symbol"] <- "SYMBOL"

# 再合并（此时by参数为"SYMBOL"，两个数据框都有该列）
gene_merged <- merge(GSEAanalyse, gene_entrez, by = "SYMBOL")

# 过滤掉未匹配到ENTREZ ID的基因（NA行）
gene_merged <- gene_merged[!is.na(gene_merged$ENTREZID), ]
GSEAanalyse<- gene_merged
# 构建命名向量
gene_list <- gene_merged$logfc
names(gene_list) <- gene_merged$ENTREZID

# 按logFC从大到小排序（GSEA会基于此排序分析通路分布）
gene_list <- sort(gene_list, decreasing = TRUE)
#
library(ReactomePA)

# 获取人类的 Reactome 通路（返回一个列表，key 是通路ID，value 是通路包含的基因ENTREZID）
library(ReactomePA)

# 获取人类的 Reactome 通路（species 参数指定物种，返回通路-基因列表）
gsea_result_more <- gsePathway(
  geneList = gene_list,       # 仅保留核心参数：命名向量（名字=ENTREZID，值=logFC，已排序）
  pvalueCutoff = 1,        # 显著性阈值
  verbose = FALSE             # 关闭冗余输出
)
gseaplot(gsea_result, geneSetID = gsea_result@result$ID[1])
print(gsea_result)





# 1. 提取包含Description的完整GSEA结果表格
gsea_table <- gsea_result_more@result

# 2. 仅提取Description列（若只需要这一列）
desc_table <- data.frame(
  通路ID = gsea_result_more@result$ID,       # 可选：补充通路ID，方便对应
  通路名称 = gsea_result_more@result$Description # 核心的Description列
)

# 3. 输出完整结果表格到CSV文件（保存到当前工作目录）
write.csv(
  gsea_table, 
  file = "GSEA_Reactome_完整结果表格.csv", # 文件名自定义
  row.names = FALSE # 不保留R的行号，表格更整洁
)

# 4. 仅输出Description列到CSV
write.csv(
  desc_table, 
  file = "GSEA_通路名称列表.csv", 
  row.names = FALSE
)
# 1. 加载必需包
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

# 2. 提取前3+后3通路ID（共6个，确保是长度为6的向量）
all_pathway_ids <- gsea_result@result$ID
target_ids <- c(head(all_pathway_ids, 3), tail(all_pathway_ids, 3))

# 3. 绘制6条曲线同图的GSEA图（匹配你的示例）
p_final <- gseaplot2(
  x = gsea_result,                  # 原始GSEA结果（不替换pvalue列，保留p和adj.p）
  geneSetID = target_ids,           # 一次性传入6个通路ID → 同图显示6条曲线
  title = "GSEA Enriched Pathways",
  pvalue_table = TRUE,              # 显示P值表格（默认会包含pvalue和p.adjust）
  color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FF5F99"), # 6条曲线的区分色
  subplots = c(1, 2, 3),            # 显示完整组件：排序曲线+富集分数+基因位置
  base_size = 6                    # 缩小字体，适配多曲线+表格
)

# 调整rel_heights参数，减小表格高度
p_final <- gseaplot2(
  x = gsea_result,
  geneSetID = target_ids,
  title = "GSEA Enriched Pathways",
  pvalue_table = TRUE,
  color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FF5F99"),
  subplots = c(1, 2, 3),  # 显示三部分：排序曲线+富集分数+基因位置
  rel_heights = c(2, 0.7, 0.3),  # 调整三部分的高度比例
  base_size = 6
)
# 5. 显示+保存（与示例图一致的高清格式）
print(p_final)

ggsave(
  filename = "GSEA_6Pathways_SinglePlot.pdf",
  plot = p_final,
  width = 14,  # 宽度适配6条曲线+右侧表格
  height = 8,
  dpi = 300
)

#lasso分析，以及svm-RFE分析
library(glmnet)
library(e1071)
library(kernlab)
library(caret)
# 假设out_exp_counts1是包含核心基因名称的向量
# 筛选行名属于total_common_genes的数据
hub_gene1 <- out_exp_counts[rownames(out_exp_counts) %in% common_genes, ]
hub_gene <- res_DEG[rownames(res_DEG) %in% common_genes, ]
group1 <- group_info3[1:24, ]
# 根据实际分组名称（"Non-diabetic"和"Diabetic"）进行0/1编码
group1$group <- ifelse(
  grepl("Non-diabetic", group1$group),  # 匹配含"Non-diabetic"的分组（非糖尿病组）
  0,  # 非糖尿病组编码为0
  1   # 其余为糖尿病组，编码为1
)
group1 = as.matrix(group1)
hub_gene1 = as.matrix(hub_gene1)
hub_gene1 <- t(hub_gene1)  # 原hub_gene1是基因×样本，转置后为样本×基因
group1 <- as.data.frame(group1)
group1 <- group1[order(group1$group), ]
# 提取group1的Id列作为排序依据（假设Id列的值与hub_gene1的行名对应）
sort_order <- group1$Id

# 按sort_order的顺序重新排列hub_gene1的行
# 仅保留hub_gene1中存在于sort_order中的行（自动忽略不匹配的行）
hub_gene1 <- hub_gene1[sort_order[sort_order %in% rownames(hub_gene1)], ]
group1 = as.matrix(group1)

y <- as.numeric(group1[, "group"])  # 确保y是数值型向量（0=疾病，1=健康）
# 检查数据维度是否匹配（样本数需一致）
cat("样本数：", nrow(hub_gene1), "\n")  # 应输出18
cat("基因数：", ncol(hub_gene1), "\n")
cat("分组维度：", length(y), "\n")  # 应输出18
# 2. 交叉验证选择最优lambda（LASSO核心步骤）

set.seed(123)  # 保证结果可重复
cv_lasso <- cv.glmnet(
  x = hub_gene1, 
  y = y, 
  family = "binomial",  # 二分类问题（疾病/健康）
  alpha = 1,  # alpha=1表示LASSO；0<alpha<1为弹性网络
  nfolds = 5,  # 10折交叉验证（样本量小，可适当减少折数如5折）
  type.measure = "deviance"  # 二分类常用偏差作为评价指标
)

# 可视化交叉验证结果（查看lambda与误差的关系）
plot(cv_lasso)
title("LASSO交叉验证曲线")

# 3. 提取最优参数及模型结果
# 最优lambda（两种常用选择）
lambda_min <- cv_lasso$lambda.min  # 最小交叉验证误差对应的lambda
lambda_1se <- cv_lasso$lambda.1se  # 1个标准差内的最小误差对应的lambda（更简洁的模型）

# 用最优lambda拟合最终LASSO模型
lasso_model <- glmnet(
  x = hub_gene1, 
  y = y, 
  family = "binomial", 
  alpha = 1, 
  lambda = lambda_1se  # 或使用lambda_1se
)


# 提取非零系数（LASSO筛选出的关键基因）
coef_lasso <- coef(lasso_model)
selected_genes <- rownames(coef_lasso)[coef_lasso@i + 1]  # 非零系数对应的基因名
coef_values <- coef_lasso@x  # 对应的系数值

# 输出筛选结果
cat("LASSO筛选出的关键基因（共", length(selected_genes), "个）：\n")
print(data.frame(基因名 = selected_genes, 系数 = coef_values))
selected_genes
best_features
common_elements <- intersect(selected_genes, best_features)
print(common_elements)

library(glmnet)


# x：样本×特征矩阵（行=样本，列=基因/特征）
# y：因变量（如二分类的0/1，连续型的数值）
lasso_fit <- glmnet(
  x = hub_gene1, 
  y = y, 
  family = "binomial",  # 二分类（疾病/健康），连续型用"gaussian"
  alpha = 1  # alpha=1 → LASSO；0<alpha<1 → 弹性网络
)

plot(lasso_fit, xvar = "lambda", label = TRUE)  



library(glmnet)
library(e1071)
library(kernlab)
library(caret)
library(dplyr)

# 数据准备

group1 <- as.data.frame(group1)
group1 <- group1[order(group1$group), ]
# 提取group1的Id列作为排序依据（假设Id列的值与hub_gene1的行名对应）
sort_order <- group1$Id

# 按sort_order的顺序重新排列hub_gene1的行
# 仅保留hub_gene1中存在于sort_order中的行（自动忽略不匹配的行）
hub_gene1 <- hub_gene1[sort_order[sort_order %in% rownames(hub_gene1)], ]
y <- group1$group

# 关键：正确转换为因子（保持原始0/1编码）
y <- as.factor(y)
levels(y) <- c("Class0", "Class1")  # 明确指定水平名称

set.seed(123)
x = hub_gene1



set.seed(123)
ctrl <- rfeControl(
  functions = caretFuncs,
  method = "LOOCV",
  verbose = FALSE
)

# 定义自定义的SVM参数网格
svm_grid <- expand.grid(C = c( 0.01))  # 更小的C值，更强的正则化

svm_rfe_strict <- rfe(
  x = x,
  y = y,
  sizes <- seq(from = 1, by = 2, length.out = 31),
  rfeControl = ctrl,
  method = "svmLinear",
  tuneGrid = svm_grid,
  preProcess = c("center", "scale")
)
print(svm_rfe_strict)

plot(svm_rfe_strict, type = c("g", "o"), 
     main = "SVM-RFE: Feature Number vs Accuracy")

# 提取最佳特征
best_features <- predictors(svm_rfe_strict)
cat("最佳特征数量:", length(best_features), "\n")
cat("最佳特征:", best_features, "\n")




library(pROC)

common_elements





hub_gene1 <- out_exp_counts[rownames(out_exp_counts) %in% common_genes, ]
hub_gene <- res_DEG[rownames(res_DEG) %in% common_genes, ]
group1 <- group_info3[1:24, ]
# 根据实际分组名称（"Non-diabetic"和"Diabetic"）进行0/1编码
group1$group <- ifelse(
  grepl("Non-diabetic", group1$group),  # 匹配含"Non-diabetic"的分组（非糖尿病组）
  0,  # 非糖尿病组编码为0
  1   # 其余为糖尿病组，编码为1
)
group1 = as.matrix(group1)
hub_gene1 = as.matrix(hub_gene1)
hub_gene1 <- t(hub_gene1)  # 原hub_gene1是基因×样本，转置后为样本×基因
group1 <- as.data.frame(group1)
group1 <- group1[order(group1$group), ]
# 提取group1的Id列作为排序依据（假设Id列的值与hub_gene1的行名对应）
sort_order <- group1$Id

# 按sort_order的顺序重新排列hub_gene1的行
# 仅保留hub_gene1中存在于sort_order中的行（自动忽略不匹配的行）
hub_gene1 <- hub_gene1[sort_order[sort_order %in% rownames(hub_gene1)], ]
group1 = as.matrix(group1)

y <- as.numeric(group1[, "group"])  # 确保y是数值型向量（0=疾病，1=健康）


# 计算每个基因的AUC
auc_values <- sapply(colnames(hub_gene1), function(gene) {
  roc_obj <- roc(y, hub_gene1[, gene])
  auc(roc_obj)
})

# 排序并查看
auc_values <- sort(auc_values, decreasing = TRUE)
print(auc_values)
# 绘制AUC最高的基因的ROC曲线
common_elements
best_gene <- names(auc_values)[1]
roc_obj <- roc(y, hub_gene1[, best_gene])
plot(roc_obj, main = paste("ROC Curve for", best_gene))

# 绘制该基因的表达分布
library(ggplot2)
df <- data.frame(Expression = x[, best_gene], Group = y)
ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  labs(title = paste("Expression of", best_gene))




#在验证集合中绘制roc曲线

hub_gene_test <- expr_anno2[rownames(expr_anno2) %in% common_genes, ]
group_test <- group_info2[1:9, ]
# 修正匹配逻辑：精准识别健康非糖尿病组
group_test$group <- ifelse(
  grepl("Diabetic Foot Ulce", group_test$group),  # 匹配含"healthy non-diabetic"的分组（非糖尿病组）
  0,  # 非糖尿病组编码为0
  1   # 其余为糖尿病组（含"diabetic"），编码为1
)
# 1. 提取group_test$id的唯一值（避免重复匹配，提升效率）
target_ids <- unique(group_test$Id)

# 2. 筛选列名属于target_ids的列（drop=FALSE确保结果仍为数据框，避免单列时变向量）
hub_gene_filtered <- hub_gene_test[, colnames(hub_gene_test) %in% target_ids, drop = FALSE]
hub_gene_test<-hub_gene_filtered
group_test = as.matrix(group_test)
hub_gene_test = as.matrix(hub_gene_test)
hub_gene_test <- t(hub_gene_test)  # 原hub_gene1是基因×样本，转置后为样本×基因
group_test <- as.data.frame(group_test)
group_test <- group_test[order(group_test$group), ]
# 提取group1的Id列作为排序依据（假设Id列的值与hub_gene1的行名对应）
sort_order <- group_test$Id

# 按sort_order的顺序重新排列hub_gene1的行
# 仅保留hub_gene1中存在于sort_order中的行（自动忽略不匹配的行）
hub_gene_test <- hub_gene_test[sort_order[sort_order %in% rownames(hub_gene_test)], ]
group_test = as.matrix(group_test)

y <- as.numeric(group_test[, "group"])  # 确保y是数值型向量（0=疾病，1=健康）




# 计算每个基因的AUC
auc_values <- sapply(colnames(hub_gene_test), function(gene) {
  roc_obj <- roc(y, hub_gene_test[, gene])
  auc(roc_obj)
})

# 排序并查看
auc_values <- sort(auc_values, decreasing = TRUE)
print(auc_values)
# 绘制AUC最高的基因的ROC曲线
common_elements
best_gene <- names(auc_values)[20]
roc_obj <- roc(y, hub_gene_test[, best_gene])
plot(roc_obj, main = paste("ROC Curve for", best_gene))

library(pROC)

# ========== 步骤1：数据筛选（同之前） ==========
target_genes <- common_elements
hub_gene_subset <- hub_gene_test[, colnames(hub_gene_test) %in% target_genes, drop = FALSE]

# ========== 步骤2：计算ROC对象和AUC（同之前） ==========
roc_list <- list()
auc_list <- c()
for (gene in colnames(hub_gene_subset)) {
  roc_obj <- roc(response = y, predictor = hub_gene_subset[, gene], quiet = TRUE)
  roc_list[[gene]] <- roc_obj
  auc_list[[gene]] <- round(auc(roc_obj), 3)
}

# ========== 步骤3：绘制ROC曲线（关闭图内AUC标注） ==========
colors <- c("#E64B35", "#4DBBD5", "#FFCC00")
ltys <- c(1, 2, 3)

# 初始化画布：关闭print.auc
first_gene <- names(roc_list)[1]
plot(roc_list[[first_gene]], 
     col = colors[1], lwd = 2,
     main = "ROC Curves for Hub Genes",
     xlab = "1 - Specificity (False Positive Rate)",
     ylab = "Sensitivity (True Positive Rate)",
     print.auc = FALSE,  # 关键修改：不显示图内AUC
     auc.polygon = FALSE
)

# 叠加其他曲线：同样关闭print.auc
for (i in 2:length(roc_list)) {
  gene <- names(roc_list)[i]
  plot(roc_list[[gene]], 
       add = TRUE,
       col = colors[i], lwd = 2, lty = ltys[i],
       print.auc = FALSE,  # 关键修改：不显示图内AUC
       auc.polygon = FALSE
  )
}

# ========== 步骤4：参考线+图例（图例保留AUC） ==========
legend("bottomright",
       legend = paste0(names(auc_list), " (AUC = ", auc_list, ")"),
       col = colors, lwd = 2, lty = ltys,
       bty = "n", cex = 0.9
)



# 1. 直接提取最佳基因的表达值（这是一个向量）
#    然后直接对这个向量进行索引，获取第7到12个元素
selected_samples <- c("GSM2114237", "GSM2114238", "GSM2114239", "GSM2114240", "GSM2114241", "GSM2114242")
expression_vector <- expr_anno2[best_gene, selected_samples]
# 3. 创建分组向量
new_groups <- rep(c("Disease", "Control"), each = 3)

# 4. 创建用于绘图的数据框
df <- data.frame(
  Expression = as.numeric(expression_vector),  # 确保转换为数值向量
  Group = factor(new_groups, levels = c("Control", "Disease"))  # 设置因子水平，Control在前
)












expression_original <- 2^as.numeric(expression_vector)  # 如果数据是log2转换过的

# 4. 创建用于绘图的数据框
df <- data.frame(
  Expression = expression_original,  # 使用转换后的表达值
  Group = factor(new_groups, levels = c("Control", "Disease"))
)

# 5. 添加统计检验
# 使用t检验
t_test_result <- t.test(Expression ~ Group, data = df)
p_value <- t_test_result$p.value

# 或者使用Wilcoxon秩和检验（对于非正态分布数据更合适）
wilcox_result <- wilcox.test(Expression ~ Group, data = df)
p_value_wilcox <- wilcox_result$p.value
# 5. 绘制箱线图
library(ggplot2)

ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("Control" = "#00BFC4", "Disease" = "#F8766D")) +
  labs(
    title = paste("Expression of", best_gene),
    x = "Patient Group",
    y = "Gene Expression Level"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )


write.csv(expr_anno2, file = "expr_anno2_table.csv", row.names = TRUE, quote = FALSE)






#随机森林2.0
# 设置随机种子，确保结果可重复
set.seed(123)
# 加载包
library(randomForest)  # 随机森林模型
library(pROC)          # 绘制ROC曲线
library(caret)         # 计算混淆矩阵等评估指标
# 构建随机森林模型（分类任务）
# 筛选行名属于total_common_genes的数据
hub_gene1 <- out_exp_counts[rownames(out_exp_counts) %in% common_genes, ]
hub_gene <- res_DEG[rownames(res_DEG) %in% common_genes, ]
group1 <- group_info3[1:24, ]
# 根据实际分组名称（"Non-diabetic"和"Diabetic"）进行0/1编码
group1$group <- ifelse(
  grepl("Non-diabetic", group1$group),  # 匹配含"Non-diabetic"的分组（非糖尿病组）
  0,  # 非糖尿病组编码为0
  1   # 其余为糖尿病组，编码为1
)
group1 = as.matrix(group1)
hub_gene1 = as.matrix(hub_gene1)
hub_gene1 <- t(hub_gene1)  # 原hub_gene1是基因×样本，转置后为样本×基因
group1 <- as.data.frame(group1)
group1 <- group1[order(group1$group), ]
# 提取group1的Id列作为排序依据（假设Id列的值与hub_gene1的行名对应）
sort_order <- group1$Id

# 按sort_order的顺序重新排列hub_gene1的行
# 仅保留hub_gene1中存在于sort_order中的行（自动忽略不匹配的行）
hub_gene1 <- hub_gene1[sort_order[sort_order %in% rownames(hub_gene1)], ]
group1 = as.matrix(group1)

y <- as.numeric(group1[, "group"])  # 确保y是数值型向量（0=疾病，1=健康）
# 检查数据维度是否匹配（样本数需一致）
y<- as.factor(y)

rf_model <- randomForest(
  x = hub_gene1,               # 特征矩阵（样本×基因）
  y = y,             # 响应变量（因子类型）
  ntree = 500,              # 树的数量（可调整，默认500）
  mtry = sqrt(ncol(hub_gene1)),# 每次分裂随机选择的特征数（分类问题默认sqrt(p)）
  importance = TRUE,        # 计算特征重要性
  proximity = FALSE         # 不计算样本相似度（节省内存）
)

# 查看模型基本结果
print(rf_model)
# 输出中包含混淆矩阵（训练集的预测结果）、OOB（袋外数据）误差等
# 1. 计算OOB预测概率（袋外数据的预测概率，更稳健）
# 提取完整的预测概率矩阵
pred_prob <- predict(rf_model, type = "prob")
# 查看列名（这一步能帮你找到正确的下标名称）
colnames(pred_prob)
oob_prob <- pred_prob[, "1"]

# 2. 绘制ROC曲线并计算AUC（用OOB概率）
roc_obj <- roc(
  response = y,      # 真实标签
  predictor = oob_prob,     # 预测概率
  levels = levels(y) # 指定标签顺序
)

# 输出AUC值
cat("OOB数据的AUC值：", auc(roc_obj), "\n")

# 3. 可视化ROC曲线
plot(roc_obj,
     main = "随机森林模型ROC曲线",
     col = "red", lwd = 2,
     xlab = "1 - 特异性", ylab = "敏感性")
abline(a = 0, b = 1, lty = 2, col = "gray")  # 随机猜测的参考线
# 提取特征重要性（两种指标：节点不纯度减少度和准确率减少度）
gene_importance <- importance(rf_model)

# 按"MeanDecreaseAccuracy"（准确率减少度）排序，查看最重要的基因
gene_importance_sorted <- gene_importance[order(-gene_importance[, "MeanDecreaseGini"]), ]

# 打印前10个最重要的基因
cat("核心基因重要性（前10）：\n")
print(head(gene_importance_sorted, 10))

# 可视化特征重要性（前15个基因）
varImpPlot(rf_model,
           main = "核心基因重要性",
           n.var = 8)  # 显示前15个


plot(rf_model, 
     main = "随机森林: 树数量 vs 错误率",
     lwd = 2)

# 添加图例
legend("topright", 
       legend = colnames(rf_model$err.rate),
       col = 1:ncol(rf_model$err.rate),
       lty = 1,
       lwd = 2)





##韦恩图
# 加载必要的包
library(VennDiagram)

# 定义三个基因列表
lasso_genes <- c("ATF3", "IRAK2", "MTHFR", "THBS1", "VEGFA")
best_features_genes <- c("THBS1", "VEGFA", "COPS4", "CDK9", "JAK1", "USP10", "H1-2", "MTHFR", "IFNAR2", "IRAK1", "GBP2", "IFRD1", "NOP53")
rf_genes <- c("VEGFA", "THBS1", "COPS4", "IRAK1", "GBP2", "USP10", "IFNAR2", "CDK9", "JAK1", "MTHFR")



library(VennDiagram)

# 创建自定义布局的正三角形韦恩图
venn.plot <- draw.triple.venn(
  area1 = length(lasso_genes),
  area2 = length(best_features_genes),
  area3 = length(rf_genes),
  n12 = length(intersect(lasso_genes, best_features_genes)),
  n23 = length(intersect(best_features_genes, rf_genes)),
  n13 = length(intersect(lasso_genes, rf_genes)),
  n123 = length(intersect(intersect(lasso_genes, best_features_genes), rf_genes)),
  category = c("LASSO", "Best Features", "Random Forest"),
  
  # 自定义填充颜色和透明度
  fill = c("lightcoral", "lightgreen", "lightblue"),
  alpha = c(0.5, 0.5, 0.5),
  
  # 设置标签
  cex = 1.5,
  cat.cex = 1.3,
  cat.col = c("darkred", "darkgreen", "darkblue"),
  
  # 设置正三角形布局
  rotation.degree = 0,  # 不旋转
  ind = TRUE,  # 显示独立标签
  
  # 调整图形布局参数以获得正三角形
  # 手动调整位置（可能需要根据实际效果微调）
  direct.area = FALSE,
  euler.d = TRUE,
  scaled = TRUE
)

# 显示图形
grid.newpage()
grid.draw(venn.plot)

# 或者保存
# ggsave("venn_triangle.png", plot = venn.plot, width = 8, height = 8, dpi = 300)
# 计算交集
n12 <- length(intersect(lasso_genes, best_features_genes))
n23 <- length(intersect(best_features_genes, rf_genes))
n13 <- length(intersect(lasso_genes, rf_genes))
n123 <- length(intersect(intersect(lasso_genes, best_features_genes), rf_genes))

# 方法1: 使用grid包手动绘制正圆形
# 这种方法可以完全控制，但更复杂

# 方法2: 使用venn.diagram函数，并确保正圆形
venn.plot <- venn.diagram(
  x = list(
    "LASSO" = lasso_genes,
    "Best Features" = best_features_genes,
    "Random Forest" = rf_genes
  ),
  filename = NULL,
  category.names = c("LASSO", "Best Features", "Random Forest"),
  
  # 设置输出为圆形
  force.unique = TRUE,
  
  # 设置颜色
  fill = c("lightcoral", "lightgreen", "lightblue"),
  alpha = c(0.5, 0.5, 0.5),
  
  # 设置标签
  cex = 1.5,
  cat.cex = 1.3,
  cat.col = c("darkred", "darkgreen", "darkblue"),
  
  # 设置类别标签位置（正三角形排列）
  cat.pos = c(0, 180, 180),  # 0=顶部, 180=底部
  cat.dist = c(0.1, 0.1, 0.1),
  
  # 使用默认布局（正三角形）
  rotation.degree = 0,
  
  # 不缩放，保持正圆形
  scaled = FALSE,
  
  # 设置输出尺寸，确保图形是正方形，这样圆就不会被拉伸
  height = 10,
  width = 10,
  units = "in",
  resolution = 300
)

# 显示图形
grid::grid.draw(venn.plot)

library(ggVennDiagram)
library(ggplot2)

# 创建列表
gene_lists <- list(
  "LASSO" = lasso_genes,
  "Best Features" = best_features_genes,
  "Random Forest" = rf_genes
)

# 创建韦恩图
ggVennDiagram(gene_lists, 
              label_alpha = 0,
              edge_size = 0.8) +
  scale_fill_gradient(low = "white", high = "lightblue") +
  scale_color_manual(values = c("LASSO" = "darkred", 
                                "Best Features" = "darkgreen", 
                                "Random Forest" = "darkblue")) +
  theme(legend.position = "none") +
  labs(title = "Venn Diagram of Gene Lists") +
  theme_void() +
  # 确保图形是正方形，这样圆形就不会被拉伸
  coord_equal()






library(ggVennDiagram)
library(ggplot2)

# 创建列表
gene_lists <- list(
  "LASSO" = lasso_genes,
  "Best Features" = best_features_genes,
  "Random Forest" = rf_genes
)

# 首先获取Venn对象
venn_obj <- Venn(gene_lists)
venn_data <- process_data(venn_obj)

# 绘制韦恩图并自定义每个区域的填充颜色
# 注意：ggVennDiagram 默认基于交集大小填充颜色，我们需要手动覆盖
ggplot() +
  # 绘制集合区域
  geom_sf(aes(fill = name), data = venn_region(venn_data), alpha = 0.5, show.legend = FALSE) +
  # 绘制集合边界
  geom_sf(aes(color = name), data = venn_setedge(venn_data), size = 1, show.legend = FALSE) +
  # 绘制集合标签
  geom_sf_text(aes(label = name), data = venn_setlabel(venn_data), size = 5) +
  # 绘制区域标签（显示数量）
  geom_sf_label(aes(label = count), data = venn_region(venn_data), alpha = 0.8, size = 4.5) +
  # 设置填充颜色
  scale_fill_manual(values = c(
    "LASSO" = "lightcoral",
    "Best Features" = "lightgreen", 
    "Random Forest" = "lightblue"
  )) +
  # 设置边界颜色
  scale_color_manual(values = c(
    "LASSO" = "darkred",
    "Best Features" = "darkgreen", 
    "Random Forest" = "darkblue"
  )) +
  theme_void() +
  coord_equal() +
  labs(title = "Venn Diagram of Gene Lists") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))




library(caret)
set.seed(123)
train_control <- trainControl(method = "cv", number = 10)
rf_cv <- train(
  x = hub_gene1,
  y = y,
  method = "rf",
  trControl = train_control,
  tuneLength = 3
)
print(rf_cv)



rf_model$err.rate  # 查看每棵树的OOB误差
set.seed(123)
train_control_loocv <- trainControl(
  method = "LOOCV",  # 留一法交叉验证
  classProbs = TRUE,
  summaryFunction = twoClassSummary  # 如果需要更详细的指标
)

# 注意：需要将y转换为因子且有意义的水平名
y_factor <- factor(y, levels = c(0, 1), labels = c("Class0", "Class1"))

rf_loocv <- train(
  x = hub_gene1,
  y = y_factor,
  method = "rf",
  trControl = train_control_loocv,
  metric = "Accuracy"  # 或"ROC"、"Sensitivity"等
)

print(rf_loocv)





















logFoldChange <-0.5


gene_up_sig1 <- subset(res_DEG, adj.P.Val < adjustP & logFC > logFoldChange)
#是用subset函数从res_DEG数据框中筛选出显著上调的基因，
#也就是校正后的p值（adj.P.Val）小于adjustP变量的值，
#且对数倍数变化（logFC）大于logFoldChange变量的值的基因，
#然后赋值给gene_up_sig变量。
gene_down_sig1 <- subset(res_DEG, adj.P.Val < adjustP & logFC < (-logFoldChange))
#用subset函数从res_DEG数据框中筛选出显著下调的基因，
#也就是校正后的p值（adj.P.Val）小于adjustP变量的值，
#且对数倍数变化（logFC）小于负的logFoldChange变量的值的基因，
#然后赋值给gene_down_sig变量。
write.table(gene_up_sig1, file = "GSE29221_sig_up.txt", sep = "\t", quote = F, col.names = F)
write.table(gene_down_sig1, file = "GSE29221_sig_down.txt", sep = "\t", quote = F, col.names = F)


#绘制差异表达基因火山图
xMax1 <- max(abs(res_DEG$logFC))
#用max函数计算res_DEG数据框中的对数倍数变化（logFC）的绝对值的最大值，
#然后赋值给xMax变量。
yMax1 <- max(-log10(res_DEG$adj.P.Val))
#用max函数计算res_DEG数据框中的校正后的p值（adj.P.Val）的负对数的最大值，
#然后赋值给yMax变量。
plot(res_DEG$logFC, -log10(res_DEG$adj.P.Val), xlab = expression(log[2]("Fold Change")), ylab = expression(-log[10]("Pvalue")),
     main = "GSE29221", ylim = c(0, yMax1), xlim = c(-xMax1-1, xMax1+1), yaxs = "i", pch = 20, cex = 0.8)

#用plot函数画出差异表达基因（DEG）的火山图，
#横坐标是对数倍数变化（logFC），
#纵坐标是校正后的p值（adj.P.Val）的负对数。
#横坐标标签是"logFC"，纵坐标标签是"-log10(adj.P.Val)“，图标题是"Volcano”。
#纵坐标范围是从0到yMax变量的值，横坐标范围是从负的xMax变量的值到正的xMax变量的值。
#纵坐标轴使用内部样式，点符号使用圆点，点大小使用0.8倍。
diffSub_up1 <- subset(res_DEG, adj.P.Val < adjustP & logFC > logFoldChange)
#用subset函数从res_DEG数据框中筛选出显著上调的基因，
#也就是校正后的p值（adj.P.Val）小于adjustP变量的值，
#且对数倍数变化（logFC）大于logFoldChange变量的值的基因，
#然后赋值给diffSub_up变量。
points(diffSub_up1$logFC, -log10(diffSub_up1$adj.P.Val), pch = 20, col = "#D73027", cex= 0.8)
#用points函数在火山图上添加显著上调的基因的点，
#用红色圆点表示，点大小使用0.8倍。
diffSub_down1 <- subset(res_DEG, adj.P.Val < adjustP & logFC < (-logFoldChange))
#用subset函数从res_DEG数据框中筛选出显著下调的基因，
#也就是校正后的p值（adj.P.Val）小于adjustP变量的值，
#且对数倍数变化（logFC）小于负的logFoldChange变量的值的基因，
#然后赋值给diffSub_down变量。
points(diffSub_down1$logFC, -log10(diffSub_down1$adj.P.Val), pch = 20, col = "#4575B4", cex = 0.8)
#用points函数在火山图上添加显著下调的基因的点，
#用绿色圆点表示，点大小使用0.8倍。
abline(v = 0, lty = 2, lwd = 3)
abline(v = c(-logFoldChange, logFoldChange), lty = 2, col = "gray")
abline(h = -log10(adjustP), lty = 2, col = "gray")
# 添加图例（字体从0.8增至1.0，点大小同步放大）
legend("topright", 
       legend = c("Down-regulated", "Up-regulated", "Non-significant"),
       col = c("#4575B4", "#D73027", "black"),
       pch = c(20, 20, 20),
       cex = 1.0,  # 字体放大
       pt.cex = 2.0,  # 图例中点的大小放大
       bty = "n",
       inset = c(-0.01, 0) )















# 假设res_DEG_sig1.1包含log2FoldChange列，按绝对值降序排序
#绘制差异表达基因热图
class(res_DEG_sig1$logFC)  # 输出可能是 "character" 或 "factor"（非 "numeric"）
res_DEG_sig1$logFC <- as.numeric(as.character(res_DEG_sig1$logFC))
res_DEG_sig1.1 <- res_DEG_sig1[order(-abs(res_DEG_sig1$logFC)), ]

# 2. 提取前100个基因的名称（若不足100个则取全部）
top100_genes <- head(rownames(res_DEG_sig1.1), 100)

hmExp1 <- out_exp_counts_log2_clean[top100_genes, ]
#用out_exp_counts矩阵中的显著差异表达基因（res_DEG_sig）的行名作为索引
#提取出相应的表达量数据，然后赋值给hmExp变量。
Type1 <- c(rep("CONTROL", conNum1), rep("DFU", treatNum1))
#用c函数和rep函数创建一个向量，包含两种类型的样本，
#分别是正常（NOR）和肿瘤（tumor）。
#正常样本的个数由conNum变量指定，
#肿瘤样本的个数由treatNum变量指定。然后赋值给Type变量。
names(Type1) <- colnames(out_exp_counts_log2_clean)
#用colnames函数获取out_exp_counts矩阵的列名，也就是样本的名字，
#然后赋值给Type向量的名字。
Type1 <- as.data.frame(Type1)
#用as.data.frame函数把Type向量转换为一个数据框，然后重新赋值给Type变量。
# 绘制热图
pheatmap(hmExp1, 
         annotation = Type1, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols = F,
         show_colnames = F,
         scale = "row",
         fontsize = 12,
         fontsize_row = 3,
         fontsize_col = 10)
#用pheatmap函数画出hmExp数据框中的表达量数据的热图，
#用Type数据框作为样本的注释，
#用从绿色到黑色再到红色的颜色渐变作为颜色方案，
#不对列进行聚类，不显示列名，对行进行标准化，
#设置字体大小为12，行名字体大小为3，列名字体大小为10。
# 修改颜色参数 + 优化可视化细节
library(pheatmap)

# 定义颜色（建议使用 ggplot 风格的蓝-白-红）
color_heatmap <- colorRampPalette(c("#4575B4", "white", "#D73027"))(50)  # 深蓝到浅红


#试试ggplot配色
pheatmap(
  hmExp1,
  annotation = Type1,
  color = color_heatmap, 
  cluster_cols = FALSE,
  show_colnames = FALSE,
  scale = "row",
  border_color = NA,          # 去除单元格边框线（更简洁）
  fontsize_row = 5 ,           # 行名字体大小（避免重叠）
  fontsize_col = 10,
  cellwidth = 20,  # 固定每个样本的宽度（关键参数）
  annotation_colors = list(   # 同步修改注释颜色（Type）
    Type = c("NOR" = "#619CFF", "tumor" = "#F8766D")  # ggplot默认蓝/红
  )
)

# 1. 打开PNG设备，设置文件名、尺寸和分辨率
png("高分辨率热图？.png", width = 3600, height = 2400, res = 300) 

# 2. 运行你的pheatmap代码
pheatmap(
  hmExp1,
  annotation = Type1,
  color = color_heatmap, 
  cluster_cols = FALSE,
  show_colnames = FALSE,
  scale = "row",
  border_color = NA,
  fontsize_row = 5,
  fontsize_col = 10,
  cellwidth = 20,
  annotation_colors = list(
    Type = c("NOR" = "#619CFF", "tumor" = "#F8766D")
  )
)

# 3. 关闭设备
dev.off()






# 8.Id 转换----
#富集分析是用于确定在特定生物学过程、功能或通路中显著富集的基因集合
#通过富集分析，可以了解这些基因在特定生物学过程中的重要性
#，并推断它们可能扮演的角色。


# 安装加载包和数据
load("genelist.Rdata")
inter_gene <- total_common


library(clusterProfiler)

# ENTREZID：NCBI 提供的编号系统，用来标识染色体上的基因位点。

# clusterProfiler 转换Id
library(org.Hs.eg.db)

id_info3 <- bitr(inter_gene,
                 fromType = "SYMBOL",
                 toType = c("ENTREZID", "GENENAME", "SYMBOL"),
                 OrgDb = org.Hs.eg.db)
#自动去除NA缺失值，因此有部分会匹配失败

#保存id转换数据
save(id_info3, file = "enrich.Rdata")
write.csv(id_info3, file = "id_info3.csv")



# Load packages
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

library(clusterProfiler)  # Core package for enrichment analysis
library(org.Hs.eg.db)    # Human gene annotation package (for mouse use org.Mm.eg.db, for rat use org.Rn.eg.db)
library(dplyr)           # Data processing

# Import table (first column contains gene names)
deg_data <- read.csv("res_DEG_sig.csv", header = TRUE)

# Extract gene names from first column (Symbol format)
gene_symbol <- deg_data[, 1]
# Convert gene Symbol to ENTREZID (note: replace OrgDb with appropriate species)
gene_entrez <- bitr(
  geneID = gene_symbol,
  fromType = "SYMBOL",   # Original ID type
  toType = "ENTREZID",   # Target ID type
  OrgDb = org.Hs.eg.db   # Species annotation package
) %>% na.omit()  # Remove genes that failed conversion

# Extract converted ENTREZID list
gene_list <- gene_entrez$ENTREZID
# Perform BP, CC, MF enrichment analysis simultaneously

mf_result <- enrichGO(
  gene = all_gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "MF",            # Analyze only Molecular Function
  pAdjustMethod = "fdr", # Adjustment method
  qvalueCutoff = 0.05,   # Significance threshold
  readable = TRUE        # Display gene Symbol in results (easier to read)
)
bp_result <- enrichGO(
  gene = all_gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "bp",            # Analyze only Molecular Function
  pAdjustMethod = "fdr", # Adjustment method
  qvalueCutoff = 0.05,   # Significance threshold
  readable = TRUE        # Display gene Symbol in results (easier to read)
)
cc_result <- enrichGO(
  gene = all_gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "cc",            # Analyze only Molecular Function
  pAdjustMethod = "fdr", # Adjustment method
  qvalueCutoff = 0.05,   # Significance threshold
  readable = TRUE        # Display gene Symbol in results (easier to read)
)
# Extract BP, CC, MF results separately
go_bp <- go_result %>% filter(ONTOLOGY == "BP")
go_cc <- go_result %>% filter(ONTOLOGY == "CC")
go_mf <- go_result %>% filter(ONTOLOGY == "mf")
kegg_result <- enrichKEGG(
  gene = all_gene_list,
  organism = "hsa",      # KEGG code for species (human=hsa, mouse=mmu, rat=rn)
  pvalueCutoff = 0.05,   # Significance threshold
  qvalueCutoff = 0.05
)
# Load visualization dependency packages
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Categorized bar plot of GO enrichment results
barplot(go_result, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scale = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Dot plot of KEGG enrichment results
dotplot(kegg_result)
# Display only BP bubble plot
bp_dotplot <- dotplot(
  bp_result, # Previously extracted BP results
  showCategory = 10, # Display top 20 significant BP terms
  font.size = 10,
  title = "GO Biological Process (bp) Enrichment Dot Plot",
  color = "p.adjust"
) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  ) +
  labs(x = "Gene Ratio", y = "BP Term", color = "Q-value")

print(bp_dotplot)

cc_dotplot <- dotplot(
  cc_result, # Previously extracted BP results
  showCategory = 10, # Display top 20 significant BP terms
  font.size = 10,
  title = "GO Cellular Component (cc) Enrichment Dot Plot",
  color = "p.adjust"
) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  ) +
  labs(x = "Gene Ratio", y = "BP Term", color = "Q-value")

print(cc_dotplot)
mf_dotplot <- dotplot(
  mf_result, # Previously extracted BP results
  showCategory = 10, # Display top 20 significant BP terms
  font.size = 10,
  title = "GO Molecular Function (mf) Enrichment Dot Plot",
  color = "p.adjust"
) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  ) +
  labs(x = "Gene Ratio", y = "BP Term", color = "Q-value")

print(mf_dotplot)
# Similarly can plot CC/MF: replace go_bp with go_cc/go_mf
kegg_dotplot <- dotplot(
  kegg_result, # Previous KEGG enrichment results
  showCategory = 20, # Display top 20 significant pathways
  font.size = 10,
  title = "KEGG Pathway Enrichment Dot Plot",
  color = "qvalue" # Color by qvalue, can also change to pvalue
) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.title = element_text(size = 10, face = "bold")
  ) +
  labs(
    x = "Gene Ratio",
    y = "KEGG Pathway",
    color = "Q-value",
    size = "Gene Count" # Legend label for bubble size
  )

# Display KEGG bubble plot
print(kegg_dotplot)

kegg_dotplot_count <- dotplot(
  kegg_result,
  showCategory = 20,
  font.size = 10,
  title = "KEGG Enrichment (Sorted by Gene Count)",
  color = "qvalue",
  orderBy = "pvalue",       # Specify sorting by 'Count'
  decreasing = F       # Descending order (from most to least)
) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.title = element_text(size = 10, face = "bold")
  ) +
  labs(
    x = "Gene Ratio",
    y = "KEGG Pathway",
    color = "Q-value",
    size = "Gene Count"
  )

print(kegg_dotplot_count)















# KEGG enrichment analysis (combine UP/DOWN genes)
kegg_result <- enrichKEGG(
  gene = c(up_list, down_list),
  organism = "hsa",  # human=hsa, mouse=mmu, rat=rn
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.05
)
# Process KEGG results and plot (dotplot version)
if (is.null(kegg_result) || nrow(kegg_result) == 0) {
  message("No significant results in KEGG enrichment!")
} else {
  # Core: directly plot using clusterProfiler's dotplot (no need to convert to data frame)
  kegg_dotplot <- dotplot(
    kegg_result,
    showCategory = 10,        # Only display top 10 significant KEGG pathways
    font.size = 10,           # Base font size
    title = "KEGG Pathway Enrichment Dot Plot (Top 10)",
    color = "qvalue",         # Color by qvalue (redder means more significant)
    orderBy = "GeneRatio"     # Sort by enrichment factor (GeneRatio)
  ) +
    # Beautification: completely unified style with DO/GO bubble plots
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Title centered and bold
      axis.text.y = element_text(size = 9),                             # Y-axis pathway name font
      axis.text.x = element_text(angle = 30, hjust = 1, size = 9),      # X-axis rotated 30 degrees to avoid overlap
      legend.title = element_text(size = 10, face = "bold"),            # Legend title bold
      panel.grid = element_blank()                                      # Hide grid lines for clarity
    ) +
    # Custom axis/legend labels (enhance readability)
    labs(
      x = "Gene Ratio (Enriched / Total)",
      y = "KEGG Pathway",
      color = "Q-value",
      size = "Gene Count"
    )
  
  print(kegg_dotplot)
  
  # Optional: save high-resolution image
  ggsave(
    "aaaKEGG_Top10_Dotplot.png",
    plot = kegg_dotplot,
    width = 10, height = 7, dpi = 300
  )
}

# Merge UP and DOWN ENTREZID (consistent with gene range of previous KEGG enrichment)
all_gene_list <- c(up_list, down_list)  # up_list/down_list are previously split gene IDs






# Immune Infiltration Correlation Analysis
# Data Transformation
# 1. Directly convert to data frame
immune_infiltration2.1 <- as.data.frame(immune_infiltration2)
# 1. Assign values of first column as row names
rownames(immune_data1) <- immune_data1$Mixture  
immune_data1 <- immune_data1[, -1]  

# Data transformation: row names as sample names, column names as gene names

immune_infiltration1.1 <- t(immune_infiltration1)
immune_infiltration1.2 <- t(immu_int1)
identical(rownames(immune_data1), rownames(immune_infiltration1.1))
# Extract top 3 hub genes
ppi_top5_genes <- c("MTHFR", "THBS1", "VEGFA")  # Multiple gene names
immune_infiltration_ppitop5.1 <- immune_infiltration1.2[, ppi_top5_genes]  # Returns sub-matrix (samples × genes)

# Convert back to data frame format
immune_infiltration_ppitop5.1 <- as.data.frame(immune_infiltration_ppitop5.1)
# Only take patient expression values
immune_infiltration_ppitop5.1 <- immune_infiltration_ppitop5.1[13:24, ]
immune_data1.1 <- immune_data1[13:24, ]
# Step 1: Check column types of immune_data1 (identify non-numeric columns)

# Check if samples are consistent
cat("Are samples consistent:", identical(rownames(immune_data1.1), rownames(immune_infiltration_ppitop5.1)), "\n")
# ============ Part 2: Correlation Calculation (multiple methods) ============
# 1. Define list of row names to delete (replace with actual row names to delete, note case sensitivity)
delete_rownames <- c("GSM722680", "GSM722683", "GSM722686", "GSM722687")  

# 2. Replace your_df with your actual data frame name (e.g., immune_data1/immune_infiltration_ppitop5.1)
immune_data1.2 <- immune_data1.1[!rownames(immune_data1.1) %in% delete_rownames, ]  
immune_infiltration_ppitop5.2 <- immune_infiltration_ppitop5.1[!rownames(immune_infiltration_ppitop5.1) %in% delete_rownames, ]  
cat("Are samples consistent:", identical(rownames(immune_data1.2), rownames(immune_infiltration_ppitop5.2)), "\n")

### Method 1: Use base R cor.test() function (pairwise calculation, includes p-value)
calculate_correlations <- function(immune_data, gene_data, method = "pearson") {
  # Matrix to store results
  cor_matrix <- matrix(NA, nrow = ncol(immune_data), ncol = ncol(gene_data),
                       dimnames = list(colnames(immune_data), colnames(gene_data)))
  pval_matrix <- matrix(NA, nrow = ncol(immune_data), ncol = ncol(gene_data),
                        dimnames = list(colnames(immune_data), colnames(gene_data)))
  
  # Pairwise calculation of correlation and p-value
  for (i in 1:ncol(immune_data)) {
    for (j in 1:ncol(gene_data)) {
      # Check for missing values
      if (sum(!is.na(immune_data[, i]) & !is.na(gene_data[, j])) > 3) {
        cor_test <- cor.test(immune_data[, i], gene_data[, j], method = method)
        cor_matrix[i, j] <- cor_test$estimate  # Correlation coefficient
        pval_matrix[i, j] <- cor_test$p.value  # p-value
      }
    }
  }
  
  return(list(correlation = cor_matrix, pvalue = pval_matrix))
}

# Calculate correlation (default Pearson, options: "spearman" or "kendall")
result <- calculate_correlations(immune_data1.1, immune_infiltration_ppitop5.1, method = "kendall")

# View results
cat("\n========== Correlation Coefficient Matrix ==========\n")
print(round(result$correlation, 3))

cat("\n========== P-value Matrix ==========\n")
print(round(result$pvalue, 4))

cor_matrix <- as.matrix(result$correlation)
pval_matrix <- as.matrix(result$pvalue)
cor_matrix <- cor_matrix[, colSums(is.na(cor_matrix)) != nrow(cor_matrix)]  
# Remove rows/columns with all NA → then remove all rows/columns containing NA (completely NA-free)
cor_matrix <- cor_matrix[, colSums(is.na(cor_matrix)) != nrow(cor_matrix)]  
cor_matrix <- cor_matrix[rowSums(is.na(cor_matrix)) != ncol(cor_matrix), ]  
cor_matrix <- na.omit(cor_matrix)  # Remove all rows containing NA (na.omit automatically deletes corresponding columns to keep matrix complete)
pval_matrix <- pval_matrix[, colSums(is.na(pval_matrix)) != nrow(pval_matrix)]  
pval_matrix <- pval_matrix[rowSums(is.na(pval_matrix)) != ncol(pval_matrix), ]  
pval_matrix <- na.omit(pval_matrix)  # Remove all rows containing NA (na.omit automatically deletes corresponding columns to keep matrix complete)

# ============ Part 3: Multiple Heatmap Drawing Methods ============

### Method 1: Basic pheatmap heatmap (with values and significance)
# Create significance marker matrix
signif_stars <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
signif_stars[pval_matrix > 0.05] <- " "
signif_stars[pval_matrix < 0.05] <- "*"
signif_stars[pval_matrix < 0.01] <- "**"
signif_stars[pval_matrix < 0.001] <- "***"
library(pheatmap)
library(corrplot)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(ggpubr)

library(pheatmap)

# ========== Step 1: Clean NA (simplified version, ensure cor/pval matrix rows and columns match completely) ==========
# Extract correlation and p-value matrices
cor_matrix <- as.matrix(result$correlation)
pval_matrix <- as.matrix(result$pvalue)

# 1. First delete all-NA rows/columns from cor_matrix, then synchronously filter pval_matrix (ensure consistent dimensions)
cor_matrix <- cor_matrix[, colSums(is.na(cor_matrix)) != nrow(cor_matrix)]  # Delete all-NA columns
cor_matrix <- cor_matrix[rowSums(is.na(cor_matrix)) != ncol(cor_matrix), ]  # Delete all-NA rows
cor_matrix <- na.omit(cor_matrix)  # Delete all rows/columns containing NA

# 2. Filter pval_matrix according to row/column names of cor_matrix (key! avoid dimension mismatch)
pval_matrix <- pval_matrix[rownames(cor_matrix), colnames(cor_matrix)]  
pval_matrix[is.na(pval_matrix)] <- 1  # Set remaining NA p-values to 1 (non-significant)

# ========== Step 2: Generate significance star matrix (*/**/***) ==========
signif_stars <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix),
                       dimnames = dimnames(pval_matrix))  # Preserve row/column names
# Assign stars based on p-value
signif_stars[pval_matrix < 0.05 & pval_matrix >= 0.01] <- "*"    # p<0.05
signif_stars[pval_matrix < 0.01 & pval_matrix >= 0.001] <- "**" # p<0.01
signif_stars[pval_matrix < 0.001] <- "***"                      # p<0.001
signif_stars[pval_matrix >= 0.05] <- ""                         # p≥0.05 (no star)

# ========== Step 3: Draw heatmap with values + stars ==========
# Define color palette (correlation range -1~1, blue=negative correlation, red=positive correlation)
color_palette <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)

# Concatenate "correlation value + stars" (keep 2 decimal places)
display_matrix <- matrix(
  paste0(round(cor_matrix, 2), signif_stars),  # e.g., "0.78**" "0.35" "-0.61*"
  nrow = nrow(cor_matrix),
  ncol = ncol(cor_matrix),
  dimnames = dimnames(cor_matrix)
)

# Draw heatmap (core code)
p1 <- pheatmap(
  cor_matrix,  # Heatmap color based on correlation values
  color = color_palette,
  main = "Immune Infiltration and Core Gene Correlation Heatmap (with significance)",
  display_numbers = display_matrix,  # Display "value + stars"
  number_format = "s",  # Disable default number formatting (avoid star truncation)
  number_color = "black",  # Color of values + stars
  fontsize_row = 9,        # Row name font size
  fontsize_col = 10,       # Column name font size
  fontsize_number = 7,     # Cell value font size
  cluster_rows = TRUE,     # Row clustering
  cluster_cols = TRUE,     # Column clustering
  treeheight_row = 30,     # Row clustering tree height
  treeheight_col = 30,     # Column clustering tree height
  border_color = "grey60", # Cell border color
  cellwidth = 40,          # Cell width
  cellheight = 30,         # Cell height
  breaks = seq(-1, 1, length.out = 101),  # Color segmentation (adapted for correlation -1~1)
  silent = FALSE
)

# Display heatmap
print(p1)

library(tidyr)
library(dplyr)

library(ggplot2)
library(dplyr)
library(tidyr)

# ========== Step 1: Matrix to long format + data preprocessing ==========
# 1. Correlation matrix to long format (rows=immune cells, columns=genes)
cor_long <- as.data.frame(cor_matrix) %>%
  rownames_to_column(var = "Immune_Cell") %>%
  pivot_longer(
    cols = -Immune_Cell,
    names_to = "Gene",
    values_to = "Correlation"
  )

# 2. P-value matrix to long format synchronously
pval_long <- as.data.frame(pval_matrix) %>%
  rownames_to_column(var = "Immune_Cell") %>%
  pivot_longer(
    cols = -Immune_Cell,
    names_to = "Gene",
    values_to = "P_value"
  )

# 3. Merge data + calculate key columns
plot_data <- merge(cor_long, pval_long, by = c("Immune_Cell", "Gene")) %>%
  mutate(
    Cor_Abs = abs(Correlation),  # Absolute correlation (bubble size)
    lower = Correlation - 0.1,   # Error bar lower limit (compatible with negative numbers)
    upper = Correlation + 0.1,   # Error bar upper limit (compatible with negative numbers)
    # P-value grouping (color mapping, matches example)
    P_Group = factor(
      case_when(
        P_value < 0.001 ~ "<0.001",
        P_value < 0.01 ~ "<0.01",
        P_value < 0.05 ~ "<0.05",
        TRUE ~ ">0.05"
      ),
      levels = c("<0.001", "<0.01", "<0.05", ">0.05")
    )
  )

# ========== Step 2: Sort by correlation (core!) ==========
# For each gene group, sort by correlation descending (high positive to low, negatives later)
# If ascending order desired (negatives first, positives later), change to arrange(Correlation, .by_group = TRUE)
plot_data_sorted <- plot_data %>%
  group_by(Gene) %>%  # Group by gene
  arrange(desc(Correlation), .by_group = TRUE) %>%  # Correlation descending (can change to ascending)
  mutate(
    # Reset factor levels of Immune_Cell to sorted order (Y-axis fixed to this order)
    Immune_Cell = factor(Immune_Cell, levels = unique(Immune_Cell))
  ) %>%
  ungroup()

# ========== Step 3: Define single gene plotting function (sorted) ==========
plot_single_gene <- function(gene_data) {
  gene_name <- unique(gene_data$Gene)
  
  ggplot(gene_data, aes(x = Correlation, y = Immune_Cell)) +
    # 1. Horizontal error bars (compatible with negative numbers, example gray bars)
    geom_errorbarh(
      aes(xmin = lower, xmax = upper),
      color = "gray",
      height = 0.2
    ) +
    # 2. Bubbles (size=absolute correlation (continuous), color=P-value group)
    geom_point(
      aes(size = Cor_Abs, color = P_Group),
      alpha = 0.8
    ) +
    # 3. X-axis range (adapted for positive/negative correlations)
    xlim(-1, 1) +
    # 4. Bubble size: continuous value mapping (no grouping, based on actual values)
    scale_size_continuous(
      range = c(2, 10),  # Minimum/maximum bubble size (adjustable)
      name = "|Cor|"     # Legend title: absolute correlation
    ) +
    # 5. Bubble colors (match example color scheme)
    scale_color_manual(
      values = c(
        "<0.001" = "#E69F00",  # Orange (p<0.001)
        "<0.01" = "#56B4E9",   # Blue (p<0.01)
        "<0.05" = "#009E73",   # Green (p<0.05)
        ">0.05" = "#999999"    # Gray (p≥0.05)
      ),
      name = "pvalue"
    ) +
    # 6. Title + axes (match example style)
    ggtitle(gene_name) +
    labs(x = "cor", y = "Immune") +
    # 7. Theme beautification
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 9),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

plot_single_gene_bar <- function(gene_data) {
  gene_name <- unique(gene_data$Gene)
  
  ggplot(gene_data, aes(x = Immune_Cell, y = Correlation)) +
    # Draw bar plot, color represents p-value group
    geom_bar(
      aes(fill = P_Group),
      stat = "identity",
      width = 0.7,
      alpha = 0.8
    ) +
    # Add error bars
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.3,
      color = "black",
      linewidth = 0.5
    ) +
    # Add correlation value labels
    geom_text(
      aes(label = format(round(Correlation, 2), nsmall = 2)),
      vjust = -0.5,  # Positive value labels above
      size = 3
    ) +
    # Add significance stars
    geom_text(
      aes(label = ifelse(P_value < 0.001, "***", 
                         ifelse(P_value < 0.01, "**",
                                ifelse(P_value < 0.05, "*", ""))),
          y = ifelse(Correlation >= 0, upper + 0.05, lower - 0.05)),
      size = 4
    ) +
    # Axis adjustment
    ylim(min(gene_data$lower, -1) - 0.1, max(gene_data$upper, 1) + 0.1) +
    # Color scheme
    scale_fill_manual(
      values = c(
        "<0.001" = "#E69F00",
        "<0.01" = "#56B4E9",
        "<0.05" = "#009E73",
        ">0.05" = "#999999"
      ),
      name = "pvalue"
    ) +
    # Title and labels
    ggtitle(gene_name) +
    labs(x = "Immune Cell", y = "Correlation") +
    # Theme beautification
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

# ========== Step 4: Split by gene + batch generation/saving ==========
# 1. Split sorted data by gene
gene_list <- split(plot_data_sorted, plot_data_sorted$Gene)

# 2. Draw example for single gene (e.g., MTHFR)
plot_single_gene(gene_list[["MTHFR"]])  # Replace with your gene name
plot_single_gene_bar(gene_list[["THBS1"]])
# 3. Batch save all gene plots
for (gene in names(gene_list)) {
  p <- plot_single_gene(gene_list[[gene]])
  ggsave(
    filename = paste0("Gene_", gene, "_Immune_Infiltration_Correlation(Correlation_Sorted).png"),
    plot = p,
    width = 6,
    height = 8  # Can adjust height based on number of immune cells
  )
}

# Result:
immune_result1 <- cor(immune_infiltration_ppitop5.1, immune_data1.1, method = "spearman")
# Draw heatmap
library(pheatmap)
if (!require("ggcorrplot")) install.packages("ggcorrplot")
library(ggcorrplot)
library(ggplot2)  # Dependency package, needs to be loaded simultaneously
# Delete all-NA columns (columns where all values are NA)
immune_result1 <- immune_result1[, colSums(is.na(immune_result1)) != nrow(immune_result1)]  
# Unstandardized data
# Unstandardized data

# Basic heatmap
mianyijinrunretu <- ggcorrplot(
  immune_result1,               # Correlation matrix
  method = "square",           # Heatmap shape: "square" or "circle"
  colors = c("#2166ac", "white", "#B2182B"),
  type = "full",               # Display range: "full" (full matrix), "upper" (upper triangle), "lower" (lower triangle)
  title = "",  # Title
  show.legend = TRUE,          # Show color legend
  show.diag = TRUE,            # Show diagonal (self-correlation, usually 1)
  lab = TRUE,                  # Show correlation coefficient values
  lab_size = 3,                # Coefficient font size
  tl.cex = 10,                 # Axis label font size
  tl.srt = 45                  # Column label rotation angle (avoid overlap)
)

print(mianyijinrunretu)  # Display graph






# Define three gene lists
lasso_genes <- c("ATF3", "IRAK2", "MTHFR", "THBS1", "VEGFA")
best_features_genes <- c("THBS1", "VEGFA", "COPS4", "CDK9", "JAK1", "USP10", "H1-2", "MTHFR", "IFNAR2", "IRAK1", "GBP2", "IFRD1", "NOP53")
rf_genes <- c("VEGFA", "THBS1", "COPS4", "IRAK1", "GBP2", "USP10", "IFNAR2", "CDK9", "JAK1", "MTHFR")

# 1. Install and load required packages (if not already installed)
# install.packages("VennDiagram")
# install.packages("grid")
# install.packages("png")
library(VennDiagram)
library(grid)
library(png)


venn.diagram(
  x = list(
    "LASSO" = lasso_genes,
    "SVM RFE" = best_features_genes,
    "Random Forest" = rf_genes
  ),
  filename = "venn_diagram_final.png",
  
  # Use default Euler diagram layout
  scaled = TRUE,
  
  # Appearance settings (change col to transparent to remove border lines)
  fill = c("#C03B2D", "#62B69E", "#EBB81A"),
  alpha = 0.5,
  col = "transparent",  # Key modification: set outline to transparent
  
  # Adjust label positions to look like an inverted triangle
  cat.pos = c(-30, 30, 180),  # Label angles
  cat.dist = c(0.05, 0.05, 0.02),  # Label distances
  cat.cex = 1.2,
  cat.fontface = "bold",
  
  # Other settings
  cex = 1.5,
  fontface = "bold",
  height = 1000,
  width = 1000,
  resolution = 300
)

