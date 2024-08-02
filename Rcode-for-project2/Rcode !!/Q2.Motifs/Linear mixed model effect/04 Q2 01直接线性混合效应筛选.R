# 加载必要的库
library(readr)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(MASS)
library(caret)
library(DESeq2)

# 读取数据文件
counts_all <- read_csv("C://Users//Summer//Desktop//PH data//project 2//Data//counts_all.csv")
samplesheet <- read_csv("C://Users//Summer//Desktop//PH data//project 2//Data//Gat201_samplesheet.csv")

# 将tibble转换为data frame
counts_all <- as.data.frame(counts_all)
samplesheet <- as.data.frame(samplesheet)

# 设置行名
rownames(counts_all) <- counts_all$Geneid

# 删除初始数据列
# 确保数据框中只有基因表达数据和k-mer频率数据
filtered_counts_all <- counts_all[rowSums(counts_all == 0) == 0, ]
counts_all_without_names <- filtered_counts_all[,c(-1,-2,-3,-4,-5,-6)]

rownames(samplesheet) <- samplesheet$Title
samplesheet_without_ids <- samplesheet[,c(-1,-2,-3,-4,-6,-8)]
# 将必要的列转换为因子类型
samplesheet_without_ids$GAT201 <- as.factor(samplesheet_without_ids$GAT201)
samplesheet_without_ids$Condition <- as.factor(samplesheet_without_ids$Condition)
samplesheet_without_ids$BioRep <- as.factor(samplesheet_without_ids$BioRep)

# 创建DESeq2对象
dds <- DESeqDataSetFromMatrix(countData = counts_all_without_names,
                              colData = samplesheet_without_ids,
                              design = ~ GAT201 + Condition + BioRep)


# 过滤掉低表达基因
dds <- dds[rowSums(counts(dds)) > 10,]

# 运行DESeq2进行标准化
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# 提取标准化后的表达数据
expression_data_standardized <- assay(vsd)
Gene <- rownames(expression_data_standardized)
expression_data_standardized01 <- cbind(Gene, expression_data_standardized)
# 将基因表达数据转换为长格式
expression_data_standardized01 <- as.data.frame(expression_data_standardized01)
#把初始数据加入模型

initial_samples_cols <- colnames(filtered_counts_all)[grepl("_Y_0_", colnames(filtered_counts_all))]
initial_samples <- expression_data_standardized01[,c("Gene",initial_samples_cols)]
A_1 <- initial_samples$A_Y_0_1
A_2 <- initial_samples$A_Y_0_2
a_1 <- initial_samples$a_Y_0_1
a_2 <- initial_samples$a_Y_0_2
B_1 <- initial_samples$B_Y_0_1
B_2 <- initial_samples$B_Y_0_2
M_1 <- initial_samples$M_Y_0_1
M_2 <- initial_samples$M_Y_0_2
New_initial_sample_values <- cbind(A_1,A_1,A_1,A_2,A_2,A_2,
                                   A_1,A_1,A_1,A_2,A_2,A_2,
                                   a_1,a_1,a_1,a_2,a_2,a_2,
                                   a_1,a_1,a_1,a_2,a_2,a_2,
                                   B_1,B_1,B_1,B_2,B_2,B_2,
                                   B_1,B_1,B_1,B_2,B_2,B_2,
                                   M_1,M_1,M_1,M_2,M_2,M_2,
                                   M_1,M_1,M_1,M_2,M_2,M_2
)

# 将每一行转换成列并拼接成一个长列
initial_sample_repeated_values <- as.vector(t(New_initial_sample_values))

# 创建结果数据框
initial_sample_repeated_values <- data.frame(initial_value = initial_sample_repeated_values)
names(initial_sample_repeated_values) <- "initial_value"

initial_sample_repeated_geneids <- rep(initial_samples$Gene,each=48)
initial_sample_repeated_geneids <- as.data.frame(initial_sample_repeated_geneids)
names(initial_sample_repeated_geneids) <- "Gene"

initial_sample_for_modeling <- cbind(initial_sample_repeated_geneids,initial_sample_repeated_values)
write.csv(initial_sample_for_modeling,"initial_sample_for_modeling.csv",row.names = FALSE)

expression_data_standardized011 <- expression_data_standardized01[,c(-2,-6,-16,-20,-30,-34,-44,-48)]
write.csv(expression_data_standardized011,"expression_data_standardized011.csv",row.names = FALSE)

long_expression_data <- expression_data_standardized011 %>%
  pivot_longer(cols = c(-Gene), names_to = "SampleID", values_to = "Expression")
# 确保表达数据是数值型
long_expression_data$Expression <- as.numeric(long_expression_data$Expression)
write.csv(long_expression_data,"long_expression_data.csv",row.names = FALSE)



merged_initial_expression <- long_expression_data
merged_initial_expression$initial_value <- initial_sample_for_modeling$initial_value
write.csv(merged_initial_expression,"merged_initial_expression.csv",row.names = FALSE)

merged_data_2mer <- merge(merged_initial_expression, data_2mer, by = "Gene")


#筛选显著的motif
# 构建线性混合效应模型
response_2mer <- merged_data_2mer$Expression
predictors_2mer <- merged_data_2mer[, colnames(data_2mer)[-1]]  # 除去基因名列

# 构建线性混合效应模型并筛选显著motif
anova_results_2mer <- list()
i <- 1
for (kmer in colnames(predictors_2mer)) {
  print(paste("Fitting model for k-mer:", kmer, " - ", i))   # 调试输出
  
  formula <- as.formula(paste("Expression ~", kmer, "+ (1 | Gene) + (1 | SampleID)"))
  model <- lmer(formula, data = merged_data_2mer)
  anova_results_2mer[[kmer]] <- model
  
  i <- i + 1
}

# 定义函数来判断分析结果并输出结论
analyze_results <- function(anova_results, p_value_threshold = 0.05, output_prefix) {
  # 初始化结果列表和显著的2-mer名称列表
  results_list <- list()
  significant_k_mer <- data.frame(kmer = character())
  
  # 遍历每个k-mer的结果
  for (kmer in names(anova_results)) {
    model_summary <- summary(anova_results[[kmer]])
    coef_summary <- coef(model_summary)
    
    if (nrow(coef_summary) > 1) {
      p_value <- coef_summary[2, "Pr(>|t|)"]
      t_value <- coef_summary[2, "t value"]
      
      # 调试输出
      print(paste("Checking k-mer:", kmer))
      print(paste("P-value:", p_value))
      print(paste("T-value:", t_value))
      
      if (!is.na(p_value) && p_value < p_value_threshold) {
        # 提取显著性结果
        estimate <- coef_summary[2, "Estimate"]
        std_error <- coef_summary[2, "Std. Error"]
        
        # 构建结果数据框
        result_text <- paste("The k-mer", kmer, "has an estimate of", estimate, 
                             "with a standard error of", std_error, 
                             "a t value of", t_value, 
                             "and a p value of", p_value, 
                             ". This suggests that the effect of", kmer, 
                             "on gene expression is significant.")
        
        results_list[[kmer]] <- data.frame(
          kmer = kmer,
          Estimate = estimate,
          StdError = std_error,
          TValue = t_value,
          PValue = p_value,
          ResultText = result_text
        )
        
        # 将显著的2-mer名称存储到数据框中
        significant_k_mer <- rbind(significant_k_mer, data.frame(kmer = kmer))
      }
    }
  }
  
  # 将结果列表转换为数据框
  results_df <- bind_rows(results_list)
  
  # 多重检验校正
  if (nrow(results_df) > 0) {
    results_df$AdjustedPValue <- p.adjust(results_df$PValue, method = "BH")
    significant_results <- results_df[results_df$AdjustedPValue < 0.05, ]
    
    # 构建结果文本
    significant_results$ResultText <- paste("The k-mer", significant_results$kmer, 
                                            "has an estimate of", significant_results$Estimate, 
                                            "with a standard error of", significant_results$StdError, 
                                            "a t value of", significant_results$TValue, 
                                            "and a p value of", significant_results$PValue, 
                                            ". This suggests that the effect of", significant_results$kmer, 
                                            "on gene expression is significant after adjustment.")
    
    write.csv(significant_results, paste0(output_prefix, "_motifs.csv"), row.names = FALSE)
    write.csv(significant_k_mer, paste0(output_prefix, ".csv"), row.names = FALSE)
    return(significant_results)
  } else {
    return("No significant motifs found.")
  }
}


# 示例调用
significant_motifs_df_2mer <- analyze_results(anova_results_2mer,,"significant_2_mer")

# 打印结果
print(significant_motifs_df_2mer)

###########################################################################################
#筛选3mer的数据
data_3mer <- read_csv("C://Users//Summer//Desktop//PH data//project 2//Data//H99_all_genes_promoter_500nt_3mer_counts.csv")
data_3mer <- as.data.frame(data_3mer)
rownames(data_3mer) <- data_3mer$Gene

merged_data_3mer <- merge(merged_initial_expression, data_3mer, by = "Gene")

#筛选显著的motif
# 构建线性混合效应模型
response_3mer <- merged_data_3mer$Expression
predictors_3mer <- merged_data_3mer[, colnames(data_3mer)[-1]]  # 除去基因名列

# 构建线性混合效应模型并筛选显著motif
anova_results_3mer <- list()
i <- 1
for (kmer in colnames(predictors_3mer)) {
  print(paste("Fitting model for k-mer:", kmer, " - ", i))   # 调试输出
  
  formula <- as.formula(paste("Expression ~", kmer, "+ (1 | Gene) + (1 | SampleID)"))
  model <- lmer(formula, data = merged_data_3mer)
  anova_results_3mer[[kmer]] <- model
  i <- i + 1
}

# 定义函数来判断分析结果并输出结论
# 示例调用
significant_motifs_df_3mer <- analyze_results(anova_results_3mer,,"significant_3_mer")

# 打印结果
print(significant_motifs_df_3mer)

##################################################################################
#筛选4mer的数据
data_4mer <- read_csv("C://Users//Summer//Desktop//PH data//project 2//Data//H99_all_genes_promoter_500nt_4mer_counts.csv")
data_4mer <- as.data.frame(data_4mer)
rownames(data_4mer) <- data_4mer$Gene

merged_data_4mer <- merge(merged_initial_expression data_4mer, by = "Gene")

#筛选显著的motif
# 构建线性混合效应模型
response_4mer <- merged_data_4mer$Expression
predictors_4mer <- merged_data_4mer[, colnames(data_4mer)[-1]]  # 除去基因名列

# 构建线性混合效应模型并筛选显著motif
anova_results_4mer <- list()
i <- 1
for (kmer in colnames(predictors_4mer)) {
  print(paste("Fitting model for k-mer:", kmer, " - ", i))  # 调试输出
  
  formula <- as.formula(paste("Expression ~", kmer, "+ (1 | Gene) + (1 | SampleID)"))
  model <- lmer(formula, data = merged_data_4mer)
  anova_results_4mer[[kmer]] <- model
  i <- i + 1
}

# 示例调用
significant_motifs_df_4mer <- analyze_results(anova_results_4mer,,"significant_4_mer")

# 打印结果
print(significant_motifs_df_4mer)
#############################################################################
#筛选5mer的数据
data_5mer <- read_csv("C://Users//Summer//Desktop//PH data//project 2//Data//H99_all_genes_promoter_500nt_5mer_counts.csv")
data_5mer <- as.data.frame(data_5mer)
rownames(data_5mer) <- data_5mer$Gene
long_expression_data <- read.csv("C://Users//Summer//Desktop//PH data//project 2//240719//long_expression_data.csv")
long_expression_data <- data.frame(long_expression_data)

data_5mer_01 <- data_5mer[,1:257]
data_5mer_02 <- data_5mer[,c(1,258:513)]
data_5mer_03 <- data_5mer[,c(1,514:769)]
data_5mer_04 <- data_5mer[,c(1,770:1025)]
write.csv(data_5mer_01,"data_5mer_01.csv",row.names = FALSE)
write.csv(data_5mer_02,"data_5mer_02.csv",row.names = FALSE)
write.csv(data_5mer_03,"data_5mer_03.csv",row.names = FALSE)
write.csv(data_5mer_04,"data_5mer_04.csv",row.names = FALSE)

merged_data_5mer_01 <- merge(merged_initial_expression, data_5mer_01, by = "Gene")

#筛选显著的motif
# 构建线性混合效应模型
response_5mer_01 <- merged_data_5mer_01$Expression
predictors_5mer_01 <- merged_data_5mer_01[, colnames(data_5mer_01)[-1]]  # 除去基因名列

# 构建线性混合效应模型并筛选显著motif
anova_results_5mer_01 <- list()
i <- 1
batch_size <- 50
for (batch_start in seq(1, ncol(predictors_5mer_01), by = batch_size)) {
  batch_end <- min(batch_start + batch_size - 1, ncol(predictors_5mer_01))
  kmer_batch <- colnames(predictors_5mer_01)[batch_start:batch_end]
  
  for (kmer in kmer_batch) {
    print(paste("Fitting model for k-mer:", kmer, " - ", i))  # 调试输出
    
    
    formula <- as.formula(paste("Expression ~", kmer, "+ (1 | Gene) + (1 | SampleID)"))
    model <- lmer(formula, data = merged_data_5mer_01)
    anova_results_5mer_01[[kmer]] <- model
    i <- i + 1
    
  }
  
  # 手动释放内存
  
  gc()
}



# 示例调用
significant_motifs_df_5mer_01 <- analyze_results(anova_results_5mer_01,,"significant_5_mer_01")

# 打印结果
print(significant_motifs_df_5mer_01)

merged_data_5mer_02 <- merge(merged_initial_expression, data_5mer_02, by = "Gene")

#筛选显著的motif
# 构建线性混合效应模型
response_5mer_02 <- merged_data_5mer_02$Expression
predictors_5mer_02 <- merged_data_5mer_02[, colnames(data_5mer_02)[-1]]  # 除去基因名列

# 构建线性混合效应模型并筛选显著motif
anova_results_5mer_02 <- list()
i <- 1
batch_size <- 50
for (batch_start in seq(1, ncol(predictors_5mer_02), by = batch_size)) {
  batch_end <- min(batch_start + batch_size - 1, ncol(predictors_5mer_02))
  kmer_batch <- colnames(predictors_5mer_02)[batch_start:batch_end]
  
  for (kmer in kmer_batch) {
    print(paste("Fitting model for k-mer:", kmer, " - ", i))  # 调试输出
    
    
    formula <- as.formula(paste("Expression ~", kmer, "+ (1 | Gene) + (1 | SampleID)"))
    model <- lmer(formula, data = merged_data_5mer_02)
    anova_results_5mer_02[[kmer]] <- model
    i <- i + 1
    
  }
  
  # 手动释放内存
  
  gc()
}
# 示例调用
significant_motifs_df_5mer_02 <- analyze_results(anova_results_5mer_02,,"significant_5_mer_02")

# 打印结果
print(significant_motifs_df_5mer_02)


merged_data_5mer_03 <- merge(merged_initial_expression, data_5mer_03, by = "Gene")

#筛选显著的motif
# 构建线性混合效应模型
response_5mer_03 <- merged_data_5mer_03$Expression
predictors_5mer_03 <- merged_data_5mer_03[, colnames(data_5mer_03)[-1]]  # 除去基因名列

# 构建线性混合效应模型并筛选显著motif
anova_results_5mer_03 <- list()
i <- 1
batch_size <- 50
for (batch_start in seq(1, ncol(predictors_5mer_03), by = batch_size)) {
  batch_end <- min(batch_start + batch_size - 1, ncol(predictors_5mer_03))
  kmer_batch <- colnames(predictors_5mer_03)[batch_start:batch_end]
  
  for (kmer in kmer_batch) {
    print(paste("Fitting model for k-mer:", kmer, " - ", i))  # 调试输出
    
    
    formula <- as.formula(paste("Expression ~", kmer, "+ (1 | Gene) + (1 | SampleID)"))
    model <- lmer(formula, data = merged_data_5mer_03)
    anova_results_5mer_03[[kmer]] <- model
    i <- i + 1
    
  }
  
  # 手动释放内存
  
  gc()
}



# 示例调用
significant_motifs_df_5mer_03 <- analyze_results(anova_results_5mer_03,,"significant_5_mer_03")

# 打印结果
print(significant_motifs_df_5mer_03)

merged_data_5mer_04 <- merge(merged_initial_expression, data_5mer_04, by = "Gene")

#筛选显著的motif
# 构建线性混合效应模型
response_5mer_04 <- merged_data_5mer_04$Expression
predictors_5mer_04 <- merged_data_5mer_04[, colnames(data_5mer_04)[-1]]  # 除去基因名列

# 构建线性混合效应模型并筛选显著motif
anova_results_5mer_04 <- list()
i <- 1
batch_size <- 50
for (batch_start in seq(1, ncol(predictors_5mer_04), by = batch_size)) {
  batch_end <- min(batch_start + batch_size - 1, ncol(predictors_5mer_04))
  kmer_batch <- colnames(predictors_5mer_04)[batch_start:batch_end]
  
  for (kmer in kmer_batch) {
    print(paste("Fitting model for k-mer:", kmer, " - ", i))  # 调试输出
    
    
    formula <- as.formula(paste("Expression ~", kmer, "+ (1 | Gene) + (1 | SampleID)"))
    model <- lmer(formula, data = merged_data_5mer_04)
    anova_results_5mer_04[[kmer]] <- model
    i <- i + 1
    
  }
  
  # 手动释放内存
  
  gc()
}



# 示例调用
significant_motifs_df_5mer_04 <- analyze_results(anova_results_5mer_04,,"significant_5_mer_04")

# 打印结果
print(significant_motifs_df_5mer_04)




#########################################################################
#将筛选出来的显著的K-mer数据放到一个模型里面去
sdata_2mer <- read_csv("C://Users//Summer//Desktop//PH data//project 2//significant_2_mer.csv")
sdata_2mer <- as.data.frame(sdata_2mer)


# 函数：生成反向互补序列
reverse_complement <- function(sequence) {
  complement <- chartr("ACGT", "TGCA", sequence)
  return(paste(rev(strsplit(complement, NULL)[[1]]), collapse = ""))
}

# 函数：筛选不存在反向互补序列的K-mer
filter_reverse_complements <- function(kmer_df) {
  # 检查是否存在名为 "kmer" 的列
  if (!"kmer" %in% colnames(kmer_df)) {
    stop("数据框中不存在名为 'kmer' 的列")
  }
  
  unique_kmers <- c()
  seen_kmers <- c()
  
  for (kmer in kmer_df$kmer) {
    rev_comp <- reverse_complement(kmer)
    if (!(kmer %in% seen_kmers) && !(rev_comp %in% seen_kmers)) {
      unique_kmers <- c(unique_kmers, kmer)
      seen_kmers <- c(seen_kmers, kmer, rev_comp)
    }
  }
  
  # 返回筛选后的数据框
  return(kmer_df %>% filter(kmer %in% unique_kmers))
}

filtered_significant_2_mer <- filter_reverse_complements(sdata_2mer)
write.csv(filtered_significant_2_mer,"filtered_significant_2_mer.csv",row.names = FALSE)

#将筛选出来的显著的K-mer数据放到一个模型里面去
sdata_3mer <- read_csv("C://Users//Summer//Desktop//PH data//project 2//significant_3_mer.csv")
sdata_3mer <- as.data.frame(sdata_3mer)
filtered_significant_3_mer <- filter_reverse_complements(sdata_3mer)
write.csv(filtered_significant_3_mer,"filtered_significant_3_mer.csv",row.names = FALSE)

#将筛选出来的显著的K-mer数据放到一个模型里面去
sdata_4mer <- read_csv("C://Users//Summer//Desktop//PH data//project 2//significant_4_mer.csv")
sdata_4mer <- as.data.frame(sdata_4mer)
filtered_significant_4_mer <- filter_reverse_complements(sdata_4mer)
write.csv(filtered_significant_4_mer,"filtered_significant_4_mer.csv",row.names = FALSE)

#算出5-mer的数据再合并在一起
将筛选出来的显著的K-mer数据放到一个模型里面去
sdata_5mer_01 <- read_csv("C://Users//Summer//Desktop//PH data//project 2//240721//significant_5_mer_01.csv")
sdata_5mer_01 <- as.data.frame(sdata_5mer_01)
sdata_5mer_02 <- read_csv("C://Users//Summer//Desktop//PH data//project 2//240721//significant_5_mer_02.csv")
sdata_5mer_02 <- as.data.frame(sdata_5mer_02)
sdata_5mer_03 <- read_csv("C://Users//Summer//Desktop//PH data//project 2//240721//significant_5_mer_03.csv")
sdata_5mer_03 <- as.data.frame(sdata_5mer_03)
sdata_5mer_04 <- read_csv("C://Users//Summer//Desktop//PH data//project 2//240721//significant_5_mer_04.csv")
sdata_5mer_04 <- as.data.frame(sdata_5mer_04)
sdata_5mer <- rbind(sdata_5mer_01,sdata_5mer_02,sdata_5mer_03,sdata_5mer_04)
filtered_significant_5_mer <- filter_reverse_complements(sdata_5mer)
write.csv(filtered_significant_5_mer,"filtered_significant_5_mer.csv",row.names = FALSE)

