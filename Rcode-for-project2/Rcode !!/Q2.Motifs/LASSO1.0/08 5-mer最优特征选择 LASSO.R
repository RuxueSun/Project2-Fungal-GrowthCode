data_5mer <- read_csv("C://Users//Summer//Desktop//PH data//project 2//Data//H99_all_genes_promoter_500nt_5mer_counts.csv")
data_5mer <- as.data.frame(data_5mer)
rownames(data_5mer) <- data_5mer$Gene
data_5mer01 <- as.data.frame(scale(data_5mer[, -1])) 
data_5mer01 <- cbind("Gene" = data_5mer$Gene,data_5mer01)

merged_initial_expression <- read.csv("C://Users//Summer//Desktop//PH data//project 2//240719//merged_initial_expression.csv")
merged_initial_expression <- as.data.frame(merged_initial_expression)

merged_data_5mer <- merge(merged_initial_expression, data_5mer01, by = "Gene")


#筛选显著的motif
# 构建线性混合效应模型
response_5mer <- merged_data_5mer$Expression
predictors_5mer <- merged_data_5mer[, c("initial_value",colnames(data_5mer)[-1])]  # 除去基因名列

#使用LASSO回归选择合适的K-mer数据
# 准备 LASSO 分析的数据
library(glmnet)
X <- as.matrix(merged_data_5mer[, c("initial_value",colnames(data_5mer)[-1])])
y <- merged_data_5mer$Expression

# 使用交叉验证选择最优的alpha值和lambda值
# 定义alpha值的序列
alpha_values <- seq(0, 1, by = 0.01)

# 初始化存储交叉验证结果的列表
cv_results <- list()

# 对每个alpha值进行交叉验证
for (alpha_val in alpha_values) {
  cv_model <- cv.glmnet(X, y, nfolds = 10, alpha = alpha_val)
  cv_results[[as.character(alpha_val)]] <- cv_model
}

# 提取最优的alpha和对应的lambda值
best_alpha <- NULL
best_lambda <- NULL
best_mse <- Inf

for (alpha_val in alpha_values) {
  cv_model <- cv_results[[as.character(alpha_val)]]
  min_mse <- min(cv_model$cvm)
  
  if (min_mse < best_mse) {
    best_mse <- min_mse
    best_alpha <- alpha_val
    best_lambda <- cv_model$lambda.min
  }
}
lasso_result_5mer <- data.frame(name=c("alpha","lambda","mse"),value=c(best_alpha,best_lambda,best_mse))
write.csv(lasso_result_5mer,"lasso_result_5mer.csv",row.names = FALSE)
# 使用最佳alpha值和lambda值重新拟合模型
elastic_net_model <- glmnet(X, y, alpha = best_alpha, lambda = best_lambda)

# 输出模型系数
coef(elastic_net_model)

# 打印最佳的alpha和lambda值
print(paste("Best alpha:", best_alpha))
print(paste("Best lambda:", best_lambda))

# 可视化结果
plot(cv_results[[as.character(best_alpha)]])


