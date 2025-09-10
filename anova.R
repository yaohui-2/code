####插补缺失值
library(impute)
d <- read.table("E:\\西湖大学\\workto_cai\\IGA_exp_log2.txt",sep = "\t",header = TRUE,row.names = 1,check.names = F)
d <- as.matrix(d)
f<- impute.knn(d,k=10,rowmax=0.5,colmax = 0.8)
newdata <- f[["data"]]
write.csv(newdata,"IGA_exp_log2_impute.csv")

####组间蛋白方差分析
rm(list = ls())
library(tidyverse)
setwd("E:\\西湖大学\\workto_xiao\\MBPA\\Mufzz\\Muffz_input\\W\\anova")
t_test <- read.csv('ace_anova_input.csv',header = T, row.names = 1)
t_t_test <- t(t_test)###转置

####构建循环matrix
protein <- matrix(0, ncol(t_test), 2)
colnames(protein) <- c('sample','intensity')
#protein[,1] <- rep(c("0h","30min","60min"), 1, each = 27)
#protein[,1] <-  c(rep("A",3),rep("B",37),rep(c("C","D"),1,each=27))
protein[,1] <-  c(rep("A",3),rep("B",2),rep("C",3),rep("D",3),rep("E",3))
protein <- data.frame(protein)
t_test <- t_test %>% mutate(p_value = 0)###

for (i in 1:ncol(t_t_test)){
  # i = 1
  
  protein[,2] <- t_t_test[,i]
  data_aov <- aov(intensity ~ sample, data = protein)
  p_value <- summary(data_aov)[[1]][1,5]###提取p值
  t_test[i,15] <- p_value[1]
}

####统计P<0.05
#time_p <- HDP[HDP$time.p < 0.05,]###
#inflame_p <- t_test[t_test$p_value < 0.05,]
#joint_p <- HDP[HDP$joint.p < 0.05,]###
#write.csv(t_test, file = 'ace_anova_result.csv')

####计算FDR
#test_adj <- read.csv("ace_anova_result.csv",header = T, row.names = 1)
adjusted_p <- p.adjust(t_test$p_value, method = "BH")
t_test$FDR <- adjusted_p 
write.csv(t_test, file = 'ace_anova_result_adj.csv')
