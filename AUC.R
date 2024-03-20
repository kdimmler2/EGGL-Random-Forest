library(pROC)
library(tidyverse)

outfile <- file("training_auc_values.txt", "w")

pred_train <- read.table("../EGBV.txt")
pred_test <- read.table("../Predictions.txt")

ytrain<-read.table(file="../data/training.txt",header=F)
ytest<-read.table(file="../data/testing.txt",header=F)

roc_curve_train <- roc(response = ytrain$V1, predictor = pred_train$V2)

auc_value_train <- auc(roc_curve_train)
auc_value_train <- as.character(auc_value_train[1])

writeLines(auc_value_train, outfile)

close(outfile)

outfile <- file("testing_auc_values.txt", "w")

roc_curve_test <- roc(response = ytest$V1, predictor = pred_test$V2)

auc_value_test <- auc(roc_curve_test)
auc_value_test <- as.character(auc_value_test[1])

writeLines(auc_value_test, outfile)

close(outfile)

roc_data <- data.frame(
  FPR = roc_curve_test$specificities,
  TPR = roc_curve_test$sensitivities
)

plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "ROC Curve",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()

ggsave("ROC_curves/subset1/ROC_curve.png", plot = plot, width = 8, height = 5)
