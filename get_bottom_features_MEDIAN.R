library(tidyverse)

aucs <- read.table("results/iteration1/training/testing_auc.list")
aucs <- aucs$V1

median_value <- median(aucs)

median_subset <- aucs[which.min(abs(aucs - median_value))]

final_subset <- as.character(which(aucs == median_subset))

# Check if the length of median_indices is greater than or equal to 2
if (length(final_subset) >= 2) {
  # Randomly select one index from median_indices
  chosen_index <- sample(final_subset, 1)
  print(chosen_index)

  # Save the subset corresponding to the chosen index
  final_subset <- chosen_index
}

print(final_subset)

paste("subset", final_subset, sep="")

vi_path <- paste("results/iteration1/training/subsets/subset", final_subset, "/Variable_Importance.txt", sep="")
ts_path <- paste("results/iteration1/training/subsets/subset", final_subset, "/TimesSelected.txt", sep="")
vi <- read.table(vi_path)
ts <- read.table(ts_path)

vi <- merge(vi,ts,by=1)

names(vi) <- c("FEATURE_number", "Variable_Importance", "Times Selected")
vi$percentage_VI <- 100*vi$Variable_Importance/max(abs(vi$Variable_Importance))
vi_subset <- vi[order(vi$"Times Selected", decreasing = F),]

str(vi_subset)

outfile <- write.csv(vi_subset, "results/iteration1/training/feature_importance_scores.txt", row.names=FALSE)

