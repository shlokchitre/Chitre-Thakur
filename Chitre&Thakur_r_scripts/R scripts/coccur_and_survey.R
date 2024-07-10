#for total sites and species

library(cooccur)
library("writexl")



# for seasons combo-----------------------------------------------------------------

library(cooccur)

pathoctjan <- "/data_winter.csv"

dataoctjan <- read.csv(pathoctjan)

cooccur_octjan <- cooccur(mat=dataoctjan, type = "spp_site", thresh = TRUE, spp_names = FALSE,true_rand_classifier = 0.1, prob = "hyper",site_mask = NULL, only_effects = FALSE,eff_standard = TRUE, eff_matrix = FALSE)
resultsoctjan <- cooccur_octjan$results
summary(cooccur_octjan)
View(resultsoctjan)


library("writexl")
write_xlsx(resultsoctjan,"~/resultsoctjan.xlsx")

#--------------------------------------------------------------------
June-Sept
pathjunsep <- "~/data_monsoon.csv"

datajunsep <- read.csv(pathjunsep)

cooccur_junsep <- cooccur(mat=datajunsep, type = "spp_site", thresh = TRUE, spp_names = FALSE,true_rand_classifier = 0.1, prob = "hyper",site_mask = NULL, only_effects = FALSE,eff_standard = TRUE, eff_matrix = FALSE)
resultsjunsep <- cooccur_junsep$results
summary(cooccur_junsep)
View(resultsjunsep)
#-----------------------------------------------



library(cooccur)

pathfebmay <- "~/data_summer.csv"

datafebmay <- read.csv(pathfebmay)

cooccur_febmay <- cooccur(mat=datafebmay, type = "spp_site", thresh = TRUE, spp_names = FALSE,true_rand_classifier = 0.1, prob = "hyper",site_mask = NULL, only_effects = FALSE,eff_standard = TRUE, eff_matrix = FALSE)
resultsfebmay <- cooccur_febmay$results
summary(cooccur_febmay)
View(resultsfebmay)


library("writexl")
write_xlsx(resultsfebmay,"~/resultsfebmay.xlsx")



#all seasons--------------------------------------------------------------------------------------

library(cooccur)

pathall3 <- "~/all3.csv"

dataall3 <- read.csv(pathall3)

cooccur_all3 <- cooccur(mat=dataall3, type = "spp_site", thresh = TRUE, spp_names = FALSE,true_rand_classifier = 0.1, prob = "hyper",site_mask = NULL, only_effects = FALSE,eff_standard = TRUE, eff_matrix = FALSE)
resultsall3 <- cooccur_all3$results
summary(cooccur_all3)
View(resultsall3)

library("writexl")
write_xlsx(resultsall3,"~/resultsYear.xlsx")



#Fischer Exact Test or Proportion test 
# Create a matrix for the contingency table
# Define the counts
group1_with_feature <- 17       #winter
group1_without_feature <- 3.75
group2_with_feature <- 7.75       #summer
group2_without_feature <- 15.75
group3_with_feature <- 7.25       #monsoon
group3_without_feature <- 7


# Create a contingency table
contingency_table <- matrix(c(group1_with_feature, group1_without_feature,
                              group2_with_feature, group2_without_feature,
                              group3_with_feature, group3_without_feature),
                            nrow = 3, byrow = TRUE)

colnames(contingency_table) <- c("Interacting", "Non-interacting")
rownames(contingency_table) <- c("Group 1", "Group 2", "Group 3")

# Perform Fisher's Exact Test
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)


#PAIRWISE
group1_with_feature <- 17       #winter
group1_without_feature <- 3.75
group2_with_feature <- 7.75       #summer
group2_without_feature <- 15.75
group3_with_feature <- 7.25       #monsoon
group3_without_feature <- 7


perform_fisher_test <- function(group1_with, group1_without, group2_with, group2_without) {
  contingency_table <- matrix(c(group1_with, group1_without,
                                group2_with, group2_without),
                              nrow = 2, byrow = TRUE)
  
  # Perform Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Return the result
  return(fisher_test_result)
}

# Perform pairwise Fisher's Exact Tests
fisher_test_1_2 <- perform_fisher_test(group1_with_feature, group1_without_feature, group2_with_feature, group2_without_feature)
fisher_test_1_3 <- perform_fisher_test(group1_with_feature, group1_without_feature, group3_with_feature, group3_without_feature)
fisher_test_2_3 <- perform_fisher_test(group2_with_feature, group2_without_feature, group3_with_feature, group3_without_feature)

# Print the results
cat("Fisher's Exact Test for Group 1 vs Group 2\n")
print(fisher_test_1_2)
cat("\nFisher's Exact Test for Group 1 vs Group 3\n")
print(fisher_test_1_3)
cat("\nFisher's Exact Test for Group 2 vs Group 3\n")
print(fisher_test_2_3)




#Stacked barchart
#combined chart
library(ggplot2)
library(tidyr)  # Make sure tidyr package is loaded for gather()

# Define months with desired order
months <- c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", 
            "Dec", "Jan", "Feb", "Mar", "Apr", "May")

# Counts for Group A (Show and Not Show Feature)
counts_group_a_show <- c(4,5,8,12,15,17,16,20,14,5,6,6)       #anemone
counts_group_a_not_show <- c(5,5,8,10,6,3,3,3,9,19,24,11)

# Counts for Group B (Show and Not Show Feature)
counts_group_b_show <- c(2,1,5,6,10,12,17,15,13,9,6,3)      #sponge
counts_group_b_not_show <- c(14,5,6,6,6,2,3,3,3,7,18,18)

# Create data frames for Group A
data_group_a <- data.frame(
  Month = factor(months, levels = months),
  Show_Feature = counts_group_a_show,
  Not_Show_Feature = counts_group_a_not_show,
  Group = "Anemone"
)

# Create data frames for Group B
data_group_b <- data.frame(
  Month = factor(months, levels = months),
  Show_Feature = counts_group_b_show,
  Not_Show_Feature = counts_group_b_not_show,
  Group = "Sponge"
)

# Combine data frames
combined_data <- bind_rows(data_group_a, data_group_b)

# Reshape data to long format for ggplot using gather()
data_long <- combined_data %>%
  gather(key = "Type", value = "Count", Show_Feature, Not_Show_Feature)

# Plot combined bar charts
ggplot(data_long, aes(x = Month, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  facet_wrap(~ Group, scales = "free", nrow = 1, strip.position = "bottom") +
  labs(x = NULL, y = "Count", fill = NULL) +
  scale_fill_manual(values = c("Show_Feature" = "black", 
                               "Not_Show_Feature" = "grey"),
                    labels = c("Non-interacting individuals", 
                               "Interacting individuals")) +
  theme_classic() +
  theme(
    legend.position = "bottom",  
    strip.background = element_blank(),  
    strip.placement = "outside",  
    legend.text = element_text(size = 1.2 * rel(1)), strip.text = element_text(size = 1.2 * rel(1))
  ) 

