# CALLING SCRIPT
#
# by Alica Guzmán 
# (17/07/2024)

# REQUIRES: 
#   1. packages to load 
#   2. ./roiData/roitoroi/AvgCond_3DparticipantStack.mat file containing matrix of dim [#roi] x [1] x [#participants]
#   3. ./iglu.csv containing regression residuals from iglu 


# Load packages
library(R.matlab)
library(readxl)
library(grid)
library(stringr)
library(corrplot)
library(openxlsx)
library(rlang)
library(tidyr)
library(dplyr)
library(reshape2)

# PCA
library(ggplot2)
library(ggfortify)
library(factoextra)

# CV LM model
library(caret)
library(lme4)
library(nlme)
library(dplyr)
library(car)
library(glmnet)

# Faktor Analysis
library(stats)

source("R_helpers_lme.R", local = TRUE)
source("R_plotting_lme.R", local = TRUE)


# Palette for corr_plots
standard_palette <- colorRampPalette(c("darkblue", "white", "darkred"))(200)


# -- OBTAIN AND PREPROCESS AVG HYPOTHALAMUS FC DATA dim [#roi] x [1] x [#par] --

# Load and preprocess FC data
fc_data_conn <- preprocess_data("roiData/roitoroi") # CONN data
# fc_data_xcpd <- preprocess_data("roiData/roi_corr") # new pipeline data

sprintf("Preprocessed Conn Data: missing participants: %s", fc_data_conn$con_par_nan)
# sprintf("Preprocessed XCP-D Data: missing participants: %s", fc_data_xcpd$con_par_nan)

# predictors_xcpd <- as.data.frame(fc_data_xcpd$predictors)
predictors_conn <- as.data.frame(fc_data_conn$predictors)

sprintf("Dimension Conn Data: %s", dim(predictors_conn))
# sprintf("Dimension XCP-D Data: %s", dim(predictors_xcpd))

if (pipeline == "CONN"){
  complete_idx_conn <- fc_data_conn$complete_idx
  con_par_nan <- fc_data_conn$con_par_nan
  predictors <- fc_data_conn$predictors
} else if (pipeline == "XCP-D"){
  complete_idx_conn <- fc_data_xcpd$complete_idx
  con_par_nan <- fc_data_xcpd$con_par_nan
  predictors <- fc_data_xcpd$predictors
} else{
  print("This pipeline is not implemented (l69 in lme_FChypothalamus_iglu.R).")
}


# ---------------------- OBTAIN RELEVANT IGLU DATA -----------------------------
print("Obtain iglu data \n")

# Load IGLU data
iglu_path_prefix = script_dir # str_extract(script_dir, ".*(?=/Thesis_work)")
iglu_file_path <- paste(iglu_path_prefix, "iglu.xlsx", sep = "/")

iglu_data <- read.xlsx(iglu_file_path)

print(colnames(iglu_data))

# Generate random BMI
iglu_data$BMI <- apply(iglu_data, 1, function(x) sample(18:35, 1))

# Get BMI
bmi_iglu <- iglu_data$BMI
bmi <- as.data.frame(bmi_iglu)
colnames(bmi) <- "BMI"
bmi$Glucose_ID <- rownames(bmi)

# Preprocess IGLU data
iglu <- iglu_data[, c(1:ncol(iglu_data)-1)]
iglu[] <- lapply(iglu, function(x) as.numeric(x))
iglu_matrix <- as.matrix(iglu)

cat("Sum of NaN values in IGLU data:", sum(is.na(as.data.frame(iglu))), "\n")
cat("sum of inf values in IGLU data:", sum(is.infinite(iglu_matrix)), "\n")


# --------------------------- FILTER USED DATA ---------------------------------
print("Filter used data \n")

# Load participant data
match_data_path <- paste(script_dir, "participant_map_100par.xlsx", sep = "/") 
participant_data <- read_excel(match_data_path, 
                                 col_names = c("Glucose_ID","CONN_ID"), 
                                 skip = 1
                               )

# Save copies
participant_data_copy <- participant_data
iglu_copy <- iglu 

# Match IGLU rownames with processing pipeline IDs
row.names(iglu) <- c(participant_data$CONN_ID)

# Exclude missing participant 
for (i in con_par_nan){
  participant_data$CONN_ID[participant_data$CONN_ID == i] <- NaN
  cat("Participant with CONN ID ", 
        i, 
        " and Glucose ID ", 
        participant_data$Glucose_ID[participant_data$CONN_ID == i], 
        "has missing values in FC data \n"
     )
}

participant_data <- na.omit(participant_data) 

# Add index to participant table
participant_data$Index <- 1:dim(predictors)[1]

# IDs to exclude
par_exclude_glucose <- c() #c(48,21,29,30,53,85)
par_exclude_conn <- participant_data[participant_data$Glucose_ID 
                                       %in% par_exclude_glucose,]$CONN_ID
# ausschluss par 48: conn id 40
# ausschluss pars that are not identifiable:

# Filter data
cidx <- complete_idx_conn
participant_data_filtered <- participant_data[!participant_data$Glucose_ID 
                                                %in% par_exclude_glucose, ]
cidx <- cidx[!(unlist(cidx) %in% par_exclude_conn)]
cidx <- cidx[!(unlist(cidx) %in% con_par_nan)]

complete_idx_conn <- cidx
complete_idx_index <- participant_data_filtered$Index

predictors <- predictors[complete_idx_index,]

iglu <- iglu[complete_idx_conn,] 

row.names(iglu) <- complete_idx_conn

if(dim(iglu)[1] != dim(predictors)[1]){
  stop("First dim of iglu and FC after nan filtering should be the same (l94).")   
}

# Filter BMI
bmi <- merge(bmi, participant_data, by = "Glucose_ID") %>% # Merge dataframes take only smaller dim
         mutate(Glucose_ID = as.numeric(Glucose_ID)) %>% #Change Glucose_ID col to numeric
         arrange(Glucose_ID) # Sort by Glucose_ID


# -----------------RUN PCA ON IGLU DATA TO OBTAIN OUTCOME----------------------------------------
print("Run PCA on IGLU data")

# Exclude columns from iglu with too low variance
col_exclude <- c()

var_cutoff <- 0.00001

col_vars <- sapply(iglu[,c(1:ncol(iglu))], function(x) var(as.numeric(x)))

for (col in col_vars){
  if (col < var_cutoff){
    col_ <- which(col_vars == col)
    col_exclude <- c(col_exclude, col_)
  }
  else{
    sprintf("No column has var smaller than %f", var_cutoff)
  }
}
sprintf("Excluded columns from iglu data, which have smaller than %f variance: 
          %s", var_cutoff, colnames(iglu)[col_exclude]
        )

if (!is.null(col_exclude)){
  iglu <- iglu[, -col_exclude]
}

# Run PCA / Faktor Analysis
if (var_pred == "PCA"){
  
  # Extract and scale iglu predictors, run pca
  iglu.pca <- prcomp(iglu[, c(1:ncol(iglu))], scale = TRUE)
  ncomp = length(colnames(iglu.pca$x))
  
  # Extract the first 3 principal components and store as dataframe
  pc_iglu <- as.data.frame(iglu.pca$x[, 1:3])

}

if (pca_full_save == TRUE){
  
  # Exclude columns
  iglu_copy <- iglu_copy[,-col_exclude]
  
  # Extract and scale iglu predictors, run pca
  iglu_copy.pca <- prcomp(iglu_copy[, c(1:ncol(iglu_copy))], scale = TRUE)
  
  # Extract the first 3 principal components and store as dataframe
  pc_iglu_copy <- as.data.frame(iglu_copy.pca$x[, 1:3])
  
  pc_iglu_all_full <- as.data.frame(iglu_copy.pca$x)
  
  # Add participant identification
  pc_iglu_all_full$Glucose_ID <- row.names(pc_iglu_all_full)
  pc_iglu_all_full <- merge(participant_data[c("CONN_ID", "Glucose_ID")], 
                              pc_iglu_all_full, 
                              by = "Glucose_ID", 
                              all = TRUE
                            )
  pc_iglu_all_full[, "Glucose_ID"] <- sapply(pc_iglu_all_full[, "Glucose_ID"], 
                                               as.numeric
                                             )
  pc_iglu_all_full <- pc_iglu_all_full[order(pc_iglu_all_full$Glucose_ID),]
  
  # Save
  write.xlsx(pc_iglu_all_full, "pca_iglu.xlsx")
}


# ---------------- RUN FACTOR ANALYSIS TO OBTAIN OUTCOME -----------------------
print("Run Factor Analysis on IGLU data")

# Compute IGLU correlation matrix
iglu_corr <- cor(iglu)

# Highly correlated predictors
sprintf("IGLU/Outcome correlated over 0.9 and under -0.9: %s", 
          find_correlations(iglu_corr, cutoff = 0.9, how = "high")
        ) 
  
# Exclude too correlated columns for Factor Analysis
if (var_pred == "FA"){

  # Set cutoff
  corr_cutoff = 0.95
  
  # Columns to exclude
  col_ex <- unique(find_correlations(iglu_corr, cutoff = corr_cutoff, how = "abs")$col) 
  sprintf("excluded columns from iglu data for factor analysis, 
            which have higher correlation than %f:", corr_cutoff
          )
  print(col_ex)
  
  col_ex <- match(col_ex, colnames(iglu))
  
  # Filter IGLU data
  iglu_forfa <- iglu[,-col_ex]
  
  # Run Factor Analysis
  iglu_fa <- factanal(as.matrix(iglu_forfa), factors = 1, scores = "regression", rotation = "none")
  print(iglu_fa)
  

  # Extract the loadings for the first (and only) factor
  loadings <- iglu_fa$loadings[,1] 
  
  # Plot Factor Loadings
  png(file=paste0(script_dir, "/plots/", "factor_loadings_iglu.png"))
   
  # Create a bar plot
  barplot(sort(loadings, decreasing = TRUE), 
            las = 2,  # Make the labels perpendicular to the axis
            col = "skyblue", 
            main = "Factor Loadings for Factor 1",
            xlab = "Variables", 
            ylim = c(-1,1),
            cex.axis = 0.75,
            cex.names = 0.75,
            ylab = "Loadings"
          )
   graphics.off()
}


# ----------------------------- DEFINE OUTCOME ---------------------------------
print("Define outcome")

# Define outcome
if (var_pred == "PCA"){
  outcome <- data.frame("outcome" = c(pc_iglu$PC1))
  } else if (var_pred == "FA"){
  outcome <- data.frame("outcome" = c(iglu_fa$scores))
  } else if (var_pred == "BMI"){
  outcome <- data.frame("outcome" = c(bmi$BMI))
    
  }

#print(which(is.na(outcome), arr.ind = TRUE)) # par 31 & 34 missing FC data


# ------ RUN LM PREDICTING PCA COMPONENTS WITH FC WITH HYPOTHALAMUS SEED ------
# --------------------- WITH 10-FOLD CROSS VALIDATION --------------------------
print("Run LM")

# Compute correlation matrix of predictors
pred_corr <- cor(predictors)

# Specify Model

# Regions of Interest (ROI)
seed_names <- colnames(predictors) 

# Add each ROI as predictor
fm <- paste(sapply(colnames(predictors), function(x) paste0("`", x, "`")), 
              collapse = " + "
            )
model_formula <- as.formula(paste("outcome ~", fm))

# Merge outcome with predictors
processed_data <- as.data.frame(cbind(predictors, outcome = outcome))

# Fit the model
if (fit_model == TRUE){
  
  # Initialize performance lists
  perf_list <- numeric(n_seeds)
  perf_list2 <- numeric(n_seeds)
  perf_list3 <- numeric(n_seeds)
  
  # Predictor names
  columns <- colnames(processed_data) 
  columns <- c("(Intercept)",  columns[columns != "outcome"])
  
  # Initialize coefficients dataframe and fill it with 0
  nzcoefs_df <- as.data.frame(matrix(ncol = length(columns), nrow = n_seeds))
  nzcoefs_df[is.na(nzcoefs_df)] = 0
  colnames(nzcoefs_df) <- columns
  nzcoefs_df$bestLambda <- NaN
  
  # Initialize performances dataframe and fill it with 0
  performances_models_df <- as.data.frame(matrix(ncol = 0, nrow = n_seeds))

  # Initialize list to store predictions
  predictions_list <- list()
  
  # For each seed
  for (s in 1:n_seeds){
    set.seed(s)
    
    # Define 10-fold elastic net cross-validation and fit model to data
    result <- run_cv_elasticnet(processed_data, 
                                  model_formula, 
                                  nfolds = 10, 
                                  validation_metric = "RMSE", 
                                  fix_lambda = FALSE
                                )
    
    cv_model <- result$cv_model
    
    # Store predictions
    predictions_list[[s]] <- result$predictions$Predictions  # Extract the predictions for the current fit
    
    # Extract the best lambda value
    best_lambda <- cv_model$bestTune$lambda
    nzcoefs_df$bestLambda[s] <- best_lambda
  
    # Access non zero coefficients
    non_zero_coef <- result$non_zero_coef
     
    # Store them in nz_coef_df
    i = 1
    for (var in non_zero_coef$Predictor) {
      var <- gsub("[`´]", "", var)
      if (var %in% colnames(nzcoefs_df)) {
        nzcoefs_df[s, var] <- non_zero_coef$s1[i]
      }
      else{
        if (n_seeds > 5){
        print("Column not in df.l332.")
        print(var)
        }
      }
      i = i+1
    }
  
    # Store performances
    print(paste0("done computing seed ", s))
    performance2 <- mean(cv_model$resample$Rsquared)
    perf_list2[s] <- performance2
  
    performance3 <- mean(cv_model$resample$MAE)
    perf_list3[s] <- performance3
  
    performance <- mean(cv_model$resample$RMSE)
    perf_list[s] <- performance
    
  }
  
  # Combine predictions into a data frame and save
  combined_predictions <- as.data.frame(do.call(cbind, predictions_list))
  combined_predictions_path <- paste0(output_dir, 
                                        "10fold/model_predictions_",
                                        n_seeds,
                                        "seeds.xlsx"
                                      )
  write.xlsx(combined_predictions, file.path(getwd(), combined_predictions_path))
  
  # Calculate average predictions for each observation
  average_predictions <- rowMeans(combined_predictions)
  
  # Dataframe to store outcome and predictions
  df_comp <- data.frame(processed_data$outcome)
  colnames(df_comp) <- c("Outcome")
  
  # Add average predictions to the original data frame for comparison
  df_comp$Average_Predictions <- average_predictions
  
  # Evaluate performance (calculating RMSE)
  rmse <- sqrt(mean((df_comp$Outcome - df_comp$Average_Predictions)^2))  
  print(paste("Average RMSE over 1000 iterations:", round(rmse, 2)))
  
  print(nzcoefs_df)
  print(cor(t(nzcoefs_df)))
  
  # Store performances
  nzcoefs_df$Performance <- perf_list
  
  performances_models_df$MAE <- perf_list3
  performances_models_df$RMSE <- perf_list
  performances_models_df$RSquared <- perf_list2
  
  # Save performances
  model_performances_path <- paste0(output_dir, 
                                      "10fold/model_performances_",
                                      n_seeds,
                                      "seeds.xlsx")
  write.xlsx(performances_models_df, file.path(getwd(), model_performances_path))
  
  # Save non-zero coefficients
  coef_path <- paste0(script_dir, output_dir, "10fold/coefficients_",n_seeds,"seeds.xlsx")
  write.xlsx(nzcoefs_df, file.path(coef_path))
  
  # Print coefficients that are never 0 across seeds
  print(colnames(nzcoefs_df %>% select_if(colSums(.) != 0)))
  
  # Store them in column
  non_zero_columns <- nzcoefs_df[, apply(nzcoefs_df, 2, function(col) all(col != 0))]
  
  # For single model use
  #predictions <- predict(cv_model, processed_data, s = cv_model$lambda)
  #predictions <- as.data.frame(predictions)
  # Calculate performance metrics (e.g., RMSE, R-squared)
  #performance <- postResample(pred = predictions, obs = outcome)
  #print(performance)
} 

# Read in data 
coef_path_fit <- paste0(script_dir, output_dir, "10fold/coefficients_",n_seeds,"seeds.xlsx")
nzcoefs_df <- read.xlsx(file.path(coef_path_fit))
  
# Get variables
perf_list <- nzcoefs_df$Performance
bestLambdas <- nzcoefs_df$bestLambda
a <- nzcoefs_df[,c("Performance", "bestLambda")]
  
# Extract the bestLambda by picking the one which is fitted most often
best_lambda <- Mode(bestLambdas)#a[median_index, "bestLambda"]
print(best_lambda)
  
nzcoefs_df <- select(nzcoefs_df, -c(Performance, bestLambda)) 

print('b')

# Extract predictions
fn_pred <- paste0(script_dir, output_dir, "10fold/model_predictions_", n_seeds,"seeds.xlsx")
combined_predictions <- read.xlsx(file.path(fn_pred))                              
print('after')
  
# Calculate average predictions for each observation
average_predictions <- rowMeans(combined_predictions)

df_comp <- data.frame(processed_data$outcome)
  
colnames(df_comp) <- c("Outcome")
  
# Add average predictions to the original data frame for comparison
df_comp$Average_Predictions <- average_predictions



if (plot_sub_models == TRUE){
  
  colnames(df_comp) <- c("Outcome", "Average_Predictions")
  
  # Create the dot plot using ggplot2
  plot <- ggplot(df_comp, aes(x = Outcome, y = Average_Predictions)) +
    geom_point(size = 1.5) +
    geom_text(aes(label = participant_data$Glucose_ID), 
                    vjust = -0.5, 
                    hjust = 0.5, 
                    size = 2.5, 
                    color = "blue"
              ) +  # Add labels with subject IDs
    geom_smooth(method = "lm", color = "red", se = FALSE) +  # Regression line
    labs(title = sprintf("Average Model Predictions from %s Seeds", n_seeds), 
           x = "Iglu PC1", 
           y = "Average Prediction"
         ) +
    theme_classic() +
    theme(
      plot.background = element_rect(fill = "white"),  # Set the plot background to white
      panel.background = element_rect(fill = "white"),  # Set the panel background to white
      legend.position = "none"  # Remove legend
    ) +  
    scale_x_continuous(limits = c(min(c(df_comp$Outcome, 
                                          df_comp$Average_Predictions)
                                      ), 
                                  max(c(df_comp$Outcome, 
                                          df_comp$Average_Predictions)
                                      )
                                  )
                       ) +  # Set x-axis limits
    scale_y_continuous(limits = c(min(c(df_comp$Outcome, 
                                          df_comp$Average_Predictions)
                                      ), 
                                  max(c(df_comp$Outcome, 
                                          df_comp$Average_Predictions)
                                      )
                                  )
                       ) +  # Set y-axis limits
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +  # Add a diagonal line
    coord_fixed()  # Ensure equal scaling
  
  # Save plot
  ggsave(paste0(script_dir,output_dir,"10fold/model_predictions/final_model_",n_seeds,".png"), 
           plot, 
           width = 6, 
           height = 6, 
           bg = "white"
         )
  
  if (plot_singles == TRUE){
    
    for (i in 1:dim(combined_predictions)[2]){
      
      df_i <- as.data.frame(cbind(processed_data$outcome, combined_predictions[, i]))
      colnames(df_i) <- c("Outcome", "Average_Predictions")
      
      # Create the dot plot using ggplot2
      plot <- ggplot(df_i, aes(x = Outcome, y = Average_Predictions)) +
        geom_point(size = 1.5) +
        geom_text(aes(label = participant_data$Glucose_ID), 
                        vjust = -0.5, 
                        hjust = 0.5, 
                        size = 2.5, 
                        color = "blue") +  # Add labels with subject IDs
        geom_smooth(method = "lm", color = "red", se = FALSE) +  # regression line
        labs(title = sprintf("Average Model Predictions from %s Seeds", n_seeds), 
               x = "Iglu PC1", y = "Average Prediction"
             ) +
        theme_classic() +
        theme(
          plot.background = element_rect(fill = "white"),  # Set the plot background to white
          panel.background = element_rect(fill = "white"),  # Set the panel background to white
          legend.position = "none"  # Remove legend
        ) +  # Make sure to add a comma here
        scale_x_continuous(limits = c(min(c(df_i$Outcome, df_i$Average_Predictions)), 
                                        max(c(df_i$Outcome, df_i$Average_Predictions)))) +  # Set x-axis limits
        scale_y_continuous(limits = c(min(c(df_i$Outcome, df_i$Average_Predictions)), 
                                        max(c(df_i$Outcome, df_i$Average_Predictions)))) +  # Set y-axis limits
        geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +  # Add a diagonal line
        coord_fixed()  # Ensure equal scaling
    
      ggsave(paste0(script_dir,output_dir,"10fold/model_predictions/sub_models/model_s",i,".png"), 
               plot, 
               width = 6, 
               height = 6, 
               bg = "white"
             )
    }
  }
}

if (fit_final_model == TRUE){
  set.seed(1)
  
  # Fit model with final model (best_lambda)
  full_model <- run_cv_elasticnet(processed_data, 
                                    model_formula, 
                                    nfolds = 10, 
                                    validation_metric = "RMSE", 
                                    fix_lambda = best_lambda
                                  )
  
  # Predict
  single_model_pred <- predict(full_model$cv_model, newdata = processed_data)
  single_model_pred <- as.data.frame(single_model_pred)
  print(single_model_pred)
  
  # Bind with outcomes
  merged_data <- as.data.frame(cbind(outcome$outcome, single_model_pred$single_model_pred))
  colnames(merged_data) <- c("Outcome", "Predictions")
  
  # Create the dot plot using ggplot2
  plot <- ggplot(merged_data, aes(x = Outcome, y = Predictions)) +
    geom_point(size = 1.5) +
    geom_text(aes(label = participant_data$Glucose_ID), 
                vjust = -0.5, 
                hjust = 0.5, 
                size = 2.5, 
                color = "blue"
              ) +  # Add labels with subject IDs
    geom_smooth(method = "lm", color = "red", se = FALSE) +  # regression line
    labs(title = sprintf("Model fit with \n alpha = %.2f, lambda = %f", 0.5, best_lambda), 
           x = "Iglu PC1", 
           y = "Prediction"
         ) +
    theme_classic() +
    theme(
      plot.background = element_rect(fill = "white"),  # Set the plot background to white
      panel.background = element_rect(fill = "white"),  # Set the panel background to white
      legend.position = "none"  # Remove legend
    ) +  # Make sure to add a comma here
    scale_x_continuous(limits = c(min(c(merged_data$Outcome, merged_data$Predictions)), 
                                    max(c(merged_data$Outcome, merged_data$Predictions)))) +  # Set x-axis limits
    scale_y_continuous(limits = c(min(c(merged_data$Outcome, merged_data$Predictions)), 
                                    max(c(merged_data$Outcome, merged_data$Predictions)))) +  # Set y-axis limits
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +  # Add a diagonal line
    coord_fixed()  # Ensure equal scaling
  
  # Show the plot
  print(plot)
  ggsave(paste0(script_dir,output_dir,"10fold/model_predictions/fix_lambda_model_based_on_mode.png"), 
           plot, 
           width = 6, 
           height = 6, 
           bg = "white"
         )

}


 # -------------------------- PERMUTATION TEST ---------------------------------
print("Run permutation test")

if(perform_permutation_test == TRUE){
  set.seed(1)
  
  # Compute permutation test
  permutation_test_results <- permutation_test(processed_data, 
                                                 model_formula = model_formula, 
                                                 nfolds = 10, 
                                                 validation_metric = "RMSE", 
                                                 n_permutations = n_seeds)
  
  # Save
  write.xlsx(as.data.frame(permutation_test_results), 
               file.path(getwd(), paste0(output_dir,"10fold/",n_seeds,"permutation_results.xlsx"))
             )
  
  permutation_test_results <- permutation_test_results$permuted_performance
  # save pca with conn identification
  #pc_iglu_all <- as.data.frame(iglu.pca$x)
  #pc_iglu_all$CONN_ID <- row.names(pc_iglu_all)
  #pc_df <-merge(participant_data, pc_iglu_all, by = "CONN_ID")
  #write.xlsx(pc_df, "pca_iglu_old.xlsx")
  
} else if (perform_permutation_test == FALSE){
  
  # Read in data
  file_name <- file.path(paste0(script_dir, output_dir, "10fold/",n_seeds,"permutation_results.xlsx"))
  permutation_test_results <- read.xlsx(file_name)
  permutation_test_results <- permutation_test_results$permuted_performance 
}

print(mean((mean(perf_list) > permutation_test_results)))
print(sprintf("Observed Performance Mean: %f", mean(perf_list))) 
print(sprintf("Permuted Performance Mean: %f", mean(permutation_test_results))) 


# -------------------------------- PLOTTING ------------------------------------
print("Plotting")

#data_corr <- sapply(1:nrow(as.matrix(predictors_xcpd[1:79,])), function(i) cor(as.matrix(predictors_xcpd[1:79,])[i,], as.matrix(predictors_conn)[i,]))

if (plot_corrplots == TRUE){
  # remove the one par which is not yet in XCP pipeline
  #predictors_conn <- predictors_conn[-c(participant_data[participant_data$Glucose_ID == 42,]$Index),]
  
  #m <- create_corr_matrix(predictors_xcpd, predictors_conn, threshold = 0.6, write = TRUE)
  #print(m)
  
  # plot heatmop of mse between pipelines
  # h <- plot_mse(predictors_xcpd, predictors_conn, output_dir)

  # Compare the two pipelines
  # results <- compare_pipelines(predictors_xcpd, predictors_conn)

  # Print the results
  # print("Correlation values for each ROI:")
  # print(results$correlation)
  # print("Mean Squared Error for each ROI:")
  # print(results$mse)
  
  # Ensure the number of rows in both matrices is the same
  # if (nrow(predictors_xcpd) != nrow(predictors_conn)) {
    # stop("The number of participants (rows) in the two matrices is not the same.")
  # }
}
if(plot_all == TRUE){
    # plot_roi_corr_between_pipeline <- function(predictors_xcpd, predictors_conn, low_corr_df)
  plot_pca_results(iglu.pca, pc_iglu, outcome, output_dir)
    
  plot_correlation_outcome_predictors(output_dir, iglu, predictors)
  
  plot_coef_consistency(output_dir, nzcoefs_df, n_seeds)
    
  # plot_each_coef_distr(output_dir, nzcoefs_df)
    
  plot_all_coef_hist(output_dir, nzcoefs_df, nseeds=n_seeds)
    
  #plot_prediction_final_model(output_dir, combined_predictions, outcome, model_formula, cv_model, s)
    
  plot_performance_across_seeds(output_dir, perf_list, n_seeds)
    
  plot_permuted_performance(output_dir, permutation_test_results, n_seeds)
  
  plot_best_lambda(output_dir, bestLambdas)
} 

# corr bmi
#for (i in 1:dim(processed_data)[2]){
 # print(colnames(processed_data)[i])
  #print(cor(processed_data$outcome, subset(processed_data, select = - c(outcome))[, i]))
#}

