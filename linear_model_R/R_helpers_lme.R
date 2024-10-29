# HELPER FUNCTIONS FOR LM IGLU PREDICITIONS FROM FC
#
# by Alica Guzmán 
# (13/08/2024)

library(corrplot)
library(stringr)


# Preprocess FC data
preprocess_data <- function(folder){
  
  # Get avg cond .mat file
  mat_file_path <- paste(script_dir, folder, "AvgCond_3DparticipantStack.mat", sep = "/")
  mat_data <- readMat(mat_file_path)
  
  # Extract the matrix (variable name is matrix3D)
  matrix_data <- mat_data$matrix3D
  matrix_data[is.infinite(matrix_data)] <- NaN #inf was nan! in new data
  
  arr_idx <- as.data.frame(which(is.nan(matrix_data), arr.ind = TRUE))
  
  nan_counts_per_dim3 <- arr_idx %>%
    group_by(dim3) %>%
    summarise(nan_count = n()) 
  print(nan_counts_per_dim3, n = 100)
  
  # Extract hypothalamus data (left right conn seed 147 168 respectively)
  fc_hypothalamus_left_right <- matrix_data[, c(147,148), ]
  
  # Avg FC of left and right hypothalamus
  fc_hypothalamus_avg <- array(apply(fc_hypothalamus_left_right, MARGIN = c(1, 3), 
                                       FUN = function(x) mean(x, na.rm = TRUE)), 
                                 dim = c(dim(fc_hypothalamus_left_right)[1], 
                                           1, 
                                           dim(fc_hypothalamus_left_right)[3])
                               )
  
  # Define predictors (reshape fc_hypothalamus_avg), drop singleton dimension
  predictors_all <- drop(aperm(fc_hypothalamus_avg, c(3, 1, 2)))
  
  # Import the roi labels from excel table
  roi_labels <- read_excel(paste(script_dir, "CONN_labels_Atlas.xlsx", sep = "/"), 
                             col_names = c("Index", "Region"), skip=1
                           )
  
  # Process the matrix by averaging "l" and "r" columns and renaming them
  processed_matrix <- predictors_all
  
  # Get columnames by region
  new_colnames <- roi_labels$Region
  new_colnames <- gsub("\\s*\\(.*?\\)", "", new_colnames) # Remove parentheses and their content
  new_colnames <- new_colnames[new_colnames != "Hypothalamus l"]  # Remove left hypothalamus
  new_colnames[147] <- "Hypothalamus avg" # rename right hypothalamus to avg since this was taken before
  
  # Initialize counters, FC matrix and new columnames
  c_avg = 0
  c_no = 0
  fc_avg <- NULL
  col_names <- c()
  
  # Iterate through unique regions without "l" or "r"
  for (region in unique(gsub(" [lr]$", "", new_colnames))) {
    left_col <- which(new_colnames == paste0(region, " l"))
    right_col <- which(new_colnames == paste0(region, " r"))
    
    # When there are left and right columns
    if (length(left_col) > 0 & length(right_col) > 0) {
      c_avg <- c_avg+1
      # Compute avg and append it to new matrix: fc_avg
      avg_col <- rowMeans(processed_matrix[, c(left_col, right_col)], na.rm = TRUE)
      fc_avg <- cbind(fc_avg, avg_col)
      col_names <- c(col_names, paste0(region, " avg"))
    }
    # Otherwise just append column
    else{
      c_no <- c_no+1
      col_name <- which(new_colnames == region)
      fc_avg <- cbind(fc_avg, processed_matrix[,col_name])   
      col_names <- c(col_names, region)
    }
    
  }
  
  # Update the column names of predictors
  colnames(fc_avg) <- col_names
  predictors_all <- fc_avg
  
  # Remove nan rows of predictors from outcome and predictors
  complete_idx <- which(complete.cases(predictors_all)) # get non-nan rows
  con_par_nan <- setdiff(1:dim(predictors_all)[1], complete_idx)
  predictors <- predictors_all[complete_idx, ] # remove corresponding entries for FC nan values
  
  # Create result list
  result_list <- list(predictors = predictors,
                        complete_idx = complete_idx,
                        con_par_nan = con_par_nan,
                        arr_idx = arr_idx
                      )
  return (result_list)
}


# Give correlated pairs over/under certain cutoff
find_correlations <- function(corr_matrix, cutoff = 0.8, how = "high"){
  
  # Find the indices where the correlation exceeds the cutoff
  if (how == "high"){ # greater than cutoff
    corr_indices <- which(corr_matrix > cutoff, arr.ind = TRUE)
  }
  else if (how == "low"){ # smaller than cutoff
    corr_indices <- which(corr_matrix < cutoff, arr.ind = TRUE)
  }
  else if (how == "abs"){ # absolute cutoff (smaller than and greater than abs number)
    corr_indices <- which(abs(corr_matrix) > abs(cutoff), arr.ind = TRUE)
  }
  
  # Remove the diagonal elements (where correlation is always 1)
  corr_indices <- corr_indices[corr_indices[, 1] != corr_indices[, 2], ]
  
  # Get the row and column names corresponding to those indices
  corr_pairs <- data.frame(
    row = rownames(corr_matrix)[corr_indices[, 1]],
    col = colnames(corr_matrix)[corr_indices[, 2]],
    corr_value = corr_matrix[corr_indices]
  )
  
  return (corr_pairs)
}


# Define a function to split the title
split_title <- function(title, max_length = 50) {
  # Split the title into words
  words <- str_split(title, " ")[[1]]
  lines <- c()
  current_line <- ""
  
  for (word in words) {
    # Check if adding the word exceeds the max_length
    if (nchar(current_line) + nchar(word) + 1 > max_length) {
      # Add current line to lines and start a new line
      lines <- c(lines, current_line)
      current_line <- word
    } else {
      # Add word to current line
      if (nchar(current_line) > 0) {
        current_line <- paste(current_line, word)
      } else {
        current_line <- word
      }
    }
  }
  # Add the last line
  lines <- c(lines, current_line)
  
  # Combine lines with newline characters
  return(paste(lines, collapse = "\n"))
}


# Define Grid for hyperparameter in elastic net
lambda_grid <- 10^seq(-4, 2, length = 100)

# Define a function to do a permumation test on elastic net (alpha=0.5) cross validation regression model
permutation_test <- function(data, model_formula, nfolds, validation_metric, n_permutations) {
  # data (for prediciton) - type dataframe
  # outcome (to be predicted) - type dataframe (single column)
  # validation_metric - str
  # Function to perform 10-fold CV and calculate performance metric
  
  # Calculates performance
  cv_performance <- function(df, model_formula, nfolds, validation_metric) {
    
    cv_model <- train(model_formula, data = df,
                        method = "glmnet",
                        trControl = trainControl(method = "cv", number = nfolds, 
                                      preProc = c("center", "scale"),  # Center and scale within each fold
                                      savePredictions = TRUE
                                                ),
                        tuneGrid = expand.grid(alpha = 0.5, lambda = lambda_grid), #seq(0.001, 5, length = 500)),  # Fixed alpha, tune lambda
                        preProcess = c("center", "scale"),  # Center and scale each fold
                        metric = validation_metric
                      )
    
    return(mean(cv_model$resample[[validation_metric]]))
  }
  
  # 1. Train the model and calculate observed performance
  observed_performance <- cv_performance(data, model_formula, nfolds, validation_metric)
  
  # 2. Permutation test
  permuted_performance <- numeric(n_permutations)
  
  # Store perumuted performances in list
  for (i in 1:n_permutations) {
    cat("permutation ", i, "\n")
    permuted_df <- data
    permuted_df$outcome <- sample(permuted_df$outcome)  # Shuffle the target variable
    permuted_performance[i] <- cv_performance(permuted_df, 
                                                model_formula, 
                                                nfolds, 
                                                validation_metric
                                              )
  }
  
  # 3. Calculate p-value
  p_value <- mean(permuted_performance >= observed_performance)
  
  # Print results
  cat("Observed Performance (",validation_metric, "):", 
        observed_performance, "\n")
  cat("Mean Permuted Performance  (",validation_metric, "):", 
        mean(permuted_performance), "\n")
  cat("p-value:", p_value, "\n")
  
  # Create result list
  result_list <- list(permuted_performance = permuted_performance,
                        observed_performance = observed_performance,
                        p_value = p_value
                      )
  
  return(result_list)
}


# Run elastic net cross validation 
run_cv_elasticnet <- function(df, model_formula, nfolds, validation_metric, fix_lambda) {
  # df is concatenated data and outcome dataframe

  # define folds
  cv_folds <- trainControl(method = "cv", 
                             number = nfolds, 
                             preProc = c("center", "scale"),  # Center and scale within each fold
                             verboseIter = TRUE,
                             savePredictions = "final")

  # Initialize an empty data frame to store predictions
  predictions_df <- data.frame()
    
  # Train model using glmnet
  if (fix_lambda == FALSE){
    cv_model <- train(model_formula, data = df,
                        method = "glmnet",
                        trControl = cv_folds,
                        tuneGrid = expand.grid(alpha = 0.5, lambda = lambda_grid),  # Fixed alpha, tune lambda
                        preProcess = c("center", "scale"),  # Center and scale each fold
                        metric = validation_metric  # Metric to optimize
                      )  
  }
  else{  
    cv_model <- train(model_formula, data = df,
                        method = "glmnet",
                        trControl = cv_folds,
                        tuneGrid = expand.grid(alpha = 0.5, lambda = fix_lambda),  # Fixed alpha, tune lambda
                        preProcess = c("center", "scale"),  # Center and scale each fold
                        metric = validation_metric  # Metric to optimize
                      ) 
  }

  # Extract the best model from the `train` object
  best_model <- cv_model$finalModel

  # Extract the best lambda value
  best_lambda <- cv_model$bestTune$lambda

  # Get coefficients for the best lambda
  coefficients <- coef(best_model, s = best_lambda)

  # Convert to a data frame for easier handling
  coef_df <- as.data.frame(as.matrix(coefficients))
  coef_df$Predictor <- rownames(coef_df)

  # Filter non-zero coefficients
  non_zero_coef <- coef_df[coef_df$s1 != 0, ]
  print(cv_model$results[cv_model$results$lambda == best_lambda, ])
  print(cv_model$resample[])

  # Access predictions
  all_predictions <- cv_model$pred
  
  all_predictions$Rescaled_Predictions <- NA
  
  # Loop through each fold and compute the mean and standard deviation of outcomes
  for (i in unique(all_predictions$Resample)) {
    fold_pred <- all_predictions[all_predictions$Resample == i, ]  # Get predictions for the current fold
  
    # Extract the indices for the test set in this fold
    test_indices <- as.numeric(row.names(fold_pred))
    
    # Calculate mean and standard deviation of the actual outcomes on the test set
    actual_outcomes <- df[test_indices, c("outcome")]    
    mean_outcome <- mean(actual_outcomes)
    sd_outcome <- sd(actual_outcomes)
    
    rescaled_predictions <- (fold_pred$pred * sd_outcome) + mean_outcome
    
    all_predictions$Rescaled_Predictions[all_predictions$Resample == i] <- rescaled_predictions
  }

  # Combine with actual values
  predictions_df <- data.frame(
    Actual = df[all_predictions$rowIndex, "outcome"],  # Ensure you"re accessing the correct rows
    Predictions = all_predictions$Rescaled_Predictions 
  )
  
  # Save results in list
  result_list <- list(
    cv_model = cv_model,
    non_zero_coef = non_zero_coef,
    predictions = predictions_df
    )

  return(result_list)
}


# Compute Mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}