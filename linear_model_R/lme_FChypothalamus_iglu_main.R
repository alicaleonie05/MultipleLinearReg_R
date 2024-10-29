# CALLING SCRIPT
#
# Author: Alica Guzmán 
# (24/10/2024)
#
# Fits multiple linear regression model specified in model_formula to data using 10 fold cross-validation and returns model comparison metrices
#
# REQUIRES: 
#   1. packages to load 
#   2. .roiData/roitoroi/ containing CONN 3D FC Data for all participants and .roiData/roi_corr/ for XCP-D pipeline
#   3. ./iglu.csv containing regression residuals from iglu
#   4. participant_map_100par.xlsx in working directory matching Glucose and Conn ID 
#   5. lme_FChypothalamus_iglu_main.R, R_helpers_lme.R and R_plotting_lme.R in working directory for processing


# Clear all variables
rm(list = ls())

# Set the maximum number of rows to print
options(max.print = 10000)

# Get the directory of the current script file
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

# Set output directory
output_dir <- "/plots/data_CONN_100par_predictPCA/"
output_path <- paste0(script_dir, output_dir, "10fold/coefs")
output_path_mp <- paste0(script_dir, output_dir, "10fold/model_predictions/sub_models")

# Create the directory only if it doesn"t exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
if (!dir.exists(output_path_mp)) {
  dir.create(output_path_mp, recursive = TRUE)
}

# For Reproducibility
s <- 1
set.seed(s) 

#------------------------------ SET VARIABLES ----------------------------------

# SCRIPT STEPS
# 1. obtain and preprocess data from rs-fmri with hypothalamus seed (matrix: [#roi] x [1] x [#participants]) 
# 2. obtain and preprocess iglu data
# 3. filter both data
# 4. run pca or fa on iglu data (outcome)
# 5. run lme with elastic net 10-fold cv to predict outcome from rs-fmri data
# 6. permutation test
# plotting


# Set variables

# Decide on variable to predict 
var_pred <- "PCA" # "FA" for Factor Analysis or "PCA" for Principal Component Analysis

pipeline <- "CONN" # Choose pipeline, either 'CONN' or "XCP-D"
pca_full_save = FALSE # Choose whether to save PCA on full data or not
n_seeds <- 100 # Number of Seeds
fit_model <- FALSE # Fit the model
plot_sub_models <- FALSE # Plot average predictions
plot_singles <- FALSE # Plot Single seed predictions (Only possible if plot_sub_model == TRUE)
fit_final_model <- FALSE # Plot predictions of model fitted to full data with lambda of mean RMSE from #n_seeds seeds
perform_permutation_test <- FALSE # Perform permutation test
plot_corrplots <- FALSE # Plot correlation figures (independent of number of folds in Cross Validation)
plot_all <- TRUE # Plot all figures

# Create String of variable that is predicted
if (var_pred == "PCA"){
  str_var <- "PC1" 
} else if (var_pred == "FA"){
  str_var <- "Factor1"
} else if (var_pred == "BMI"){
  str_var <- "BMI"
}

source("lme_FChypothalamus_iglu.R", local = TRUE)



