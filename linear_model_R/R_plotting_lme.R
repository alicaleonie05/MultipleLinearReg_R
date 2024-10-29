# FUNCTIONS FOR PLOTTING
#
# by Alica Guzmán 
# (17/07/2024)


# Load the R packages
library(readxl)
library(corrplot)
library(openxlsx)
library(tidyr)
library(dplyr)
library(reshape2)
library(stats)  
library(ggplot2)  
library(ggfortify)
library(factoextra)

# Get the directory of the current script file
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

create_corr_matrix <- function(predictors_all_new, predictors_all_old, threshold, write) {
  
  # Create correlation matrix between predictors
  data_corr <- sapply(1:nrow(as.matrix(predictors_all_new)), 
                        function(i) cor(as.matrix(predictors_all_new)[i, ], 
                                          as.matrix(predictors_all_old)[i,]
                                        )
                      )
  # Define cutoff
  cutoff_corr <- threshold
  low_corr_par <-  which(data_corr <= cutoff_corr)
  
  # Filter participants (if necessary)
  participant_data_old <- participant_data#[1:79,]
  low_corr_df <- participant_data_old[low_corr_par,]
  low_corr_df$corr <- data_corr[data_corr <= cutoff_corr]
  low_corr_df$index <- low_corr_par
  
  # Save
  if (write == TRUE){
    write.xlsx(low_corr_df, paste0(getwd(), output_dir, "ParticipantMap_Corr.xlsx"))}
  return(low_corr_df)  
}


compute_mse <- function(actual, predicted) {
  # Calculate the Mean Squared Error for each row
  mse <- rowMeans((actual - predicted)^2)
  return(mse)
}


plot_mse <- function(predictors_all_new, predictors_all_old, output_dir){
  
  # Create correlation matrix
  low_corr_df <- create_corr_matrix(predictors_all_new, predictors_all_old, 1, write = FALSE)
  low_corr_par <- low_corr_df$index
  
  # Compute MSE for each participant
  mse_per_participant <- compute_mse(as.matrix(predictors_all_new), as.matrix(predictors_all_old))
  
  # Create dataframe with Glucose ID
  df2 <- data.frame(mse_per_participant)
  df2$Glucose_ID <- low_corr_df$Glucose_ID
  
  # Save dataframe
  write.xlsx(df2, paste0(script_dir, output_dir, "MSE_oldVSnewPipeline.xlsx"))
  
  # Compute MSE
  df <- (predictors_all_new[low_corr_par,] - predictors_all_old[low_corr_par,])**2
  row.names(df) <- NULL
  
  # Reshape the DataFrame into a long format for ggplot2
  df_melt <- melt(as.matrix(df))
  
  # Add Glucose_ID to the data for labeling (create a new column)
  df_melt$Glucose_ID <- low_corr_df$Glucose_ID
  graphics.off()
  filtered_df_melt <- df_melt[!is.na(df_melt$Glucose_ID), ]
  
  # Create a specific mapping for Glucose_IDs (assuming there are 4 unique Glucose_IDs)
  filtered_df_melt$Glucose_ID <- factor(filtered_df_melt$Glucose_ID, levels = low_corr_df$Glucose_ID)
  
  # Create the heatmap using ggplot2
  heatmap_plot <- ggplot(filtered_df_melt, aes(x = Var2, y = Glucose_ID, fill = value)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "red") +  # Customize color gradient
   labs(title = "MSE between CONN and new pipeline for low corr participants",
          x = "Columns", 
          y = "Glucose ID", 
          fill = "MSE"
        ) +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),  # Rotate x-axis labels
           axis.text.y = element_text(size = 6),  # Adjust y-axis (Glucose ID) label size
           plot.title = element_text(hjust = 0.5),  # Center the title
           plot.margin = margin(10, 10, 10, 20)    # Adjust margins for labels
         ) +
   coord_fixed(ratio = 1)  # Ensure tiles remain square
  
  # Save
  filename <- paste0(".",output_dir, "heatmap_low_corr_allpar.png")
  ggsave(filename, heatmap_plot, width = 12, height = 6, dpi = 300, bg="white")
  return(heatmap_plot)
}


# Function to compute correlation and MSE for each ROI
compare_pipelines <- function(pipe1, pipe2) {
  if (!all(dim(pipe1) == dim(pipe2))) {
    stop("Both pipelines must have the same dimensions.")
  }
  
  # Calculate Pearson correlation for each ROI
  correlation_values <- apply(pipe1, 
                                2, 
                                function(x) cor(x, 
                                                  pipe2[, which(colnames(pipe2) == colnames(pipe1))]
                                                )
                              )
  
  # Calculate MSE for each ROI
  mse_values <- colMeans((pipe1 - pipe2)^2)
  
  return(list(correlation = correlation_values, mse = mse_values))
}


# Function to plot roi data for both pipelines in ./pipeline_comparisons/roiData/
# Scatter plot for a specific ROI (e.g., ROI 1)
plot_roi_corr_between_pipeline <- function(pipeline1_avg, pipeline2_avg, low_corr_df){
  
  for (roi_index in 1:length(colnames(pipeline1_avg))){
    
    # Create a data frame for plotting with Glucose_ID labels
    plot_data <- data.frame(Pipeline1 = pipeline1_avg[, roi_index],
                              Pipeline2 = pipeline2_avg[, roi_index],
                              Glucose_ID = low_corr_df$Glucose_ID  # Add the Glucose_ID for labeling
                            )
    # Write in PNG
    fn <- paste0(script_dir, 
                 "pipeline_comparison/ROIplots/",  
                 paste("ScatterPlot_, colnames(pipeline1_avg[roi_index])") ,".png")
    png(file = fn, res = 200, height = 1000, width = 1000)
        
    # Define limits    
    min_value <- min(c(plot_data$Pipeline1, plot_data$Pipeline2), na.rm = TRUE)
    max_value <- max(c(plot_data$Pipeline1, plot_data$Pipeline2), na.rm = TRUE)

    # Plot with labels
    sc_plot <- ggplot(plot_data, aes(x = Pipeline1, y = Pipeline2)) +
    geom_point(color = "blue", size = 3) +  # Scatter plot points
    geom_text(aes(label = Glucose_ID), vjust = -0.5, hjust = 0.5, size = 3) +  # Add participant Glucose_ID labels
    labs(title = paste("Scatter Plot for", colnames(pipeline1_avg[roi_index])),
           x = "New Pipeline ROI Values",
           y = "CONN Pipeline ROI Values"
         ) +
   theme_minimal() +
   geom_smooth(method = "lm", col = "red") + # Add regression line for trend
   coord_fixed(ratio = 1) +  # Ensure equal scaling of x and y axes
   xlim(min_value, max_value) +  # Set x-axis limits
   ylim(min_value, max_value)    # Set y-axis limits
 
 # Save
 sfn <- paste0("pipeline_comparison/ROIplots/ScatterPlot_", 
                 colnames(pipeline1_avg[roi_index]), 
                 ".png")
 ggsave(sfn, sc_plot, width = 10, height = 8, dpi = 300, bg = "white")
   
 print(sprintf("Plot for ROI %s saved.", roi_index))
  }

}


plot_pca_results <- function(iglu.pca, pc_iglu, outcome, output_dir){
  
  graphics.off()
  
  # Plot explained variance from PCA
  png(file = paste0(script_dir, output_dir, "PCA_iglu_variance_expl.png"))
  
  # Increase the top margin if needed
  print(fviz_eig(iglu.pca))  # Explicitly print the plot
  graphics.off()
  
  # Plot PCA results (PCA space)
  png(file = paste0(script_dir, output_dir, "PCA_iglu_space.png"))
  print(autoplot(iglu.pca))  # Explicitly print the plot
  graphics.off()
  
  # Plot correlation of principal components
  png(file = paste0(script_dir, output_dir, "PCA_iglu_corr.png"))
  corrplot(cor(pc_iglu), 
             method = "color", 
             type = "lower", 
             tl.cex = 0.8, 
             number.cex = 0.5, 
             col = standard_palette
           )
  graphics.off()
  
  
  # Plot histogram of PC1
  file_name <- paste0(script_dir, output_dir, "histogram_PC1.png")
  png(file_name)
  hist(outcome$outcome, breaks = 50, xlab = "PC1", main = "")
  graphics.off()
}

# Plot the correlation between outcome and predictors
plot_correlation_outcome_predictors <- function(output_dir, iglu, predictors){
  
  # Plot iglu/outcome correlations
  png(file=paste0(script_dir, output_dir, "iglu_corr.png"), 
        height = 800, 
        width = 800, 
        res = 200)
  out_corr <- corrplot(cor(iglu), 
                         method="color", 
                         type="lower", 
                         tl.cex=0.3, 
                         number.cex=0.5, 
                         col=standard_palette
                       )
  graphics.off()
  
  
  # Plot correlations of predictors
  png(file=paste0(script_dir, output_dir, "predictors_corr.png"), 
        height = 1000, 
        width = 1000,
        res = 200)
  pred_corr <- corrplot(cor(predictors), 
                          method = "color", 
                          type = "lower", 
                          tl.cex = 0.3,
                          number.cex = 0.5,
                          col = standard_palette
                        )
  
  graphics.off()
  
}


# Plot correlation between fitted coefficients and distribution of the never zero fitted coefficients
plot_coef_consistency <- function(output_dir, nzcoefs_df, s){
  
  # Get nonzero coefficients
  non_zero_columns <- nzcoefs_df[, apply(nzcoefs_df, 2, function(col) all(col != 0))]
  
  # Plot consistency of coefficients across seeds
  fn <- paste0(script_dir, output_dir, "10fold/", "coefficient_consistency.png")
  png(file = fn, height = 1000, width = 1000, res = 200)
  
  par(mar = c(5, 5, 5, 5))
  
  coef_consistency <- corrplot(cor(t(nzcoefs_df)), 
                                 method = "color", 
                                 type = "lower", 
                                 tl.cex = 0.2, 
                                 number.cex = 0.5, 
                                 col = standard_palette, 
                                 add = FALSE
                               )
  
  grid.text(paste0("Correlation of coefficients \n between", s, "seeds"), # Add title
            x = unit(0.95, "npc"), y = unit(0.95, "npc"), 
            just = c("right", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
  
  graphics.off()
  
  
  # Create dataframe to plot non-zero coefficients
  df_long <- pivot_longer(non_zero_columns, 
                            cols = everything(), 
                            names_to = "Coefficient",  
                            values_to = "Value"
                          )
  
  #Determine dynamic x-axis limits
  x_limits <- range(df_long$Value, na.rm = TRUE)
  x_range <- max(abs(x_limits))  # Extend the range equally in both directions
  dynamic_limits <- c(-x_range, x_range)
  
  # Determine the number of coefficients
  num_coefficients <- length(unique(df_long$Coefficient))
  
  # Calculate the height of the output image based on the number of coefficients
  # Adjust height factor as needed to control spacing
  height_factor <- 1.5  # Height per coefficient
  plot_height <- num_coefficients * height_factor
  # Define output directory
  address_dir <- paste0(script_dir, output_dir, "10fold/coefs")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  graphics.off()
  # Create the density plot with facets and dynamic x-axis limits
  p <- ggplot(df_long, aes(x = Value, fill = Coefficient, color = Coefficient)) +
    geom_density(alpha = 0.5, adjust = 1.5) +  # Overlay density plots with transparency
    facet_wrap(~ Coefficient, scales = "free_y", ncol = 1) +  # Stack plots vertically
    scale_x_continuous(name = "Coefficient Value", limits = dynamic_limits) +  # Center around zero with dynamic limits
    labs(title = "Distribution of Non-Zero Coefficients Across Seeds", y = "Density") +
    theme_minimal() +
    theme(
           plot.title = element_text(size = 22, face = "bold", hjust = 0.5),  # Title settings
           axis.title.x = element_text(size = 20),  # X-axis title size
           axis.title.y = element_text(size = 20),  # Y-axis title size
           axis.text.x = element_text(size = 16),  # X-axis title size
           axis.text.y = element_text(size = 16),  # Y-axis title size
           strip.text = element_text(size = 22),    # Facet label size
           legend.position="none",
          ) +
    coord_cartesian(xlim = dynamic_limits)  # Center around zero with dynamic limits
  
  # Define filename and path
  file_path <- file.path(address_dir, "coefficient_distributions_vertical_stacked.png")
  # Save the plot as a PNG file with dynamic height
  ggsave(filename = file_path, plot = p, width = 20, height = plot_height, dpi = 100)
}


# Plot each coefficient distribution (also non-zero)
plot_each_coef_distr <- function(output_dir, nzcoefs_df){
  
  # Loop over each non-zero column and create histograms for coefficients
  for (col_name in colnames(nzcoefs_df)) {
    print(col_name)
    
    # Calculate density and dynamically determine xlim
    density_data <- density(nzcoefs_df[[col_name]], na.rm = TRUE)
    xlim_range <- range(density_data$x)
    
    # Plot density with dynamically adjusted xlim
    p <- ggplot(nzcoefs_df, aes(x = !!sym(col_name))) +
      geom_density(fill = "blue", color = "black") +
      labs(title = paste("Density Plot of", col_name), x = col_name, y = "Density") +
      theme_minimal() +
      xlim(xlim_range) +  # Dynamically set xlim based on the density range
      theme(
              plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Title settings
              axis.title.x = element_text(size = 12),  # X-axis title size
              axis.title.y = element_text(size = 12)   # Y-axis title size
            )
    
    # Define filename and path, replacing special characters in column names
    filep <- paste0(script_dir, output_dir, "10fold/coefs")
    file_path <- file.path(filep, paste0(gsub("[^[:alnum:]_]", "_", col_name), "_histogram.png"))
    
    # Save the histogram as a PNG file
    ggsave(filename = file_path, plot = p, width = 8, height = 6, dpi = 100, bg = "white")
  }
}


# Plot histogram of all coefficients
plot_all_coef_hist <- function(output_dir, nzcoefs_df, nseeds){
  
  # Select all columns
  non_zero_columns <- nzcoefs_df
  
  # Obtain columns that are never zero
  never_zero_columns <- nzcoefs_df[, apply(nzcoefs_df, 2, function(col) all(col != 0))]
  
  # plot non zero coefficients for one seed
  file_name <- paste0(script_dir, output_dir, "/10fold/coefs/nonzero_model_coefs_",nseeds,"seeds.png")
  
  # Calculate column means
  mean_values <- colMeans(non_zero_columns)
  
  # Create a data frame for ggplot
  plot_data <- data.frame(Predictor = names(mean_values),
                            CoefficientValue = mean_values
                          )
  
  # Add a column to indicate if the predictor is "never zero" (TRUE/FALSE)
  plot_data$is_never_zero <- plot_data$Predictor %in% colnames(never_zero_columns)
  plot_data <- plot_data[order(plot_data$Predictor), ]
  print(plot_data)
  
  # Create the bar plot
  pl <- ggplot(plot_data, aes(x = Predictor, y = CoefficientValue, fill = is_never_zero)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste0("Mean Coefficients from Elastic Net Reg \n#Seeds=", nseeds),
         x = "Predictor",
         y = "Mean Coefficient Value"
         ) +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "darkgray"), 
                      labels = c("TRUE" = "Never Zero", "FALSE" = "Sometimes Zero")) +  # Color for never zero (red) and sometimes zero (blue)
    theme_minimal() +
    theme(
      # Change the font size of the y-axis labels
      axis.text.y = element_text(size = 8, color = ifelse(plot_data$is_never_zero, "darkred", "darkgray")),# Adjust the size as needed
      axis.text.x = element_text(size = 8)
         ) 
                                
  
  # Save
  ggsave(file_name, pl, width = 6, height = 8, dpi = 200, bg = "white")
  graphics.off()

}

# Plot prediction of final model
plot_prediction_final_model <- function(output_dir, predictions, outcome, model_formula, cv_model, s){
  
  file_name <- paste0(script_dir, output_dir, "10fold/predictions_seed", s, "_",str_var,".jpg")
  
  jpeg(file_name, width = 6, height = 6, units = "in", res = 600)
  
  par(mar = c(5, 4, 12, 2))  # Increase the top margin to 7 lines
  
  x <- c(predictions[,1])
  y <- c(outcome[,1])
  
  p <- plot(x, y, 
            main=split_title(paste0(paste0(str_var, "~", model_formula)[3],
                                      "\n\n", 
                                      "with alpha =0.05, lambda=",
                                      sprintf(cv_model$bestTune$lambda, fmt = "%#.2f"), 
                                      ", R_squared=", sprintf(performance[2], fmt = "%#.2f")), 
                              max_length=100), 
            xlab=paste0("Predicted ",str_var," from Iglu"),
            ylab=paste0(str_var, "from Iglu"),
            cex.main = 0.6,       # Make main title smaller
            xlim = range(x, y),  # Align x and y axes
            ylim = range(x, y))  # Align x and y axes)
  graphics.off()
  
}


# Plot permuted performance
plot_permuted_performance <- function(output_dir, permutation_test_results, n_seeds){
  
  fn <- paste0(output_dir, "10fold/histogram_permuted_performance_",n_seeds,"permuations.png")
  file_name <- file.path(getwd(), fn)
  
  png(file_name)
  
  hist(permutation_test_results, 
         breaks=50, 
         xlab = "Permuted Performance Avg. RMSE", 
         main=sprintf("mean = %#.2f", mean(permutation_test_results))
       )
  
  graphics.off()
}


# Plot performance across seeds
plot_performance_across_seeds <- function(output_dir, perf_list, n_seeds){
  
  # Histogram
  file_name <- file.path(getwd(), paste0(output_dir, 
                                           "/10fold/histogram_performance_", 
                                           n_seeds, 
                                           "seeds.png")
                         )
  png(file_name)
  
  hist(perf_list, 
         breaks=50, 
         xlab = sprintf("Performance Avg. RMSE"), 
         main=sprintf("mean = %#.2f", mean(perf_list))
       )
  
  graphics.off()
}


# Plot best lambda
plot_best_lambda <- function(output_dir, bestLambdas){
  
  # Plot histogram of bestLambdas
  file_name <- paste0(script_dir, output_dir, "10fold/bestLambdas.png")
  png(file_name)
  
  hist(bestLambdas, 
         breaks = 50, 
         xlab = "Best Lambda", 
         main = "Distribution of fitted Lambdas"
       )
  
  graphics.off()
}

