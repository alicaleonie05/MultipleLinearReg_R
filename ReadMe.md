**Multiple Linear Regression**

Fits a simple Multiple Linear Regression model without interaction effects to predict the first principal component (PC1) of generated data (iglu.xlsx) using generated resting state functional Magnetic Resonance Imaging (rs-FMRI) data.


Test data was generate using ./linear_model_R/generate_test_data.m to control the fitting procedure of the model.

./linear_model_R/roiData contains generated rs-fMRI data.

./linear_model_R/CONN_labels_Atlas.xlsx contains column names of Regions of Interest (ROI) from the used/assumed fMRI Atlas.


The predictors are set to main effects of each column, whereas I created effects for the columns Mammilary body r, Mammilary body l, Subthalamic nucleus r, Subthalamic nucleus l, NTS and LC which can be observed wiht significant non-zero coefficients throughout all 100 model fits.
Plots of the fitting procedure and data can be found in ./linear_model_R/plots/.


