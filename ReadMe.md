<br> Multiple Linear Regression </br>

Fits a simple Multiple Linear Regression model without interactions but to generated data (iglu.xlsx) using generated resting state functional Magnetic Resonance Imaging (rs-FMRI) data.
Test data was generate using ./linear_model_R/generate_test_data.m to control the fitting procedure of the model.
The predictors are set to main effects of each column, whereas we created effects for the columns Mammilary body r, Mammilary body l, Subthalamic nucleus r, Subthalamic nucleus l, NTS and LC which can be observed wiht significant non-zero coefficients throughout all 100 model fits.
Plots of the fitting procedure and data can be found in ./linear_model_R/plots/.

./linear_model_R/roiData contains generated rs-fMRI data.
./CONN_labels_Atlas.xlsx contains column names of Regions of Interest (ROI) from the used/assumed fMRI Atlas. 

