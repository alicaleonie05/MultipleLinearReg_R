% generate sample of rs-FC data

% Initialize a 154 x 154 x 100 matrix with Gaussian noise (z-transformed values)
% Initialize dimensions
num_rows = 154;
num_cols = 154;
num_participants = 100;

% Create the base 3D matrix with random noise centered around 0
matrix3D = 0.4 * randn(num_rows, num_cols, num_participants);

% Define participant-wise increment/decrement factors for correlation strengths
pos_increment = linspace(0.5, 2.5, num_participants);  % positive correlations increase over participants
neg_increment = linspace(-2.5, -0.5, num_participants); % negative correlations decrease over participants

% Apply correlations to hypothalamus connections
for p = 1:num_participants
    % Strong positive correlation with columns 149 and 150
    matrix3D(147, 149, p) = pos_increment(p);
    matrix3D(148, 149, p) = pos_increment(p);
    matrix3D(147, 150, p) = pos_increment(p);
    matrix3D(148, 150, p) = pos_increment(p);

    % Strong negative correlation with columns 151 and 152
    matrix3D(147, 151, p) = neg_increment(p);
    matrix3D(148, 151, p) = neg_increment(p);
    matrix3D(147, 152, p) = neg_increment(p);
    matrix3D(148, 152, p) = neg_increment(p);

    % Weak positive correlation with column 153
    matrix3D(147, 153, p) = 0.5 * pos_increment(p);
    matrix3D(148, 153, p) = 0.5 * pos_increment(p);

    % Weaker positive correlation with column 154
    matrix3D(147, 154, p) = 0.25 * pos_increment(p);
    matrix3D(148, 154, p) = 0.25 * pos_increment(p);
end

% Symmetrize the matrix for each participant by mirroring the upper triangle
for p = 1:num_participants
    matrix3D(:,:,p) = triu(matrix3D(:,:,p),1) + triu(matrix3D(:,:,p),1)';
end

% Add more noise
noise = 0.25 * randn(matrix_dim);
matrix3D = matrix3D + noise;

% Plot one of the slices to verify the correlations visually
slice_number = 1;  % Choose a slice to display
imagesc(matrix3D(:, :, slice_number));
colorbar;
title(['Functional Connectivity Matrix Slice ', num2str(slice_number)]);
xlabel('Regions');
ylabel('Regions');
% Save matrix3D to a .mat file
save('./roiData/roitoroi/AvgCond_3DparticipantStack.mat', 'matrix3D');

% Initialize the iglu matrix with 5 columns and 153 rows
iglu = zeros(100, 5);

% Assigning values to iglu columns
iglu(:, 1) = linspace(0, 100, 100);                % Increasing values from 0 to 100 for column 1
iglu(:, 2) = linspace(60, 80, 100);                % Increasing values from 60 to 80 for column 2
iglu(:, 3) = linspace(4, 0, 100);                  % Decreasing values from 4 to 0 for column 3
iglu(:, 4:5) = -60 + 120 * rand(100, 2);           % Random noise in range [-60, 60] for columns 4 and 5

% Convert iglu to a table and rename the columns
iglu_df = array2table(iglu, 'VariableNames', {'var1', 'var2', 'var3', 'var4', 'var5'});

% Save iglu_df to an Excel file
writetable(iglu_df, 'iglu.xlsx');