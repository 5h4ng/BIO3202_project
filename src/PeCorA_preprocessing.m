function [t] = PeCorA_preprocessing(t, area_column_name, threshold_to_filter, control_name)
%PECORA_PREPROCESSING Preprocess data for PeCorA analysis
%   This function preprocesses the input data for PeCorA analysis by:
%   1. Creating peptide-charge combinations
%   2. Filtering data based on thresholds
%   3. Log2 transformation and scaling
%   4. Centering peptides relative to control group
%
%   Inputs:
%   t - Table containing the data in PeCorA format
%   area_column_name - Column number containing peak areas
%   threshold_to_filter - Threshold value for peak areas
%   control_name - Name of the control group
%
%   Outputs:
%   t - Processed table ready for PeCorA analysis

% Create peptide-charge combination column
t.modpep_z = strcat(t.Peptide_Modified_Sequence, '_all');

% Create condition-replicate combination column
t.condition_rep = strcat(t.Condition, '_', string(t.BioReplicate));

% Get total number of replicates
n_reps_total = length(unique(t.condition_rep));

% Filter data
% Remove NA values
t(isnan(t{:, area_column_name}), :) = [];

% Filter by threshold
t(t{:, area_column_name} <= threshold_to_filter, :) = [];

% Keep only peptides with all measurements
peptide_counts = countcats(categorical(t.modpep_z));
valid_peptides = peptide_counts == n_reps_total;
t = t(ismember(t.modpep_z, categories(categorical(t.modpep_z))(valid_peptides)), :);

% Log2 transform and scale
t.ms1log2 = log2(t{:, area_column_name});

% Scale the data
% Note: MATLAB doesn't have a direct equivalent to R's scale_by
% We'll implement a custom scaling function
t.ms1scaled = scale_data(t.ms1log2, t.Condition, t.BioReplicate);

% Center peptides relative to control
fprintf('Scaling peptides to control == 0\n');
t.ms1adj = zeros(height(t), 1);

% Get control data
t_cntrl = t(strcmp(t.Condition, control_name), :);
ms1scaled_cntrl = t_cntrl.ms1scaled;
ms1scaled_full = t.ms1scaled;

% Process each unique peptide
unique_peptides = unique(t.modpep_z);
h = waitbar(0, 'Processing peptides...');
for i = 1:length(unique_peptides)
    x = unique_peptides{i};
    peptide_idx = strcmp(t.modpep_z, x);
    control_idx = strcmp(t_cntrl.modpep_z, x);
    
    t.ms1adj(peptide_idx) = ms1scaled_full(peptide_idx) - mean(ms1scaled_cntrl(control_idx));
    
    waitbar(i/length(unique_peptides), h);
end
close(h);

fprintf('Peptide scaling completed\n');
end

function scaled_data = scale_data(data, condition, replicate)
% Custom scaling function to replace R's scale_by
% Centers the data by condition and replicate
    scaled_data = zeros(size(data));
    unique_conditions = unique(condition);
    unique_replicates = unique(replicate);
    
    for i = 1:length(unique_conditions)
        for j = 1:length(unique_replicates)
            idx = strcmp(condition, unique_conditions{i}) & ...
                  replicate == unique_replicates(j);
            if any(idx)
                scaled_data(idx) = data(idx) - mean(data(idx));
            end
        end
    end
end 