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
%   area_column_name - Column name containing peak areas
%   threshold_to_filter - Threshold value for peak areas
%   control_name - Name of the control group
%
%   Outputs:
%   t - Processed table ready for PeCorA analysis

% Input validation
if ~istable(t)
    error('Input must be a table');
end

if ~ischar(area_column_name) && ~isstring(area_column_name)
    error('Area column name must be a string');
end

if ~isnumeric(threshold_to_filter) || threshold_to_filter <= 0
    error('Threshold must be a positive number');
end

if ~ischar(control_name) && ~isstring(control_name)
    error('Control name must be a string');
end

fprintf('Initial table size: %d rows\n', height(t));

% Verify column names exist
required_columns = {'Peptide Modified Sequence', 'Condition', 'BioReplicate', area_column_name};
missing_columns = setdiff(required_columns, t.Properties.VariableNames);
if ~isempty(missing_columns)
    error('Missing required columns: %s', strjoin(missing_columns, ', '));
end

% Verify data types
area_data = table2array(t(:, area_column_name));
if ~isnumeric(area_data)
    error('Area column must contain numeric values');
end

% Print summary of area data
fprintf('Area data summary:\n');
fprintf('Min: %f\n', min(area_data));
fprintf('Max: %f\n', max(area_data));
fprintf('Mean: %f\n', mean(area_data, 'omitnan'));
fprintf('Number of NaN: %d\n', sum(isnan(area_data)));
fprintf('Number of Inf: %d\n', sum(isinf(area_data)));

% Create peptide-charge combination column
t.modpep_z = strcat(t.('Peptide Modified Sequence'), '_all');

% Create condition-replicate combination column
t.condition_rep = strcat(t.Condition, '_', string(t.BioReplicate));

% Get total number of replicates
n_reps_total = length(unique(t.condition_rep));
fprintf('Total number of replicates: %d\n', n_reps_total);

% Filter data
% Remove NA values
valid_rows = ~isnan(area_data) & ~isinf(area_data);
fprintf('Number of valid rows (non-NA, non-Inf): %d\n', sum(valid_rows));
t = t(valid_rows, :);
fprintf('After removing NA and Inf values: %d rows\n', height(t));

% Filter by threshold
valid_rows = area_data(valid_rows) > threshold_to_filter;
fprintf('Number of rows above threshold: %d\n', sum(valid_rows));
t = t(valid_rows, :);
fprintf('After threshold filtering: %d rows\n', height(t));

if height(t) == 0
    error('No data remains after filtering. Please check your threshold value and data.');
end

% Keep only peptides with all measurements
fprintf('Filtering peptides with complete measurements...\n');

% Get unique peptides and their counts
[unique_peptides, ~, peptide_idx] = unique(t.modpep_z);
peptide_counts = accumarray(peptide_idx, 1);
fprintf('Number of unique peptides before filtering: %d\n', length(unique_peptides));

% Find peptides with complete measurements
valid_peptides = peptide_counts == n_reps_total;
fprintf('Number of peptides with complete measurements: %d\n', sum(valid_peptides));

if sum(valid_peptides) == 0
    error('No peptides have complete measurements across all replicates');
end

% Filter table
valid_peptide_names = unique_peptides(valid_peptides);
t = t(ismember(t.modpep_z, valid_peptide_names), :);
fprintf('After keeping complete peptides: %d rows\n', height(t));

% Log2 transform and scale
area_values = table2array(t(:, area_column_name));
if any(area_values <= 0)
    error('Cannot perform log2 transformation: data contains zero or negative values');
end

t.ms1log2 = log2(area_values);

% Scale the data
t.ms1scaled = scale_data(t.ms1log2, t.Condition, t.BioReplicate);

% Center peptides relative to control
fprintf('Scaling peptides to control == 0\n');
t.ms1adj = zeros(height(t), 1);

% Get control data
control_idx = strcmp(t.Condition, control_name);
if ~any(control_idx)
    error('Control group "%s" not found in the data', control_name);
end

t_cntrl = t(control_idx, :);
ms1scaled_cntrl = t_cntrl.ms1scaled;
ms1scaled_full = t.ms1scaled;

% Process each unique peptide
unique_peptides = unique(t.modpep_z);
h = waitbar(0, 'Processing peptides...');
for i = 1:length(unique_peptides)
    x = unique_peptides{i};
    peptide_idx = strcmp(t.modpep_z, x);
    control_idx = strcmp(t_cntrl.modpep_z, x);
    
    if any(control_idx)
        t.ms1adj(peptide_idx) = ms1scaled_full(peptide_idx) - mean(ms1scaled_cntrl(control_idx));
    end
    
    waitbar(i/length(unique_peptides), h);
end
close(h);

% Add summary statistics
fprintf('\nFinal data summary:\n');
fprintf('Number of unique proteins: %d\n', length(unique(t.Protein)));
fprintf('Number of unique peptides: %d\n', length(unique(t.modpep_z)));
fprintf('Number of conditions: %d\n', length(unique(t.Condition)));
fprintf('Number of replicates: %d\n', n_reps_total);

fprintf('Peptide scaling completed\n');
end

function scaled_data = scale_data(data, condition, replicate)
% Custom scaling function to replace R's scale_by
% Centers the data by condition and replicate combination (like R's scale_by)
% This scales within each combination of condition and replicate
    scaled_data = zeros(size(data));
    
    % Create condition-replicate combinations
    unique_conditions = unique(condition);
    unique_replicates = unique(replicate);
    
    % Process each unique condition-replicate combination
    for i = 1:length(unique_conditions)
        for j = 1:length(unique_replicates)
            % Find indices for this specific condition-replicate combination
            idx = strcmp(condition, unique_conditions{i}) & ...
                  replicate == unique_replicates(j);
            
            if any(idx)
                group_data = data(idx);
                % Scale within this group: subtract mean and divide by std
                if length(group_data) > 1 && std(group_data) > 0
                    scaled_data(idx) = (group_data - mean(group_data)) / std(group_data);
                else
                    % If only one value or no variation, just center
                    scaled_data(idx) = group_data - mean(group_data);
                end
            end
        end
    end
end 