function disagree_peptides = PeCorA(scaled_peptides)
%PECORA Peptide Correlation Analysis
%   Performs statistical analysis to identify peptides that show changes
%   inconsistent with their parent protein's overall changes
%   
%   Based on the original R implementation logic:
%   For each protein, tests if individual peptides have different 
%   condition-response patterns compared to all other peptides in the protein

% Input validation
assert(istable(scaled_peptides), 'Input must be a table');
required_cols = {'Protein', 'modpep_z', 'ms1adj', 'Condition'};
missing_cols = setdiff(required_cols, scaled_peptides.Properties.VariableNames);
assert(isempty(missing_cols), 'Missing required columns: %s', strjoin(missing_cols, ', '));

fprintf('Scaled peptides are ready for PeCorA\n');

% Get all unique protein groups
fprintf('Checking which proteins still have at least 2 peptides\n');
unique_proteins = unique(scaled_peptides.Protein);
pgs_morethan2 = {};

for i = 1:length(unique_proteins)
    protein = unique_proteins{i};
    protein_data = scaled_peptides(strcmp(scaled_peptides.Protein, protein), :);
    unique_peptides = unique(protein_data.modpep_z);
    
    if length(unique_peptides) > 1
        pgs_morethan2{end+1} = protein;
    end
end

fprintf('Found %d proteins with at least 2 peptides\n', length(pgs_morethan2));

% Initialize storage for p-values
allp = containers.Map('KeyType', 'char', 'ValueType', 'any');

fprintf('Computing the interaction p-values\n');

% Loop through proteins with >1 peptide
for i = 1:length(pgs_morethan2)
    protein = pgs_morethan2{i};
    
    % Get subset dataframe for this protein
    tmpdf = scaled_peptides(strcmp(scaled_peptides.Protein, protein), :);
    
    % Get unique peptides for this protein
    unique_peptides = unique(tmpdf.modpep_z);
    pvalues = zeros(length(unique_peptides), 1);
    
    % Loop through each peptide in this protein
    for j = 1:length(unique_peptides)
        target_peptide = unique_peptides{j};
        
        % Create allothers column: 'allothers' for all, then change target peptide
        allothers = repmat({'allothers'}, height(tmpdf), 1);
        target_idx = strcmp(tmpdf.modpep_z, target_peptide);
        allothers(target_idx) = {target_peptide};
        
        % Convert to categorical for ANOVA
        condition_cat = categorical(tmpdf.Condition);
        allothers_cat = categorical(allothers);
        
        % Check if we have enough variation for ANOVA
        if length(unique(condition_cat)) < 2 || length(unique(allothers_cat)) < 2
            pvalues(j) = NaN;
            continue;
        end
        
        try
            % Perform ANOVA: ms1adj ~ Condition * allothers
            [p, tbl, stats] = anovan(tmpdf.ms1adj, ...
                {condition_cat, allothers_cat}, ...
                'model', 'interaction', 'display', 'off');
            
            % Get interaction p-value (should be 3rd element)
            if length(p) >= 3 && ~isnan(p(3))
                pvalues(j) = p(3);
            else
                pvalues(j) = NaN;
            end
        catch ME
            warning('Error in ANOVA for protein %s, peptide %s: %s', ...
                protein, target_peptide, ME.message);
            pvalues(j) = NaN;
        end
    end
    
    % Store p-values for this protein
    allp(protein) = pvalues;
end

% Post testing analysis
total_proteins = length(allp.keys);
all_pvals = [];
protein_keys = allp.keys;
for i = 1:length(protein_keys)
    all_pvals = [all_pvals; allp(protein_keys{i})];
end

fprintf('Number of proteins tested = %d\n', total_proteins);
fprintf('Number of peptides tested = %d\n', length(all_pvals));

% Make table of results
fprintf('Started making data table\n');
protein_list = {};
peptide_list = {};
pval_list = [];

for i = 1:length(protein_keys)
    protein = protein_keys{i};
    tmpdf = scaled_peptides(strcmp(scaled_peptides.Protein, protein), :);
    tmp_peps = unique(tmpdf.modpep_z);
    tmp_pval = allp(protein);
    
    if length(tmp_peps) > 0 && length(tmp_pval) == length(tmp_peps)
        protein_list = [protein_list; repmat({protein}, length(tmp_peps), 1)];
        peptide_list = [peptide_list; tmp_peps];
        pval_list = [pval_list; tmp_pval];
    end
end

% Create results table
disagree_peptides = table(protein_list, peptide_list, pval_list, ...
    'VariableNames', {'protein', 'peptide', 'pvalue'});

% Remove NaN p-values
valid_idx = ~isnan(disagree_peptides.pvalue);
disagree_peptides = disagree_peptides(valid_idx, :);

if height(disagree_peptides) > 0
    % Correct p-values using Benjamini-Hochberg
    fprintf('Correcting p-values\n');
    disagree_peptides.adj_pval = mafdr(disagree_peptides.pvalue, 'BHFDR', true);
    
    % Sort by adjusted p-value
    disagree_peptides = sortrows(disagree_peptides, 'adj_pval');
    
    % Calculate fold changes for each peptide
    log2FC_list = zeros(height(disagree_peptides), 1);
    CI_lower_list = zeros(height(disagree_peptides), 1);
    CI_upper_list = zeros(height(disagree_peptides), 1);
    
    for i = 1:height(disagree_peptides)
        protein = disagree_peptides.protein{i};
        peptide = disagree_peptides.peptide{i};
        
        % Get peptide data
        peptide_data = scaled_peptides(strcmp(scaled_peptides.Protein, protein) & ...
                                      strcmp(scaled_peptides.modpep_z, peptide), :);
        
        if height(peptide_data) > 0
            % Calculate fold change (treatment vs control)
            conditions = unique(peptide_data.Condition);
            if length(conditions) >= 2
                control_idx = strcmp(peptide_data.Condition, 'cntrl');
                if any(control_idx)
                    control_mean = mean(peptide_data.ms1adj(control_idx));
                    treatment_data = peptide_data.ms1adj(~control_idx);
                    if ~isempty(treatment_data)
                        treatment_mean = mean(treatment_data);
                        
                        % Calculate log2 fold change
                        if control_mean ~= 0 && treatment_mean ~= 0
                            log2FC_list(i) = treatment_mean - control_mean; % Already log2 scaled
                        else
                            log2FC_list(i) = 0;
                        end
                        
                        % Calculate confidence intervals
                        SE = std(peptide_data.ms1adj) / sqrt(height(peptide_data));
                        if ~isnan(SE) && SE > 0
                            CI_lower_list(i) = log2FC_list(i) - 1.96 * SE;
                            CI_upper_list(i) = log2FC_list(i) + 1.96 * SE;
                        end
                    end
                end
            end
        end
    end
    
    % Add fold change and CI columns
    disagree_peptides.log2FC = log2FC_list;
    disagree_peptides.CI_lower = CI_lower_list;
    disagree_peptides.CI_upper = CI_upper_list;
    
    % Rename pvalue column to match expected format
    disagree_peptides.pval = disagree_peptides.pvalue;
    disagree_peptides.pvalue = [];
    
    % Print summary results
    significant_peptides = sum(disagree_peptides.adj_pval <= 0.01);
    significant_proteins = length(unique(disagree_peptides.protein(disagree_peptides.adj_pval <= 0.01)));
    
    fprintf('Number of uncorrelated peptides = %d\n', significant_peptides);
    fprintf('Number of proteins with uncorrelated peptides = %d\n', significant_proteins);
else
    fprintf('No valid results found\n');
end

end 