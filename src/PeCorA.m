function [alldf_ordered] = PeCorA(t)
%PECORA Core PeCorA analysis function
%   This function performs the main PeCorA analysis to identify peptides
%   that show discordant behavior compared to other peptides in the same protein
%
%   Inputs:
%   t - Preprocessed table from PeCorA_preprocessing
%
%   Outputs:
%   alldf_ordered - Table containing results ordered by adjusted p-value

fprintf('Checking which proteins still have at least 2 peptides\n');

% Get unique proteins and filter those with at least 2 peptides
pgs = unique(t.Protein);
pgs_morethan2 = {};
for i = 1:length(pgs)
    if length(unique(t{strcmp(t.Protein, pgs{i}), 'modpep_z'})) > 1
        pgs_morethan2{end+1} = pgs{i};
    end
end

% Initialize storage for p-values
allp = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Process each protein
fprintf('Computing the interaction p-values\n');
h = waitbar(0, 'Processing proteins...');
for i = 1:length(pgs_morethan2)
    x = pgs_morethan2{i};
    tmpdf = t(strcmp(t.Protein, x), :);
    tmpdf.allothers = repmat({'allothers'}, height(tmpdf), 1);
    
    % Get unique peptides for this protein
    unique_peps = unique(tmpdf.modpep_z);
    pvalues = zeros(length(unique_peps), 1);
    
    % Process each peptide
    for j = 1:length(unique_peps)
        y = unique_peps{j};
        subtmpdf = tmpdf;
        subtmpdf.allothers(strcmp(tmpdf.modpep_z, y)) = {y};
        
        % Fit linear model
        X = [ones(height(subtmpdf), 1), ...
             categorical(subtmpdf.Condition), ...
             categorical(subtmpdf.allothers), ...
             categorical(subtmpdf.Condition).*categorical(subtmpdf.allothers)];
        
        [~, ~, stats] = anovan(subtmpdf.ms1adj, ...
            {categorical(subtmpdf.Condition), categorical(subtmpdf.allothers)}, ...
            'model', 'interaction', ...
            'display', 'off');
        
        pvalues(j) = stats.p(3); % Get p-value for interaction term
    end
    
    allp(x) = pvalues;
    waitbar(i/length(pgs_morethan2), h);
end
close(h);

fprintf('PeCorA finished\n');

% Create results table
fprintf('Creating results table\n');
alldf = table();
for i = 1:length(pgs_morethan2)
    x = pgs_morethan2{i};
    tmpdf = t(strcmp(t.Protein, x), :);
    tmp_peps = unique(tmpdf.modpep_z);
    
    if ~isempty(tmp_peps)
        tmp_pval = allp(x);
        tmpout = table(repmat({x}, length(tmp_pval), 1), ...
                      tmp_peps, ...
                      tmp_pval, ...
                      'VariableNames', {'protein', 'peptide', 'pvalue'});
        alldf = [alldf; tmpout];
    end
end

% Adjust p-values
fprintf('Correcting p-values\n');
alldf.adj_pval = mafdr(alldf.pvalue, 'BHFDR', true);

% Sort by adjusted p-value
alldf_ordered = sortrows(alldf, 'adj_pval');

% Print summary
fprintf('Number of uncorrelated peptides = %d\n', sum(alldf_ordered.adj_pval <= 0.01));
fprintf('Number of proteins with uncorrelated peptides = %d\n', ...
    length(unique(alldf_ordered.protein(alldf_ordered.adj_pval <= 0.01))));

end 