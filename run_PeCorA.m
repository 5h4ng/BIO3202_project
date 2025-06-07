addpath('src');

% 1. Import data
opts = detectImportOptions('data/PeCorA_noZ.csv');
opts.VariableNamingRule = 'preserve';
t = readtable('data/PeCorA_noZ.csv', opts);

% 2. Data preprocessing
scaled_peptides = PeCorA_preprocessing(t, 'Normalized Area', 100, 'cntrl');

% 3. Run PeCorA analysis
disagree_peptides = PeCorA(scaled_peptides);

% 4. Save results
if ~exist('results', 'dir')
    mkdir('results');
end

writetable(disagree_peptides, 'results/PeCorA_results.csv');
writetable(scaled_peptides, 'results/PeCorA_scaled_peptides.csv');

% 5. Display summary 
fprintf('\nPeCorA Analysis Complete\n');
fprintf('========================\n');
if height(disagree_peptides) > 0
    sig_05 = sum(disagree_peptides.adj_pval < 0.05);
    sig_01 = sum(disagree_peptides.adj_pval < 0.01);
    sig_proteins = length(unique(disagree_peptides.protein(disagree_peptides.adj_pval < 0.05)));
    
    fprintf('Total peptides tested: %d\n', height(disagree_peptides));
    fprintf('Significant peptides (adj.p < 0.05): %d\n', sig_05);
    fprintf('Significant peptides (adj.p < 0.01): %d\n', sig_01);
    fprintf('Proteins with significant peptides: %d\n', sig_proteins);
    
    if sig_05 > 0
        fprintf('\nTop 5 most significant peptides:\n');
        top5 = disagree_peptides(1:min(5, height(disagree_peptides)), :);
        for i = 1:height(top5)
            fprintf('  %s (adj.p = %.2e)\n', top5.peptide{i}, top5.adj_pval(i));
        end
    end
else
    fprintf('No significant peptides found.\n');
end

fprintf('\nResults saved to results/ directory\n');
fprintf('Run plot_PeCorA_results to generate visualizations\n');
    
