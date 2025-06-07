% generate_summary_statistics.m
% Script to generate comprehensive summary statistics from PeCorA results
% Outputs to CLI and text file only, no figures

fprintf('=== PeCorA Summary Statistics Generation ===\n');

% Add src directory to MATLAB path
addpath('src');

% Check if results directory exists
if ~exist('results', 'dir')
    error('Results directory not found. Please run PeCorA analysis first.');
end

% Load PeCorA results
results_file = 'results/PeCorA_results.csv';
if ~exist(results_file, 'file')
    error('PeCorA results file not found: %s', results_file);
end

fprintf('Loading PeCorA analysis results...\n');
disagree_peptides = readtable(results_file);

%% ============================================================================
%% 1. PROTEIN-LEVEL STATISTICS
%% ============================================================================

fprintf('\n=== Protein-Level Analysis ===\n');

% Calculate protein-level statistics
unique_proteins = unique(disagree_peptides.protein);
total_proteins = length(unique_proteins);

% Proteins with significant peptides at different thresholds
sig_001_proteins = unique(disagree_peptides.protein(disagree_peptides.adj_pval < 0.01));
sig_005_proteins = unique(disagree_peptides.protein(disagree_peptides.adj_pval < 0.05));
sig_01_proteins = unique(disagree_peptides.protein(disagree_peptides.adj_pval < 0.1));

% Count peptides per protein
protein_peptide_counts = table();
protein_peptide_counts.protein = unique_proteins;
protein_peptide_counts.total_peptides = zeros(total_proteins, 1);
protein_peptide_counts.sig_001_peptides = zeros(total_proteins, 1);
protein_peptide_counts.sig_005_peptides = zeros(total_proteins, 1);

for i = 1:total_proteins
    protein = unique_proteins{i};
    protein_data = disagree_peptides(strcmp(disagree_peptides.protein, protein), :);
    
    protein_peptide_counts.total_peptides(i) = height(protein_data);
    protein_peptide_counts.sig_001_peptides(i) = sum(protein_data.adj_pval < 0.01);
    protein_peptide_counts.sig_005_peptides(i) = sum(protein_data.adj_pval < 0.05);
end

% Print summary statistics
fprintf('Total proteins analyzed: %d\n', total_proteins);
fprintf('Proteins with significant peptides (p < 0.01): %d (%.1f%%)\n', ...
    length(sig_001_proteins), 100*length(sig_001_proteins)/total_proteins);
fprintf('Proteins with significant peptides (p < 0.05): %d (%.1f%%)\n', ...
    length(sig_005_proteins), 100*length(sig_005_proteins)/total_proteins);
fprintf('Proteins with significant peptides (p < 0.1): %d (%.1f%%)\n', ...
    length(sig_01_proteins), 100*length(sig_01_proteins)/total_proteins);

%% ============================================================================
%% 2. PEPTIDE-LEVEL STATISTICS
%% ============================================================================

fprintf('\n=== Peptide-Level Analysis ===\n');

total_peptides = height(disagree_peptides);
sig_001_peptides = sum(disagree_peptides.adj_pval < 0.01);
sig_005_peptides = sum(disagree_peptides.adj_pval < 0.05);
sig_01_peptides = sum(disagree_peptides.adj_pval < 0.1);

fprintf('Total peptides analyzed: %d\n', total_peptides);
fprintf('Significant peptides (FDR < 0.01): %d (%.1f%%)\n', ...
    sig_001_peptides, 100*sig_001_peptides/total_peptides);
fprintf('Significant peptides (FDR < 0.05): %d (%.1f%%)\n', ...
    sig_005_peptides, 100*sig_005_peptides/total_peptides);
fprintf('Significant peptides (FDR < 0.1): %d (%.1f%%)\n', ...
    sig_01_peptides, 100*sig_01_peptides/total_peptides);

%% ============================================================================
%% 3. EFFECT SIZE ANALYSIS
%% ============================================================================

fprintf('\n=== Effect Size Analysis ===\n');

% Effect sizes for significant peptides
sig_peptides = disagree_peptides(disagree_peptides.adj_pval < 0.05, :);
if height(sig_peptides) > 0
    median_abs_fc = median(abs(sig_peptides.log2FC));
    mean_abs_fc = mean(abs(sig_peptides.log2FC));
    min_abs_fc = min(abs(disagree_peptides.log2FC));
    max_abs_fc = max(abs(disagree_peptides.log2FC));
    
    fprintf('Median |log2FC| of significant peptides (FDR < 0.05): %.3f\n', median_abs_fc);
    fprintf('Mean |log2FC| of significant peptides (FDR < 0.05): %.3f\n', mean_abs_fc);
    fprintf('Range of |log2FC| (all peptides): %.3f - %.3f\n', min_abs_fc, max_abs_fc);
    
    % Effect size categories
    abs_fc = abs(sig_peptides.log2FC);
    small_effect = sum(abs_fc >= 0.5 & abs_fc < 1);
    medium_effect = sum(abs_fc >= 1 & abs_fc < 2);
    large_effect = sum(abs_fc >= 2);
    
    fprintf('Effect size distribution (FDR < 0.05):\n');
    fprintf('  Small effects (0.5-1): %d peptides\n', small_effect);
    fprintf('  Medium effects (1-2): %d peptides\n', medium_effect);
    fprintf('  Large effects (≥2): %d peptides\n', large_effect);
end

%% ============================================================================
%% 4. PEPTIDES PER PROTEIN ANALYSIS
%% ============================================================================

fprintf('\n=== Peptides per Protein Analysis ===\n');

mean_total_peptides = mean(protein_peptide_counts.total_peptides);
median_total_peptides = median(protein_peptide_counts.total_peptides);

sig_protein_counts = protein_peptide_counts.sig_005_peptides(protein_peptide_counts.sig_005_peptides > 0);
if ~isempty(sig_protein_counts)
    mean_sig_peptides = mean(sig_protein_counts);
    median_sig_peptides = median(sig_protein_counts);
    max_sig_peptides = max(sig_protein_counts);
    
    fprintf('Mean total peptides per protein: %.1f\n', mean_total_peptides);
    fprintf('Median total peptides per protein: %.0f\n', median_total_peptides);
    fprintf('Mean significant peptides per protein (affected proteins only): %.1f\n', mean_sig_peptides);
    fprintf('Median significant peptides per protein (affected proteins only): %.0f\n', median_sig_peptides);
    fprintf('Maximum significant peptides in a single protein: %d\n', max_sig_peptides);
end

%% ============================================================================
%% 5. GENERATE COMPREHENSIVE SUMMARY TABLE
%% ============================================================================

fprintf('\n=== Generating Summary Statistics Table ===\n');

% Create comprehensive summary table
summary_stats = create_summary_statistics_table(disagree_peptides, protein_peptide_counts);
writetable(summary_stats, 'results/Summary_Statistics.txt', 'Delimiter', '\t');
fprintf('✓ Summary statistics table saved to results/Summary_Statistics.txt\n');

% Also save as CSV
writetable(summary_stats, 'results/Summary_Statistics.csv');
fprintf('✓ Summary statistics table saved to results/Summary_Statistics.csv\n');

%% ============================================================================
%% 6. TOP SIGNIFICANT RESULTS
%% ============================================================================

fprintf('\n=== Top 10 Most Significant Results ===\n');

% Sort by adjusted p-value
sorted_results = sortrows(disagree_peptides, 'adj_pval');
top_results = sorted_results(1:min(10, height(sorted_results)), :);

fprintf('Rank\tProtein\tPeptide\tAdj.P-val\tLog2FC\n');
for i = 1:height(top_results)
    fprintf('%d\t%s\t%s\t%.2e\t%.3f\n', i, ...
        top_results.protein{i}, top_results.peptide{i}, ...
        top_results.adj_pval(i), top_results.log2FC(i));
end

fprintf('\n=== Summary Statistics Generation Complete ===\n');
fprintf('All statistics saved in results/ directory\n');

%% ============================================================================
%% SUPPORTING FUNCTIONS
%% ============================================================================

function summary_table = create_summary_statistics_table(disagree_peptides, protein_counts)
    % Create comprehensive summary statistics table
    
    % Calculate statistics
    total_peptides = height(disagree_peptides);
    total_proteins = height(protein_counts);
    
    sig_001_peptides = sum(disagree_peptides.adj_pval < 0.01);
    sig_005_peptides = sum(disagree_peptides.adj_pval < 0.05);
    sig_01_peptides = sum(disagree_peptides.adj_pval < 0.1);
    
    sig_001_proteins = sum(protein_counts.sig_001_peptides > 0);
    sig_005_proteins = sum(protein_counts.sig_005_peptides > 0);
    
    % Effect sizes for significant peptides
    sig_peptides = disagree_peptides(disagree_peptides.adj_pval < 0.05, :);
    median_abs_fc = median(abs(sig_peptides.log2FC));
    mean_abs_fc = mean(abs(sig_peptides.log2FC));
    
    % Create summary table
    categories = {
        'Analysis Overview';
        'Total Peptides Analyzed';
        'Total Proteins Analyzed';
        '';
        'Significance Thresholds';
        'Significant Peptides (FDR < 0.01)';
        'Significant Peptides (FDR < 0.05)';
        'Significant Peptides (FDR < 0.1)';
        '';
        'Protein-Level Results';
        'Proteins with Significant Peptides (FDR < 0.01)';
        'Proteins with Significant Peptides (FDR < 0.05)';
        '';
        'Effect Size Analysis';
        'Median |Log2FC| (FDR < 0.05)';
        'Mean |Log2FC| (FDR < 0.05)';
        '';
        'Analysis Date';
        'Timestamp'
    };
    
    values = {
        '';
        sprintf('%d', total_peptides);
        sprintf('%d', total_proteins);
        '';
        '';
        sprintf('%d (%.1f%%)', sig_001_peptides, 100*sig_001_peptides/total_peptides);
        sprintf('%d (%.1f%%)', sig_005_peptides, 100*sig_005_peptides/total_peptides);
        sprintf('%d (%.1f%%)', sig_01_peptides, 100*sig_01_peptides/total_peptides);
        '';
        '';
        sprintf('%d (%.1f%%)', sig_001_proteins, 100*sig_001_proteins/total_proteins);
        sprintf('%d (%.1f%%)', sig_005_proteins, 100*sig_005_proteins/total_proteins);
        '';
        '';
        sprintf('%.3f', median_abs_fc);
        sprintf('%.3f', mean_abs_fc);
        '';
        '';
        datestr(now)
    };
    
    summary_table = table(categories, values, 'VariableNames', {'Category', 'Value'});
end 