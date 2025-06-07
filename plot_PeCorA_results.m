% plot_PeCorA_results.m
% Script to generate plots from saved PeCorA analysis results
% 
% This script reads the saved results from PeCorA analysis and generates
% visualization plots for significant peptides

fprintf('=== PeCorA Results Plotting Script ===\n');

% Add src directory to MATLAB path
addpath('src');

% Check if results directory exists
if ~exist('results', 'dir')
    error('Results directory not found. Please run run_PeCorA first to generate results.');
end

% Check if required files exist
results_file = 'results/PeCorA_results.csv';
scaled_peptides_file = 'results/PeCorA_scaled_peptides.csv';

if ~exist(results_file, 'file')
    error('Results file not found: %s\nPlease run run_PeCorA first.', results_file);
end

if ~exist(scaled_peptides_file, 'file')
    error('Scaled peptides file not found: %s\nPlease run run_PeCorA first.', scaled_peptides_file);
end

% Load saved results
fprintf('Loading saved results...\n');
try
    % Load disagree peptides results
    disagree_peptides = readtable(results_file);
    fprintf('Loaded results: %d peptides\n', height(disagree_peptides));
    
    % Load scaled peptides data
    scaled_peptides = readtable(scaled_peptides_file);
    fprintf('Loaded scaled data: %d rows\n', height(scaled_peptides));
    
catch ME
    fprintf('Error loading results: %s\n', ME.message);
    rethrow(ME);
end

% Display analysis summary
fprintf('\nAnalysis Summary:\n');
fprintf('Total peptides analyzed: %d\n', height(disagree_peptides));
if height(disagree_peptides) > 0
    sig_05 = sum(disagree_peptides.adj_pval < 0.05);
    sig_01 = sum(disagree_peptides.adj_pval < 0.01);
    fprintf('Significant peptides (p < 0.05): %d\n', sig_05);
    fprintf('Significant peptides (p < 0.01): %d\n', sig_01);
    fprintf('Proteins with significant peptides: %d\n', ...
        length(unique(disagree_peptides.protein(disagree_peptides.adj_pval < 0.05))));
end

% Generate visualizations
fprintf('\nGenerating visualizations...\n');
try
    % Get significant peptides (p < 0.05)
    sig_peptides = disagree_peptides(disagree_peptides.adj_pval < 0.05, :);
    
    if height(sig_peptides) > 0
        num_plots = min(10, height(sig_peptides));  % Limit to 10 plots max
        fprintf('Creating %d plots for significant peptides...\n', num_plots);
        
        for i = 1:num_plots
            fprintf('  Creating plot %d/%d: %s from %s\n', i, num_plots, ...
                sig_peptides.peptide{i}, sig_peptides.protein{i});
            
            try
                % Create R-style single boxplot
                fig = PeCorA_single_boxplot(scaled_peptides, ...
                    sig_peptides.protein{i}, sig_peptides.peptide{i});
                
                % Save figure
                filename = sprintf('results/PeCorA_plot_%d.png', i);
                saveas(fig, filename);
                close(fig);
                fprintf('    Saved: %s\n', filename);
                
            catch plot_error
                fprintf('    Error creating plot %d: %s\n', i, plot_error.message);
                % Continue with next plot
            end
        end
    else
        % No significant peptides, plot top 5 by p-value
        fprintf('No significant peptides found. Creating plots for top 5 peptides by p-value...\n');
        
        if height(disagree_peptides) > 0
            num_plots = min(5, height(disagree_peptides));
            
            for i = 1:num_plots
                fprintf('  Creating plot %d/%d: %s from %s (p=%.2e)\n', i, num_plots, ...
                    disagree_peptides.peptide{i}, disagree_peptides.protein{i}, ...
                    disagree_peptides.adj_pval(i));
                
                try
                    % Create R-style single boxplot
                    fig = PeCorA_single_boxplot(scaled_peptides, ...
                        disagree_peptides.protein{i}, disagree_peptides.peptide{i});
                    
                    % Save figure
                    filename = sprintf('results/PeCorA_plot_%d.png', i);
                    saveas(fig, filename);
                    close(fig);
                    fprintf('    Saved: %s\n', filename);
                    
                catch plot_error
                    fprintf('    Error creating plot %d: %s\n', i, plot_error.message);
                    % Continue with next plot
                end
            end
        else
            fprintf('No peptides to plot\n');
        end
    end
    
    fprintf('\nPlotting completed! Check results directory for PNG files.\n');
    
catch ME
    fprintf('Error during visualization: %s\n', ME.message);
    rethrow(ME);
end

fprintf('\n=== Plotting Complete ===\n'); 