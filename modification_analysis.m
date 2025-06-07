% modification_analysis.m
% 分析显著肽段中的修饰情况并与其他肽段进行比较

fprintf('=== Peptide Modification Analysis ===\n');

% Load data
if ~exist('results/PeCorA_results.csv', 'file')
    error('PeCorA results file not found. Please run PeCorA analysis first.');
end

if ~exist('data/PeCorA_noZ.csv', 'file')
    error('Input data file not found: data/PeCorA_noZ.csv');
end

% Read PeCorA results
results = readtable('results/PeCorA_results.csv');

% Read original data
opts = detectImportOptions('data/PeCorA_noZ.csv');
opts.VariableNamingRule = 'preserve';
raw_data = readtable('data/PeCorA_noZ.csv', opts);

% Get unique peptides from original data
unique_peptides_all = unique(raw_data.('Peptide Modified Sequence'));
total_peptides = length(unique_peptides_all);

% Count modified peptides in all data (those containing '[')
modified_all = sum(contains(unique_peptides_all, '['));
unmodified_all = total_peptides - modified_all;
modification_rate_all = modified_all / total_peptides * 100;

fprintf('\n=== Overall Data Statistics ===\n');
fprintf('Total Unique Peptides: %d\n', total_peptides);
fprintf('Modified Peptides: %d\n', modified_all);
fprintf('Unmodified Peptides: %d\n', unmodified_all);
fprintf('Modification Rate: %.1f%%\n', modification_rate_all);

% Define significance threshold (FDR < 0.05)
significant_peptides = results(results.adj_pval < 0.05, :);

% Function to analyze modification in a subset
function stats = analyze_modifications(peptide_list, subset_name)
    n_peptides = height(peptide_list);
    if n_peptides == 0
        stats = struct('name', subset_name, 'total', 0, 'modified', 0, ...
                      'unmodified', 0, 'mod_rate', 0);
        return;
    end
    
    % Extract peptide sequences for modification analysis
    if ismember('peptide', peptide_list.Properties.VariableNames)
        sequences = peptide_list.peptide;
    elseif ismember('Peptide', peptide_list.Properties.VariableNames)
        sequences = peptide_list.Peptide;
    else
        error('Cannot find peptide sequence column');
    end
    
    modified = sum(contains(sequences, '['));
    unmodified = n_peptides - modified;
    mod_rate = modified / n_peptides * 100;
    
    stats = struct('name', subset_name, 'total', n_peptides, ...
                   'modified', modified, 'unmodified', unmodified, ...
                   'mod_rate', mod_rate);
    
    fprintf('\n=== %s ===\n', subset_name);
    fprintf('Total Peptides: %d\n', n_peptides);
    fprintf('Modified Peptides: %d\n', modified);
    fprintf('Unmodified Peptides: %d\n', unmodified);
    fprintf('Modification Rate: %.1f%%\n', mod_rate);
end

% Analyze significant peptides
stats_significant = analyze_modifications(significant_peptides, 'Significant Peptides (FDR < 0.05)');

% Calculate other peptides (all - significant)
other_modified = modified_all - stats_significant.modified;
other_unmodified = unmodified_all - stats_significant.unmodified;
other_total = other_modified + other_unmodified;
other_mod_rate = other_modified / other_total * 100;

fprintf('\n=== Other Peptides (Non-significant) ===\n');
fprintf('Total Peptides: %d\n', other_total);
fprintf('Modified Peptides: %d\n', other_modified);
fprintf('Unmodified Peptides: %d\n', other_unmodified);
fprintf('Modification Rate: %.1f%%\n', other_mod_rate);

% Statistical test (Chi-square test for 2x2 contingency table)
function [p_val, chi2_stat, enrichment] = test_modification_enrichment(sig_mod, sig_total, other_mod, other_total)
    % Create contingency table
    % Rows: Modified/Unmodified, Columns: Significant/Other
    sig_unmod = sig_total - sig_mod;
    other_unmod = other_total - other_mod;
    
    % 2x2 contingency table
    %                Modified    Unmodified
    % Significant      a           b
    % Other            c           d
    a = sig_mod;
    b = sig_unmod;
    c = other_mod;
    d = other_unmod;
    
    % Manual Chi-square calculation for 2x2 table
    n = a + b + c + d;  % Total observations
    
    % Expected frequencies
    expected_a = (a + b) * (a + c) / n;
    expected_b = (a + b) * (b + d) / n;
    expected_c = (c + d) * (a + c) / n;
    expected_d = (c + d) * (b + d) / n;
    
    % Chi-square statistic with Yates' continuity correction
    chi2_stat = n * (abs(a*d - b*c) - n/2)^2 / ((a+b)*(c+d)*(a+c)*(b+d));
    
    % Degrees of freedom = 1 for 2x2 table
    df = 1;
    
    % Calculate p-value
    p_val = 1 - chi2cdf(chi2_stat, df);
    
    % Calculate enrichment (odds ratio)
    if b > 0 && d > 0
        odds_sig = a / b;
        odds_other = c / d;
        enrichment = odds_sig / odds_other;
    else
        enrichment = NaN;  % Cannot calculate if any denominator is 0
    end
    
    fprintf('\n=== Statistical Test Results ===\n');
    fprintf('Contingency Table:\n');
    fprintf('               Modified    Unmodified\n');
    fprintf('Significant      %6d    %6d\n', a, b);
    fprintf('Other            %6d    %6d\n', c, d);
    fprintf('\nExpected Frequencies:\n');
    fprintf('               Modified    Unmodified\n');
    fprintf('Significant      %6.1f    %6.1f\n', expected_a, expected_b);
    fprintf('Other            %6.1f    %6.1f\n', expected_c, expected_d);
    fprintf('\nChi-square = %.3f, df = %d, p-value = %.3e\n', chi2_stat, df, p_val);
    fprintf('Enrichment (Odds Ratio) = %.2f\n', enrichment);
end

% Perform statistical test
[p_val, chi2_stat, enrichment] = test_modification_enrichment(...
    stats_significant.modified, stats_significant.total, other_modified, other_total);

% Create visualization with two pie charts
fig = figure('Position', [100, 100, 900, 450], 'Color', 'white');

% Set default text color to black
set(0, 'DefaultTextColor', 'k');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');

% Subplot 1: Significant peptides pie chart
subplot(1, 2, 1);
pie_data_sig = [stats_significant.modified, stats_significant.unmodified];
labels_sig = {sprintf('Modified\n(%d, %.1f%%)', stats_significant.modified, stats_significant.mod_rate), ...
              sprintf('Unmodified\n(%d, %.1f%%)', stats_significant.unmodified, 100-stats_significant.mod_rate)};

pie_handle_sig = pie(pie_data_sig, labels_sig);

% Set colors and text properties
pie_colors = [0.8, 0.3, 0.3; 0.3, 0.6, 0.8];  % Red for modified, Blue for unmodified
for i = 1:2:length(pie_handle_sig)
    pie_handle_sig(i).FaceColor = pie_colors((i+1)/2, :);
    pie_handle_sig(i).EdgeColor = 'k';
    pie_handle_sig(i).LineWidth = 1;
end

% Set text properties to black
for i = 2:2:length(pie_handle_sig)
    pie_handle_sig(i).Color = 'k';
    pie_handle_sig(i).FontWeight = 'bold';
    pie_handle_sig(i).FontSize = 11;
end

title('Significant Peptides (FDR < 0.05)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

% Subplot 2: Other peptides pie chart
subplot(1, 2, 2);
pie_data_other = [other_modified, other_unmodified];
labels_other = {sprintf('Modified\n(%d, %.1f%%)', other_modified, other_mod_rate), ...
                sprintf('Unmodified\n(%d, %.1f%%)', other_unmodified, 100-other_mod_rate)};

pie_handle_other = pie(pie_data_other, labels_other);

% Set colors and text properties
for i = 1:2:length(pie_handle_other)
    pie_handle_other(i).FaceColor = pie_colors((i+1)/2, :);
    pie_handle_other(i).EdgeColor = 'k';
    pie_handle_other(i).LineWidth = 1;
end

% Set text properties to black
for i = 2:2:length(pie_handle_other)
    pie_handle_other(i).Color = 'k';
    pie_handle_other(i).FontWeight = 'bold';
    pie_handle_other(i).FontSize = 11;
end

title('Other Peptides (Non-significant)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

% Add overall title with statistical results
if p_val < 0.001
    sig_text = '***';
    sig_level = 'p < 0.001';
elseif p_val < 0.01
    sig_text = '**';
    sig_level = 'p < 0.01';
elseif p_val < 0.05
    sig_text = '*';
    sig_level = 'p < 0.05';
else
    sig_text = 'ns';
    sig_level = sprintf('p = %.3f', p_val);
end

main_title = sprintf('Peptide Modification Status Comparison (Enrichment: %.2fx, %s %s)', ...
    enrichment, sig_level, sig_text);
sgtitle(main_title, 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

% Save figure
if ~exist('results/modification_analysis', 'dir')
    mkdir('results/modification_analysis');
end

saveas(fig, 'results/modification_analysis/peptide_modification_comparison.png');
close(fig);

% Create detailed summary table
summary_table = table();
summary_table.Category = {'Significant Peptides (FDR<0.05)'; 'Other Peptides'};
summary_table.Total_Peptides = [stats_significant.total; other_total];
summary_table.Modified_Peptides = [stats_significant.modified; other_modified];
summary_table.Unmodified_Peptides = [stats_significant.unmodified; other_unmodified];
summary_table.Modification_Rate_Percent = [stats_significant.mod_rate; other_mod_rate];

% Save summary table
writetable(summary_table, 'results/modification_analysis/peptide_modification_comparison.csv');

% Print final conclusions
fprintf('\n=== Conclusions ===\n');
fprintf('1. Modification Rate in Significant Peptides: %.1f%% (%d/%d)\n', ...
    stats_significant.mod_rate, stats_significant.modified, stats_significant.total);
fprintf('2. Modification Rate in Other Peptides: %.1f%% (%d/%d)\n', ...
    other_mod_rate, other_modified, other_total);
fprintf('3. Enrichment (Odds Ratio) in Significant Peptides: %.2fx\n', enrichment);

if p_val < 0.001
    fprintf('4. Statistical Test Result: Highly Significant Difference (p < 0.001, Chi² = %.3f)\n', chi2_stat);
    conclusion = 'Modified peptides are significantly more abundant in significant peptides compared to other peptides';
elseif p_val < 0.01
    fprintf('4. Statistical Test Result: Highly Significant Difference (p < 0.01, Chi² = %.3f)\n', chi2_stat);
    conclusion = 'Modified peptides are significantly more abundant in significant peptides compared to other peptides';
elseif p_val < 0.05
    fprintf('4. Statistical Test Result: Significant Difference (p < 0.05, Chi² = %.3f)\n', chi2_stat);
    conclusion = 'Modified peptides are significantly more abundant in significant peptides compared to other peptides';
else
    fprintf('4. Statistical Test Result: No Significant Difference (p = %.3f, Chi² = %.3f)\n', p_val, chi2_stat);
    conclusion = 'Modified peptides are not significantly more abundant in significant peptides compared to other peptides';
end

fprintf('5. %s\n', conclusion);

if p_val < 0.05
    fprintf('6. Biological Significance: Post-translational modifications may be an important factor in peptide expression differences\n');
    fprintf('   - Modified peptides are more likely to show significant expression changes\n');
    fprintf('   - This may reflect the regulatory effect of modifications on protein function and stability\n');
    fprintf('   - PeCorA method is effective in capturing expression inconsistencies related to modifications\n');
else
    fprintf('6. Biological Significance: Post-translational modifications are not significantly associated with peptide expression differences\n');
    fprintf('   - Significant expression changes are mainly driven by other factors\n');
end

fprintf('\n✓ Peptide Modification Analysis Completed, Results Saved to results/modification_analysis/\n'); 