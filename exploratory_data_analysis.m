% exploratory_data_analysis.m
% Comprehensive Exploratory Data Analysis for PeCorA input data
% Generates statistical visualizations and summary reports

fprintf('=== 开始探索性数据分析 (EDA) ===\n');

% Add src directory to MATLAB path
addpath('src');

% Load raw data
if ~exist('data/PeCorA_noZ.csv', 'file')
    error('Input data file not found: data/PeCorA_noZ.csv');
end

% Read data
opts = detectImportOptions('data/PeCorA_noZ.csv');
opts.VariableNamingRule = 'preserve';
raw_data = readtable('data/PeCorA_noZ.csv', opts);

fprintf('数据加载完成，开始分析...\n');

%% Basic Data Overview
fprintf('\n=== 基础数据概览 ===\n');
fprintf('总记录数: %d\n', height(raw_data));
fprintf('变量数: %d\n', width(raw_data));
fprintf('数据时间范围: %s\n', datestr(now));

% Create results directory
if ~exist('results/eda', 'dir')
    mkdir('results/eda');
end

%% 1. Generate Overview Statistics
generate_overview_statistics(raw_data);

%% 2. Generate Distribution Plots
generate_distribution_plots(raw_data);

%% 3. Generate Data Quality Assessment
generate_quality_assessment(raw_data);

%% 4. Generate Biological Insights
generate_biological_insights(raw_data);

%% 5. Generate Correlation Analysis
generate_correlation_analysis(raw_data);

fprintf('\n=== EDA分析完成 ===\n');
fprintf('所有图表已保存到: results/eda/\n');

%% ============================================================================
%% SUPPORTING FUNCTIONS
%% ============================================================================

function generate_overview_statistics(data)
    fprintf('\n--- 生成概览统计 ---\n');
    
    % Basic statistics
    n_records = height(data);
    n_peptides = length(unique(data.Peptide));
    n_proteins = length(unique(data.Protein));
    n_conditions = length(unique(data.Condition));
    n_replicates = length(unique(data.BioReplicate));
    
    conditions = unique(data.Condition);
    replicates = unique(data.BioReplicate);
    
    % Calculate per-condition statistics
    cond_stats = table();
    for i = 1:length(conditions)
        cond = conditions{i};
        cond_data = data(strcmp(data.Condition, cond), :);
        cond_stats.Condition{i} = cond;
        cond_stats.Records(i) = height(cond_data);
        cond_stats.Peptides(i) = length(unique(cond_data.Peptide));
        cond_stats.Proteins(i) = length(unique(cond_data.Protein));
        cond_stats.Mean_Area(i) = mean(cond_data.('Normalized Area'));
        cond_stats.Median_Area(i) = median(cond_data.('Normalized Area'));
        cond_stats.Std_Area(i) = std(cond_data.('Normalized Area'));
    end
    
    % Save statistics to file
    writetable(cond_stats, 'results/eda/condition_statistics.csv');
    
    % Print summary
    fprintf('独特肽段数: %d\n', n_peptides);
    fprintf('独特蛋白质数: %d\n', n_proteins);
    fprintf('实验条件数: %d (%s)\n', n_conditions, strjoin(conditions, ', '));
    fprintf('生物学重复数: %d\n', n_replicates);
    
    % Create overview figure
    fig = figure('Position', [100, 100, 1200, 800], 'Color', 'white');
    
    % Subplot 1: Records per condition
    subplot(2, 3, 1);
    bar(cond_stats.Records, 'FaceColor', [0.3, 0.6, 0.8]);
    set(gca, 'XTickLabel', cond_stats.Condition);
    title('Records per Condition', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Number of Records', 'FontSize', 12);
    grid on; set(gca, 'GridAlpha', 0.3);
    
    % Subplot 2: Peptides per condition
    subplot(2, 3, 2);
    bar(cond_stats.Peptides, 'FaceColor', [0.8, 0.4, 0.2]);
    set(gca, 'XTickLabel', cond_stats.Condition);
    title('Unique Peptides per Condition', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Number of Peptides', 'FontSize', 12);
    grid on; set(gca, 'GridAlpha', 0.3);
    
    % Subplot 3: Proteins per condition
    subplot(2, 3, 3);
    bar(cond_stats.Proteins, 'FaceColor', [0.2, 0.7, 0.3]);
    set(gca, 'XTickLabel', cond_stats.Condition);
    title('Unique Proteins per Condition', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Number of Proteins', 'FontSize', 12);
    grid on; set(gca, 'GridAlpha', 0.3);
    
    % Subplot 4: Mean intensity per condition
    subplot(2, 3, 4);
    bar(cond_stats.Mean_Area, 'FaceColor', [0.7, 0.2, 0.6]);
    set(gca, 'XTickLabel', cond_stats.Condition);
    title('Mean Peak Area per Condition', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Mean Peak Area', 'FontSize', 12);
    grid on; set(gca, 'GridAlpha', 0.3);
    
    % Subplot 5: CV per condition
    subplot(2, 3, 5);
    cv_values = cond_stats.Std_Area ./ cond_stats.Mean_Area * 100;
    bar(cv_values, 'FaceColor', [0.9, 0.6, 0.1]);
    set(gca, 'XTickLabel', cond_stats.Condition);
    title('Coefficient of Variation (%)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('CV (%)', 'FontSize', 12);
    grid on; set(gca, 'GridAlpha', 0.3);
    
    % Subplot 6: Summary table
    subplot(2, 3, 6);
    axis off;
    summary_text = {
        sprintf('总数据记录: %d', n_records),
        sprintf('独特肽段: %d', n_peptides),
        sprintf('独特蛋白质: %d', n_proteins),
        sprintf('实验条件: %d', n_conditions),
        sprintf('生物学重复: %d', n_replicates),
        '',
        '峰面积范围:',
        sprintf('最小值: %.2e', min(data.('Normalized Area'))),
        sprintf('最大值: %.2e', max(data.('Normalized Area'))),
        sprintf('中位数: %.2e', median(data.('Normalized Area')))
    };
    text(0.1, 0.9, summary_text, 'FontSize', 11, 'FontName', 'Arial', ...
        'VerticalAlignment', 'top', 'Units', 'normalized');
    
    sgtitle('数据概览统计', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save figure
    saveas(fig, 'results/eda/overview_statistics.png');
    close(fig);
    
    fprintf('✓ 概览统计图表已保存\n');
end

function generate_distribution_plots(data)
    fprintf('\n--- 生成分布图 ---\n');
    
    % Figure 1: Peak Area Distributions
    fig1 = figure('Position', [100, 100, 1000, 600], 'Color', 'white');
    
    conditions = unique(data.Condition, 'stable');
    colors = [0.3, 0.6, 0.8; 0.8, 0.4, 0.2; 0.2, 0.7, 0.3];
    
    % Raw peak area histogram
    subplot(2, 2, 1);
    hold on;
    for i = 1:length(conditions)
        cond = conditions{i};
        cond_data = data(strcmp(data.Condition, cond), :);
        histogram(cond_data.('Normalized Area'), 50, 'FaceColor', colors(i,:), ...
            'FaceAlpha', 0.6, 'EdgeColor', 'black');
    end
    hold off;
    xlabel('Peak Area', 'FontSize', 12);
    ylabel('Frequency', 'FontSize', 12);
    title('Raw Peak Area Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    legend(conditions, 'Location', 'best');
    set(gca, 'XScale', 'log');
    grid on;
    
    % Log2 transformed histogram
    subplot(2, 2, 2);
    hold on;
    for i = 1:length(conditions)
        cond = conditions{i};
        cond_data = data(strcmp(data.Condition, cond), :);
        log_data = log2(cond_data.('Normalized Area'));
        histogram(log_data, 50, 'FaceColor', colors(i,:), ...
            'FaceAlpha', 0.6, 'EdgeColor', 'black');
    end
    hold off;
    xlabel('Log₂(Peak Area)', 'FontSize', 12);
    ylabel('Frequency', 'FontSize', 12);
    title('Log₂ Transformed Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    legend(conditions, 'Location', 'best');
    grid on;
    
    % Box plot by condition
    subplot(2, 2, 3);
    all_data = [];
    all_groups = [];
    for i = 1:length(conditions)
        cond = conditions{i};
        cond_data = data(strcmp(data.Condition, cond), :);
        log_data = log2(cond_data.('Normalized Area'));
        all_data = [all_data; log_data];
        all_groups = [all_groups; repmat(i, length(log_data), 1)];
    end
    boxplot(all_data, all_groups, 'Labels', conditions);
    ylabel('Log₂(Peak Area)', 'FontSize', 12);
    title('Distribution by Condition', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Violin plot by replicate
    subplot(2, 2, 4);
    replicates = unique(data.BioReplicate);
    all_rep_data = [];
    all_rep_groups = [];
    for i = 1:length(replicates)
        rep = replicates(i);
        rep_data = data(data.BioReplicate == rep, :);
        log_data = log2(rep_data.('Normalized Area'));
        all_rep_data = [all_rep_data; log_data];
        all_rep_groups = [all_rep_groups; repmat(i, length(log_data), 1)];
    end
    boxplot(all_rep_data, all_rep_groups, 'Labels', arrayfun(@(x) sprintf('Rep %d', x), replicates, 'UniformOutput', false));
    ylabel('Log₂(Peak Area)', 'FontSize', 12);
    title('Distribution by Replicate', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    sgtitle('Peak Area Distributions', 'FontSize', 16, 'FontWeight', 'bold');
    
    saveas(fig1, 'results/eda/peak_area_distributions.png');
    close(fig1);
    
    fprintf('✓ 分布图已保存\n');
end

function generate_quality_assessment(data)
    fprintf('\n--- 生成数据质量评估 ---\n');
    
    fig = figure('Position', [100, 100, 1200, 800], 'Color', 'white');
    
    conditions = unique(data.Condition, 'stable');
    replicates = unique(data.BioReplicate);
    
    % 1. Missing value analysis
    subplot(2, 3, 1);
    missing_by_condition = [];
    condition_labels = {};
    for i = 1:length(conditions)
        cond = conditions{i};
        cond_data = data(strcmp(data.Condition, cond), :);
        missing_count = sum(isnan(cond_data.('Normalized Area')) | cond_data.('Normalized Area') <= 0);
        missing_pct = missing_count / height(cond_data) * 100;
        missing_by_condition = [missing_by_condition, missing_pct];
        condition_labels{i} = cond;
    end
    bar(missing_by_condition, 'FaceColor', [0.8, 0.3, 0.3]);
    set(gca, 'XTickLabel', condition_labels);
    title('Missing/Zero Values (%)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Percentage (%)', 'FontSize', 12);
    grid on;
    
    % 2. Dynamic range analysis
    subplot(2, 3, 2);
    dynamic_ranges = [];
    for i = 1:length(conditions)
        cond = conditions{i};
        cond_data = data(strcmp(data.Condition, cond), :);
        valid_data = cond_data.('Normalized Area')(cond_data.('Normalized Area') > 0);
        if ~isempty(valid_data)
            dynamic_range = log10(max(valid_data) / min(valid_data));
            dynamic_ranges = [dynamic_ranges, dynamic_range];
        else
            dynamic_ranges = [dynamic_ranges, 0];
        end
    end
    bar(dynamic_ranges, 'FaceColor', [0.3, 0.6, 0.8]);
    set(gca, 'XTickLabel', condition_labels);
    title('Dynamic Range (log₁₀)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Log₁₀(Max/Min)', 'FontSize', 12);
    grid on;
    
    % 3. Data completeness heatmap
    subplot(2, 3, 3);
    % Create peptide-condition matrix
    unique_peptides = unique(data.Peptide);
    completeness_matrix = zeros(length(conditions), length(replicates));
    for i = 1:length(conditions)
        for j = 1:length(replicates)
            subset = data(strcmp(data.Condition, conditions{i}) & data.BioReplicate == replicates(j), :);
            completeness = height(subset) / length(unique_peptides) * 100;
            completeness_matrix(i, j) = completeness;
        end
    end
    imagesc(completeness_matrix);
    colorbar;
    set(gca, 'XTick', 1:length(replicates), 'XTickLabel', arrayfun(@(x) sprintf('Rep%d', x), replicates, 'UniformOutput', false));
    set(gca, 'YTick', 1:length(conditions), 'YTickLabel', conditions);
    title('Data Completeness (%)', 'FontSize', 14, 'FontWeight', 'bold');
    colormap('hot');
    
    % 4. Coefficient of variation
    subplot(2, 3, 4);
    cv_by_peptide = [];
    for i = 1:min(100, length(unique_peptides))  % Sample first 100 peptides
        peptide = unique_peptides{i};
        peptide_data = data(strcmp(data.Peptide, peptide), :);
        if height(peptide_data) > 1
            cv = std(peptide_data.('Normalized Area')) / mean(peptide_data.('Normalized Area')) * 100;
            cv_by_peptide = [cv_by_peptide, cv];
        end
    end
    histogram(cv_by_peptide, 30, 'FaceColor', [0.6, 0.4, 0.8], 'EdgeColor', 'black');
    xlabel('Coefficient of Variation (%)', 'FontSize', 12);
    ylabel('Frequency', 'FontSize', 12);
    title('CV Distribution (Sample)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % 5. Replicate correlation
    subplot(2, 3, 5);
    if length(replicates) >= 2
        rep1_data = data(data.BioReplicate == replicates(1), :);
        rep2_data = data(data.BioReplicate == replicates(2), :);
        
        % Match peptides
        [common_peptides, ia, ib] = intersect(rep1_data.Peptide, rep2_data.Peptide);
        if length(common_peptides) > 10
            x_data = log2(rep1_data.('Normalized Area')(ia));
            y_data = log2(rep2_data.('Normalized Area')(ib));
            scatter(x_data, y_data, 20, 'filled', 'MarkerFaceAlpha', 0.6);
            xlabel(sprintf('Rep %d Log₂(Area)', replicates(1)), 'FontSize', 12);
            ylabel(sprintf('Rep %d Log₂(Area)', replicates(2)), 'FontSize', 12);
            title('Replicate Correlation', 'FontSize', 14, 'FontWeight', 'bold');
            
            % Add correlation coefficient
            r = corrcoef(x_data, y_data);
            text(0.05, 0.95, sprintf('r = %.3f', r(1,2)), 'Units', 'normalized', ...
                'FontSize', 12, 'BackgroundColor', 'white');
            grid on;
        end
    end
    
    % 6. Quality metrics summary
    subplot(2, 3, 6);
    axis off;
    
    % Calculate quality metrics
    total_missing = sum(isnan(data.('Normalized Area')) | data.('Normalized Area') <= 0);
    missing_pct = total_missing / height(data) * 100;
    overall_cv = std(data.('Normalized Area')) / mean(data.('Normalized Area')) * 100;
    
    quality_text = {
        '数据质量指标:',
        '',
        sprintf('缺失值比例: %.1f%%', missing_pct),
        sprintf('总体CV: %.1f%%', overall_cv),
        sprintf('动态范围: %.1f log₁₀', log10(max(data.('Normalized Area')) / min(data.('Normalized Area')(data.('Normalized Area') > 0)))),
        '',
        '推荐阈值:',
        '• 缺失值 < 20%',
        '• CV < 30%',
        '• 动态范围 > 3 log₁₀'
    };
    
    text(0.1, 0.9, quality_text, 'FontSize', 11, 'FontName', 'Arial', ...
        'VerticalAlignment', 'top', 'Units', 'normalized');
    
    sgtitle('Data Quality Assessment', 'FontSize', 16, 'FontWeight', 'bold');
    
    saveas(fig, 'results/eda/quality_assessment.png');
    close(fig);
    
    fprintf('✓ 数据质量评估图表已保存\n');
end

function generate_biological_insights(data)
    fprintf('\n--- 生成生物学洞察 ---\n');
    
    fig = figure('Position', [100, 100, 1200, 800], 'Color', 'white');
    
    % 1. Peptides per protein distribution
    subplot(2, 3, 1);
    protein_peptide_counts = [];
    unique_proteins = unique(data.Protein);
    for i = 1:length(unique_proteins)
        protein = unique_proteins{i};
        peptides_for_protein = unique(data.Peptide(strcmp(data.Protein, protein)));
        protein_peptide_counts = [protein_peptide_counts, length(peptides_for_protein)];
    end
    histogram(protein_peptide_counts, 'BinEdges', 0.5:1:max(protein_peptide_counts)+0.5, ...
        'FaceColor', [0.4, 0.7, 0.4], 'EdgeColor', 'black');
    xlabel('Peptides per Protein', 'FontSize', 12);
    ylabel('Number of Proteins', 'FontSize', 12);
    title('Peptide Coverage Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % 2. Protein length analysis (based on peptide positions)
    subplot(2, 3, 2);
    protein_lengths = [];
    for i = 1:length(unique_proteins)
        protein = unique_proteins{i};
        protein_data = data(strcmp(data.Protein, protein), :);
        if ~isempty(protein_data)
            max_pos = max(protein_data.('End Pos'));
            protein_lengths = [protein_lengths, max_pos];
        end
    end
    histogram(protein_lengths, 30, 'FaceColor', [0.7, 0.4, 0.6], 'EdgeColor', 'black');
    xlabel('Protein Length (aa)', 'FontSize', 12);
    ylabel('Number of Proteins', 'FontSize', 12);
    title('Protein Length Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % 3. Peptide length distribution
    subplot(2, 3, 3);
    peptide_lengths = [];
    unique_peptides = unique(data.Peptide);
    for i = 1:length(unique_peptides)
        peptide_lengths = [peptide_lengths, length(unique_peptides{i})];
    end
    histogram(peptide_lengths, 'BinEdges', 4.5:1:max(peptide_lengths)+0.5, ...
        'FaceColor', [0.6, 0.6, 0.3], 'EdgeColor', 'black');
    xlabel('Peptide Length (aa)', 'FontSize', 12);
    ylabel('Number of Peptides', 'FontSize', 12);
    title('Peptide Length Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % 4. Modification frequency
    subplot(2, 3, 4);
    modified_peptides = sum(contains(data.('Peptide Modified Sequence'), '['));
    unmodified_peptides = height(data) - modified_peptides;
    pie([modified_peptides, unmodified_peptides], {'Modified', 'Unmodified'});
    title('Peptide Modification Status', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 5. Top proteins by peptide count
    subplot(2, 3, 5);
    [sorted_counts, sort_idx] = sort(protein_peptide_counts, 'descend');
    top_proteins = unique_proteins(sort_idx(1:min(10, length(sort_idx))));
    top_counts = sorted_counts(1:min(10, length(sorted_counts)));
    
    barh(top_counts, 'FaceColor', [0.8, 0.5, 0.2]);
    xlabel('Number of Peptides', 'FontSize', 12);
    title('Top 10 Proteins by Peptide Count', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Clean protein names for display
    clean_names = cell(size(top_proteins));
    for i = 1:length(top_proteins)
        name = top_proteins{i};
        if contains(name, '|')
            parts = split(name, '|');
            if length(parts) >= 3
                clean_names{i} = parts{3};
            else
                clean_names{i} = parts{end};
            end
        else
            clean_names{i} = name;
        end
        if length(clean_names{i}) > 15
            clean_names{i} = [clean_names{i}(1:12) '...'];
        end
    end
    set(gca, 'YTick', 1:length(top_proteins), 'YTickLabel', clean_names);
    grid on;
    
    % 6. Summary statistics
    subplot(2, 3, 6);
    axis off;
    
    summary_text = {
        '生物学统计摘要:',
        '',
        sprintf('独特蛋白质: %d', length(unique_proteins)),
        sprintf('独特肽段: %d', length(unique_peptides)),
        sprintf('平均每蛋白肽段数: %.1f', mean(protein_peptide_counts)),
        sprintf('修饰肽段比例: %.1f%%', modified_peptides/height(data)*100),
        '',
        sprintf('肽段长度范围: %d-%d aa', min(peptide_lengths), max(peptide_lengths)),
        sprintf('平均肽段长度: %.1f aa', mean(peptide_lengths)),
        '',
        sprintf('蛋白质长度范围: %d-%d aa', min(protein_lengths), max(protein_lengths))
    };
    
    text(0.1, 0.9, summary_text, 'FontSize', 11, 'FontName', 'Arial', ...
        'VerticalAlignment', 'top', 'Units', 'normalized');
    
    sgtitle('Biological Insights', 'FontSize', 16, 'FontWeight', 'bold');
    
    saveas(fig, 'results/eda/biological_insights.png');
    close(fig);
    
    fprintf('✓ 生物学洞察图表已保存\n');
end

function generate_correlation_analysis(data)
    fprintf('\n--- 生成相关性分析 ---\n');
    
    conditions = unique(data.Condition, 'stable');
    replicates = unique(data.BioReplicate);
    
    % Create condition-replicate correlation matrix
    n_samples = length(conditions) * length(replicates);
    sample_names = {};
    sample_data = [];
    
    % Get common peptides across all samples
    common_peptides = unique(data.Peptide);
    for i = 1:length(conditions)
        for j = 1:length(replicates)
            subset = data(strcmp(data.Condition, conditions{i}) & data.BioReplicate == replicates(j), :);
            available_peptides = unique(subset.Peptide);
            common_peptides = intersect(common_peptides, available_peptides);
        end
    end
    
    fprintf('发现 %d 个共同肽段用于相关性分析\n', length(common_peptides));
    
    if length(common_peptides) < 10
        fprintf('警告: 共同肽段数量太少，跳过相关性分析\n');
        return;
    end
    
    % Build data matrix
    sample_counter = 1;
    for i = 1:length(conditions)
        for j = 1:length(replicates)
            sample_names{sample_counter} = sprintf('%s_Rep%d', conditions{i}, replicates(j));
            
            subset = data(strcmp(data.Condition, conditions{i}) & data.BioReplicate == replicates(j), :);
            sample_values = zeros(length(common_peptides), 1);
            
            for k = 1:length(common_peptides)
                peptide_row = subset(strcmp(subset.Peptide, common_peptides{k}), :);
                if height(peptide_row) > 0
                    sample_values(k) = log2(peptide_row.('Normalized Area')(1));
                else
                    sample_values(k) = NaN;
                end
            end
            
            sample_data = [sample_data, sample_values];
            sample_counter = sample_counter + 1;
        end
    end
    
    % Remove peptides with too many NaN values
    valid_peptides = sum(~isnan(sample_data), 2) >= size(sample_data, 2) * 0.7;
    sample_data = sample_data(valid_peptides, :);
    common_peptides = common_peptides(valid_peptides);
    
    fprintf('使用 %d 个有效肽段进行相关性分析\n', length(common_peptides));
    
    if size(sample_data, 1) < 5
        fprintf('警告: 有效肽段数量太少，跳过相关性分析\n');
        return;
    end
    
    % Calculate correlation matrix
    corr_matrix = corrcoef(sample_data', 'Rows', 'complete');
    
    % Create correlation heatmap
    fig = figure('Position', [100, 100, 1000, 800], 'Color', 'white');
    
    % Main correlation heatmap
    subplot(2, 2, [1, 3]);
    imagesc(corr_matrix);
    colorbar;
    colormap('RdYlBu');
    caxis([-1, 1]);
    
    set(gca, 'XTick', 1:length(sample_names), 'XTickLabel', sample_names, 'XTickLabelRotation', 45);
    set(gca, 'YTick', 1:length(sample_names), 'YTickLabel', sample_names);
    title('Sample Correlation Matrix', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add correlation values as text
    for i = 1:size(corr_matrix, 1)
        for j = 1:size(corr_matrix, 2)
            if abs(corr_matrix(i,j)) > 0.5  % Only show strong correlations
                text(j, i, sprintf('%.2f', corr_matrix(i,j)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8, ...
                    'Color', 'white', 'FontWeight', 'bold');
            end
        end
    end
    
    % Correlation distribution
    subplot(2, 2, 2);
    upper_tri = triu(corr_matrix, 1);
    upper_tri_vals = upper_tri(upper_tri ~= 0);
    histogram(upper_tri_vals, 20, 'FaceColor', [0.6, 0.4, 0.8], 'EdgeColor', 'black');
    xlabel('Correlation Coefficient', 'FontSize', 12);
    ylabel('Frequency', 'FontSize', 12);
    title('Correlation Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Average correlation by condition
    subplot(2, 2, 4);
    condition_corrs = [];
    condition_labels = {};
    
    for i = 1:length(conditions)
        cond_indices = contains(sample_names, conditions{i});
        cond_corr_matrix = corr_matrix(cond_indices, cond_indices);
        upper_tri = triu(cond_corr_matrix, 1);
        avg_corr = mean(upper_tri(upper_tri ~= 0));
        condition_corrs = [condition_corrs, avg_corr];
        condition_labels{i} = conditions{i};
    end
    
    bar(condition_corrs, 'FaceColor', [0.3, 0.7, 0.5]);
    set(gca, 'XTickLabel', condition_labels);
    ylabel('Average Correlation', 'FontSize', 12);
    title('Within-Condition Correlation', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    ylim([0, 1]);
    
    sgtitle('Sample Correlation Analysis', 'FontSize', 16, 'FontWeight', 'bold');
    
    saveas(fig, 'results/eda/correlation_analysis.png');
    close(fig);
    
    % Save correlation matrix
    corr_table = array2table(corr_matrix, 'VariableNames', sample_names, 'RowNames', sample_names);
    writetable(corr_table, 'results/eda/correlation_matrix.csv', 'WriteRowNames', true);
    
    fprintf('✓ 相关性分析图表已保存\n');
end 