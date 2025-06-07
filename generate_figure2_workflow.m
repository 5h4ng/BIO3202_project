% generate_figure2_workflow.m
% Script to generate Figure 2: PeCorA workflow demonstration
% Generates 3 simple figures showing the PeCorA workflow steps

fprintf('=== Generating Figure 2: PeCorA Workflow (3 figures) ===\n');

% Add src directory to MATLAB path
addpath('src');

% Load data
if ~exist('results/PeCorA_scaled_peptides.csv', 'file')
    error('Scaled peptides file not found. Please run PeCorA analysis first.');
end

scaled_peptides = readtable('results/PeCorA_scaled_peptides.csv');

% Select a protein with significant peptides for demonstration
if exist('results/PeCorA_results.csv', 'file')
    results = readtable('results/PeCorA_results.csv');
    demo_protein = results.protein{1};  % Most significant protein
    demo_peptide = results.peptide{1};  % Most significant peptide
else
    % Fallback to a known protein
    demo_protein = 'sp|Q9DBC7|KAP0_MOUSE';
    demo_peptide = 'LGPSDYFGEIALLM[+16]NRPR_all';
end

fprintf('Using demonstration protein: %s\n', demo_protein);
fprintf('Using demonstration peptide: %s\n', demo_peptide);

% Get data for demonstration protein
protein_data = scaled_peptides(strcmp(scaled_peptides.Protein, demo_protein), :);

if height(protein_data) == 0
    warning('Demo protein not found, using first available protein');
    available_proteins = unique(scaled_peptides.Protein);
    demo_protein = available_proteins{1};
    protein_data = scaled_peptides(strcmp(scaled_peptides.Protein, demo_protein), :);
    demo_peptide = protein_data.modpep_z{1};
end

% Generate 3 figures (removed workflow diagram)
fig2A = create_figure_2A_raw_data(protein_data);
saveas(fig2A, 'results/Figure_2A_Raw_Data.png');
close(fig2A);
fprintf('✓ Figure 2A: Raw Data saved\n');

fig2B = create_figure_2B_global_scaling(protein_data);
saveas(fig2B, 'results/Figure_2B_Global_Scaling.png');
close(fig2B);
fprintf('✓ Figure 2B: Global Scaling saved\n');

fig2C = create_figure_2C_scaling_comparison(protein_data, demo_peptide);
saveas(fig2C, 'results/Figure_2C_Scaling_Comparison.png');
close(fig2C);
fprintf('✓ Figure 2C: Scaling Comparison saved\n');

fprintf('=== Figure 2 Generation Complete (3 figures) ===\n');

%% ============================================================================
%% SUPPORTING FUNCTIONS
%% ============================================================================

function fig = create_figure_2A_raw_data(protein_data)
    % Figure 2A: Raw data distributions
    fig = figure('Position', [100, 100, 600, 400], 'Color', 'white');
    
    conditions = unique(protein_data.Condition, 'stable');
    n_conditions = length(conditions);
    
    % Create box plots for raw data
    raw_data = [];
    group_labels = [];
    
    for i = 1:n_conditions
        condition = conditions{i};
        condition_data = protein_data(strcmp(protein_data.Condition, condition), :);
        log2_areas = log2(condition_data.NormalizedArea);
        raw_data = [raw_data; log2_areas];
        group_labels = [group_labels; repmat(i, length(log2_areas), 1)];
    end
    
    boxplot(raw_data, group_labels, 'Labels', conditions, 'Colors', [0.7, 0.7, 0.7]);
    ylabel_handle = ylabel('Log₂(Peak Area)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    xlabel_handle = xlabel('Treatment Groups', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    title_handle = title('A. Raw Data Distributions', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    set(ylabel_handle, 'Color', 'black');
    set(xlabel_handle, 'Color', 'black');
    set(title_handle, 'Color', 'black');
    
    set(gca, 'Color', 'white', 'Box', 'on', 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Arial');
    set(gca, 'XColor', 'black', 'YColor', 'black');
    grid on; set(gca, 'GridAlpha', 0.15, 'GridColor', [0.5 0.5 0.5]);
    
    set(gcf, 'InvertHardcopy', 'off');
end

function fig = create_figure_2B_global_scaling(protein_data)
    % Figure 2B: After global scaling
    fig = figure('Position', [100, 100, 600, 400], 'Color', 'white');
    
    conditions = unique(protein_data.Condition, 'stable');
    n_conditions = length(conditions);
    
    % Create box plots for scaled data
    scaled_data = [];
    group_labels = [];
    
    for i = 1:n_conditions
        condition = conditions{i};
        condition_data = protein_data(strcmp(protein_data.Condition, condition), :);
        scaled_values = condition_data.ms1scaled;
        scaled_data = [scaled_data; scaled_values];
        group_labels = [group_labels; repmat(i, length(scaled_values), 1)];
    end
    
    boxplot(scaled_data, group_labels, 'Labels', conditions, 'Colors', [0.5, 0.5, 0.5]);
    ylabel_handle = ylabel('Log₂(Scaled Peak Area)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    xlabel_handle = xlabel('Treatment Groups', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    title_handle = title('B. Global Scaling', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    set(ylabel_handle, 'Color', 'black');
    set(xlabel_handle, 'Color', 'black');
    set(title_handle, 'Color', 'black');
    
    % Add horizontal line at zero
    hold on;
    yline(0, 'k--', 'LineWidth', 2, 'Alpha', 0.8);
    hold off;
    
    set(gca, 'Color', 'white', 'Box', 'on', 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Arial');
    set(gca, 'XColor', 'black', 'YColor', 'black');
    grid on; set(gca, 'GridAlpha', 0.15, 'GridColor', [0.5 0.5 0.5]);
    
    % Print the actual mean to check
    fprintf('Global scaling mean: %.3f\n', mean(scaled_data));
    
    set(gcf, 'InvertHardcopy', 'off');
end

function fig = create_figure_2C_scaling_comparison(protein_data, demo_peptide)
    % Figure 2C: Before and after scaling comparison with single-track style
    fig = figure('Position', [100, 100, 1000, 400], 'Color', 'white');
    
    conditions = unique(protein_data.Condition, 'stable');
    n_conditions = length(conditions);
    
    % Colors
    gray_color = [0.5, 0.5, 0.5];
    green_color = [0.2, 0.7, 0.3];
    
    % Get protein accession for legend
    protein_accession = protein_data.Protein{1};
    clean_peptide_seq = regexprep(demo_peptide, '_all$', '');
    
    % Before scaling (left panel)
    subplot(1, 2, 1);
    hold on;
    
    for i = 1:n_conditions
        condition = conditions{i};
        condition_data = protein_data(strcmp(protein_data.Condition, condition), :);
        
        % Other peptides data (gray)
        other_data = condition_data(~strcmp(condition_data.modpep_z, demo_peptide), :);
        if height(other_data) > 0
            x_pos = i - 0.2;  % Left side of the condition
            boxplot_data = other_data.ms1scaled;
            create_single_boxplot(x_pos, boxplot_data, gray_color);
            
            % Add scatter points on top
            y_vals = boxplot_data;
            x_vals = x_pos + 0.1 * (rand(length(y_vals), 1) - 0.5);
            scatter(x_vals, y_vals, 15, gray_color, 'filled', 'MarkerFaceAlpha', 0.6, ...
                'MarkerEdgeColor', 'black', 'LineWidth', 0.5);
        end
        
        % Target peptide data (green)
        target_data = condition_data(strcmp(condition_data.modpep_z, demo_peptide), :);
        if height(target_data) > 0
            x_pos = i + 0.2;  % Right side of the condition
            boxplot_data = target_data.ms1scaled;
            create_single_boxplot(x_pos, boxplot_data, green_color);
            
            % Add scatter points on top
            y_vals = boxplot_data;
            x_vals = x_pos + 0.05 * (rand(length(y_vals), 1) - 0.5);
            scatter(x_vals, y_vals, 15, green_color, 'filled', 'MarkerFaceAlpha', 0.8, ...
                'MarkerEdgeColor', 'black', 'LineWidth', 0.5);
        end
    end
    
    set(gca, 'XTick', 1:n_conditions, 'XTickLabel', conditions);
    ylabel_handle = ylabel('Log₂(Intensity)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    title_handle = title('Before Peptide Scaling', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    % Create legend with protein and peptide info
    h1 = scatter(NaN, NaN, 15, gray_color, 'filled', 'MarkerFaceAlpha', 0.6);
    h2 = scatter(NaN, NaN, 15, green_color, 'filled', 'MarkerFaceAlpha', 0.8);
    
    legend_labels = {
        sprintf('Other peptides'), ...
        sprintf(clean_peptide_seq)
    };
    
    leg = legend([h1, h2], legend_labels, ...
        'Location', 'best', 'FontSize', 10, 'FontName', 'Arial');
    set(leg, 'TextColor', 'black', 'Color', 'white', 'EdgeColor', 'black');
    
    set(ylabel_handle, 'Color', 'black');
    set(title_handle, 'Color', 'black');
    set(gca, 'Color', 'white', 'Box', 'on', 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Arial');
    set(gca, 'XColor', 'black', 'YColor', 'black');
    grid on; set(gca, 'GridAlpha', 0.15, 'GridColor', [0.5 0.5 0.5]);
    xlim([0.5, n_conditions + 0.5]);
    
    hold off;
    
    % After scaling (right panel)
    subplot(1, 2, 2);
    hold on;
    
    for i = 1:n_conditions
        condition = conditions{i};
        condition_data = protein_data(strcmp(protein_data.Condition, condition), :);
        
        % Other peptides data (gray)
        other_data = condition_data(~strcmp(condition_data.modpep_z, demo_peptide), :);
        if height(other_data) > 0
            x_pos = i - 0.2;  % Left side of the condition
            boxplot_data = other_data.ms1adj;
            create_single_boxplot(x_pos, boxplot_data, gray_color);
            
            % Add scatter points on top
            y_vals = boxplot_data;
            x_vals = x_pos + 0.1 * (rand(length(y_vals), 1) - 0.5);
            scatter(x_vals, y_vals, 15, gray_color, 'filled', 'MarkerFaceAlpha', 0.6, ...
                'MarkerEdgeColor', 'black', 'LineWidth', 0.5);
        end
        
        % Target peptide data (green)
        target_data = condition_data(strcmp(condition_data.modpep_z, demo_peptide), :);
        if height(target_data) > 0
            x_pos = i + 0.2;  % Right side of the condition
            boxplot_data = target_data.ms1adj;
            create_single_boxplot(x_pos, boxplot_data, green_color);
            
            % Add scatter points on top
            y_vals = boxplot_data;
            x_vals = x_pos + 0.05 * (rand(length(y_vals), 1) - 0.5);
            scatter(x_vals, y_vals, 15, green_color, 'filled', 'MarkerFaceAlpha', 0.8, ...
                'MarkerEdgeColor', 'black', 'LineWidth', 0.5);
        end
    end
    
    set(gca, 'XTick', 1:n_conditions, 'XTickLabel', conditions);
    ylabel_handle = ylabel('Log₂(Intensity)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    title_handle = title('After Peptide Scaling', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    % Create legend with protein and peptide info
    h1 = scatter(NaN, NaN, 15, gray_color, 'filled', 'MarkerFaceAlpha', 0.6);
    h2 = scatter(NaN, NaN, 15, green_color, 'filled', 'MarkerFaceAlpha', 0.8);
    
    leg = legend([h1, h2], legend_labels, ...
        'Location', 'best', 'FontSize', 10, 'FontName', 'Arial');
    set(leg, 'TextColor', 'black', 'Color', 'white', 'EdgeColor', 'black');
    
    set(ylabel_handle, 'Color', 'black');
    set(title_handle, 'Color', 'black');
    set(gca, 'Color', 'white', 'Box', 'on', 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Arial');
    set(gca, 'XColor', 'black', 'YColor', 'black');
    grid on; set(gca, 'GridAlpha', 0.15, 'GridColor', [0.5 0.5 0.5]);
    xlim([0.5, n_conditions + 0.5]);
    
    hold off;
    
    set(gcf, 'InvertHardcopy', 'off');
end

function create_single_boxplot(x_pos, data, color)
    % Create a single boxplot at specified position (similar to PeCorA_single_boxplot)
    if isempty(data) || length(data) < 2
        return;
    end
    
    % Calculate statistics
    q1 = prctile(data, 25);
    q2 = median(data);
    q3 = prctile(data, 75);
    iqr = q3 - q1;
    lower_whisker = max(min(data), q1 - 1.5*iqr);
    upper_whisker = min(max(data), q3 + 1.5*iqr);
    
    % Box width
    width = 0.3;
    
    % Draw box
    rectangle('Position', [x_pos - width/2, q1, width, q3-q1], ...
        'FaceColor', color, 'EdgeColor', 'black', 'LineWidth', 1);
    
    % Draw median line
    line([x_pos - width/2, x_pos + width/2], [q2, q2], ...
        'Color', 'black', 'LineWidth', 2);
    
    % Draw whiskers
    line([x_pos, x_pos], [q3, upper_whisker], 'Color', 'black', 'LineWidth', 1);
    line([x_pos, x_pos], [q1, lower_whisker], 'Color', 'black', 'LineWidth', 1);
    line([x_pos - width/4, x_pos + width/4], [upper_whisker, upper_whisker], 'Color', 'black', 'LineWidth', 1);
    line([x_pos - width/4, x_pos + width/4], [lower_whisker, lower_whisker], 'Color', 'black', 'LineWidth', 1);
end 