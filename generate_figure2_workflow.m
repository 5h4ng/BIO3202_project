% generate_figure2_workflow.m
% Script to generate Figure 2: PeCorA workflow demonstration
% Generates 3 figures showing the PeCorA workflow steps with academic visualization

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

% Generate 3 figures
fig2A = create_figure_2A_raw_data(protein_data);
saveas(fig2A, 'results/Figure_2A_Raw_Data.png');
close(fig2A);
fprintf('✓ Figure 2A: Raw Data saved\n');

fig2B = create_figure_2B_scaling_comparison(protein_data, demo_peptide);
saveas(fig2B, 'results/Figure_2B_Global_Scaling.png');
close(fig2B);
fprintf('✓ Figure 2B: Scaling Comparison saved\n');

fig2C = create_figure_2C_peptide_adjustment_control(protein_data);
saveas(fig2C, 'results/Figure_2C_Scaling_Comparison.png');
close(fig2C);
fprintf('✓ Figure 2C: Peptide Adjustment saved\n');

fprintf('=== Figure 2 Generation Complete (3 figures) ===\n');

%% ============================================================================
%% SUPPORTING FUNCTIONS
%% ============================================================================

function fig = create_figure_2A_raw_data(protein_data)
    % Figure 2A: Raw data distributions across all treatment groups
    fig = figure('Position', [100, 100, 800, 500], 'Color', 'white');
    
    conditions = unique(protein_data.Condition, 'stable');
    n_conditions = length(conditions);
    
    % Create violin plots for raw data
    raw_data = [];
    group_labels = [];
    
    for i = 1:n_conditions
        condition = conditions{i};
        condition_data = protein_data(strcmp(protein_data.Condition, condition), :);
        log2_areas = log2(condition_data.NormalizedArea);
        raw_data = [raw_data; log2_areas];
        group_labels = [group_labels; repmat(i, length(log2_areas), 1)];
    end
    
    % Create violin plots
    create_violin_plots(raw_data, group_labels, conditions, [0.7, 0.7, 0.7]);
    
    ylabel_handle = ylabel('Log₂(Peak Area)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    xlabel_handle = xlabel('Treatment Groups', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    title_handle = title('Raw Data Distributions', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    set(ylabel_handle, 'Color', 'black');
    set(xlabel_handle, 'Color', 'black');
    set(title_handle, 'Color', 'black');
    
    set(gca, 'Color', 'white', 'Box', 'on', 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Arial');
    set(gca, 'XColor', 'black', 'YColor', 'black');
    grid on; set(gca, 'GridAlpha', 0.15, 'GridColor', [0.5 0.5 0.5]);
    
    set(gcf, 'InvertHardcopy', 'off');
end

function fig = create_figure_2B_scaling_comparison(protein_data, demo_peptide)
    % Figure 2B: Before and after global scaling comparison
    fig = figure('Position', [100, 100, 1000, 400], 'Color', 'white');
    
    conditions = unique(protein_data.Condition, 'stable');
    n_conditions = length(conditions);
    
    % Colors
    gray_color = [0.5, 0.5, 0.5];
    green_color = [0.2, 0.7, 0.3];
    
    % Get protein accession for legend
    clean_peptide_seq = regexprep(demo_peptide, '_all$', '');
    
    % Before scaling (left panel) - show log2 transformed data
    subplot(1, 2, 1);
    hold on;
    
    for i = 1:n_conditions
        condition = conditions{i};
        condition_data = protein_data(strcmp(protein_data.Condition, condition), :);
        
        % Other peptides data (gray)
        other_data = condition_data(~strcmp(condition_data.modpep_z, demo_peptide), :);
        if height(other_data) > 0
            x_pos = i - 0.2;  % Left side of the condition
            boxplot_data = other_data.ms1log2;  % Use log2 transformed data
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
            boxplot_data = target_data.ms1log2;  % Use log2 transformed data
            create_single_boxplot(x_pos, boxplot_data, green_color);
            
            % Add scatter points on top
            y_vals = boxplot_data;
            x_vals = x_pos + 0.05 * (rand(length(y_vals), 1) - 0.5);
            scatter(x_vals, y_vals, 15, green_color, 'filled', 'MarkerFaceAlpha', 0.8, ...
                'MarkerEdgeColor', 'black', 'LineWidth', 0.5);
        end
    end
    
    set(gca, 'XTick', 1:n_conditions, 'XTickLabel', conditions);
    ylabel_handle = ylabel('Log₂ Transformed Intensity', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    title_handle = title('Before Global Scaling', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

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
    ylabel_handle = ylabel('Global Scaled Intensity', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    title_handle = title('After Global Scaling', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    % Add horizontal line at zero
    yline(0, 'k--', 'LineWidth', 2, 'Alpha', 0.7);
    
    % Create legend
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

function fig = create_figure_2C_peptide_adjustment_control(protein_data)
    % Figure 2C: Peptide-level adjustment effect - Control group only, by replicate
    fig = figure('Position', [100, 100, 800, 500], 'Color', 'white');
    
    % Filter for control group only
    control_data = protein_data(strcmp(protein_data.Condition, 'cntrl'), :);
    
    if height(control_data) == 0
        error('No control data found. Please check your control group name.');
    end
    
    replicates = unique(control_data.BioReplicate, 'stable');
    n_replicates = length(replicates);
    
    % Prepare data for violin plots
    adjusted_data = [];
    group_labels = [];
    replicate_names = {};
    
    for i = 1:n_replicates
        rep = replicates(i);
        rep_data = control_data(control_data.BioReplicate == rep, :);
        adjusted_values = rep_data.ms1adj;
        adjusted_data = [adjusted_data; adjusted_values];
        group_labels = [group_labels; repmat(i, length(adjusted_values), 1)];
        replicate_names{i} = sprintf('Rep %d', rep);
    end
    
    % Create violin plots
    create_violin_plots(adjusted_data, group_labels, replicate_names, [0.2, 0.7, 0.3]);
    
    % Add horizontal line at zero
    hold on;
    yline(0, 'k--', 'LineWidth', 2, 'Alpha', 0.7);
    hold off;
    
    ylabel_handle = ylabel('Peptide-Adjusted Intensity', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    xlabel_handle = xlabel('Biological Replicates (Control Group)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    title_handle = title('Peptide-Level Adjustment Effect', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    set(ylabel_handle, 'Color', 'black');
    set(xlabel_handle, 'Color', 'black');
    set(title_handle, 'Color', 'black');
    
    set(gca, 'Color', 'white', 'Box', 'on', 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Arial');
    set(gca, 'XColor', 'black', 'YColor', 'black');
    grid on; set(gca, 'GridAlpha', 0.15, 'GridColor', [0.5 0.5 0.5]);
    
    % Add text annotation
    mean_val = mean(adjusted_data);
    text(0.02, 0.98, sprintf('Mean: %.3f', mean_val), 'Units', 'normalized', ...
        'FontSize', 11, 'FontName', 'Arial', 'BackgroundColor', 'white', ...
        'EdgeColor', 'black', 'VerticalAlignment', 'top');
    
    set(gcf, 'InvertHardcopy', 'off');
end

function create_violin_plots(data, groups, labels, color)
    % Create violin plots for the given data
    % data: vector of values
    % groups: group assignments for each data point
    % labels: cell array of group labels
    % color: RGB color for the violins
    
    unique_groups = unique(groups);
    n_groups = length(unique_groups);
    
    for i = 1:n_groups
        group_data = data(groups == unique_groups(i));
        
        if length(group_data) < 2
            % If too few points, just plot scatter
            scatter(i * ones(size(group_data)), group_data, 50, color, 'filled', ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'black', 'LineWidth', 1);
            continue;
        end
        
        % Calculate kernel density estimation
        [density, values] = ksdensity(group_data, 'NumPoints', 100);
        
        % Normalize density for width
        max_density = max(density);
        if max_density > 0
            density = density / max_density * 0.4; % Control violin width
        end
        
        % Create violin shape
        x_left = i - density;
        x_right = i + density;
        
        % Plot the violin
        fill([x_left, fliplr(x_right)], [values, fliplr(values)], color, ...
            'FaceAlpha', 0.6, 'EdgeColor', 'black', 'LineWidth', 1);
        hold on;
        
        % Add median line
        median_val = median(group_data);
        plot([i-0.3, i+0.3], [median_val, median_val], 'k-', 'LineWidth', 2);
        
        % Add quartile lines
        q25 = prctile(group_data, 25);
        q75 = prctile(group_data, 75);
        plot([i-0.15, i+0.15], [q25, q25], 'k-', 'LineWidth', 1);
        plot([i-0.15, i+0.15], [q75, q75], 'k-', 'LineWidth', 1);
        
        % Add scatter points for individual data
        jitter = 0.1 * (rand(size(group_data)) - 0.5);
        scatter(i + jitter, group_data, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
    end
    
    % Set x-axis
    set(gca, 'XTick', 1:n_groups, 'XTickLabel', labels);
    xlim([0.5, n_groups + 0.5]);
    
    hold off;
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