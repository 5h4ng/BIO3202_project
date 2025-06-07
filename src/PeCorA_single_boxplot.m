function fig = PeCorA_single_boxplot(scaled_peptides, target_protein, target_peptide)
%PECORA_SINGLE_BOXPLOT Generate R-style single peptide boxplot
%   Creates a publication-ready plot showing all peptides from a protein with the target peptide highlighted
%   Similar to the original R PeCorA output style
%
%   Inputs:
%   scaled_peptides - Table containing preprocessed data
%   target_protein - Protein identifier
%   target_peptide - Target peptide to highlight
%
%   Output:
%   fig - Figure handle

% Create figure
fig = figure('Position', [100, 100, 800, 600], 'Color', 'white');

% Set figure properties for publication quality
set(fig, 'Color', 'white');
set(fig, 'Units', 'inches');
pos = get(fig, 'Position');
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'inches', 'PaperSize', [pos(3), pos(4)]);

% Create main axes
ax = axes('Parent', fig);
hold(ax, 'on');
box(ax, 'on');

% Set axes properties
set(ax, 'Color', 'white');
set(ax, 'LineWidth', 1.5);
set(ax, 'FontSize', 12);
set(ax, 'FontName', 'Arial');
set(ax, 'XColor', 'black');
set(ax, 'YColor', 'black');
set(ax, 'TickDir', 'out');
set(ax, 'TickLength', [0.02 0.02]);

% Get data for this protein
protein_data = scaled_peptides(strcmp(scaled_peptides.Protein, target_protein), :);

if height(protein_data) == 0
    error('No data found for protein %s', target_protein);
end

% Get unique conditions
conditions = unique(protein_data.Condition, 'stable');
num_conditions = length(conditions);

% Get all peptides for this protein
all_peptides = unique(protein_data.modpep_z);

% Colors: gray for others, green for target
gray_color = [0.5, 0.5, 0.5];  % Darker gray for better contrast
green_color = [0.2, 0.7, 0.3];  % Slightly muted green for professional look

% Initialize data storage for boxplots and statistics
boxplot_data_others = cell(num_conditions, 1);
boxplot_data_target = cell(num_conditions, 1);
p_values = zeros(num_conditions, 1);

% First pass: collect data and calculate statistics
for i = 1:num_conditions
    condition = conditions{i};
    condition_data = protein_data(strcmp(protein_data.Condition, condition), :);
    
    % Separate target peptide from others
    target_data = condition_data(strcmp(condition_data.modpep_z, target_peptide), :);
    other_data = condition_data(~strcmp(condition_data.modpep_z, target_peptide), :);
    
    % Store data for statistics and plotting
    boxplot_data_others{i} = other_data.ms1adj;
    boxplot_data_target{i} = target_data.ms1adj;
    
    % Perform statistical test between target and others for this condition
    if ~isempty(boxplot_data_others{i}) && ~isempty(boxplot_data_target{i})
        [~, p_values(i)] = ttest2(boxplot_data_others{i}, boxplot_data_target{i});
    end
end

% Create boxplots first (so they appear in background)
create_enhanced_boxplot(boxplot_data_others, boxplot_data_target, num_conditions, gray_color, green_color);

% Then plot scatter points on top
for i = 1:num_conditions
    condition = conditions{i};
    condition_data = protein_data(strcmp(protein_data.Condition, condition), :);
    
    % Separate target peptide from others
    target_data = condition_data(strcmp(condition_data.modpep_z, target_peptide), :);
    other_data = condition_data(~strcmp(condition_data.modpep_z, target_peptide), :);
    
    % Plot scatter points for other peptides (gray) - on top of boxes
    if height(other_data) > 0
        x_positions = repmat(i - 0.1, height(other_data), 1) + 0.2 * rand(height(other_data), 1) - 0.1;
        scatter(x_positions, other_data.ms1adj, 40, gray_color, 'filled', ...
            'MarkerFaceAlpha', 0.8, 'LineWidth', 0.5, 'MarkerEdgeColor', 'black');
    end
    
    % Plot scatter points for target peptide (green) - on top of boxes
    if height(target_data) > 0
        x_positions = repmat(i + 0.1, height(target_data), 1) + 0.1 * rand(height(target_data), 1) - 0.05;
        scatter(x_positions, target_data.ms1adj, 50, green_color, 'filled', ...
            'MarkerFaceAlpha', 0.9, 'LineWidth', 0.5, 'MarkerEdgeColor', 'black');
    end
end

% Add statistical significance annotations
add_statistical_annotations(p_values, num_conditions);

% Set up axes with publication-ready styling
set(ax, 'XTick', 1:num_conditions, 'XTickLabel', conditions);

% Labels
xlabel_handle = xlabel('Treatment Group', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel_handle = ylabel('Log₂(Intensity) - Log₂(Control)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
set(xlabel_handle, 'Color', 'black');
set(ylabel_handle, 'Color', 'black');

% Clean up protein and peptide names for display
title_str = strrep(target_protein, '_MOUSE', '');
title_str = regexprep(title_str, '_all$', '');
target_peptide_display = regexprep(target_peptide, '_all$', '');

% Set title
title_handle = title(title_str, 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial', 'Interpreter', 'none');
set(title_handle, 'Color', 'black');

% Create enhanced legend with black border and no vertical lines
h1 = scatter(NaN, NaN, 40, gray_color, 'filled', 'MarkerFaceAlpha', 0.6);
h2 = scatter(NaN, NaN, 50, green_color, 'filled', 'MarkerFaceAlpha', 0.8);
leg = legend([h1, h2], {'Other peptides', target_peptide_display}, ...
    'Location', 'northoutside', ...
    'Orientation', 'horizontal', ...
    'FontSize', 12, ...
    'FontName', 'Arial', ...
    'EdgeColor', 'black', ...
    'Color', 'white', ...
    'Box', 'off');

% Adjust legend position to remove gaps
legPos = get(leg, 'Position');
set(leg, 'Position', [legPos(1), legPos(2)-0.02, legPos(3), legPos(4)]);

% Create black border around the entire plot
set(ax, 'Box', 'on', 'LineWidth', 1.5);

% Set all text elements to use Arial font and black color
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
try
    % Only set TextColor if the property exists
    textObjects = findall(gcf, '-property', 'TextColor');
    if ~isempty(textObjects)
        set(textObjects, 'TextColor', 'black');
    end
catch
    % Ignore if TextColor property doesn't exist
end

% Final styling
grid on;
set(ax, 'GridAlpha', 0.15, 'GridColor', [0.5 0.5 0.5]);
set(ax, 'Color', 'white');
set(gcf, 'InvertHardcopy', 'off');

% Adjust layout
xlim([0.5, num_conditions + 0.5]);

hold off;

% Create a black border around the entire figure
annotation('rectangle', [0.01 0.01 0.98 0.98], 'Color', 'black', 'LineWidth', 1.5);

end

function create_enhanced_boxplot(boxplot_data_others, boxplot_data_target, num_conditions, gray_color, green_color)
    % Create publication-quality boxplots
    for type = 1:2
        if type == 1
            data = boxplot_data_others;
            color = gray_color;
            offset = -0.15;
        else
            data = boxplot_data_target;
            color = green_color;
            offset = 0.15;
        end
        
        for i = 1:num_conditions
            if ~isempty(data{i})
                pos = i + offset;
                
                % Calculate statistics
                q1 = prctile(data{i}, 25);
                q2 = median(data{i});
                q3 = prctile(data{i}, 75);
                iqr = q3 - q1;
                lower_whisker = max(min(data{i}), q1 - 1.5*iqr);
                upper_whisker = min(max(data{i}), q3 + 1.5*iqr);
                
                % Draw enhanced box
                box_width = 0.2;
                rectangle('Position', [pos - box_width/2, q1, box_width, q3-q1], ...
                    'FaceColor', color, 'EdgeColor', 'black', 'LineWidth', 1.5);
                
                % Draw median line
                line([pos - box_width/2, pos + box_width/2], [q2, q2], ...
                    'Color', 'black', 'LineWidth', 2);
                
                % Draw whiskers with enhanced style
                line([pos, pos], [q3, upper_whisker], 'Color', 'black', 'LineWidth', 1.5);
                line([pos, pos], [q1, lower_whisker], 'Color', 'black', 'LineWidth', 1.5);
                line([pos - box_width/4, pos + box_width/4], [upper_whisker, upper_whisker], ...
                    'Color', 'black', 'LineWidth', 1.5);
                line([pos - box_width/4, pos + box_width/4], [lower_whisker, lower_whisker], ...
                    'Color', 'black', 'LineWidth', 1.5);
            end
        end
    end
end

function add_statistical_annotations(p_values, num_conditions)
    % Add statistical significance annotations
    y_range = range(ylim);
    y_top = max(ylim) + 0.05 * y_range;
    
    for i = 1:num_conditions
        if p_values(i) > 0
            if p_values(i) < 0.001
                sig_str = '***';
            elseif p_values(i) < 0.01
                sig_str = '**';
            elseif p_values(i) < 0.05
                sig_str = '*';
            else
                sig_str = 'ns';
            end
            
            % Add significance marker with black color
            text_handle = text(i, y_top, sig_str, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 14, ...
                'FontWeight', 'bold', ...
                'FontName', 'Arial');
            set(text_handle, 'Color', 'black');
        end
    end
    
    % Adjust y-axis to accommodate annotations
    ylim([min(ylim), y_top + 0.1 * y_range]);
end 