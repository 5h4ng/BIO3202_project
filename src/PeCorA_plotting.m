function [fig] = PeCorA_plotting(w, u, v)
%PECORA_PLOTTING Generate plots for PeCorA results
%   This function creates boxplots showing the distribution of peptide
%   intensities for proteins with discordant peptides
%
%   Inputs:
%   w - Table of disagree peptides (PeCorA output)
%   u - Selected row from disagree peptides
%   v - Scaled peptides table
%
%   Outputs:
%   fig - Figure handle containing the plot

% Get significant proteins
sign_prots = unique(u.protein(u.adj_pval <= 0.01));

% Create figure
fig = figure('Position', [100, 100, 800, 600]);

% Process each significant protein
for i = 1:length(sign_prots)
    x = sign_prots{i};
    tmpdf = v(strcmp(v.Protein, x), :);
    tmpdf.allothers = repmat({'allothers'}, height(tmpdf), 1);
    
    % Get significant peptides for this protein
    tmpalldf = w(strcmp(w.protein, x), :);
    tmp_peps = unique(tmpalldf.peptide(tmpalldf.adj_pval < 0.01));
    
    if ~isempty(tmp_peps)
        for j = 1:length(tmp_peps)
            y = tmp_peps{j};
            subtmpdf = tmpdf;
            subtmpdf.allothers(strcmp(tmpdf.modpep_z, y)) = {y};
            
            % Create boxplot
            clf(fig);
            boxplot(subtmpdf.ms1adj, {subtmpdf.Condition, subtmpdf.allothers}, ...
                'ColorGroup', subtmpdf.allothers, ...
                'Colors', [0.7 0.7 0.7; 0 0.73 0.22], ...
                'Symbol', '');
            
            % Add individual points
            hold on;
            conditions = unique(subtmpdf.Condition);
            for k = 1:length(conditions)
                cond_idx = strcmp(subtmpdf.Condition, conditions{k});
                scatter(repmat(k, sum(cond_idx), 1) + ...
                    (strcmp(subtmpdf.allothers(cond_idx), y) - 0.5) * 0.3, ...
                    subtmpdf.ms1adj(cond_idx), ...
                    'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
            end
            
            % Customize plot
            title(x, 'Interpreter', 'none');
            xlabel('Group');
            ylabel('Log2(intensity)-Log2(control)');
            set(gca, 'FontSize', 12);
            
            % Add legend
            legend('Other peptides', 'Discordant peptide', 'Location', 'southoutside');
            
            % Adjust layout
            set(gcf, 'Color', 'white');
            set(gca, 'Box', 'on');
            
            % Save figure if needed
            % saveas(fig, sprintf('%s_%s.png', x, y));
        end
    end
end

end 