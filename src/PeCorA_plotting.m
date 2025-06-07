function fig = PeCorA_plotting(w, u, v)
%PECORA_PLOTTING Generate plots for PeCorA results
%   w: Table of disagree peptides
%   u: Selected row from w
%   v: Scaled peptides table

% Input validation
assert(istable(w), 'First input must be a table');
assert(istable(u), 'Second input must be a table');
assert(istable(v), 'Third input must be a table');

% Create figure
fig = figure('Position', [100, 100, 1200, 800], 'Color', 'white');

% Get significant proteins
sig_proteins = unique(w.protein(w.adj_pval < 0.05));

% Process each significant protein
for i = 1:length(sig_proteins)
    protein = sig_proteins{i};
    
    % Get peptides for this protein
    protein_peptides = w(strcmp(w.protein, protein), :);
    
    % Create subplots
    subplot(2, 3, 1); % Peptide Abundance Boxplot
    hold on;
    set(gca, 'Color', 'white');  % Set white background
    % Convert Condition to categorical for plotting
    conditions = categorical(v.Condition(strcmp(v.Protein, protein)));
    boxplot(v.ms1adj(strcmp(v.Protein, protein)), conditions);
    scatter(conditions, v.ms1adj(strcmp(v.Protein, protein)), 'filled');
    title('Peptide Abundance Boxplot');
    xlabel('Condition');
    ylabel('MS1 Adjusted Abundance');
    legend('Boxplot', 'Individual Measurements');
    hold off;
    
    % Statistical Summary
    subplot(2, 3, 2);
    hold on;
    set(gca, 'Color', 'white');  % Set white background
    text(0.1, 0.9, sprintf('Protein: %s', protein), 'FontSize', 12);
    text(0.1, 0.8, sprintf('Number of Peptides: %d', height(protein_peptides)), 'FontSize', 12);
    text(0.1, 0.7, sprintf('Significant Peptides: %d', sum(protein_peptides.adj_pval < 0.05)), 'FontSize', 12);
    text(0.1, 0.6, sprintf('Average Log2FC: %.2f', mean(protein_peptides.log2FC)), 'FontSize', 12);
    text(0.1, 0.5, sprintf('Min p-value: %.2e', min(protein_peptides.pval)), 'FontSize', 12);
    axis off;
    title('Statistical Summary');
    hold off;
    
    % Protein Level Changes
    subplot(2, 3, 3);
    hold on;
    set(gca, 'Color', 'white');  % Set white background
    boxplot(v.ms1adj(strcmp(v.Protein, protein)), conditions);
    plot(1:length(unique(conditions)), mean(v.ms1adj(strcmp(v.Protein, protein))), 'r-', 'LineWidth', 2);
    title('Protein Level Changes');
    xlabel('Condition');
    ylabel('MS1 Adjusted Abundance');
    legend('Boxplot', 'Mean');
    hold off;
    
    % Peptide-Protein Relationship
    subplot(2, 3, 4);
    hold on;
    set(gca, 'Color', 'white');  % Set white background
    scatter(v.ms1adj(strcmp(v.Protein, protein)), v.ms1adj(strcmp(v.Protein, protein)), 'filled');
    xlabel('Peptide Abundance');
    ylabel('Protein Abundance');
    title('Peptide-Protein Relationship');
    % Add regression line
    p = polyfit(v.ms1adj(strcmp(v.Protein, protein)), v.ms1adj(strcmp(v.Protein, protein)), 1);
    x = linspace(min(v.ms1adj(strcmp(v.Protein, protein))), max(v.ms1adj(strcmp(v.Protein, protein))), 100);
    plot(x, polyval(p, x), 'r-', 'LineWidth', 2);
    legend('Data Points', 'Regression Line');
    hold off;
    
    % Volcano Plot
    subplot(2, 3, 5);
    hold on;
    set(gca, 'Color', 'white');  % Set white background
    scatter(protein_peptides.log2FC, -log10(protein_peptides.adj_pval), 'filled');
    xline(0, 'k--');
    yline(-log10(0.05), 'k--');
    xlabel('Log2 Fold Change');
    ylabel('-log10(Adjusted p-value)');
    title('Volcano Plot');
    % Add significance labels
    sig_idx = protein_peptides.adj_pval < 0.05;
    if any(sig_idx)
        text(protein_peptides.log2FC(sig_idx), -log10(protein_peptides.adj_pval(sig_idx)), ...
            protein_peptides.peptide(sig_idx), 'FontSize', 8);
    end
    hold off;
    
    % Add main title
    sgtitle(sprintf('PeCorA Analysis Results for %s', protein), 'FontSize', 14);
end

% Ensure white background is preserved when saving
set(gcf, 'InvertHardcopy', 'off');

end 