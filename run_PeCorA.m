% PeCorA Main Script
% This script demonstrates how to use the PeCorA MATLAB implementation

% Add src directory to MATLAB path
addpath('src');

try
    t = readtable('data/PeCorA_noZ.csv');
catch
    error('请确保数据文件存在于data目录下，并且文件名正确');
end

% 2. 数据预处理
fprintf('开始数据预处理...\n');
scaled_peptides = PeCorA_preprocessing(t, ...
    8, ...    % 峰面积数据所在的列号
    100, ...  % 过滤阈值
    'cntrl'); % 对照组名称
fprintf('数据预处理完成\n');

% 3. 运行PeCorA分析
fprintf('开始PeCorA分析...\n');
disagree_peptides = PeCorA(scaled_peptides);
fprintf('PeCorA分析完成\n');

% 4. 保存结果
% 创建results目录（如果不存在）
if ~exist('results', 'dir')
    mkdir('results');
end

% 保存分析结果
writetable(disagree_peptides, 'results/PeCorA_results.csv');

% 5. 可视化结果
fprintf('生成可视化结果...\n');
% 选择前5个显著差异的肽段进行可视化
num_plots = min(5, height(disagree_peptides));
for i = 1:num_plots
    fig = PeCorA_plotting(disagree_peptides, disagree_peptides(i,:), scaled_peptides);
    % 保存图形
    saveas(fig, sprintf('results/PeCorA_plot_%d.png', i));
end
fprintf('分析完成！结果已保存到results目录\n');

% 显示结果摘要
fprintf('\n结果摘要:\n');
fprintf('总分析肽段数: %d\n', height(disagree_peptides));
fprintf('显著差异肽段数 (p < 0.01): %d\n', sum(disagree_peptides.adj_pval <= 0.01));
fprintf('显著差异蛋白质数: %d\n', length(unique(disagree_peptides.protein(disagree_peptides.adj_pval <= 0.01)))); 