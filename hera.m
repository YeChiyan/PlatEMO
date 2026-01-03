%% 通用算法演化过程展示脚本 (图片与数据同步保存)
clear; clc;

%% 1. 参数设置区域
algName  = 'HREA'; % 算法名称
probName = 'MMF6';            % 目标函数名称
popSize  =200;
maxFE    = popSize * 100;
savePts  = 20;

% 动画速度
baseDelay = 0.4;
slowDelay = 0.8;

%% 2. 运行算法并检索最新数据路径
fprintf('正在使用 %s 运行 %s...\n', algName, probName);
platemo('algorithm', str2func(algName), 'problem', str2func(probName), ...
    'N', popSize, 'maxFE', maxFE, 'save', savePts);

% 自动寻找最新生成的 .mat 文件
searchPath = fullfile('Data', algName);
filePattern = sprintf('%s_%s_*.mat', algName, probName);
files = dir(fullfile(searchPath, filePattern));
if isempty(files); error('未找到结果文件。'); end
[~, idx] = sort([files.datenum], 'descend');

% 记录数据文件的完整路径和文件夹路径
dataFileName = files(idx(1)).name;
dataFolder   = files(idx(1)).folder;
fullDataPath = fullfile(dataFolder, dataFileName);

fprintf('加载数据：%s\n', dataFileName);
load(fullDataPath);

%% 3. 获取问题对象与参考数据 (PS 和 PF)
% 确保路径已添加 (通过运行一次 platemo)
addpath(genpath(fileparts(mfilename('fullpath'))));

proObject = [];
try
    % 模仿 platemo.m 的实例化方式
    PRO_func = str2func(probName);
    proObject = PRO_func('N', popSize, 'maxFE', maxFE);
    % 确保 optimum 和 POS 属性已初始化
    proObject.GetOptimum(1000);
catch ME
    fprintf('警告: 无法完整实例化问题对象 (%s)，尝试备选方案。\n', ME.message);
    try
        proObject = feval(probName);
        proObject.GetOptimum(1000);
    catch
        fprintf('错误: 无法获取参考数据。\n');
    end
end

lb = []; ub = [];
if ~isempty(proObject)
    lb = proObject.lower; ub = proObject.upper;
else
    if ~isempty(result)
        sampleDecs = result{end, 2}.decs;
        lb = min(sampleDecs, [], 1); ub = max(sampleDecs, [], 1);
    end
end

%% 4. 动态展示演化动画
fig = figure('Color', 'w', 'Position', [100 100 1100 450]);

numStages = size(result, 1);
for i = 1 : numStages
    currentFE = result{i, 1};
    popObj    = result{i, 2};
    progress  = i / numStages;
    
    objs = popObj.objs;
    decs = popObj.decs;
    
    % --- 左：目标空间 (PF) ---
    subplot(1, 2, 1); cla; hold on;
    % 绘制参考 PF (黑点/黑线)
    if ~isempty(proObject)
        if ~isempty(proObject.PF) && ismatrix(proObject.PF) && size(proObject.PF, 2) == 2
            plot(proObject.PF(:, 1), proObject.PF(:, 2), 'k-', 'LineWidth', 1.5);
        elseif ~isempty(proObject.optimum) && ismatrix(proObject.optimum) && size(proObject.optimum, 2) == 2
            plot(proObject.optimum(:, 1), proObject.optimum(:, 2), 'k.', 'MarkerSize', 5);
        end
    end
    % 绘制当前种群 (统一绿色散点)
    scatter(objs(:, 1), objs(:, 2), 15, 'filled', 'MarkerFaceColor', [0.1 0.6 0.3]);
    xlabel('f1'); ylabel('f2');
    title(sprintf('%s (PF) | 进度: %.0f%%', probName, progress*100));
    grid on; axis square; hold off;
    
    % --- 右：决策空间 (PS) ---
    subplot(1, 2, 2); cla; hold on;
    if ~isempty(proObject) && isprop(proObject, 'POS') && ~isempty(proObject.POS)
        % 绘制参考 PS (Pareto Optimal Set)
        % MMF 问题通常在 POS 属性中存储参考点
        scatter(proObject.POS(:, 1), proObject.POS(:, 2), 2, [0.8 0.8 0.8], 'filled');
    end
    % 绘制当前种群
    scatter(decs(:, 1), decs(:, 2), 15, 'filled', 'MarkerFaceColor', [0.6 0.3 0.1]);
    
    if ~isempty(lb) && length(lb) >= 2
        xlim([lb(1), ub(1)]); ylim([lb(2), ub(2)]);
    else
        axis tight;
    end
    xlabel('x1'); ylabel('x2');
    title(sprintf('%s (PS) | FE: %d', probName, currentFE));
    grid on; axis square;
    drawnow;
    
    if progress > 0.7; pause(slowDelay); else; pause(baseDelay); end
end

%% 5. 保存图像至数据文件同级目录
% 将 .mat 后缀替换为 .png
imgFileName = strrep(dataFileName, '.mat', '_Plot.png');
imgSavePath = fullfile(dataFolder, imgFileName);

% 保存截图
saveas(gcf, imgSavePath);

fprintf('\n【保存成功】\n');
fprintf('数据文件: %s\n', fullDataPath);
fprintf('配套图像: %s\n', imgSavePath);