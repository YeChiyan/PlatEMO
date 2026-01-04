%% FPITSEA 算法演化过程展示脚本
clear; clc;

%% 1. 参数设置区域
algName  = 'FPITSEADPC'; % 算法名称 (指向 Algorithms 文件夹中的类)
probName = 'MMF8';    % 目标函数名称
popSize  = 200;
maxFE    = popSize * 200;
savePts  = 20;

% 动画速度
baseDelay = 0.2;
slowDelay = 0.5;

%% 2. 运行算法并检索最新数据路径
fprintf('正在使用 %s 运行 %s...\n', algName, probName);
% 确保路径已添加
addpath(genpath(fileparts(mfilename('fullpath'))));

platemo('algorithm', str2func(algName), 'problem', str2func(probName), ...
    'N', popSize, 'maxFE', maxFE, 'save', savePts);

% 自动寻找最新生成的 .mat 文件
searchPath = fullfile('Data', algName);
filePattern = sprintf('%s_%s_*.mat', algName, probName);
files = dir(fullfile(searchPath, filePattern));
if isempty(files); error('未找到结果文件。'); end
[~, idx] = sort([files.datenum], 'descend');

dataFileName = files(idx(1)).name;
dataFolder   = files(idx(1)).folder;
fullDataPath = fullfile(dataFolder, dataFileName);
fprintf('加载数据：%s\n', dataFileName);
load(fullDataPath);

%% 3. 获取问题对象与参考数据 (PS 和 PF)
proObject = [];
try
    PRO_func = str2func(probName);
    proObject = PRO_func('N', popSize, 'maxFE', maxFE);
    proObject.GetOptimum(1000);
catch
    fprintf('警告: 无法获取参考数据。\n');
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
fig = figure('Color', 'w', 'Name', ['Testing ' algName], 'Position', [100 100 1100 450]);

numStages = size(result, 1);
for i = 1 : numStages
    currentFE = result{i, 1};
    popObj    = result{i, 2};
    progress  = i / numStages;
    
    decs = popObj.decs;
    
    objs = popObj.objs;
    
    % --- 计算指标 (IGD & IGDX) ---
    currentIGD  = NaN;
    currentIGDX = NaN;
    if ~isempty(proObject)
        % 计算 IGD (目标空间)
        PF_ref = [];
        if ~isempty(proObject.PF); PF_ref = proObject.PF;
        elseif ~isempty(proObject.optimum); PF_ref = proObject.optimum; end
        
        if ~isempty(PF_ref) && ismatrix(PF_ref) && size(PF_ref, 2) == size(objs, 2)
            dists = pdist2(PF_ref, objs);
            currentIGD = mean(min(dists, [], 2));
        end
        
        % 计算 IGDX (决策空间)
        POS_ref = [];
        if isprop(proObject, 'POS') && ~isempty(proObject.POS); POS_ref = proObject.POS; end
        
        if ~isempty(POS_ref) && ismatrix(POS_ref) && size(POS_ref, 2) == size(decs, 2)
            distsX = pdist2(POS_ref, decs);
            currentIGDX = mean(min(distsX, [], 2));
        end
    end
    
    % --- 左：目标空间 (PF) ---
    subplot(1, 2, 1); cla; hold on;
    % 绘制参考 PF (蓝色空心点 - True PF)
    if ~isempty(proObject)
        if ~isempty(proObject.PF) && ismatrix(proObject.PF) && size(proObject.PF, 2) == 2
            scatter(proObject.PF(:, 1), proObject.PF(:, 2), 10, 'b');
        elseif ~isempty(proObject.optimum) && ismatrix(proObject.optimum) && size(proObject.optimum, 2) == 2
            scatter(proObject.optimum(:, 1), proObject.optimum(:, 2), 10, 'b');
        end
    end
    % 绘制当前种群 (红色空心点 - Obtained PF)
    scatter(objs(:, 1), objs(:, 2), 15, 'r');
    xlabel('f1'); ylabel('f2');
    titleStr = sprintf('%s (PF) | 进度: %.0f%%', probName, progress*100);
    if ~isnan(currentIGD); titleStr = sprintf('%s\nIGD: %.4f', titleStr, currentIGD); end
    title(titleStr);
    grid on; axis square; hold off;
    
    % --- 右：决策空间 (PS) ---
    subplot(1, 2, 2); cla; hold on;
    if ~isempty(proObject) && isprop(proObject, 'POS') && ~isempty(proObject.POS)
        % 绘制参考 PS (蓝色小点 - True PS)
        scatter(proObject.POS(:, 1), proObject.POS(:, 2), 2, 'b');
    end
    % 绘制当前种群 (红色空心点 - Obtained PS)
    scatter(decs(:, 1), decs(:, 2), 15, 'r');
    
    if ~isempty(lb) && length(lb) >= 2
        xlim([lb(1), ub(1)]); ylim([lb(2), ub(2)]);
    end
    xlabel('x1'); ylabel('x2');
    titleStrX = sprintf('%s (PS) | FE: %d', probName, currentFE);
    if ~isnan(currentIGDX); titleStrX = sprintf('%s\nIGDX: %.4f', titleStrX, currentIGDX); end
    title(titleStrX);
    grid on; axis square;
    drawnow;
    
    if progress > 0.7; pause(slowDelay); else; pause(baseDelay); end
end

%% 5. 保存图像
imgFileName = strrep(dataFileName, '.mat', '_Plot.png');
imgSavePath = fullfile(dataFolder, imgFileName);
saveas(gcf, imgSavePath);
fprintf('\n【保存成功】\n数据文件: %s\n图像文件: %s\n', fullDataPath, imgSavePath);
