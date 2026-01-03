%% 通用算法演化过程展示脚本 (图片与数据同步保存)
clear; clc;

%% 1. 参数设置区域
algName  = 'DNNSGAII'; % 算法名称
probName = 'MMF8';            % 目标函数名称
popSize  =400;               
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

%% 3. 获取问题边界
lb = []; ub = [];
try
    tempFunc = str2func(probName);
    proObject = tempFunc(); 
    lb = proObject.lower; ub = proObject.upper;
catch
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
    
    % --- 左：目标空间 ---
    subplot(1, 2, 1); cla;
    scatter(objs(:, 1), objs(:, 2), 15, 'filled', 'MarkerFaceColor', [0.1 0.6 0.3]);
    xlabel('f1'); ylabel('f2');
    title(sprintf('%s (PF) | 进度: %.0f%%', probName, progress*100));
    grid on; axis square;
    
    % --- 右：决策空间 ---
    subplot(1, 2, 2); cla;
    scatter(decs(:, 1), decs(:, 2), 15, 'filled', 'MarkerFaceColor', [0.6 0.3 0.1]);
    title(sprintf('%s (PS) | FE: %d', probName, currentFE));
    
    if ~isempty(lb) && length(lb) >= 2
        xlim([lb(1), ub(1)]); ylim([lb(2), ub(2)]); 
    else
        axis tight;
    end
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