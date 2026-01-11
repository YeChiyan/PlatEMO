%% DCMMEA-HA 历史档案 算法演化过程展示脚本
clear; clc;
%% 含义
%% 1. 参数设置区域
algName  = 'DCMMEAHA'; % 算法名称 (指向 Algorithms 文件夹中的类)
probName = 'CMMF11';    % 目标函数名称
popSize  = 200;
iters =100;
maxFE    = popSize * iters;
savePts  = 20;

str='引入档案机制,引入 niche mating小生境';


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

% 提取文件名中的序号 (例如从 ..._8.mat 提取 8)
[~, nameOnly, ~] = fileparts(dataFileName);
nameTokens = split(nameOnly, '_');
fileIdxStr = nameTokens{end};

fprintf('加载数据：%s (序号: %s)\n', dataFileName, fileIdxStr);
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
%% 4. 动态展示演化动画 (四宫格：PF, PS, IGD, IGDX)
fig = figure('Color', 'w', 'Name', ['Testing ' algName], 'Position', [100 100 1000 800]);

numStages = size(result, 1);
igdHistory  = zeros(1, numStages);
igdxHistory = zeros(1, numStages);
feHistory   = zeros(1, numStages);
iterHistory = zeros(1, numStages);

for i = 1 : numStages    % 获取当前代数据
    currentFE = result{i, 1};
    objs      = result{i, 2}.objs;
    decs      = result{i, 2}.decs;
    progress  = currentFE / maxFE;
    feHistory(i) = currentFE;
    
    % --- 计算指标 ---
    currentIGD = NaN; currentIGDX = NaN;
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
    
    igdHistory(i)  = currentIGD;
    igdxHistory(i) = currentIGDX;
    iterHistory(i) = currentFE / popSize;
    
    clf(fig);
    % --- 顶部信息行 ---
    currentTime = datestr(now, 'mm-dd-HHMM');
    headerStr = {sprintf('Algorithm: %s | Problem: %s | PopSize: %d | MaxIter: %d | Time: %s', ...
        algName, probName, popSize, iters, currentTime), str};
    annotation('textbox', [0.1, 0.92, 0.8, 0.07], 'String', headerStr, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    
    % --- (1) 左上：目标空间 (PF) ---
    ax1 = subplot(2, 2, 1);
    ax1.Position(2) = ax1.Position(2) - 0.05; % 下移一点给 header 留空间
    hold(ax1, 'on');
    
    % 绘制参考 PF (蓝色小点)
    h_pf_true = [];
    if ~isempty(proObject)
        PF_ref = [];
        if ~isempty(proObject.PF); PF_ref = proObject.PF;
        elseif ~isempty(proObject.optimum); PF_ref = proObject.optimum; end
        if ~isempty(PF_ref)
            if size(PF_ref, 2) == 2
                h_pf_true = scatter(ax1, PF_ref(:, 1), PF_ref(:, 2), 5, [0.4 0.7 1], '.');
            elseif size(PF_ref, 2) == 3
                h_pf_true = scatter3(ax1, PF_ref(:, 1), PF_ref(:, 2), PF_ref(:, 3), 5, [0.4 0.7 1], '.');
            end
        end
    end
    % 绘制获得 PF (红色空心圆)
    h_pf_obs = [];
    if size(objs, 2) == 2
        h_pf_obs = scatter(ax1, objs(:, 1), objs(:, 2), 20, 'r', 'o');
        xlabel(ax1, 'f_1'); ylabel(ax1, 'f_2');
    elseif size(objs, 2) == 3
        h_pf_obs = scatter3(ax1, objs(:, 1), objs(:, 2), objs(:, 3), 20, 'r', 'o');
        xlabel(ax1, 'f_1'); ylabel(ax1, 'f_2'); zlabel(ax1, 'f_3');
        view(ax1, [45, 20]);
    end
    grid(ax1, 'on'); axis(ax1, 'square');
    if ~isempty(h_pf_true) && ~isempty(h_pf_obs)
        legend(ax1, [h_pf_true, h_pf_obs], {'True PF', 'Obtained PF'}, 'Location', 'best');
    end
    title(ax1, sprintf('Objective Space (PF) | FE: %d', currentFE));
    
    % --- (2) 右上：决策空间 (PS) ---
    ax2 = subplot(2, 2, 2);
    ax2.Position(2) = ax2.Position(2) - 0.05;
    hold(ax2, 'on');
    
    % 绘制参考 PS (淡蓝色小点)
    h_ps_true = [];
    if ~isempty(proObject) && isprop(proObject, 'POS') && ~isempty(proObject.POS)
        if size(proObject.POS, 2) == 2
            h_ps_true = scatter(ax2, proObject.POS(:, 1), proObject.POS(:, 2), 2, [0.6 0.8 1], '.');
        elseif size(proObject.POS, 2) == 3
            h_ps_true = scatter3(ax2, proObject.POS(:, 1), proObject.POS(:, 2), proObject.POS(:, 3), 2, [0.6 0.8 1], '.');
        end
    end
    % 绘制获得 PS (红色空心圆)
    h_ps_obs = [];
    if size(decs, 2) == 2
        h_ps_obs = scatter(ax2, decs(:, 1), decs(:, 2), 15, 'r', 'o');
        xlabel(ax2, 'x_1'); ylabel(ax2, 'x_2');
    elseif size(decs, 2) == 3
        h_ps_obs = scatter3(ax2, decs(:, 1), decs(:, 2), decs(:, 3), 15, 'r', 'o');
        xlabel(ax2, 'x_1'); ylabel(ax2, 'x_2'); zlabel(ax2, 'x_3');
        view(ax2, [45, 20]);
    end
    grid(ax2, 'on'); axis(ax2, 'square');
    if ~isempty(h_ps_true) && ~isempty(h_ps_obs)
        legend(ax2, [h_ps_true, h_ps_obs], {'True PS', 'Obtained PS'}, 'Location', 'best');
    end
    title(ax2, 'Decision Space (PS)');
    
    % --- (3) 左下：IGD 收敛曲线 ---
    ax3 = subplot(2, 2, 3);
    plot(ax3, iterHistory(1:i), igdHistory(1:i), '-b.', 'LineWidth', 1.2);
    grid(ax3, 'on');
    xlabel(ax3, 'Iterations'); ylabel(ax3, 'IGD');
    title(ax3, sprintf('IGD: %.4f', currentIGD));
    
    % --- (4) 右下：IGDX 收敛曲线 ---
    ax4 = subplot(2, 2, 4);
    plot(ax4, iterHistory(1:i), igdxHistory(1:i), '-r.', 'LineWidth', 1.2);
    grid(ax4, 'on');
    xlabel(ax4, 'Iterations'); ylabel(ax4, 'IGDX');
    title(ax4, sprintf('IGDX: %.4f', currentIGDX));
    
    drawnow;
    if progress > 0.8; pause(slowDelay); else; pause(baseDelay); end
end

%% 5. 保存图像
% 动态获取 M 和 D (默认为 2)
M_val = 2; D_val = 2;
if ~isempty(proObject)
    M_val = proObject.M;
    D_val = proObject.D;
end


% 构造符合要求的文件名：算法名称小写_problem名称_M2_D2_pop 种群大小_迭代次数_序号
imgFileName = sprintf('%s_%s_M%d_D%d_pop%d_%d_%s.png', ...
    lower(algName), lower(probName), M_val, D_val, popSize, iters, fileIdxStr);

imgSavePath = fullfile(dataFolder, imgFileName);
saveas(gcf, imgSavePath);
fprintf('\n【可视化完成】\n数据: %s\n图像: %s\n', fullDataPath, imgSavePath);
if ~isnan(currentIGD); fprintf('Final IGD: %.4f\n', currentIGD); end
if ~isnan(currentIGDX); fprintf('Final IGDX: %.4f\n', currentIGDX); end