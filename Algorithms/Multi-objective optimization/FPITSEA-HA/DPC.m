function cluster = DPC(Pop, fi)
% Fitness-Aware DPC (结合了 NBC 思想的 DPC)
% Input: Pop, fi
% Output: cluster (cell array)

    Data = Pop.decs;
    Objs = Pop.objs;
    [N, ~] = size(Data);
    
    %% 1. 准备工作：计算距离与适应度
    Dist = pdist2(Data, Data);
    
    % [借鉴 NBC]: 计算适应度 (距离理想点的负距离)
    % 归一化目标空间
    f_min = min(Objs, [], 1);
    f_max = max(Objs, [], 1);
    Objs_norm = (Objs - f_min) ./ (f_max - f_min + 1e-6);
    % 理想点为原点 (0,0)
    Fitness = -sqrt(sum(Objs_norm.^2, 2)); % 值越大越好
    
    %% 2. 计算局部密度 Rho (保持物理聚集属性)
    % 使用 KNN 策略计算密度，增强鲁棒性
    K = floor(sqrt(N));
    sorted_dist = sort(Dist, 2);
    k_dist = sorted_dist(:, K+1); % 第 K 个邻居的距离
    % Rho = 1 ./ (k_dist + 1e-6); % 简单倒数密度
    Rho = exp(-k_dist); % 或者高斯核密度
    
    %% 3. [核心改进] 计算相对距离 Delta (引入 NBC 思想)
    % Delta 定义为：离"适应度比我好"的最近邻居的距离
    % 而不是原版的"密度比我高"
    
    Delta = zeros(N, 1);
    [~, fit_ord] = sort(Fitness, 'descend'); % 按适应度从好到差排序
    
    % 全局最优解 (Fitness 最高) 的 Delta 取最大距离
    Delta(fit_ord(1)) = max(max(Dist));
    
    % 对于其他人，找比自己 Fitness 高的最近邻居
    for i = 2:N
        current_idx = fit_ord(i);
        
        % 找到适应度比当前点高的所有点
        better_indices = fit_ord(1:i-1);
        
        % 在这些"强者"中，找一个物理距离最近的
        Delta(current_idx) = min(Dist(current_idx, better_indices));
    end
    
    %% 4. 确定聚类中心
    % 中心必须满足：1. 它是局部收敛的聚集区 (Rho大); 2. 它离其他更好的解很远 (Delta大)
    Gamma = Rho .* Delta;
    
    % 自动筛选阈值
    threshold = mean(Gamma) + 2.0 * std(Gamma);
    center_indices = find(Gamma > threshold);
    
    % 兜底
    if isempty(center_indices)
        [~, max_idx] = max(Gamma);
        center_indices = max_idx;
    end
    
    %% 5. 分配剩余点 (Assignment)
    % [借鉴 NBC]: 所有的点顺着"适应度梯度"归属
    % 也就是每个点归属到：离它最近的、且适应度比它好的那个邻居所在的类
    
    labels = zeros(N, 1);
    num_clusters = length(center_indices);
    
    for k = 1:num_clusters
        labels(center_indices(k)) = k;
    end
    
    % 按适应度从高到低遍历 (除了 Global Best，它已经是中心)
    for i = 2:N
        current_idx = fit_ord(i);
        
        % 如果它是中心，跳过
        if labels(current_idx) > 0
            continue;
        end
        
        % 找适应度比自己高、且距离最近的点 (Nearest-Better Neighbor)
        better_indices = fit_ord(1:i-1);
        [~, relative_idx] = min(Dist(current_idx, better_indices));
        nearest_better_neighbor = better_indices(relative_idx);
        
        % 继承大哥的标签
        labels(current_idx) = labels(nearest_better_neighbor);
    end
    
    %% 6. 输出
    cluster = cell(1, num_clusters);
    for k = 1:num_clusters
        cluster{k} = find(labels == k);
    end
end