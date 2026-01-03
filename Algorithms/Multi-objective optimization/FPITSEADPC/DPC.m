function cluster = DPC(Pop, fi)
% Density Peak Clustering (DPC) for FPITSEA
% Input: 
%   Pop: Population object (PlatEMO)
%   fi:  Parameter for cutoff distance percentile (e.g., 0.1 for 10%)
% Output:
%   cluster: Cell array of indices for each sub-population

    Data = Pop.decs;
    [N, ~] = size(Data);
    
    %% 1. 计算距离矩阵 (欧氏距离)
    Dist = pdist2(Data, Data);
    
    %% 2. 计算截断距离 dc (Cutoff Distance)
    % 简单的自适应策略：取所有距离中第 fi% 小的值作为 dc
    % 这种方式比固定数值更适应不同的测试问题尺度
    if N < 2
        cluster = {1:N};
        return;
    end
    
    % 将距离矩阵铺平并排序
    all_dists = sort(Dist(:));
    % 排除0（自己到自己）
    all_dists = all_dists(all_dists > 0);
    
    if isempty(all_dists)
        dc = 1e-6;
    else
        % 计算位置索引
        percentile_idx = max(1, round(length(all_dists) * fi));
        dc = all_dists(percentile_idx);
    end
    
    %% 3. 计算局部密度 Rho (使用高斯核)
    % Rho_i = sum( exp( -(d_ij / dc)^2 ) )
    % 减1是为了去掉点到自身的距离贡献 (exp(0)=1)
    Rho = sum(exp(-(Dist./dc).^2), 2) - 1;
    
    %% 4. 计算相对距离 Delta
    % Delta_i: 距离“比我密度高的点”中最近的那个点的距离
    Delta = zeros(N, 1);
    [~, ord] = sort(Rho, 'descend'); % 密度从大到小排序的索引
    
    % 处理密度最大的点 (它没有比它密度更高的点，Delta 取离它最远的点的距离)
    Delta(ord(1)) = max(Dist(ord(1), :));
    
    % 处理其余点
    for i = 2:N
        current_idx = ord(i);
        % 找到比当前点密度高的所有点的索引
        higher_density_indices = ord(1:i-1);
        % 在这些点中找距离最近的
        Delta(current_idx) = min(Dist(current_idx, higher_density_indices));
    end
    
    %% 5. 确定聚类中心 (Automatic Center Selection)
    % 计算决策值 Gamma
    Gamma = Rho .* Delta;
    
    % 自动筛选策略：利用统计阈值 (Mean + 2*Std) 识别离群点作为中心
    % 这避免了人工观察决策图
    threshold = mean(Gamma) + 2.0 * std(Gamma);
    center_indices = find(Gamma > threshold);
    
    % 兜底策略：如果没选出任何中心（极少见），选 Gamma 最大的点
    if isempty(center_indices)
        [~, max_idx] = max(Gamma);
        center_indices = max_idx;
    end
    
    %% 6. 分配剩余点 (Assignment)
    labels = zeros(N, 1);
    num_clusters = length(center_indices);
    
    % 先给中心点打标签
    for k = 1:num_clusters
        labels(center_indices(k)) = k;
    end
    
    % 按照密度从大到小的顺序，将非中心点分配给其“最近的、密度更高的邻居”所在的类
    for i = 1:N
        current_idx = ord(i);
        if labels(current_idx) == 0
            higher_density_indices = ord(1:i-1);
            [~, relative_neighbor_idx] = min(Dist(current_idx, higher_density_indices));
            nearest_neighbor = higher_density_indices(relative_neighbor_idx);
            
            % 继承邻居的标签
            labels(current_idx) = labels(nearest_neighbor);
        end
    end
    
    %% 7. 格式化输出
    cluster = cell(1, num_clusters);
    for k = 1:num_clusters
        cluster{k} = find(labels == k);
    end
end