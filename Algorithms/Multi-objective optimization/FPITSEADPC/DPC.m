function [cluster, noise_idx] = DPC(Pop, fi)
% Density Peak Clustering with Noise Detection
% Input: Pop, fi
% Output: 
%   cluster: valid clusters (cell array)
%   noise_idx: indices of noise points (vector)

    Data = Pop.decs;
    [N, ~] = size(Data);
    
    %% 1. 计算距离矩阵
    Dist = pdist2(Data, Data);
    
    %% 2. 计算截断距离 dc (自适应)
    if N < 2
        cluster = {1:N};
        noise_idx = [];
        return;
    end
    
    all_dists = sort(Dist(:));
    all_dists = all_dists(all_dists > 0);
    
    if isempty(all_dists)
        dc = 1e-6;
    else
        percentile_idx = max(1, round(length(all_dists) * fi));
        dc = all_dists(percentile_idx);
    end
    
    %% 3. 计算局部密度 Rho
    Rho = sum(exp(-(Dist./dc).^2), 2) - 1;
    
    %% 4. 计算相对距离 Delta
    Delta = zeros(N, 1);
    [~, ord] = sort(Rho, 'descend');
    
    Delta(ord(1)) = max(Dist(ord(1), :));
    
    for i = 2:N
        current_idx = ord(i);
        higher_density_indices = ord(1:i-1);
        Delta(current_idx) = min(Dist(current_idx, higher_density_indices));
    end
    
    %% 5. 确定聚类中心
    Gamma = Rho .* Delta;
    threshold = mean(Gamma) + 2.0 * std(Gamma);
    center_indices = find(Gamma > threshold);
    
    if isempty(center_indices)
        [~, max_idx] = max(Gamma);
        center_indices = max_idx;
    end
    
    %% 6. 分配剩余点 (增加噪声过滤机制)
    labels = zeros(N, 1);
    num_clusters = length(center_indices);
    
    % [新增] 定义噪声阈值：平均密度的 10%
    % 只有非中心点才会被判定为噪声
    noise_rho_threshold = mean(Rho) * 0.1; 
    
    % 先给中心点打标签
    for k = 1:num_clusters
        labels(center_indices(k)) = k;
    end
    
    % 分配其余点
    for i = 1:N
        current_idx = ord(i);
        
        % 如果已经作为中心处理过，跳过
        if labels(current_idx) > 0
            continue; 
        end
        
        % [核心修改] 噪声检测
        if Rho(current_idx) < noise_rho_threshold
            labels(current_idx) = -1; % 标记为噪声
            continue;
        end
        
        % 正常的归属逻辑：找密度比自己高的最近邻居
        higher_density_indices = ord(1:i-1);
        % 注意：这里要在“非噪声”的高密度点里找，或者直接找所有高密度点
        % 为了简化，我们找所有高密度点，哪怕它是噪声，也会因为传递性最终指向中心
        % 但为了稳健，最好找已分类的点
        valid_higher = higher_density_indices(labels(higher_density_indices) > 0);
        
        if isempty(valid_higher)
            % 极罕见情况：比它密度高的都是噪声，那它也暂定为噪声
            labels(current_idx) = -1;
        else
            [~, nearest_idx_in_higher] = min(Dist(current_idx, valid_higher));
            nearest_neighbor = valid_higher(nearest_idx_in_higher);
            labels(current_idx) = labels(nearest_neighbor);
        end
    end
    
    %% 7. 格式化输出
    % 输出有效聚类
    cluster = cell(1, num_clusters);
    for k = 1:num_clusters
        cluster{k} = find(labels == k);
    end
    
    % 输出噪声索引
    noise_idx = find(labels == -1);
end