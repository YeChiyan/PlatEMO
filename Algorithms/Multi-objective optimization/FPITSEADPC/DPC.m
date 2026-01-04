function cluster = DPC(Pop, fi)
% Standard DPC: No noise removal, assign every point.
% Input: Pop, fi
% Output: cluster (cell array)

    Data = Pop.decs;
    [N, ~] = size(Data);
    
    %% 1. 计算距离矩阵
    Dist = pdist2(Data, Data);
    
    %% 2. 计算截断距离 dc
    if N < 2
        cluster = {1:N};
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
    
    %% 5. 确定聚类中心 (自动筛选)
    Gamma = Rho .* Delta;
    threshold = mean(Gamma) + 2.0 * std(Gamma);
    center_indices = find(Gamma > threshold);
    
    if isempty(center_indices)
        [~, max_idx] = max(Gamma);
        center_indices = max_idx;
    end
    
    %% 6. 分配剩余点 (确保所有人都有归属)
    labels = zeros(N, 1);
    num_clusters = length(center_indices);
    
    for k = 1:num_clusters
        labels(center_indices(k)) = k;
    end
    
    for i = 1:N
        current_idx = ord(i);
        if labels(current_idx) == 0
            higher_density_indices = ord(1:i-1);
            % 找密度比自己高且最近的点
            [~, nearest_idx_in_higher] = min(Dist(current_idx, higher_density_indices));
            nearest_neighbor = higher_density_indices(nearest_idx_in_higher);
            labels(current_idx) = labels(nearest_neighbor);
        end
    end
    
    %% 7. 输出
    cluster = cell(1, num_clusters);
    for k = 1:num_clusters
        cluster{k} = find(labels == k);
    end
end