function Next = ApplyClearing(Population, Fitness, Next, N, Sigma)
% Population: 当前候选种群
% Fitness: 适应度值 (越小越好)
% Next: 当前被选中的个体逻辑索引 (true/false)
% N: 目标保留数量
% Sigma: 小生境半径 (决策空间距离阈值)

% 1. 获取当前存活的个体索引
Idx = find(Next);
if isempty(Idx); return; end
Idx = Idx(:); % 强制为列向量

PopDec = Population(Idx).decs;
FitSubset = Fitness(Idx);

% 2. 对存活个体按适应度排序 (最好的排前面)
[~, SortIdx] = sort(FitSubset);
SortedOriginalIdx = Idx(SortIdx);
SortedOriginalIdx = SortedOriginalIdx(:); % 强制为列向量

% 3. 执行清除
Keep = true(size(SortedOriginalIdx)); % 临时标记谁在清除中存活

for i = 1 : length(SortedOriginalIdx)
    if Keep(i) % 如果当前这个“赢家”还没被清除
        WinnerDec = Population(SortedOriginalIdx(i)).decs;
        
        % 检查排在它后面的所有人
        for j = i+1 : length(SortedOriginalIdx)
            if Keep(j) % 如果后面这个人还活着
                % 计算距离
                dist = norm(WinnerDec - Population(SortedOriginalIdx(j)).decs);
                
                % 如果距离小于半径，说明在同一个峰，杀掉较差的那个 (j)
                if dist < Sigma
                    Keep(j) = false;
                end
            end
        end
    end
end

% 4. 更新 Next
% 保护：如果清除得太狠，导致数量小于 N，我们尝试保留至少 N 个（按进度排序）
% 或者遵循标准清除：严格执行。但为了稳健性，如果存活数量不足，我们可以补齐。
FinalIdx = SortedOriginalIdx(Keep);

if length(FinalIdx) < N
    % 如果清除后不足 N 个，补充原始 Next 中最好的（且不在 FinalIdx 中的）
    RemainingIdx = setdiff(SortedOriginalIdx, FinalIdx, 'stable');
    needed = N - length(FinalIdx);
    FinalIdx = [FinalIdx; RemainingIdx(1:min(needed, length(RemainingIdx)))];
end

% 重置 Next
Next(:) = false;
Next(FinalIdx) = true;
end
