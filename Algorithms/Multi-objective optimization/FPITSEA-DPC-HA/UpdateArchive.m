function Archive = UpdateArchive(Archive, NewPopulation, N)
% Update the external archive
% Logic: Merge -> Non-dominated Sort -> IDSS Pruning

%% 1. 合并
% 将旧档案和新种群合并
Union = [Archive, NewPopulation];

%% 2. 简单去重 (可选，防止完全一样的解占位)
% 这里利用 PlatEMO 的 SOLUTION 类特性，暂时不做硬性去重，依赖 IDSS 筛选

%% 3. 非支配排序预筛选
% 我们希望档案里存的尽量是收敛性好的解
[FrontNo, MaxFNo] = NDSort(Union.objs, N);

% 保留前几层的解 (不仅仅是第一层，防止漏掉稍弱的局部最优)
Next = FrontNo <= MaxFNo;
Candidate = Union(Next);

%% 4. 使用 IDSS 维护档案规模
% 如果候选解超过 N，用 IDSS 挑选分布最好的 N 个
% 这一步保证了档案既收敛，又在决策空间分布均匀
if length(Candidate) > N
    Archive = IDSS(Candidate, N);
else
    Archive = Candidate;
end
end
