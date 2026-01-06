function Population = IndependentEvolution(Problem,Population,N,fi)
% Independent evolution strategy
% Improvement:
% 1. Small Cluster Protection (Reproduction phase)
% 2. Global Competition (Selection phase) - [New!]

%% 1. Use DPC to divide the population
CD = Crowding(Population.decs);

% 调用 DPC (请确保 DPC.m 是无噪声过滤的标准版)
HD = DPC(Population, fi);

lpop = cell(1,length(HD));
lCD = cell(1,length(HD));

for i = 1:length(HD)
    lpop{i} = Population(HD{i});
    lCD{i} = CD(HD{i});
end

%% 2. Reproduction within the subpopulation
for i = 1:length(HD)
    subpop = lpop{i};
    subCD = lCD{i};
    n_sub = length(subpop);
    
    if n_sub >= 2
        % === 正常大种群 ===
        % 内部锦标赛选择，交叉变异
        mateidx = TournamentSelection(2, n_sub, -subCD);
        Off = OperatorGA(Problem, subpop(mateidx));
        lpop{i} = [lpop{i} Off];
        
    elseif n_sub == 1
        % === 孤儿种群保护 ===
        % 原地变异，不进行交叉，防止被拉偏
        Parent = subpop(1);
        Off = OperatorGA(Problem, [Parent, Parent]); % 触发变异
        lpop{i} = [lpop{i} Off(1)];
    end
end

%% 3. Environmental selection
% 调用修改后的全局竞争选择函数
Population = EnvironmentalSelectionS2(lpop,N);
end

%% ================== 核心修改区域：全局竞争选择 ==================

function Pop = EnvironmentalSelectionS2(lpop,N)
% Modified: Local Elitism + Global IDSS
% 专门修复 MMF10-16 丢失局部 PS 的问题

ndpop = []; % 存放第一梯队（每个子种群的精英）
dpop = [];  % 存放备胎（每个子种群的普通成员）

%% 1. 组内独立排序 (Intra-Cluster Sorting)
% 核心：只要你是子种群的老大，哪怕被全局支配，也先把你保下来
for i = 1:length(lpop)
    cpop = lpop{i};
    if isempty(cpop)
        continue;
    end
    
    % [关键修改] 只在 cpop 内部进行排序，不看全局
    [FNo1, ~] = NDSort(cpop.objs, length(cpop));
    
    % 收集该子种群的 Rank 1
    local_elites = cpop(FNo1 == 1);
    ndpop = [ndpop, local_elites];
    
    % 收集该子种群的 Rank > 1
    others = cpop(FNo1 ~= 1);
    dpop = [dpop, others];
end

%% 2. 溢出处理与补充
% 情况 A：精英太多 (> N) -> 需要内卷
if length(ndpop) > N
    % 使用 IDSS 从精英中筛选。
    % IDSS 会倾向于保留决策空间距离远的点，从而保护 Local PS
    Pop = IDSS(ndpop, N);
    return;
end

% 情况 B：精英太少 (< N) -> 需要从备胎里补
% 此时，为了保证补充进来的解质量不要太差，我们在备胎池里做一次全局排序
if length(ndpop) < N && ~isempty(dpop)
    [FNo2, MaxFNo2] = NDSort(dpop.objs, length(dpop));
    k = 1;
    while length(ndpop) < N && k <= MaxFNo2
        % 这一步是从备胎里选，按全局 Rank 选比较安全
        candidates = dpop(FNo2 == k);
        
        % 如果这一层加进去不超标，全加
        if length(ndpop) + length(candidates) <= N
            ndpop = [ndpop, candidates];
        else
            % 如果超标，用 IDSS 选最好的几个填满
            needed = N - length(ndpop);
            selected = IDSS(candidates, needed);
            ndpop = [ndpop, selected];
            break;
        end
        k = k + 1;
    end
end

Pop = ndpop;
end