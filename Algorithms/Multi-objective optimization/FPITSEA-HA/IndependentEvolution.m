function Population = IndependentEvolution(Problem,Population,Archive,N,fi)
% Independent evolution strategy
% Improvement:
% 1. Small Cluster Protection (Reproduction phase)
% 2. Global Competition (Selection phase)
% 3. Pre-DPC Archive Injection (Active Archive) - [New!]

%% 0. Archive Injection Strategy (Active Injection)
% Logic: Merge current population with the archive and use NDSort + IDSS to select the best N individuals.
% This ensures that the DPC is performed on the best solutions found so far ("all-star team").

% Merge
Union = [Population, Archive];

% Selection (NDSort + IDSS)
[FrontNo, MaxFNo] = NDSort(Union.objs, N);
Next = FrontNo <= MaxFNo;
Candidate = Union(Next);

if length(Candidate) > N
    Population = IDSS(Candidate, N);
else
    Population = Candidate;
end

%% 1. Use DPC to divide the population
CD = Crowding(Population.decs);

% Call DPC
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
% Environmental Selection S2: Global Competition Strategy
% 汇聚所有子种群 -> 全局非支配排序 -> IDSS 筛选

%% 1. 汇聚所有解 (Merge)
AllPop = [];
for i = 1:length(lpop)
    AllPop = [AllPop, lpop{i}];
end

if isempty(AllPop)
    Pop = [];
    return;
end

%% 2. 全局非支配排序 (Global Non-dominated Sorting)
% 这一步是杀手锏。那些“半吊子”离群点在这里会被判定为高层级(Rank high)，
% 从而被 MaxFNo 挡在门外，直接淘汰。
[FrontNo, MaxFNo] = NDSort(AllPop.objs, N);

% 只保留前几层的优良解
Next = FrontNo <= MaxFNo;
CandidatePop = AllPop(Next);

%% 3. 基于 IDSS 的最终筛选 (IDSS Selection)
% 如果优良解数量 <= N，直接全收
if length(CandidatePop) <= N
    Pop = CandidatePop;
else
    % 如果优良解 > N，使用 IDSS 从中挑选分布最好的 N 个
    % IDSS 会平衡决策空间和目标空间的分布，确保局部最优和全局最优都能保留
    Pop = IDSS(CandidatePop, N);
end
end