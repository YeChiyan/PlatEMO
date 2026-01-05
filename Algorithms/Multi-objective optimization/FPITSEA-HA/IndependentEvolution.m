function Population = IndependentEvolution(Problem,Population,Archive,N,fi)
% Independent evolution strategy
% Strategy: Post-Clustering Injection & Lost Peak Reactivation

%% 1. Use DPC to divide the Population (ONLY)
% 这里的 Population 是纯粹的当前代，不含档案，保证拓扑结构准确
CD = Crowding(Population.decs);

% 获取聚类结果 和 截断距离 dc
[HD, dc] = DPC(Population, fi);

% 构建初始子种群
lpop = cell(1,length(HD));
lCD = cell(1,length(HD));

% 为了快速查找，建立一个映射：IndIndex -> ClusterID
pop_labels = zeros(1, length(Population));
for i = 1:length(HD)
    lpop{i} = Population(HD{i});
    lCD{i} = CD(HD{i});
    pop_labels(HD{i}) = i;
end

%% 2. [核心改进] Archive Injection logic
% 将档案中的解分配到最近的子种群，或者创建新种群（找回丢失的山峰）

if ~isempty(Archive)
    % 计算档案到当前种群的距离矩阵
    % Rows: Archive individuals, Cols: Population individuals
    DistMatrix = pdist2(Archive.decs, Population.decs);
    
    for i = 1:length(Archive)
        arc_sol = Archive(i);
        
        % 找到离这个档案解最近的 Population 个体
        [min_dist, nearest_idx] = min(DistMatrix(i, :));
        
        % 判断逻辑：
        if min_dist < dc
            % === 情况 A：锦上添花 ===
            % 档案解离当前某个子种群很近，加入它，增强优势
            target_cluster_id = pop_labels(nearest_idx);
            
            % 防御性编程：确保 ID 有效
            if target_cluster_id > 0 && target_cluster_id <= length(lpop)
                lpop{target_cluster_id} = [lpop{target_cluster_id}, arc_sol];
            end
        else
            % === 情况 B：雪中送炭 (Lost Peak Reactivation) ===
            % 档案解离所有当前个体都很远 -> 这是一个被当前种群遗忘的山峰！
            % 必须把它加回来，作为一个独立的子种群
            
            % 创建一个新的子种群
            lpop{end+1} = arc_sol;
        end
    end
end

% 此时 lpop 中包含了：
% 1. 原有的 DPC 子种群 (可能被档案加强了)
% 2. 新增的单一个体子种群 (来自档案的遗失山峰)

%% 3. Reproduction within the subpopulation
% 重新计算 HD 的长度，因为可能增加了新种群
n_clusters = length(lpop);

for i = 1:n_clusters
    subpop = lpop{i};
    
    % 由于注入了档案，我们需要重新计算一下拥挤度，或者简单处理
    n_sub = length(subpop);
    
    if n_sub >= 2
        % === 正常种群 ===
        % 需要计算当前 subpop 的拥挤度用于选择
        sub_decs = subpop.decs;
        if n_sub > 2
            current_subCD = CalCD_Internal(sub_decs);
        else
            current_subCD = ones(n_sub, 1);
        end
        
        mateidx = TournamentSelection(2, n_sub, -current_subCD);
        Off = OperatorGA(Problem, subpop(mateidx));
        lpop{i} = [lpop{i} Off];
        
    elseif n_sub == 1
        % === 孤儿种群 (通常是刚被档案找回的 Lost Peak) ===
        % 保护机制：原地变异
        Parent = subpop(1);
        Off = OperatorGA(Problem, [Parent, Parent]);
        lpop{i} = [lpop{i} Off(1)];
    end
end

%% 4. Environmental selection
Population = EnvironmentalSelectionS2(lpop,N);
end

% --- 辅助函数 ---
function CD = CalCD_Internal(PopDecs)
% 简单的内部拥挤度计算，避免调用外部依赖
N = size(PopDecs, 1);
if N <= 1, CD = 1; return; end
Dist = pdist2(PopDecs, PopDecs);
Dist(logical(eye(N))) = inf;
Dist = sort(Dist, 2);
% 取最近邻距离作为拥挤度估计
CD = Dist(:, 1);
end

function Pop = EnvironmentalSelectionS2(lpop,N)
% Environmental Selection S2: Global Competition Strategy

%% 1. 汇聚所有解
AllPop = [];
for i = 1:length(lpop)
    AllPop = [AllPop, lpop{i}];
end

if isempty(AllPop)
    Pop = [];
    return;
end

%% 2. 全局非支配排序
[FrontNo, MaxFNo] = NDSort(AllPop.objs, N);

% 只保留前几层的优良解
Next = FrontNo <= MaxFNo;
CandidatePop = AllPop(Next);

%% 3. 基于 IDSS 的最终筛选
if length(CandidatePop) <= N
    Pop = CandidatePop;
else
    Pop = IDSS(CandidatePop, N);
end
end

function ChoosePop = IDSS(Population,N)
% Improved distance-based subset selection method (IDSS)

%% 初始化
pop = [];
decdist = pdist2(Population.decs,Population.decs,'euclidean','Smallest',2);
objdist = pdist2(Population.objs,Population.objs,'euclidean','Smallest',2);

decdist = decdist(2,:);
objdist = objdist(2,:);

if max(decdist) == min(decdist)
    nordecdist = decdist;
else
    nordecdist = (decdist-min(decdist))./(max(decdist)-min(decdist));
end

if max(objdist) == min(objdist)
    norobjdist = objdist;
else
    norobjdist = (objdist-min(objdist))./(max(objdist)-min(objdist));
end

distsum = nordecdist+norobjdist;
[~,midx] = max(distsum);
pop = [pop Population(midx)];

retainpop = Population;
retainpop(midx) = [];

%% 迭代选择
while length(pop) < N && ~isempty(retainpop)
    decdist = pdist2(pop.decs,retainpop.decs,'euclidean','Smallest',1);
    objdist = pdist2(pop.objs,retainpop.objs,'euclidean','Smallest',1);
    
    if max(decdist) == min(decdist)
        nordecdist = decdist;
    else
        nordecdist = (decdist-min(decdist))./(max(decdist)-min(decdist));
    end
    
    if max(objdist) == min(objdist)
        norobjdist = objdist;
    else
        norobjdist = (objdist-min(objdist))./(max(objdist)-min(objdist));
    end
    
    distsum = nordecdist+norobjdist;
    [~,midx] = max(distsum);
    
    pop = [pop retainpop(midx)];
    retainpop(midx) = [];
end

ChoosePop = pop;
end