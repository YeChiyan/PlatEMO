function Population = IndependentEvolution(Problem,Population,N,fi)
% Independent evolution strategy
% Fix: Local Elitism + Pure IDSS (No Global NDSort)

    %% 1. Use DPC to divide the population
    CD = Crowding(Population.decs);
    % 建议 fi 设置为 0.05 以确保 MMF10 的上下两层能分开
    HD = DPC(Population, fi); 

    lpop = cell(1,length(HD));
    
    for i = 1:length(HD)
        lpop{i} = Population(HD{i});
    end

    %% 2. Reproduction (DE + GA)
    for i = 1:length(HD)
        subpop = lpop{i};
        n_sub = length(subpop);
        
        if n_sub >= 4
            % 大种群：DE 强力收敛
            % F=0.4 保持温和，防止过早收敛成单点
            Off = OperatorDE(Problem, subpop, {1, 0.4, 1, 20});
            lpop{i} = [lpop{i} Off];
            
        elseif n_sub >= 2
            % 中种群：GA
            subCD = CD(HD{i});
            mateidx = TournamentSelection(2, n_sub, -subCD);
            Off = OperatorGA(Problem, subpop(mateidx));
            lpop{i} = [lpop{i} Off];
            
        elseif n_sub == 1
            % 孤儿：变异
            Parent = subpop(1);
            Off = OperatorGA(Problem, [Parent, Parent]); 
            lpop{i} = [lpop{i} Off(1)];
        end
    end

    %% 3. Environmental selection
    Population = EnvironmentalSelectionS2(lpop, N);
end

%% ================== 核心修改：分封制选择 ==================

function Pop = EnvironmentalSelectionS2(lpop,N)
% Strategy: 
% 1. Local Filter: 每个山头自己清理门户 (Rank > 1 的内部垃圾去掉)
% 2. Global Merge: 精英会师
% 3. Direct IDSS:  不分高低贵贱，只看分布

    Candidates = [];
    
    %% 1. 组内筛选 (Local Elitism)
    for i = 1:length(lpop)
        cpop = lpop{i};
        if isempty(cpop), continue; end
        
        % [关键] 只在子种群内部排序！
        [FNo, ~] = NDSort(cpop.objs, length(cpop));
        
        % [策略] 只保留组内 Rank 1 的解 (绝对精英)
        % 这样能把那种"云团"状的宽带压成细线
        % 如果子种群还没收敛(Rank 1很少)，可以放宽到前 50%
        
        % 方案：保留组内 Rank 1。如果 Rank 1 太少(<2个)，则保留前 50%
        n_rank1 = sum(FNo == 1);
        
        if n_rank1 >= 2
            indices = (FNo == 1);
        else
            % 还没收敛好，保留前 50% 继续进化
            [~, sort_idx] = sort(FNo);
            n_keep = ceil(length(cpop) * 0.5);
            indices = false(1, length(cpop));
            indices(sort_idx(1:n_keep)) = true;
        end
        
        Candidates = [Candidates, cpop(indices)];
    end

    %% 2. IDSS 最终裁决 (Direct IDSS)
    % 此时 Candidates 里混合了：
    % - 全局 PS (Global Rank 1)
    % - 局部 PS (Global Rank 2, 但它是 Local Rank 1)
    
    if isempty(Candidates)
        Pop = [];
        return;
    end

    % [绝对禁止] 在这里做 Global NDSort！
    % 直接进 IDSS
    
    if length(Candidates) <= N
        Pop = Candidates;
    else
        % IDSS 会看到局部 PS 离全局 PS 决策距离很远，因此保留它
        Pop = IDSS(Candidates, N);
    end
end

% ... (请保留 IDSS 和 Crowding 函数在文件底部) ...
function ChoosePop = IDSS(Population,N)
    pop = [];
    decdist = pdist2(Population.decs,Population.decs,'euclidean','Smallest',2);
    objdist = pdist2(Population.objs,Population.objs,'euclidean','Smallest',2);
    decdist = decdist(2,:); objdist = objdist(2,:);
    if max(decdist) == min(decdist), nordecdist = decdist; else, nordecdist = (decdist-min(decdist))./(max(decdist)-min(decdist)); end
    if max(objdist) == min(objdist), norobjdist = objdist; else, norobjdist = (objdist-min(objdist))./(max(objdist)-min(objdist)); end
    distsum = nordecdist+norobjdist;
    [~,midx] = max(distsum);
    pop = [pop Population(midx)];
    retainpop = Population; retainpop(midx) = []; 
    while length(pop) < N && ~isempty(retainpop)
        decdist = pdist2(pop.decs,retainpop.decs,'euclidean','Smallest',1);
        objdist = pdist2(pop.objs,retainpop.objs,'euclidean','Smallest',1);
        if max(decdist) == min(decdist), nordecdist = decdist; else, nordecdist = (decdist-min(decdist))./(max(decdist)-min(decdist)); end
        if max(objdist) == min(objdist), norobjdist = objdist; else, norobjdist = (objdist-min(objdist))./(max(objdist)-min(objdist)); end
        distsum = nordecdist+norobjdist;
        [~,midx] = max(distsum);
        pop = [pop retainpop(midx)];
        retainpop(midx) = [];
    end
    ChoosePop = pop;
end

function Crowd = Crowding(Pop)
    [N,~] = size(Pop);
    if N == 0, Crowd=[]; return; end
    if N == 1, Crowd=inf; return; end
    K = N-1;
    Z = min(Pop,[],1);
    Zmax = max(Pop,[],1);
    diff = Zmax - Z; diff(diff==0) = 1;
    Norpop = (Pop-repmat(Z,N,1))./repmat(diff,N,1);
    distance = pdist2(Norpop,Norpop);
    distance(logical(eye(N))) = inf;
    distance = sort(distance,2);
    Crowd = K./sum(1./distance(:,1:min(N-1, 2)), 2); 
end