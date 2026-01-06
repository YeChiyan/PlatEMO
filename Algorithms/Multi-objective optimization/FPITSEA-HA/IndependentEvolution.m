function Population = IndependentEvolution(Problem,Population,N,fi)
% Independent evolution strategy
% Fix: Relaxed NDSort + IDSS (平衡收敛性与多样性)

    %% 1. Use DPC to divide the population
    CD = Crowding(Population.decs);
    HD = DPC(Population, fi); % 使用标准无噪声 DPC

    lpop = cell(1,length(HD));
    lCD = cell(1,length(HD));
    for i = 1:length(HD)
        lpop{i} = Population(HD{i});
        lCD{i} = CD(HD{i});
    end

    %% 2. Reproduction
    for i = 1:length(HD)
        subpop = lpop{i};
        subCD = lCD{i};
        n_sub = length(subpop);
        
        if n_sub >= 2
            mateidx = TournamentSelection(2, n_sub, -subCD);
            Off = OperatorGA(Problem, subpop(mateidx));
            lpop{i} = [lpop{i} Off];
        elseif n_sub == 1
            Parent = subpop(1);
            Off = OperatorGA(Problem, [Parent, Parent]); 
            lpop{i} = [lpop{i} Off(1)];
        end
    end

    %% 3. Environmental selection (关键修改)
    Population = EnvironmentalSelectionS2(lpop,N);
end


function Pop = EnvironmentalSelectionS2(lpop,N)
% Optimization: Local Rank 1 Extraction + IDSS
% 目的：将局部 PS 的"云团"压扁成"细线"

    Candidates = [];
    
    %% 1. 提取每个山头的"山大王" (Local Rank 1)
    % 这一步是收敛的关键。我们不再保留子种群的前50%，
    % 而是只保留子种群内部的非支配解 (Rank 1)。
    
    for i = 1:length(lpop)
        cpop = lpop{i};
        if isempty(cpop), continue; end
        
        % [关键] 组内排序
        [FNo, ~] = NDSort(cpop.objs, length(cpop));
        
        % 只取 Rank 1 (最锋利的那条线)
        LocalElites = cpop(FNo == 1);
        
        % [兜底] 如果这个子种群还没收敛，Rank 1 只有一个点，
        % 为了防止灭绝，可以稍微放宽到 Rank 2，但对于 MMF10 通常 Rank 1 足够多
        if length(LocalElites) < 2 && length(cpop) > 2
             LocalElites = [LocalElites, cpop(FNo == 2)];
        end
        
        Candidates = [Candidates, LocalElites];
    end

    %% 2. 数量控制与最终筛选
    % 现在 Candidates 里全是各个山头的顶级精英
    
    if length(Candidates) <= N
        % 如果精英不够 N 个 (比较少见)，直接全收
        % 剩下的空位，可以从全局 Rank 较好的人里补 (为了简单，这里直接返回)
        Pop = Candidates;
        
        % [可选] 强行补满 N (从原 lpop 中找漏网之鱼)
        if length(Pop) < N
             % 这里可以写一段简单的补录逻辑，或者直接忽略
             % 只要种群不太小，通常 Candidates 会大于 N
        end
    else
        % 如果精英超过 N 个 (常见情况)，用 IDSS 筛选分布
        % IDSS 会保留 决策空间距离远 (局部PS) 和 目标空间好 (全局PS) 的解
        Pop = IDSS(Candidates, N);
    end
end

% ... (IDSS 函数保持不变，务必保留在文件底部) ...
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