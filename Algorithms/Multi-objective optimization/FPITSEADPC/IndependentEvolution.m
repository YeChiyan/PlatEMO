function Population = IndependentEvolution(Problem,Population,N,fi)
% Independent evolution strategy with DPC

    %% 1. Use DPC to divide the population
    % 计算拥挤度用于后续的交配选择
    CD = Crowding(Population.decs);
    
    % [修改点]：这里不再计算 Fitness，也不再调用 NBC
    % 直接调用 DPC 进行聚类，fi 作为截断距离的百分比参数
    HD = DPC(Population, fi);

    % 将种群按聚类结果分组
    lpop = cell(1,length(HD));
    lCD  = cell(1,length(HD));
    for i = 1:length(HD)
        lpop{i} = Population(HD{i});
        lCD{i}  = CD(HD{i});
    end

    %% 2. Reproduction within the subpopulation (独立进化)
    for i = 1:length(HD)
        subpop = lpop{i};
        subCD  = lCD{i};
        % 只有当子种群个体数大于等于2时才能进行交叉变异
        if length(subpop) >= 2
            % 在子种群内部进行锦标赛选择
            mateidx = TournamentSelection(2,length(subpop),-subCD);
            Off     = OperatorGA(Problem,subpop(mateidx));
            lpop{i} = [lpop{i} Off];
        end
    end

    %% 3. Environmental selection (环境选择)
    Population = EnvironmentalSelectionS2(lpop,N);
end

function Pop = EnvironmentalSelectionS2(lpop,N)
% Environmental Selection in the second stage of evolution
% 这是一个局部函数，用于处理合并后的筛选

    ndpop = [];
    dpop  = [];
    
    % 对每个子种群分别提取非支配解
    for i = 1:length(lpop)
        cpop = lpop{i};
        [FNo1,~] = NDSort(cpop.objs,length(cpop));
        ndpop = [ndpop cpop(FNo1==1)];  % 子种群中的非支配解
        dpop  = [dpop cpop(FNo1~=1)];   % 被支配解
    end

    % 如果精英解不够，从被支配解中按层级补充
    [FNo2,~] = NDSort(dpop.objs,length(dpop));
    k = 1;
    while length(ndpop) < N && k <= max(FNo2)
        ndpop = [ndpop dpop(FNo2==k)];
        k = k+1;
    end

    % 使用 IDSS 进行最终筛选，确保双空间分布均匀
    Pop = ndpop;
    if length(Pop) > N
        Pop = IDSS(Pop,N);
    end
end