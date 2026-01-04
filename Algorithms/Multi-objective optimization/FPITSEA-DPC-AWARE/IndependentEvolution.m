function Population = IndependentEvolution(Problem,Population,N,fi)
% Independent evolution strategy
% Includes: Small Cluster Protection & DPC integration

    %% 1. Use DPC to divide the population
    % 计算拥挤度用于后续选择
    CD = Crowding(Population.decs);
    
    % 调用 DPC 进行聚类 (请确保 DPC.m 已还原为无噪声过滤的标准版)
    HD = DPC(Population, fi);

    lpop = cell(1,length(HD));
    lCD = cell(1,length(HD));
    
    % 准备子种群数据
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
            % === 情况 A：正常大种群 (>=2) ===
            % 内部锦标赛选择，正常交叉变异
            mateidx = TournamentSelection(2, n_sub, -subCD);
            Off = OperatorGA(Problem, subpop(mateidx));
            lpop{i} = [lpop{i} Off];
            
        elseif n_sub == 1
            % === 情况 B：孤儿种群 (Local Explorer) ===
            % 策略：原地变异 (Self-Mutation)
            % 避免与外部杂交导致位置偏移
            Parent = subpop(1);
            
            % 输入两个相同的父代给 OperatorGA
            % SBX 交叉无效化，多项式变异生效
            Off = OperatorGA(Problem, [Parent, Parent]);
            
            % 仅保留变异后的这一个后代
            lpop{i} = [lpop{i} Off(1)];
        end
        % n_sub == 0 的情况会自动跳过
    end

    %% 3. Environmental selection
    % 这里调用下方的子函数
    Population = EnvironmentalSelectionS2(lpop,N);
end

%% ================== 以下是之前缺失的辅助函数 ==================

function Pop = EnvironmentalSelectionS2(lpop,N)
% Environmental Selection in the second stage of evolution

    ndpop = [];
    dpop = [];
    
    % 1. 提取每个子种群的非支配解
    for i = 1:length(lpop)
        cpop = lpop{i};
        if isempty(cpop)
            continue;
        end
        [FNo1,~] = NDSort(cpop.objs,length(cpop));
        ndpop = [ndpop cpop(FNo1==1)];  % 精英保留
        dpop = [dpop cpop(FNo1~=1)];    % 淘汰候选
    end

    % 2. 如果精英不足 N，从被淘汰的解中补充
    % 先对所有备选解进行一次全局排序
    if ~isempty(dpop)
        [FNo2,~] = NDSort(dpop.objs,length(dpop));
        k = 1;
        while length(ndpop) < N && k <= max(FNo2)
            ndpop = [ndpop dpop(FNo2==k)];
            k = k+1;
        end
    end

    % 3. 如果精英超过 N，使用 IDSS 进行筛选
    Pop = ndpop;
    if length(Pop) > N
        Pop = IDSS(Pop,N);
    end
end

function ChoosePop = IDSS(Population,N)
% Improved distance-based subset selection method (IDSS)
% 同时也确保 IDSS 被包含在文件内
   
    %% 初始化
    pop = [];
    % 计算内部距离矩阵 (只取最近的几个距离)
    decdist = pdist2(Population.decs,Population.decs,'euclidean','Smallest',2);
    objdist = pdist2(Population.objs,Population.objs,'euclidean','Smallest',2);
    
    % 取最近邻距离 (排除自己到自己)
    decdist = decdist(2,:);
    objdist = objdist(2,:);
    
    % 归一化
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
    
    % 贪心选择第一个点：综合距离最大的点
    distsum = nordecdist+norobjdist;
    [~,midx] = max(distsum);
    pop = [pop Population(midx)];
    
    retainpop = Population;
    retainpop(midx) = []; % 从备选集中移除

    %% 迭代选择剩余点
    while length(pop) < N && ~isempty(retainpop)
        % 计算备选集到已选集的距离
        decdist = pdist2(pop.decs,retainpop.decs,'euclidean','Smallest',1);
        objdist = pdist2(pop.objs,retainpop.objs,'euclidean','Smallest',1);
        
        % 归一化
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
        
        % 选择离已选集最远的点
        distsum = nordecdist+norobjdist;
        [~,midx] = max(distsum);
        
        pop = [pop retainpop(midx)];
        retainpop(midx) = [];
    end
    
    ChoosePop = pop;
end