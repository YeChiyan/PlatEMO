function Population = IndependentEvolution(Problem,Population,N,fi)
% Independent evolution strategy with Noise Handling

    %% 1. Use DPC to divide the population
    CD = Crowding(Population.decs);
    
    % [修改] 获取聚类结果和噪声索引
    [HD, noise_idx] = DPC(Population, fi);

    lpop = cell(1, length(HD));
    lCD  = cell(1, length(HD));
    
    % 准备有效子种群
    for i = 1:length(HD)
        lpop{i} = Population(HD{i});
        lCD{i}  = CD(HD{i});
    end

    %% 2. Reproduction for Valid Clusters (正常子种群独立进化)
    for i = 1:length(HD)
        subpop = lpop{i};
        subCD  = lCD{i};
        if length(subpop) >= 2
            mateidx = TournamentSelection(2, length(subpop), -subCD);
            Off     = OperatorGA(Problem, subpop(mateidx));
            lpop{i} = [lpop{i} Off];
        end
    end
    
    %% 3. [核心修改] Reproduction for Noise (噪声处理机制)
    % 策略：让噪声点与全局优秀个体杂交，将其拉回正常区域
    if ~isempty(noise_idx)
        NoisePop = Population(noise_idx);
        N_noise = length(NoisePop);
        
        % 从全局种群中选择配偶 (Global Partners)
        % 使用锦标赛选择，倾向于选择拥挤距离大(稀疏)的全局个体，保持多样性
        % 注意：这里是在整个 Population 中选
        global_mate_idx = TournamentSelection(2, N_noise, -CD);
        GlobalMates = Population(global_mate_idx);
        
        % 构建配对池：[Noise_1, Global_1, Noise_2, Global_2, ...]
        % 这样 OperatorGA 在按顺序两两配对时，就是 Noise vs Global
        MatingPool = [NoisePop, GlobalMates];
        
        % 产生后代
        Off_Noise = OperatorGA(Problem, MatingPool);
        
        % 将噪声及其后代作为一个单独的组，加入到 lpop 中
        % 这样在 EnvironmentalSelectionS2 中，它们会和大家一起竞争
        % 如果变好了，IDSS 会保留；如果还差，就会被淘汰
        lpop{end+1} = [NoisePop, Off_Noise];
    end

    %% 4. Environmental selection
    Population = EnvironmentalSelectionS2(lpop, N);
end

% ... EnvironmentalSelectionS2 和 IDSS 函数保持不变 ...
% ... 请确保文件底部包含这两个函数的定义 ...
function Pop = EnvironmentalSelectionS2(lpop,N)
% Environmental Selection in the second stage of evolution

%%
ndpop = [];
dpop = [];
for i = 1:length(lpop)
    cpop = lpop{i};
    [FNo1,~] = NDSort(cpop.objs,length(cpop));
    ndpop = [ndpop cpop(FNo1==1)];  % Non-dominated solutions in subpopulation
    dpop = [dpop cpop(FNo1~=1)];    % Dominated solutions in subpopulation
end

[FNo2,~] = NDSort(dpop.objs,length(dpop));
k = 1;
while length(ndpop) < N
    ndpop = [ndpop dpop(FNo2==k)];
    k = k+1;
end

% Use IDSS to select the next generation of the population
Pop = ndpop;
if length(Pop) > N
    Pop = IDSS(Pop,N);
end
end

function ChoosePop = IDSS(Population,N)
% Improved distance-based subset selection method
   
    %%
    pop = [];
    decdist = pdist2(Population.decs,Population.decs,'euclidean','Smallest',2);
    objdist = pdist2(Population.objs,Population.objs,'euclidean','Smallest',2);
    decdist = decdist(2,:);
    objdist = objdist(2,:);
    nordecdist = (decdist-min(decdist))./(max(decdist)-min(decdist));
    norobjdist = (objdist-min(objdist))./(max(objdist)-min(objdist));    
    distsum = nordecdist+norobjdist;
    [~,midx] = max(distsum);
    pop = [pop Population(midx)];
    retainpop = Population;
    retainpop(midx) = [];

    %% Iterative selection
    while length(pop) < N
        decdist = pdist2(pop.decs,retainpop.decs,'euclidean','Smallest',1);
        objdist = pdist2(pop.objs,retainpop.objs,'euclidean','Smallest',1);
        nordecdist = (decdist-min(decdist))./(max(decdist)-min(decdist));
        norobjdist = (objdist-min(objdist))./(max(objdist)-min(objdist));
        distsum = nordecdist+norobjdist;
        [~,midx] = max(distsum);
        pop = [pop retainpop(midx)];
        retainpop(midx) = [];
    end
    
    ChoosePop = pop;
end