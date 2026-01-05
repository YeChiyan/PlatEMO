function ChoosePop = IDSS(Population,N)
% Improved distance-based subset selection method (IDSS)
% 独立文件，供全局调用

%% 初始化
pop = [];
% 计算内部距离矩阵
decdist = pdist2(Population.decs,Population.decs,'euclidean','Smallest',2);
objdist = pdist2(Population.objs,Population.objs,'euclidean','Smallest',2);

decdist = decdist(2,:);
objdist = objdist(2,:);

% 归一化处理
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
