function Offspring = OperatorDE(Problem,Parent,Parameter)
% 标准差分进化算子 (DE/rand/1/bin)

if nargin > 2
    [CR,F,ProM,DisM] = deal(Parameter{:});
else
    [CR,F,ProM,DisM] = deal(1,0.5,1,20);
end

ParentDecs = Parent.decs;
[N,D]  = size(ParentDecs);

% DE/rand/1 变异
% 随机选择三个父代 r1, r2, r3
P1 = ParentDecs(randperm(N),:);
P2 = ParentDecs(randperm(N),:);
P3 = ParentDecs(randperm(N),:);
Site = rand(N,D) < CR;
OffspringDecs = P1;
OffspringDecs(Site) = P1(Site) + F*(P2(Site)-P3(Site));

% 多项式变异 (增加一点点扰动，防止彻底重合)
Lower = repmat(Problem.lower,N,1);
Upper = repmat(Problem.upper,N,1);
Site  = rand(N,D) < ProM/D;
mu    = rand(N,D);
temp  = Site & mu<=0.5;
OffspringDecs       = min(max(OffspringDecs,Lower),Upper);
OffspringDecs(temp) = OffspringDecs(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
    (1-(OffspringDecs(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(DisM+1)).^(1/(DisM+1))-1);
temp = Site & mu>0.5;
OffspringDecs(temp) = OffspringDecs(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
    (1-(Upper(temp)-OffspringDecs(temp))./(Upper(temp)-Lower(temp))).^(DisM+1)).^(1/(DisM+1)));

% 转换回 SOLUTION 对象
Offspring = Problem.Evaluation(OffspringDecs);
end
