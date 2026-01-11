function [Population,D_Dec,R] = EnvironmentalSelection2(Population,N,epsilon_k,Problem)
% The environmental selection of SPEA2 based on objective space

%% Calculate the fitness of each solution
[D_Dec,R,Fitness] = CalFitness2(Population.objs,Population.decs,epsilon_k,Problem);

%% Environmental selection
Next = Fitness == 0;
if sum(Next) < N
    [~,Rank] = sort(Fitness);
    Next(Rank(1:N)) = true;
elseif sum(Next) > N
    % --- Improvement: Clearing Method (Explicit Niche) ---
    Sigma = 0.1;
    Next = ApplyClearing(Population, Fitness, Next, N, Sigma);
    
    while sum(Next) > N
        Population2 = Population(Next);
        K=3;
        Index=find(Next==true);
        dobj = sort(pdist2(Population2.objs,Population2.objs)); dobj = sum(dobj(1:K,:));
        ddec = sort(pdist2(Population2.decs,Population2.decs)); ddec = sum(ddec(1:K,:));
        CrowdDis = dobj./max(dobj) + ddec./max(ddec);
        [A,~] = sortrows([CrowdDis',Index'],1);
        index=A(:,2);
        Next(index(1))=false;
    end
    
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    end
end
% Population for next generation
Population = Population(Next);
D_Dec      = D_Dec(Next);
R      = R(Next);
end

function Del = Truncation(PopObj,PopDec,K)
% Select part of the solutions by truncation

%% Truncation

Distance_Dec = pdist2(PopDec,PopDec);
Distance_Dec(logical(eye(length(Distance_Dec)))) = inf;

D = Distance_Dec;

Del = false(1,size(PopObj,1));
while sum(Del) < K
    Remain   = find(~Del);
    Temp     = sort(D(Remain,Remain),2);
    [~,Rank] = sortrows(Temp);
    Del(Remain(Rank(1))) = true;
end
end