function [Population,D_Dec,D_Pop] = EnvironmentalSelection(Population,N,epsilon_k,Problem)
% The environmental selection of SPEA2 based on objective space

    %% Calculate the fitness of each solution
    [D_Dec,D_Pop,Fitness] = CalFitness(Population.objs,Population.cons,Population.decs,epsilon_k,Problem);

    %% Environmental selection
    Next = Fitness == 0;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
%         Del  = Truncation(Population(Next).objs,Population(Next).decs,sum(Next)-N);
%         Temp = find(Next);
%         Next(Temp(Del)) = false;
        while sum(Next)>N
        Population1 = Population(Next);
        K=3;
        Index=find(Next==true);
        

            dobj = sort(pdist2(Population1.objs,Population1.objs)); dobj = sum(dobj(1:K,:));%dobj./max(max(dobj));
            ddec = sort(pdist2(Population1.decs,Population1.decs)); ddec = sum(ddec(1:K,:));%ddec./max(max(ddec));
            CrowdDis = dobj./max(dobj) + ddec./max(ddec);
            [A,~] = sortrows([CrowdDis',Index'],1);
            index=A(:,2);
            Next(index(1))=false; 
%             Next(index(1:(sum(Next)-N)))=false; 
        end
    end
    % Population for next generation
%     [D_Dec,D_Pop,Fitness] = CalFitness(Population.objs,Population.cons,Population.decs,epsilon_k,Problem);
    Population = Population(Next);
    D_Dec      = D_Dec(Next);
    D_Pop      = D_Pop(Next);
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