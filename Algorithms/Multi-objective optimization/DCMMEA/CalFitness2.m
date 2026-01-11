function [D_Dec,D_Pop,Fitness] = CalFitness2(PopObj,PopDec,epsilon_k,Problem)
% Calculate the fitness of each solution based on original objective space

    N = size(PopObj,1);

    %% Detect the dominance relation between each two solutions
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(PopObj(i,:)<(1+epsilon_k).*PopObj(j,:)) - any((1+epsilon_k).*PopObj(i,:)>PopObj(j,:));
            if k == 1
                Dominate(i,j) = true;
            elseif k == -1
                Dominate(j,i) = true;
            end
        end
    end
    
    DIS=pdist([Problem.upper;Problem.lower]);
    dis=pdist2(PopDec,PopDec);
    Distance=dis./DIS;
    Dominate1=Dominate.*Distance;    
    
    %% Calculate S(i)
    S = sum(Dominate1,2);
    
    %% Calculate R(i)
    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S.*Dominate(:,i));
    end
    
    %% Calculate D(i)
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    D_Pop = 1./(Distance(:,floor(sqrt(N)))+2);
    
    Distance = pdist2(PopDec,PopDec);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    D_Dec = 1./(Distance(:,floor(sqrt(N)))+2);

%     D = D_Pop + D_Dec;
    D_Pop = D_Pop + D_Dec;
    %% Calculate the fitnesses
%     Fitness = R + D';
    Fitness = R ;
    
end