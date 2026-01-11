classdef DCMMEAHA < ALGORITHM
    % <multi/multimodal> <real/binary/permutation><multimodal> <constrained>
    % Balancing exploration and exploitation in dynamic constrained multimodal multi-objective co-evolutionary algorithm
    % Swarm and Evolutionary Computation, 89, 101652, 2024
    % Guoqing Li (li241700@126.com)
    
    methods
        function main(Algorithm,Problem)
            
            %% Generate random population
            Population1 = Problem.Initialization(Problem.N/2);
            Population2 = Problem.Initialization(Problem.N/2);
            
            % Initialize external archive
            Archive = [];
            
            %                 x1=std(Population1.cons);
            %                 y1=std(Population1.objs);
            %                 H=1;
            %                 epsilon_x1=H* sum(x1)/Problem.M;
            %                 epsilon_y1=H* sum(y1)/Problem.M;
            [epsilon_x1,epsilon_y1]=parameters(Population1,Population2,Problem);
            
            
            %% Calculate fitness of populations
            [D_Dec,D_Pop,~]    = CalFitness(Population1.objs,Population1.cons,Population1.decs,epsilon_x1,Problem);
            [D,P,~]    = CalFitness2(Population2.objs,Population2.decs,epsilon_y1,Problem);
            
            len1=length(Population1);
            len2=length(Population2);
            
            
            
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                
                % --- Improvement: Niche Mating Restriction ---
                % Using dedicated niche selection to prevent cross-peak mating
                MatingPool1 = NicheTournamentSelection(Population1, D_Dec, len1, 0.1);
                MatingPool2 = NicheTournamentSelection(Population2, D, len2, 0.1);
                
                Parents1 = Population1(MatingPool1);
                Parents2 = Population2(MatingPool2);
                
                % --- Archive Feedback Mechanism (Guidance) ---
                if ~isempty(Archive) && rand < 0.2
                    % Inject archive individuals into mating pool (50% replacement)
                    n1 = length(Parents1);
                    Parents1(1:floor(n1/2)) = Archive(randi(length(Archive), 1, floor(n1/2)));
                    n2 = length(Parents2);
                    Parents2(1:floor(n2/2)) = Archive(randi(length(Archive), 1, floor(n2/2)));
                end
                
                Offspring1  = OperatorGA(Problem,Parents1);
                Offspring2  = OperatorGA(Problem,Parents2);
                
                [Population1,D_Dec,D_Pop] = EnvironmentalSelection([Population1,Offspring1,Offspring2],len1,epsilon_x1,Problem);
                [Population2,D,P] = EnvironmentalSelection2([Population2,Offspring1,Offspring2],len2,epsilon_y1,Problem);
                
                % --- Improvement: Update External Archive ---
                % Archive keeps all feasible non-dominated solutions found so far
                Archive = UpdateArchive(Archive, [Population1, Offspring1, Offspring2], Problem.N);
                
                len1=ceil(Problem.N/2*(1+sin((Problem.FE/Problem.maxFE)*0.5*pi)));
                len2=ceil(Problem.N/2*(1-sin((Problem.FE/Problem.maxFE)*0.5*pi)));
                if len2<5
                    len2=5;
                end
                [epsilon_x1,epsilon_y1]=parameters(Population1,Population2,Problem);
                %                 num=Problem.FE/Problem.maxFE;
                %                 if num<0.5
                %                     H=1-0.5*exp((num-0.5)/0.1);
                % %                     H=0.5*exp(((Problem.FE/Problem.maxFE)-1)/0.5);
                %                 else
                %                     H=0.5*exp(-(num-0.5)/0.1);
                % %                     H=1-0.5*exp(((Problem.FE/Problem.maxFE)-1)/0.5);
                %                 end
                %                 x1=std(Population1.cons);
                %                 y1=std(Population2.objs);
                %                 epsilon_x1=0.5*H* sum(x1)/Problem.M;
                %                 epsilon_y1=0.5*H* sum(y1)/Problem.M;
                
                
                %                 epsilon_k =H* sum(min(Population1.objs,[],1),2)/Problem.M;
                [D_Dec,D_Pop,~]    = CalFitness(Population1.objs,Population1.cons,Population1.decs,epsilon_x1,Problem);
                [D,P,~]    = CalFitness2(Population2.objs,Population2.decs,epsilon_y1,Problem);
            end
        end
    end
end

function MatingPool = NicheTournamentSelection(Population, Fitness, N, Sigma)
% Population: 当前种群 (INDIVIDUAL 数组)
% Fitness: 适应度/多样性评分 (D_Dec 或 D)，这里数值越小代表多样性越好 (密度越低)
% N: 需要选出的亲本数量
% Sigma: 小生境半径比例 (例如 0.1 代表种群中最靠近的 10% 个体)

MatingPool = zeros(1, N);
PopDec = Population.decs;
PopSize = length(Population);

% 确定邻域大小 (小生境规模)
nicheSize = max(2, ceil(PopSize * Sigma));

for i = 1 : 2 : N
    % 1. 锦标赛选出第一个亲本 (P1)
    candidates = randi(PopSize, 1, 2);
    if Fitness(candidates(1)) < Fitness(candidates(2))
        p1Idx = candidates(1);
    else
        p1Idx = candidates(2);
    end
    MatingPool(i) = p1Idx;
    
    % 2. 如果还需要选第二个亲本，则在 P1 的邻域内寻找
    if i + 1 <= N
        % 计算 P1 到所有人的距离
        distances = pdist2(PopDec(p1Idx, :), PopDec);
        [~, sortedIdx] = sort(distances);
        
        % 3. 在邻域内通过锦标赛选出第二个亲本 (P2)
        nicheIndices = sortedIdx(1:nicheSize);
        c2 = nicheIndices(randi(nicheSize, 1, 2));
        
        if Fitness(c2(1)) < Fitness(c2(2))
            p2Idx = c2(1);
        else
            p2Idx = c2(2);
        end
        MatingPool(i+1) = p2Idx;
    end
end
end