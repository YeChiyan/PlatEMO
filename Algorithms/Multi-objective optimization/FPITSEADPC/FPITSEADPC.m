classdef FPITSEADPC < ALGORITHM
    % <multi/many> <real> <multimodal>
    % Two-stage evolutionary algorithm with fuzzy preference indicator
    % p  --- 0.25 --- The proportion of the first stage
    % fi --- 0.1  --- The parameter for DPC cutoff distance (percentile)

    methods
        function main(Algorithm,Problem)
            %% Parameter Setting
            [p, fi] = Algorithm.ParameterSet(0.25, 0.1);
            
            %% Generate random population
            N = Problem.N;
            Population = Problem.Initialization();
            [~,CD] = EnvironmentalSelectionS1(Population,N);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Problem.FE < p*Problem.maxFE
                    %% The first stage of evolution (Global Search)
                    Matingpool = TournamentSelection(2,N,CD);
                    Offspring  = OperatorGA(Problem,Population(Matingpool));
                    Union      = [Population Offspring];
                    [Population,CD] = EnvironmentalSelectionS1(Union,N);
                else
                    %% The second stage of evolution (Local Mining with DPC)
                    % 这里调用修改后的独立进化策略
                    Population = IndependentEvolution(Problem,Population,N,fi);
                end
            end
        end
    end
end