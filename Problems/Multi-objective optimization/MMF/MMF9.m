classdef MMF9 < PROBLEM
    % <multi> <real> <multimodal>
    % MMF9 test function
    
    %------------------------------- Reference --------------------------------
    % C. Yue, B. Qu, and J. Liang, A multi-objective particle swarm optimizer
    % using ring topology for solving multimodal multiobjective Problems, IEEE
    % Transactions on Evolutionary Computation, 2018, 22(5): 805-817.
    %--------------------------------------------------------------------------
    
    properties
        POS; % Pareto optimal set
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            obj.lower    = [0.1, 0.1];
            obj.upper    = [1.1, 1.1];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            PopObj(:,1) = X(:,1);
            g = 2 - sin(2*pi*X(:,2)).^6;
            PopObj(:,2) = g ./ X(:,1);
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj,N)
            % Generate points in Pareto optimal set
            % MMF9 has two Pareto sets in the decision space: y = 0.25 and y = 0.75
            x1 = linspace(0.1, 1.1, N)';
            obj.POS = [x1, repmat(0.25, N, 1); x1, repmat(0.75, N, 1)];
            
            % Generate points on Pareto front
            R(:,1) = linspace(0.1, 1.1, N)';
            R(:,2) = 1 ./ R(:,1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            % Generate POS for IGDX calculation
            % MMF9 has two Pareto sets in the decision space: y = 0.25 and y = 0.75
            N = 500;
            x1 = linspace(0.1, 1.1, N)';
            obj.POS = [x1, repmat(0.25, N, 1); x1, repmat(0.75, N, 1)];
            
            switch metName
                case 'IGDX'
                    score = feval(metName,Population,obj.POS);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
    end
end