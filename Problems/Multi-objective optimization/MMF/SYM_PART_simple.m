classdef SYM_PART_simple < PROBLEM
    % <multi> <real> <multimodal>
    % SYM-PART-simple test function
    
    %------------------------------- Reference --------------------------------
    % Y. Liu, G. G. Yen, and D. Gong, A multi-modal multi-objective
    % evolutionary algorithm using two-archive and recombination strategies,
    % IEEE Transactions on Evolutionary Computation, 2019, 23(4): 660-674.
    %--------------------------------------------------------------------------
    
    properties
        POS; % Pareto optimal set
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            obj.lower    = [-20, -20];
            obj.upper    = [20, 20];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
            a = 1; b = 10; c = 8;
            [N, ~] = size(X);
            PopObj = zeros(N, 2);
            for i = 1:N
                % Replicate the cell-shifting logic from original SYM-PART
                temp_t1 = sign(X(i,1))*ceil((abs(X(i,1))-(a+c/2))/(2*a+c));
                temp_t2 = sign(X(i,2))*ceil((abs(X(i,2))-b/2)/b);
                t1 = sign(temp_t1)*min(abs(temp_t1),1);
                t2 = sign(temp_t2)*min(abs(temp_t2),1);
                x1 = X(i,1)-t1*(c+2*a);
                x2 = X(i,2)-t2*b;
                
                PopObj(i,1) = (x1+a)^2 + x2^2;
                PopObj(i,2) = (x1-a)^2 + x2^2;
            end
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj,N)
            % Generate points in Pareto optimal set
            % SYM-PART-simple has 9 Pareto sets (3x3 grid)
            a = 1; b = 10; c = 8;
            npos = 100;
            x1_base = linspace(-a, a, npos)';
            x2_base = zeros(npos, 1);
            t1_vals = [-1, 0, 1] * (c+2*a);
            t2_vals = [-1, 0, 1] * b;
            pos = [];
            for t1 = t1_vals
                for t2 = t2_vals
                    pos = [pos; x1_base + t1, x2_base + t2];
                end
            end
            obj.POS = pos;
            
            % Base PF for one cell
            x = linspace(-1, 1, N)';
            R(:,1) = (x+1).^2;
            R(:,2) = (x-1).^2;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            % Generate POS for IGDX calculation
            % SYM-PART-simple has 9 Pareto sets (3x3 grid)
            if isempty(obj.POS)
                a = 1; b = 10; c = 8;
                N = 100;
                x1_base = linspace(-a, a, N)';
                x2_base = zeros(N, 1);
                
                t1_vals = [-1, 0, 1] * (c+2*a);
                t2_vals = [-1, 0, 1] * b;
                
                pos = [];
                for t1 = t1_vals
                    for t2 = t2_vals
                        pos = [pos; x1_base + t1, x2_base + t2];
                    end
                end
                obj.POS = pos;
            end
            
            switch metName
                case 'IGDX'
                    score = feval(metName,Population,obj.POS);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
    end
end