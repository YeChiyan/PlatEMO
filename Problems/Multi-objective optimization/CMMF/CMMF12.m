classdef CMMF12 < PROBLEM
    % <multi> <real> <multimodal> <constrained>
    % Constrained multi-modal multi-objective test function
    
    properties
        POS;    % Pareto optimal set for IGDX calculation
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D); obj.D = 2; end
            obj.lower    = repmat(-1, 1, obj.D);
            obj.upper    = repmat(1, 1, obj.D);
            obj.encoding = ones(1, obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, X)
            M = obj.M;
            [N, D] = size(X);
            OptX = 0.2;
            
            % Angle calculation
            THETA = 2/pi * atan(abs(X(:, 2)) ./ abs(X(:, 1)));
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            % Optimal Ellipse: X1^2 + 4*X2^2 = 1
            val = X(:, 1).^2 + 4 * X(:, 2).^2;
            T = (1 - val).^2 + h;
            G = [1-THETA, THETA];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            PopCon = zeros(N, 2);
            for i = 1 : N
                val = X(i, 1)^2 + 4 * X(i, 2)^2;
                % 1. Piecewise distance bands for quadrants
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    % Q1: Broad band [0.9, 1.1]
                    PopCon(i, 1) = max(0.9 - val, val - 1.1);
                elseif X(i, 1) < 0 && X(i, 2) >= 0
                    % Q2: Moderate band [0.95, 1.05]
                    PopCon(i, 1) = max(0.95 - val, val - 1.05);
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    % Q3: Narrow band [0.98, 1.02]
                    PopCon(i, 1) = max(0.98 - val, val - 1.02);
                else
                    % Q4: Offset band [0.95, 1.05]
                    PopCon(i, 1) = max(0.95 - val, val - 1.05);
                end
                
                % 2. Sector constraint
                THETA = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                PopCon(i, 2) = min(max(0.1 - THETA, THETA - 0.4), max(0.6 - THETA, THETA - 0.9));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            phi = linspace(-pi, pi, 15000)';
            Combined = [cos(phi), 0.5*sin(phi)];
            Combined(Combined(:,1) < -1|Combined(:,1) > 1 | Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            PopCon = obj.CalCon(Combined);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = Combined(Feasible, :);
            R = obj.CalObj(obj.POS);
            if obj.D > 2
                obj.POS = [obj.POS, repmat(repmat(0.2, 1, obj.D-2), size(obj.POS, 1), 1)];
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if isempty(obj.POS)
                obj.GetOptimum(1000);
            end
            R = obj.CalObj(obj.POS);
            [~, idx] = sort(R(:, 1));
            R = R(idx, :);
        end
        %% Calculate the metric value
        function score = CalMetric(obj, metName, Population)
            if isempty(obj.POS)
                obj.GetOptimum(2000);
            end
            switch metName
                case 'IGDX'
                    score = feval(metName, Population, obj.POS);
                otherwise
                    score = feval(metName, Population, obj.CalObj(obj.POS));
            end
        end
    end
end