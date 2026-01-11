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
            
            THETA = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) > 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ X(i, 1));
                elseif X(i, 1) == 0
                    THETA(i) = 1;
                elseif X(i, 1) < 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                end
            end
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = (1 - X(:, 1).^2 - 4 * X(:, 2).^2).^2 + h;
            G = zeros(N, M);
            for i = 1 : N
                G(i, :) = 1 - [1, cumprod(sin(pi/2 * THETA(i)), 2)] .* [cos(pi/2 * THETA(i)), 1];
            end
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            M = obj.M;
            [N, D] = size(X);
            THETA = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) > 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ X(i, 1));
                elseif X(i, 1) == 0
                    THETA(i) = 1;
                elseif X(i, 1) < 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                end
            end
            
            PopCon = zeros(N, 4);
            for i = 1 : N
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -THETA(i) + 1/2;
                    PopCon(i, 2) = THETA(i) - 2/3;
                    PopCon(i, 3) = X(i, 1)^2 + X(i, 2)^2 - 1;
                elseif X(i, 1) < 0 && X(i, 2) >= 0
                    PopCon(i, 1) = THETA(i) - 1/2;
                    PopCon(i, 2) = -THETA(i) + 1/3;
                    PopCon(i, 3) = X(i, 1)^2 + 3 * X(i, 2)^2 - 1;
                    PopCon(i, 4) = -X(i, 1)^2 - 5 * X(i, 2)^2 + 1;
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    PopCon(i, 1) = THETA(i) - 1/2;
                    PopCon(i, 2) = -THETA(i) + 1/3;
                    PopCon(i, 3) = -X(i, 1)^2 - 4 * X(i, 2)^2 + 1 - 0.04;
                    PopCon(i, 4) = X(i, 1)^2 + 4 * X(i, 2)^2 - 1 - 0.04;
                else
                    PopCon(i, 1) = -THETA(i) + 1/2;
                    PopCon(i, 2) = THETA(i) - 2/3;
                    PopCon(i, 3) = -X(i, 1)^2 - 4 * (X(i, 2) - 0.04)^2 + 1;
                    PopCon(i, 4) = X(i, 1)^2 + 4 * (X(i, 2) + 0.04)^2 - 1;
                end
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS points (Ellipse)
            % Ellipse: X1^2 + 4*X2^2 = 1 (a=1, b=1/2)
            phi = linspace(-pi, pi, 15000)';
            Combined = [cos(phi), 0.5*sin(phi)];
            
            % 2. Clip to bounds
            Combined(Combined(:,1) < -1|Combined(:,1) > 1 | Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            
            % 3. Filter through constraints
            PopCon = obj.CalCon(Combined);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = Combined(Feasible, :);
            
            % 4. Generate PF
            objs = obj.CalObj(obj.POS);
            R = objs;
            
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