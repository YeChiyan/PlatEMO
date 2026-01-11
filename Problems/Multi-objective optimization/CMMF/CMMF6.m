classdef CMMF6 < PROBLEM
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
                else
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                end
            end
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = (1/16 - 0.25 .* X(:, 1).^2 - X(:, 2).^2).^2 + h;
            G = [ones(N, 1), cumprod(sin(pi/2 * THETA), 2)] .* [cos(pi/2 * THETA), ones(N, 1)];
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
                else
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                end
            end
            
            PopCon = zeros(N, 4);
            for i = 1 : N
                PopCon(i, 1) = THETA(i, 1) - 2/3;
                PopCon(i, 2) = -THETA(i, 1) + 1/3;
                PopCon(i, 3) = 1/25 - 1/4 * X(i, 1)^2 - X(i, 2)^2;
                PopCon(i, 4) = 1/4 * X(i, 1)^2 + X(i, 2)^2 - 1/9;
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS points (Ellipse)
            % Ellipse: 0.25*X1^2 + X2^2 = 1/16 (a=1/2, b=1/4)
            phi = linspace(-pi, pi, 15000)';
            Combined = [0.5*cos(phi), 0.25*sin(phi)];
            
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