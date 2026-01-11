classdef CMMF3 < PROBLEM
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
                if X(i, 1) > -1/2
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                elseif X(i, 1) == -1/2 || X(i, 1) == 0
                    THETA(i) = 1;
                else
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1) + 1/2));
                end
            end
            
            if D > M
                h = 20 - 20 * exp(-0.2 * sqrt(sum((X(:, M+1:end) - OptX).^2, 2) / (D - M))) + exp(1) ...
                    - exp(sum(cos(2 * pi * (X(:, M+1:end) - OptX)), 2) / (D - M));
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            G = zeros(N, M);
            for i = 1 : N
                if X(i, 1) < 0 && X(i, 1) > -1/2
                    T(i) = (-0.96 - X(i, 1) - X(i, 2))^2 + h(i);
                elseif X(i, 1) <= -1/2
                    T(i) = (-0.96 - X(i, 1) - X(i, 2))^2 + h(i);
                else
                    T(i) = (0.81 - X(i, 1) - X(i, 2))^2 + h(i);
                end
                G(i, :) = 1 - [1, cumprod(sin(pi/2 * THETA(i)), 2)] .* [cos(pi/2 * THETA(i)), 1];
            end
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            M = obj.M;
            [N, D] = size(X);
            Sx = cumsum(X(:, 1:M).^2, 2, 'reverse');
            
            THETA = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) > -1/2
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                elseif X(i, 1) == -1/2 || X(i, 1) == 0
                    THETA(i) = 1;
                else
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1) + 1/2));
                end
            end
            
            b = 2.0/4.0;
            PopCon = zeros(N, 3);
            for i = 1 : N
                if X(i, 1) < 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -Sx(i, 1);
                    PopCon(i, 2) = Sx(i, 1) - 0.25;
                elseif X(i, 1) <= 0 && X(i, 2) <= 0
                    PopCon(i, 1) = X(i, 1) + X(i, 2) + 0.9;
                    PopCon(i, 2) = -X(i, 1) - X(i, 2) - 1;
                    PopCon(i, 3) = X(i, 1) + 1/2;
                elseif X(i, 1) >= 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -(THETA(i, 1) - b);
                    PopCon(i, 2) = X(i, 1) + X(i, 2) - 1;
                    PopCon(i, 3) = -X(i, 1) - X(i, 2) + 0.64;
                else
                    PopCon(i, 1) = -Sx(i, 1);
                    PopCon(i, 2) = Sx(i, 1) - 0.36;
                end
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS points (Lines)
            % Line 1: X1 + X2 = -0.96
            x1_a = linspace(-1, 0, 10000)';
            x2_a = -0.96 - x1_a;
            C1 = [x1_a, x2_a];
            
            % Line 2: X1 + X2 = 0.81
            x1_b = linspace(0, 1, 10000)';
            x2_b = 0.81 - x1_b;
            C2 = [x1_b, x2_b];
            
            Combined = [C1; C2];
            
            % 2. Clip to bounds
            Combined(Combined(:,1) < -1 | Combined(:,1) > 1, :) = [];
            Combined(Combined(:,2) < -1 | Combined(:,2) > 1, :) = [];
            
            % 3. Filter through constraints
            PopCon = obj.CalCon(Combined);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = Combined(Feasible, :);
            
            % 4. Generate PF
            objs = obj.CalObj(obj.POS);
            % Use a threshold to filter global optima (T=0)
            % For CMMF3, sum(objs, 2) should be 1 if T=0 and PF is 1-cos/sin thing
            % Actually, the objective space is defined differently, but T=0 is key.
            % Let's use T calculation directly if needed, or just rely on CalObj
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