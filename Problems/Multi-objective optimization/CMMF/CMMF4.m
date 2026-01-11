classdef CMMF4 < PROBLEM
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
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                elseif X(i, 1) == 0 || X(i, 1) == -1
                    THETA(i) = 1;
                else
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ (X(i, 1) + 1));
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
                if X(i, 1) >= 0
                    T(i) = (0.96 - X(i, 1)^2 - X(i, 2)^2)^2 + h(i);
                else
                    T(i) = (0.96 - (X(i, 1) + 1)^2 - X(i, 2)^2)^2 + h(i);
                end
                G(i, :) = 1 - [1, cumprod(sin(pi/2 * THETA(i)), 2)] .* [cos(pi/2 * THETA(i)), 1];
            end
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            PopCon = zeros(N, 3);
            PopCon(:, 1) = -X(:, 2);
            for i = 1 : N
                if X(i, 1) < 0
                    if (X(i, 1) + 1)^2 + X(i, 2)^2 <= 0.49
                        PopCon(i, 2) = (X(i, 1) + 1)^2 + X(i, 2)^2 - 0.36;
                        PopCon(i, 3) = -(X(i, 1) + 1)^2 - X(i, 2)^2 + 0.04;
                    else
                        PopCon(i, 2) = (X(i, 1) + 1)^2 + X(i, 2)^2 - 1;
                        PopCon(i, 3) = -(X(i, 1) + 1)^2 - X(i, 2)^2 + 0.64;
                    end
                elseif X(i, 1) >= 0
                    if X(i, 1)^2 + X(i, 2)^2 <= 0.49
                        PopCon(i, 2) = X(i, 1)^2 + X(i, 2)^2 - 0.36;
                        PopCon(i, 3) = -X(i, 1)^2 - X(i, 2)^2 + 0.04;
                    else
                        PopCon(i, 2) = X(i, 1)^2 + X(i, 2)^2 - 1;
                        PopCon(i, 3) = -X(i, 1)^2 - X(i, 2)^2 + 0.64;
                    end
                end
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS points (Circles)
            r = sqrt(0.96);
            phi = linspace(-pi, pi, 15000)';
            
            % Circle 1: centered at (0,0)
            C1 = [r*cos(phi), r*sin(phi)];
            % Circle 2: centered at (-1,0)
            C2 = [-1 + r*cos(phi), r*sin(phi)];
            
            Combined = [C1; C2];
            
            % 2. Strictly clip to problem bounds [-1, 1]
            Combined(Combined(:,1) < -1|Combined(:,1) > 1 | Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            
            % 3. Filter through constraints
            PopCon = obj.CalCon(Combined);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = Combined(Feasible, :);
            
            % 4. Generate PF
            objs = obj.CalObj(obj.POS);
            T_vals = sum(objs, 2) - 1; % For this problem, G1+G2 = 1? Actually verify.
            % Let's use a simpler check for global optimum if needed.
            % Since we sampled on optimal circles, we just need to verify feasibility.
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