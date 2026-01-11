classdef CMMF17 < PROBLEM
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
                if X(i, 1) == 0 || X(i, 1) == 1 || X(i, 1) == -1
                    THETA(i) = 1;
                else
                    if X(i, 1) > 0 && X(i, 2) >= 0
                        THETA(i) = 2/pi * atan(X(i, 2) ./ X(i, 1));
                    elseif X(i, 1) < 0 && X(i, 2) >= 0
                        THETA(i) = 2/pi * atan((1 - X(i, 2)) ./ (X(i, 1) + 1));
                    elseif X(i, 1) < 0 && X(i, 2) < 0
                        THETA(i) = 2/pi * atan((X(i, 2) + 1) ./ (X(i, 1) + 1));
                    elseif X(i, 1) > 0 && X(i, 2) < 0
                        THETA(i) = 2/pi * atan((X(i, 2) + 1) ./ (1 - X(i, 1)));
                    end
                end
            end
            
            if D > M
                h = 20 - 20 * exp(-0.2 * sqrt(sum((X(:, M+1:end) - OptX).^2, 2) / (D - M))) + exp(1) ...
                    - exp(sum(cos(2 * pi .* (X(:, M+1:end) - OptX)), 2) / (D - M));
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            G = zeros(N, M);
            for i = 1 : N
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    T(i) = (1 - X(i, 1)^2 - X(i, 2)^2)^2 + h(i);
                elseif X(i,1) < 0 && X(i,2) >= 0
                    T(i) = (1 - (1+X(i,1))^2 - (X(i,2)-1)^2)^2 + h(i);
                elseif X(i,1) < 0 && X(i,2) < 0
                    T(i) = (1 - (1+X(i,1))^2 - (X(i,2)+1)^2)^2 + h(i);
                elseif X(i,1) >= 0 && X(i,2) < 0
                    T(i) = (1 - (X(i,1)-1)^2 - (X(i,2)+1)^2)^2 + h(i);
                end
                G(i, :) = [ones(1, 1), cumprod(THETA(i), 2)] .* [1-THETA(i), ones(1, 1)];
            end
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            PopCon = zeros(N, 2);
            for i = 1 : N
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -X(i, 2) + sqrt(8) * X(i, 1)^4;
                    PopCon(i, 2) = X(i, 2) - sqrt(4) * X(i, 1)^3;
                elseif X(i, 1) < 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -sqrt(4) * (1 + X(i, 1))^3 + 1 - X(i, 2);
                    PopCon(i, 2) = X(i, 2) + sqrt(8) * (1 + X(i, 1))^4 - 1;
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    PopCon(i, 1) = -X(i, 2) - 1 + ((X(i, 1) + 1) ./ sqrt(4)).^(1/3);
                    PopCon(i, 2) = X(i, 2) + 1 - ((X(i, 1) + 1) ./ sqrt(8)).^(1/4);
                elseif X(i, 1) >= 0 && X(i, 2) < 0
                    PopCon(i, 1) = -X(i, 2) - 1 + ((1 - X(i, 1)) ./ sqrt(4)).^(1/3);
                    PopCon(i, 2) = X(i, 2) + 1 - ((1 - X(i, 1)) ./ sqrt(8)).^(1/4);
                end
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS points (Circles)
            phi = linspace(-pi, pi, 15000)';
            C1 = [cos(phi), sin(phi)];
            C2 = [-1 + cos(phi), 1 + sin(phi)];
            C3 = [-1 + cos(phi), -1 + sin(phi)];
            C4 = [1 + cos(phi), -1 + sin(phi)];
            Combined = [C1; C2; C3; C4];
            
            % 2. Clip to bounds
            Combined(Combined(:,1) < -1|Combined(:,1) > 1 | Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            
            % 3. Filter through constraints
            PopCon = obj.CalCon(Combined);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = Combined(Feasible, :);
            
            % 4. Filter for global optima (T=0)
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