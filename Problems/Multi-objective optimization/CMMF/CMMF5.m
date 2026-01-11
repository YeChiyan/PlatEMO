classdef CMMF5 < PROBLEM
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
            
            Sx = cumsum(X(:, 1:M).^2, 2, 'reverse');
            T = (0.64 - Sx(:, 1)).^2 + h;
            G = [ones(N, 1), cumprod(THETA, 2)] .* [1-THETA, ones(N, 1)];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            M = obj.M;
            [N, D] = size(X);
            Sx = cumsum(X(:, 1:M).^2, 2, 'reverse');
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
            
            a = 1.0/4.0;
            b = 2.0/4.0;
            c = 3.0/4.0;
            d = 4.0/4.0;
            
            PopCon = zeros(N, 5);
            for i = 1 : N
                PopCon(i, 1) = -X(i, 1) * X(i, 2) + eps;
                if X(i, 1) <= 0 && X(i, 2) <= 0
                    if THETA(i, 1) <= b
                        PopCon(i, 2) = THETA(i, 1) - a;
                        PopCon(i, 3) = -THETA(i, 1);
                        PopCon(i, 4) = Sx(i, 1) - 0.81;
                        PopCon(i, 5) = -Sx(i, 1) + 0.49;
                    else
                        PopCon(i, 2) = THETA(i, 1) - d;
                        PopCon(i, 3) = -THETA(i, 1) + c;
                        PopCon(i, 4) = Sx(i, 1) - 0.81;
                        PopCon(i, 5) = -Sx(i, 1) + 0.49;
                    end
                elseif X(i, 1) >= 0 && X(i, 2) >= 0
                    if THETA(i, 1) <= b
                        PopCon(i, 2) = -THETA(i, 1);
                        PopCon(i, 3) = THETA(i, 1) - a;
                        PopCon(i, 4) = -X(i, 1) - X(i, 2) + 0.5;
                        PopCon(i, 5) = X(i, 1) + X(i, 2) - 1.5;
                    else
                        PopCon(i, 2) = -THETA(i, 1) + c;
                        PopCon(i, 3) = THETA(i, 1) - d;
                        PopCon(i, 4) = -X(i, 1) - X(i, 2) + 0.5;
                        PopCon(i, 5) = X(i, 1) + X(i, 2) - 1.5;
                    end
                end
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS points (Circle)
            r = sqrt(0.64);
            phi = linspace(-pi, pi, 15000)';
            Combined = [r*cos(phi), r*sin(phi)];
            
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