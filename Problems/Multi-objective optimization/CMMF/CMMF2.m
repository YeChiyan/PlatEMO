classdef CMMF2 < PROBLEM
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
            
            Sx = cumsum(X(:, 1:M).^2, 2, 'reverse');
            THETA = 2/pi * atan(sqrt(Sx(:, 2:end)) ./ abs(X(:, 1:M-1)));
            
            if D > M
                h = sum(100 * ((X(:, M+1:D-1) - OptX).^2 - (X(:, M+2:D) - OptX)).^2 + (X(:, M+1:D-1) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) < 0
                    T(i) = (1/16 - X(i, 1)^2 - 0.25 * X(i, 2)^2)^2 + h(i);
                else
                    T(i) = (0.25 - X(i, 1)^2 - 0.25 * X(i, 2)^2)^2 + h(i);
                end
            end
            
            G = [ones(N, 1), cumprod(sin(pi/2 * THETA), 2)] .* [cos(pi/2 * THETA), ones(N, 1)];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            M = obj.M;
            [N, D] = size(X);
            
            Sx = cumsum(X(:, 1:M).^2, 2, 'reverse');
            THETA = 2/pi * atan(sqrt(Sx(:, 2:end)) ./ abs(X(:, 1:M-1)));
            
            b = 2.0/4.0;
            PopCon = zeros(N, 3);
            for i = 1 : N
                if X(i, 1) < 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -Sx(i, 1) - 0.04;
                    PopCon(i, 2) = Sx(i, 1) - 0.25;
                    PopCon(i, 3) = THETA(i, 1) - b;
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    PopCon(i, 1) = -Sx(i, 1) - 0.04;
                    PopCon(i, 2) = Sx(i, 1) - 0.25;
                    PopCon(i, 3) = -THETA(i, 1) + b;
                elseif X(i, 1) >= 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -(THETA(i, 1) - b);
                    PopCon(i, 2) = -Sx(i, 1) + 0.36;
                    PopCon(i, 3) = Sx(i, 1) - 1;
                else
                    PopCon(i, 1) = THETA(i, 1) - b;
                    PopCon(i, 2) = -Sx(i, 1) + 0.16;
                    PopCon(i, 3) = Sx(i, 1) - 0.49;
                end
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS points (Ellipses)
            phi = linspace(-pi, pi, 15000)';
            % PS candidates for the two ellipses
            % Ellipse 1: X1^2 + 0.25*X2^2 = 1/16 (a=1/4, b=1/2)
            C1 = [0.25*cos(phi), 0.5*sin(phi)];
            % Ellipse 2: X1^2 + 0.25*X2^2 = 0.25 (a=1/2, b=1)
            C2 = [0.5*cos(phi), 1.0*sin(phi)];
            
            Combined = [C1; C2];
            
            % 2. Clip to bounds
            Combined(Combined(:,1) < -1|Combined(:,1) > 1, :) = [];
            Combined(Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            
            % 3. Filter through constraints
            PopCon = obj.CalCon(Combined);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = Combined(Feasible, :);
            
            % 4. Filter for global optima (T=0)
            objs = obj.CalObj(obj.POS);
            T_vals = sum(objs, 2) - 1;
            GlobalIdx = abs(T_vals) < 0.05;
            obj.POS = obj.POS(GlobalIdx, :);
            R = objs(GlobalIdx, :);
            
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