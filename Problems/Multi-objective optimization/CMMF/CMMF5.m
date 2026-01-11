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
            
            % Angle calculation
            THETA = 2/pi * atan(abs(X(:, 2)) ./ abs(X(:, 1)));
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            Sx = sum(X(:, 1:M).^2, 2);
            T = (0.64 - Sx(:, 1)).^2 + h;
            G = [1-THETA, THETA];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            Sx = sum(X(:, 1:obj.M).^2, 2);
            PopCon = zeros(N, 2);
            for i = 1 : N
                % 1. Allow the circle Sx=0.64. Feasible band: [0.55, 0.75]
                PopCon(i, 1) = max(0.55 - Sx(i), Sx(i) - 0.75);
                
                % 2. Quadrant and sector constraints
                THETA = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                % Allow only Quadrants 1 and 3, and specific sectors
                if X(i, 1) * X(i, 2) < 0
                    PopCon(i, 2) = 1; % Infeasible quadrants
                else
                    PopCon(i, 2) = min(max(0.1 - THETA, THETA - 0.4), ...
                        max(0.6 - THETA, THETA - 0.9));
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