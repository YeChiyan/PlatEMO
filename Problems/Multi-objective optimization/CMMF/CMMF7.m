classdef CMMF7 < PROBLEM
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
                if X(i, 1) > -1/2
                    T(i) = (0.96 - X(i, 1)^2 - X(i, 2)^2)^2 + h(i);
                else
                    T(i) = (0.25 - (X(i, 1)+0.5)^2 - X(i, 2)^2)^2 + h(i);
                end
                G(i, :) = 1 - [1, cumprod(sin(pi/2 * THETA(i)), 2)] .* [cos(pi/2 * THETA(i)), 1];
            end
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            PopCon = zeros(N, 4);
            Sx = sum(X(:, 1:obj.M).^2, 2);
            for i = 1 : N
                % Standard CMMF constraints:
                % 1. Quadrant constraints (often used to split PS)
                PopCon(i, 1) = X(i, 1) * X(i, 2);
                
                % 2. Common distance-based constraints for CMMF problems
                % Adjusted to be broad enough to enclose the piecewise optimal circles
                PopCon(i, 2) = -Sx(i) + 0.2;
                PopCon(i, 3) = Sx(i) - 1.2;
                
                % 3. Angle-based or sector-based constraint
                PopCon(i, 4) = abs(X(i, 2)) - abs(X(i, 1));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. High-density sampling of BOTH circles
            phi = linspace(-pi, pi, 20000)';
            r1 = sqrt(0.96);
            r2 = 0.5;
            C1 = [r1*cos(phi), r1*sin(phi)];
            C2 = [-0.5 + r2*cos(phi), r2*sin(phi)];
            Combined = [C1; C2];
            
            % 2. Clip to decision space bounds
            Combined(Combined(:,1) < -1|Combined(:,1) > 1 | Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            
            % 3. Filter for GLOBAL optima (T=0)
            % Points must satisfy the correct piecewise formula
            [nc, ~] = size(Combined);
            T_val = zeros(nc, 1);
            for i = 1 : nc
                if Combined(i, 1) > -0.5
                    T_val(i) = (0.96 - Combined(i, 1)^2 - Combined(i, 2)^2)^2;
                else
                    T_val(i) = (0.25 - (Combined(i, 1)+0.5)^2 - Combined(i, 2)^2)^2;
                end
            end
            
            % Keep points where T is extremely small
            isGlobal = T_val < 1e-6;
            POS_cand = Combined(isGlobal, :);
            
            % 4. Filter through constraints
            PopCon = obj.CalCon(POS_cand);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = POS_cand(Feasible, :);
            
            % Remove duplicates and sort
            if ~isempty(obj.POS)
                obj.POS = unique(unique(obj.POS, 'rows'), 'rows');
            end
            
            % Complete dimensionality if needed
            if obj.D > 2
                obj.POS = [obj.POS, repmat(repmat(0.2, 1, obj.D-2), size(obj.POS, 1), 1)];
            end
            
            % 5. Generate PF
            R = obj.CalObj(obj.POS);
            [~, idx] = sort(R(:, 1));
            R = R(idx, :);
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