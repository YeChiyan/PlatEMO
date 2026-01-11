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
            
            % Angle calculation
            THETA = 2/pi * atan(abs(X(:, 2)) ./ abs(X(:, 1)));
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            % Ellipse-based peak
            val = 0.25 * X(:, 1).^2 + X(:, 2).^2;
            T = (0.0625 - val).^2 + h;
            G = [1-THETA, THETA];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            PopCon = zeros(N, 2);
            for i = 1 : N
                % 1. Ellipse factor
                val = 0.25 * X(i, 1)^2 + X(i, 2)^2;
                % Peak at 0.0625. Feasible band: [0.04, 0.09]
                PopCon(i, 1) = max(0.04 - val, val - 0.09);
                
                % 2. Sector constraint
                THETA = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                % Allow THETA in [0.1, 0.4] or [0.6, 0.9]
                PopCon(i, 2) = min(max(0.1 - THETA, THETA - 0.4), ...
                    max(0.6 - THETA, THETA - 0.9));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS ellipse
            % 0.25*X1^2 + X2^2 = 0.0625 => X1^2/0.25 + X2^2/0.0625 = 1 => a=0.5, b=0.25
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