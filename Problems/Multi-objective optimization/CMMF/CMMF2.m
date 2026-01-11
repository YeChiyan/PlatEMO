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
            
            % Angle calculation
            THETA = 2/pi * atan(abs(X(:, 2)) ./ abs(X(:, 1)));
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
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
            
            G = [1-THETA, THETA];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            PopCon = zeros(N, 2);
            for i = 1 : N
                % Ellipse distance from origin (squared distance proxy)
                val = X(i, 1)^2 + 0.25 * X(i, 2)^2;
                if X(i, 1) < 0
                    % Peak 1 at val = 0.0625. Feasible band: [0.04, 0.09]
                    PopCon(i, 1) = max(0.04 - val, val - 0.09);
                else
                    % Peak 2 at val = 0.25. Feasible band: [0.2, 0.3]
                    PopCon(i, 1) = max(0.2 - val, val - 0.3);
                end
                
                % Standard sector constraint
                THETA = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                PopCon(i, 2) = min(max(0.1 - THETA, THETA - 0.4), ...
                    max(0.6 - THETA, THETA - 0.9));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS points (Ellipses)
            phi = linspace(-pi, pi, 15000)';
            % Peak 1: X1^2 + 0.25*X2^2 = 1/16 (a=0.25, b=0.5)
            C1 = [0.25*cos(phi), 0.5*sin(phi)];
            % Peak 2: X1^2 + 0.25*X2^2 = 0.25 (a=0.5, b=1.0)
            C2 = [0.5*cos(phi), 1.0*sin(phi)];
            Combined = [C1; C2];
            
            % 2. Clip to bounds
            Combined(Combined(:,1) < -1|Combined(:,1) > 1 | Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            
            % 3. Filter for T = 0 (Global Optima)
            [nc, ~] = size(Combined);
            T_val = zeros(nc, 1);
            for i = 1 : nc
                if Combined(i, 1) < 0
                    T_val(i) = (0.0625 - Combined(i, 1)^2 - 0.25 * Combined(i, 2)^2)^2;
                else
                    T_val(i) = (0.25 - Combined(i, 1)^2 - 0.25 * Combined(i, 2)^2)^2;
                end
            end
            POS_cand = Combined(T_val < 1e-4, :);
            
            % 4. Filter through constraints
            PopCon = obj.CalCon(POS_cand);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = POS_cand(Feasible, :);
            
            R = obj.CalObj(obj.POS);
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