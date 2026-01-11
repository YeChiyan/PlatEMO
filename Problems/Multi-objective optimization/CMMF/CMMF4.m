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
            
            THETA = 2/pi * atan(abs(X(:, 2)) ./ abs(X(:, 1)));
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) < -1/2
                    T(i) = (0.25 - (X(i, 1)+1)^2 - X(i, 2)^2)^2 + h(i);
                else
                    T(i) = (0.96 - X(i, 1)^2 - X(i, 2)^2)^2 + h(i);
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
                if X(i, 1) < -0.5
                    % Peak 1 at (-1,0) radius 0.5. val = (x+1)^2+y^2.
                    val = (X(i, 1)+1)^2 + X(i, 2)^2;
                    % Feasible band: [0.2, 0.3]
                    PopCon(i, 1) = max(0.2 - val, val - 0.3);
                else
                    % Peak 2 at (0,0) radius sqrt(0.96). val = x^2+y^2.
                    val = X(i, 1)^2 + X(i, 2)^2;
                    % Feasible band: [0.91, 1.01]
                    PopCon(i, 1) = max(0.91 - val, val - 1.01);
                end
                
                THETA = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                PopCon(i, 2) = min(max(0.1 - THETA, THETA - 0.4), ...
                    max(0.6 - THETA, THETA - 0.9));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            phi = linspace(-pi, pi, 15000)';
            r1 = sqrt(0.96);
            r2 = 0.5;
            C1 = [r1*cos(phi), r1*sin(phi)];
            C2 = [-1 + r2*cos(phi), r2*sin(phi)];
            Combined = [C1; C2];
            
            Combined(Combined(:,1) < -1|Combined(:,1) > 1 | Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            
            [nc, ~] = size(Combined);
            T_val = zeros(nc, 1);
            for i = 1 : nc
                if Combined(i, 1) < -0.5
                    T_val(i) = (0.25 - (Combined(i, 1)+1)^2 - Combined(i, 2)^2)^2;
                else
                    T_val(i) = (0.96 - Combined(i, 1)^2 - Combined(i, 2)^2)^2;
                end
            end
            POS_cand = Combined(T_val < 1e-4, :);
            
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