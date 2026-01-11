classdef CMMF16 < PROBLEM
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
            obj.lower    = repmat(-2, 1, obj.D);
            obj.upper    = repmat(2, 1, obj.D);
            obj.encoding = ones(1, obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, X)
            M = obj.M;
            [N, D] = size(X);
            OptX = 0.2;
            
            % Angle calculation
            THETA = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) >= 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)), 1e-6));
                elseif X(i, 1) < 0 && X(i, 2) > 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)-2) / max(abs(X(i, 1)+2), 1e-6));
                else
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)), 1e-6));
                end
            end
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) >= 0
                    % Peak 1: Circle r=1 at (0,0)
                    T(i) = (1 - (X(i, 1)^2 + X(i, 2)^2))^2 + h(i);
                elseif X(i, 1) < 0 && X(i, 2) > 0
                    % Peak 2: Circle r=1 at (-2,2)
                    T(i) = (1 - ((X(i, 1)+2)^2 + (X(i, 2)-2)^2))^2 + h(i);
                else
                    % Peak 3: Circle r=2 at (0,0)
                    T(i) = (4 - (X(i, 1)^2 + X(i, 2)^2))^2 + h(i);
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
                if X(i, 1) >= 0
                    dist2 = X(i, 1)^2 + X(i, 2)^2;
                    % Peak at 1.0. Feasible: [0.9, 1.1]
                    PopCon(i, 1) = max(0.9 - dist2, dist2 - 1.1);
                elseif X(i, 1) < 0 && X(i, 2) > 0
                    dist2 = (X(i, 1)+2)^2 + (X(i, 2)-2)^2;
                    % Peak at 1.0. Feasible: [0.9, 1.1]
                    PopCon(i, 1) = max(0.9 - dist2, dist2 - 1.1);
                else
                    dist2 = X(i, 1)^2 + X(i, 2)^2;
                    % Peak at 4.0. Feasible: [3.9, 4.1]
                    PopCon(i, 1) = max(3.9 - dist2, dist2 - 4.1);
                end
                
                % Sector constraint
                if X(i, 1) >= 0
                    th = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)), 1e-6));
                elseif X(i, 1) < 0 && X(i, 2) > 0
                    th = 2/pi * atan(abs(X(i, 2)-2) / max(abs(X(i, 1)+2), 1e-6));
                else
                    th = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)), 1e-6));
                end
                PopCon(i, 2) = min(max(0.1 - th, th - 0.4), max(0.6 - th, th - 0.9));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            phi = linspace(-pi, pi, 20000)';
            C1 = [cos(phi), sin(phi)];
            C2 = [-2 + cos(phi), 2 + sin(phi)];
            C3 = [2*cos(phi), 2*sin(phi)];
            Combined = [C1; C2; C3];
            Combined(Combined(:,1) < -2|Combined(:,1) > 2 | Combined(:,2) < -2|Combined(:,2) > 2, :) = [];
            PopCon = obj.CalCon(Combined);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = Combined(Feasible, :);
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