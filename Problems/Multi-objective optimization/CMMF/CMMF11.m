classdef CMMF11 < PROBLEM
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
            
            % Angle calculation based on mirroring at X1=0
            THETA = zeros(N, 1);
            for i = 1 : N
                diff = abs(X(i, 1) + 0.5); % Shift center
                THETA(i) = 2/pi * atan(abs(X(i, 2)) / max(abs(diff), 1e-6));
            end
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) >= 0
                    % Peak 1: Line X1+X2 = 0.97
                    T(i) = (0.97 - X(i, 1) - X(i, 2))^2 + h(i);
                else
                    % Peak 2: Ellipse (X1+1)^2 + 6*X2^2 = 0.25
                    T(i) = (0.25 - (X(i, 1)+1)^2 - 6*X(i, 2)^2)^2 + h(i);
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
                    % PS Line band: [0.92, 1.02]
                    val = X(i, 1) + X(i, 2);
                    PopCon(i, 1) = max(0.92 - val, val - 1.02);
                else
                    % PS Ellipse band: [0.2, 0.3]
                    val = (X(i, 1)+1)^2 + 6*X(i, 2)^2;
                    PopCon(i, 1) = max(0.2 - val, val - 0.3);
                end
                
                % Sector constraint
                diff = abs(X(i, 1) + 0.5);
                th = 2/pi * atan(abs(X(i, 2)) / max(abs(diff), 1e-6));
                PopCon(i, 2) = min(max(0.1 - th, th - 0.4), max(0.6 - th, th - 0.9));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            phi = linspace(-pi, pi, 15000)';
            % Peak 1: Line X1+X2 = 0.97. Corrected range for mirroring?
            % Actually, if X1 >= 0, X1 in [0, 0.97]
            x1_a = linspace(0, 0.97, 10000)';
            x2_a = 0.97 - x1_a;
            C1 = [x1_a, x2_a];
            % Peak 2: (X1+1)^2 + 6*X2^2 = 0.25
            C2 = [-1 + 0.5*cos(phi), 1/sqrt(24)*sin(phi)];
            Combined = [C1; C2];
            
            Combined(Combined(:,1) < -1|Combined(:,1) > 1 | Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            
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