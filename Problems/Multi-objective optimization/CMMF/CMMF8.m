classdef CMMF8 < PROBLEM
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
            obj.lower    = repmat(0, 1, obj.D);
            obj.upper    = repmat(2, 1, obj.D);
            obj.encoding = ones(1, obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, X)
            M = obj.M;
            [N, D] = size(X);
            OptX = 0.2;
            
            % THETA calculation based on mirroring at X1=1
            THETA = zeros(N, 1);
            for i = 1 : N
                diff = abs(X(i, 1) - 1);
                if diff > 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) / diff);
                else
                    THETA(i) = 1;
                end
            end
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) > 1
                    % Centered at (2,2), b=2, a=1 => (X-2)^2 + 0.25(Y-2)^2 = 1
                    T(i) = ((X(i, 1)-2)^2 + 0.25*(X(i, 2)-2)^2 - 1)^2 + h(i);
                else
                    % Centered at (1,0), b=2, a=1 => (X-1)^2 + 0.25(Y-0)^2 = 1
                    T(i) = ((X(i, 1)-1)^2 + 0.25*(X(i, 2)-0)^2 - 1)^2 + h(i);
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
                % 1. Ellipse band constraint
                if X(i, 1) > 1
                    val = (X(i, 1)-2)^2 + 0.25*(X(i, 2)-2)^2;
                else
                    val = (X(i, 1)-1)^2 + 0.25*(X(i, 2)-0)^2;
                end
                % Peak at val = 1. Feasible band: [0.9, 1.1]
                PopCon(i, 1) = max(0.9 - val, val - 1.1);
                
                % 2. Sector constraint
                diff = abs(X(i, 1) - 1);
                THETA = 2/pi * atan(abs(X(i, 2)) / max(diff, 1e-6));
                % Allow THETA in [0.1, 0.4] or [0.6, 0.9]
                PopCon(i, 2) = min(max(0.1 - THETA, THETA - 0.4), ...
                    max(0.6 - THETA, THETA - 0.9));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample potential PS bits
            phi = linspace(-pi, pi, 15000)';
            % Peak 1: (X-2)^2 + 0.25(Y-2)^2 = 1 => X=2+cos, Y=2+2sin
            C1 = [2 + cos(phi), 2 + 2*sin(phi)];
            % Peak 2: (X-1)^2 + 0.25(Y-0)^2 = 1 => X=1+cos, Y=0+2sin
            C2 = [1 + cos(phi), 2*sin(phi)];
            Combined = [C1; C2];
            
            % 2. Clip to bounds [0, 2]
            Combined(Combined(:,1) < 0|Combined(:,1) > 2 | Combined(:,2) < 0|Combined(:,2) > 2, :) = [];
            
            % 3. Filter through constraints
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