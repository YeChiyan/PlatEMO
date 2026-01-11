classdef CMMF13 < PROBLEM
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
            THETA = zeros(N, 1);
            for i = 1 : N
                THETA(i) = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)), 1e-6));
            end
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            Sx = sum(X(:, 1:M).^2, 2);
            T = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    T(i) = (0.64 - Sx(i)).^2 + h(i);
                elseif X(i, 1) < 0 && X(i, 2) >= 0
                    T(i) = (0.36 - Sx(i)).^2 + h(i);
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    T(i) = (0.64 - Sx(i)).^2 + h(i);
                else
                    T(i) = (0.09 - Sx(i)).^2 + h(i);
                end
            end
            
            G = [1-THETA, THETA];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            Sx = sum(X(:, 1:obj.M).^2, 2);
            PopCon = zeros(N, 2);
            for i = 1 : N
                % 1. Piecewise distance bands for quadrants
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    % Peak at 0.64
                    PopCon(i, 1) = max(0.59 - Sx(i), Sx(i) - 0.69);
                elseif X(i, 1) < 0 && X(i, 2) >= 0
                    % Peak at 0.36
                    PopCon(i, 1) = max(0.31 - Sx(i), Sx(i) - 0.41);
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    % Peak at 0.64
                    PopCon(i, 1) = max(0.59 - Sx(i), Sx(i) - 0.69);
                else
                    % Peak at 0.09
                    PopCon(i, 1) = max(0.04 - Sx(i), Sx(i) - 0.14);
                end
                
                % 2. Sector constraint
                THETA = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)), 1e-6));
                PopCon(i, 2) = min(max(0.1 - THETA, THETA - 0.4), max(0.6 - THETA, THETA - 0.9));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            phi = linspace(-pi, pi, 15000)';
            r1 = sqrt(0.64);
            r2 = sqrt(0.36);
            r3 = sqrt(0.09);
            Combined = [r1*cos(phi), r1*sin(phi); r2*cos(phi), r2*sin(phi); r3*cos(phi), r3*sin(phi)];
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