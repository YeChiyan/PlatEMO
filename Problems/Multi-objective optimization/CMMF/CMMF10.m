classdef CMMF10 < PROBLEM
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
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)), 1e-6));
                elseif X(i, 1) >= 0 && X(i, 2) < 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)+1) / max(abs(X(i, 1)), 1e-6));
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)+1) / max(abs(X(i, 1)+1), 1e-6));
                else
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)+1), 1e-6));
                end
            end
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    T(i) = (0.4 - (X(i, 1)^2 + X(i, 2)^2))^2 + h(i);
                elseif X(i, 1) >= 0 && X(i, 2) < 0
                    T(i) = (0.4 - (X(i, 1)^2 + (X(i, 2)+1)^2))^2 + h(i);
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    T(i) = (0.4 - ((X(i, 1)+1)^2 + (X(i, 2)+1)^2))^2 + h(i);
                else
                    T(i) = (0.4 - ((X(i, 1)+1)^2 + X(i, 2)^2))^2 + h(i);
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
                % 1. Distance band for the 4 circles (r^2=0.4)
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    dist2 = X(i, 1)^2 + X(i, 2)^2;
                elseif X(i, 1) >= 0 && X(i, 2) < 0
                    dist2 = X(i, 1)^2 + (X(i, 2)+1)^2;
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    dist2 = (X(i, 1)+1)^2 + (X(i, 2)+1)^2;
                else
                    dist2 = (X(i, 1)+1)^2 + X(i, 2)^2;
                end
                % Feasible band: [0.35, 0.45]
                PopCon(i, 1) = max(0.35 - dist2, dist2 - 0.45);
                
                % 2. Sector constraint for each peak
                % Calculate local theta
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    th = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)), 1e-6));
                elseif X(i, 1) >= 0 && X(i, 2) < 0
                    th = 2/pi * atan(abs(X(i, 2)+1) / max(abs(X(i, 1)), 1e-6));
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    th = 2/pi * atan(abs(X(i, 2)+1) / max(abs(X(i, 1)+1), 1e-6));
                else
                    th = 2/pi * atan(abs(X(i, 2)) / max(abs(X(i, 1)+1), 1e-6));
                end
                
                PopCon(i, 2) = min(max(0.1 - th, th - 0.4), max(0.6 - th, th - 0.9));
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            phi = linspace(-pi, pi, 20000)';
            r = sqrt(0.4);
            C1 = [r*cos(phi), r*sin(phi)];
            C2 = [r*cos(phi), -1 + r*sin(phi)];
            C3 = [-1 + r*cos(phi), -1 + r*sin(phi)];
            C4 = [-1 + r*cos(phi), r*sin(phi)];
            Combined = [C1; C2; C3; C4];
            
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