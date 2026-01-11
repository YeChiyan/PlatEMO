classdef CMMF1 < PROBLEM
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
            
            % Process input similar to the original function
            tempX = X;
            for i = 1 : N
                if tempX(i,1) <= 0
                    tempX(i,1) = tempX(i,1) + 1;
                end
            end
            
            Sx = cumsum(tempX(:, 1:M).^2, 2, 'reverse');
            THETA = 2/pi * atan(sqrt(Sx(:, 2:end)) ./ tempX(:, 1:M-1));
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = (0.98 - Sx(:, 1)).^2 + h;
            G = [ones(N, 1), cumprod(THETA, 2)] .* [1-THETA, ones(N, 1)];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            M = obj.M;
            [N, D] = size(X);
            
            tempX = X;
            for i = 1 : N
                if tempX(i,1) <= 0
                    tempX(i,1) = tempX(i,1) + 1;
                end
            end
            
            Sx = cumsum(tempX(:, 1:M).^2, 2, 'reverse');
            THETA = 2/pi * atan(sqrt(Sx(:, 2:end)) ./ tempX(:, 1:M-1));
            
            a = 1.0/4.0;
            b = 2.0/4.0;
            c = 3.0/4.0;
            
            PopCon = zeros(N, 6);
            for i = 1 : N
                if Sx(i, 1) <= 0.36
                    PopCon(i, 1) = Sx(i, 1) - 0.36;
                    PopCon(i, 2) = -Sx(i, 1) + 0.04;
                elseif X(i, 1) < 0 && X(i, 2) > 0
                    if THETA(i, 1) <= b
                        PopCon(i, 3) = -(THETA(i, 1) - a);
                        PopCon(i, 4) = Sx(i, 1) - 1;
                        PopCon(i, 5) = -Sx(i, 1) + 0.96;
                    else
                        PopCon(i, 3) = -THETA(i, 1) + c;
                        PopCon(i, 4) = Sx(i, 1) - 1;
                        PopCon(i, 5) = -Sx(i, 1) + 0.96;
                    end
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    if THETA(i, 1) <= b
                        PopCon(i, 3) = THETA(i, 1) - a;
                        PopCon(i, 4) = -Sx(i, 1) + 0.96;
                        PopCon(i, 5) = Sx(i, 1) - 1;
                    else
                        PopCon(i, 3) = -THETA(i, 1) + b;
                        PopCon(i, 4) = THETA(i, 1) - c;
                        PopCon(i, 5) = -Sx(i, 1) + 0.96;
                        PopCon(i, 6) = Sx(i, 1) - 1;
                    end
                elseif X(i, 1) >= 0 && X(i, 2) >= 0
                    if THETA(i, 1) <= b
                        PopCon(i, 3) = -THETA(i, 1);
                        PopCon(i, 4) = THETA(i, 1) - a;
                        PopCon(i, 5) = -Sx(i, 1) + 0.96;
                        PopCon(i, 6) = Sx(i, 1) - 1;
                    else
                        PopCon(i, 3) = -THETA(i, 1) + b;
                        PopCon(i, 4) = THETA(i, 1) - c;
                        PopCon(i, 5) = -Sx(i, 1) + 0.96;
                        PopCon(i, 6) = Sx(i, 1) - 1;
                    end
                else
                    if THETA(i, 1) <= b
                        PopCon(i, 3) = -(THETA(i, 1) - a);
                        PopCon(i, 4) = Sx(i, 1) - 1;
                        PopCon(i, 5) = -Sx(i, 1) + 0.96;
                    else
                        PopCon(i, 3) = -THETA(i, 1) + c;
                        PopCon(i, 4) = Sx(i, 1) - 1;
                        PopCon(i, 5) = -Sx(i, 1) + 0.96;
                    end
                end
            end
            
            % PlatEMO expects individual constraints to be violated if > 0
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample high-density points on the potential optimal circles
            r = sqrt(0.98);
            phi = linspace(-pi, pi, 15000)';
            
            % Circle 1: centered at (0,0)
            C1 = [r*cos(phi), r*sin(phi)];
            % Circle 2: centered at (-1,0)
            C2 = [-1 + r*cos(phi), r*sin(phi)];
            
            Combined = [C1; C2];
            
            % 2. Strictly clip to problem bounds [-1, 1]
            Combined(Combined(:,1) < -1 | Combined(:,1) > 1, :) = [];
            Combined(Combined(:,2) < -1 | Combined(:,2) > 1, :) = [];
            
            % 3. Filter through constraints
            PopCon = obj.CalCon(Combined);
            % Use a very strict feasibility threshold for reference data
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS  = Combined(Feasible, :);
            
            % 4. Generate PF (must be f1 + f2 = 1)
            % Only keep points where T is near 0 to ensure global PF
            objs = obj.CalObj(obj.POS);
            T_vals = sum(objs, 2) - 1;
            GlobalIdx = abs(T_vals) < 0.05;
            obj.POS = obj.POS(GlobalIdx, :);
            R = objs(GlobalIdx, :);
            
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