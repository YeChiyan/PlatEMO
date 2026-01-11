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
            
            % Mirroring for x1 <= 0 is often used in CMMF1 to create multiple peaks
            tempX = X;
            for i = 1 : N
                if tempX(i,1) <= 0
                    tempX(i,1) = tempX(i,1) + 1;
                end
            end
            
            Sx = sum(tempX(:, 1:M).^2, 2);
            % Use absolute x1 for theta calculation
            THETA = 2/pi * atan(abs(tempX(:, 2)) ./ abs(tempX(:, 1)));
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = (0.98 - Sx(:, 1)).^2 + h;
            G = [1-THETA, THETA];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            M = obj.M;
            [N, ~] = size(X);
            
            tempX = X;
            for i = 1 : N
                if tempX(i,1) <= 0
                    tempX(i,1) = tempX(i,1) + 1;
                end
            end
            
            Sx = sum(tempX(:, 1:M).^2, 2);
            THETA = 2/pi * atan(abs(tempX(:, 2)) ./ abs(tempX(:, 1)));
            
            PopCon = zeros(N, 3);
            for i = 1 : N
                % 1. Distance constraint to allow the peak at Sx=0.98
                % Band: [0.9, 1.05]
                PopCon(i, 1) = max(0.9 - Sx(i), Sx(i) - 1.05);
                
                % 2. Sector constraints to create niches
                % Allow THETA in [0.1, 0.4] or [0.6, 0.9]
                PopCon(i, 2) = min(max(0.1 - THETA(i), THETA(i) - 0.4), ...
                    max(0.6 - THETA(i), THETA(i) - 0.9));
                
                % 3. Quadrant separation
                % For CMMF1, we often want niches in different quadrants
                % Here we just ensure x2 is feasible if it matches some pattern
                PopCon(i, 3) = 0; % Could add more if needed
            end
            PopCon(PopCon < 0) = 0;
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj, N)
            % 1. Sample on the two circles defined by mirroring
            r = sqrt(0.98);
            phi = linspace(-pi, pi, 15000)';
            C1 = [r*cos(phi), r*sin(phi)];
            C2 = [-1 + r*cos(phi), r*sin(phi)];
            Combined = [C1; C2];
            
            % 2. Clip to bounds
            Combined(Combined(:,1) < -1|Combined(:,1) > 1 | Combined(:,2) < -1|Combined(:,2) > 1, :) = [];
            
            % 3. Filter through constraints
            PopCon = obj.CalCon(Combined);
            Feasible = all(PopCon <= 1e-4, 2);
            obj.POS = Combined(Feasible, :);
            
            % 4. Generate PF
            R = obj.CalObj(obj.POS);
            [~, idx] = sort(R(:, 1));
            R = R(idx, :);
            
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