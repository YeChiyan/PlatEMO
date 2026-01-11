classdef CMMF8 < PROBLEM
    % <multi> <real> <multimodal> <constrained>
    % Constrained multi-modal multi-objective test function
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D); obj.D = 2; end
            % CMMF8 seems to have a different range based on Pop(i,1)>1 and Pop(i,1)<1
            % Let's use [0, 2] for X1 and [0, 2] for X2?
            % Actually, the code says Pop(i,1)>1 and Pop(i,1)<=1.
            % If X1 is in [0, 2], it works.
            obj.lower    = repmat(0, 1, obj.D);
            obj.upper    = repmat(2, 1, obj.D);
            obj.encoding = ones(1, obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, X)
            M = obj.M;
            [N, D] = size(X);
            OptX = 0.2;
            
            THETA = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) > 1
                    THETA(i) = 2/pi * atan(X(i, 2) / (X(i, 1) - 1));
                elseif X(i, 1) < 1
                    THETA(i) = 2/pi * atan(X(i, 2) / (1 - X(i, 1)));
                else
                    THETA(i) = 1;
                end
            end
            
            if D > M
                h = 20 - 20 * exp(-0.2 * sqrt(sum((X(:, M+1:end) - OptX).^2, 2) / (D - M))) + exp(1) ...
                    - exp(sum(cos(2 * pi * (X(:, M+1:end) - OptX)), 2) / (D - M));
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            G = zeros(N, M);
            for i = 1 : N
                if X(i, 1) > 1
                    T(i) = (4 - 4 * (2 - X(i, 1))^2 - (2 - X(i, 2))^2)^2 + h(i);
                else
                    T(i) = (4 - 4 * (1 - X(i, 1))^2 - (X(i, 2))^2)^2 + h(i);
                end
                G(i, :) = 1 - [1, cumprod(sin(pi/2 * THETA(i)), 2)] .* [cos(pi/2 * THETA(i)), 1];
            end
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            THETA = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) > 1
                    THETA(i) = 2/pi * atan(X(i, 2) / (X(i, 1) - 1));
                elseif X(i, 1) < 1
                    THETA(i) = 2/pi * atan(X(i, 2) / (1 - X(i, 1)));
                else
                    THETA(i) = 1;
                end
            end
            
            PopCon = zeros(N, 2);
            for i = 1 : N
                PopCon(i, 1) = THETA(i) - 2/3;
                PopCon(i, 2) = -THETA(i) + 1/3;
            end
            PopCon(PopCon < 0) = 0;
        end
    end
end