classdef CMMF14 < PROBLEM
    % <multi> <real> <multimodal> <constrained>
    % Constrained multi-modal multi-objective test function
    
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
            
            THETA = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) > -1/2
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs(X(i, 1)));
                elseif X(i, 1) == 0 || X(i, 1) == -1/2
                    THETA(i) = 1;
                else
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ abs((X(i, 1) + 1/2)));
                end
            end
            
            if D > M
                h = 20 - 20 * exp(-0.2 * sqrt(sum((X(:, M+1:end) - OptX).^2, 2) / (D - M))) + exp(1) ...
                    - exp(sum(cos(2 * pi .* (X(:, M+1:end) - OptX)), 2) / (D - M));
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            G = zeros(N, M);
            for i = 1 : N
                if X(i, 1) > -1/2
                    T(i) = (0.64 - X(i, 1)^2 - X(i, 2)^2)^2 + h(i);
                elseif X(i, 1) <= -1/2
                    T(i) = (0.25 - (X(i, 1)+0.5)^2 - X(i, 2)^2)^2 + h(i);
                end
                G(i, :) = [ones(1, 1), cumprod(sin(pi/2 * THETA(i)), 2)] .* [cos(pi/2 * THETA(i)), ones(1, 1)];
            end
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            PopCon = zeros(N, 5);
            for i = 1 : N
                PopCon(i, 1) = X(i, 1) * X(i, 2) + eps;
                if X(i, 1) <= 0
                    PopCon(i, 2) = 0.0361 - (X(i, 1) + 0.5)^2 - X(i, 2)^2;
                    PopCon(i, 3) = -0.04 + (X(i, 1) + 0.5)^2 + X(i, 2)^2;
                    PopCon(i, 4) = X(i, 1) + 0.5;
                    PopCon(i, 5) = (X(i, 1) + 0.5) + X(i, 2);
                elseif X(i, 1) > 0
                    PopCon(i, 2) = abs(X(i, 2)) - abs(X(i, 1));
                    PopCon(i, 3) = -X(i, 1)^2 - X(i, 2)^2 + 0.85;
                    PopCon(i, 4) = X(i, 1)^2 + X(i, 2)^2 - 0.86;
                end
            end
            PopCon(PopCon < 0) = 0;
        end
    end
end