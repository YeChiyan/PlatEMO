classdef CMMF17 < PROBLEM
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
                if X(i, 1) == 0 || X(i, 1) == 1 || X(i, 1) == -1
                    THETA(i) = 1;
                else
                    if X(i, 1) > 0 && X(i, 2) >= 0
                        THETA(i) = 2/pi * atan(X(i, 2) ./ X(i, 1));
                    elseif X(i, 1) < 0 && X(i, 2) >= 0
                        THETA(i) = 2/pi * atan((1 - X(i, 2)) ./ (X(i, 1) + 1));
                    elseif X(i, 1) < 0 && X(i, 2) < 0
                        THETA(i) = 2/pi * atan((X(i, 2) + 1) ./ (X(i, 1) + 1));
                    elseif X(i, 1) > 0 && X(i, 2) < 0
                        THETA(i) = 2/pi * atan((X(i, 2) + 1) ./ (1 - X(i, 1)));
                    end
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
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    T(i) = (1 - X(i, 1)^2 - X(i, 2)^2)^2 + h(i);
                elseif X(i,1) < 0 && X(i,2) >= 0
                    T(i) = (1 - (1+X(i,1))^2 - (X(i,2)-1)^2)^2 + h(i);
                elseif X(i,1) < 0 && X(i,2) < 0
                    T(i) = (1 - (1+X(i,1))^2 - (X(i,2)+1)^2)^2 + h(i);
                elseif X(i,1) >= 0 && X(i,2) < 0
                    T(i) = (1 - (X(i,1)-1)^2 - (X(i,2)+1)^2)^2 + h(i);
                end
                G(i, :) = [ones(1, 1), cumprod(THETA(i), 2)] .* [1-THETA(i), ones(1, 1)];
            end
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            [N, ~] = size(X);
            PopCon = zeros(N, 2);
            for i = 1 : N
                if X(i, 1) >= 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -X(i, 2) + sqrt(8) * X(i, 1)^4;
                    PopCon(i, 2) = X(i, 2) - sqrt(4) * X(i, 1)^3;
                elseif X(i, 1) < 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -sqrt(4) * (1 + X(i, 1))^3 + 1 - X(i, 2);
                    PopCon(i, 2) = X(i, 2) + sqrt(8) * (1 + X(i, 1))^4 - 1;
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    PopCon(i, 1) = -X(i, 2) - 1 + ((X(i, 1) + 1) ./ sqrt(4)).^(1/3);
                    PopCon(i, 2) = X(i, 2) + 1 - ((X(i, 1) + 1) ./ sqrt(8)).^(1/4);
                elseif X(i, 1) >= 0 && X(i, 2) < 0
                    PopCon(i, 1) = -X(i, 2) - 1 + ((1 - X(i, 1)) ./ sqrt(4)).^(1/3);
                    PopCon(i, 2) = X(i, 2) + 1 - ((1 - X(i, 1)) ./ sqrt(8)).^(1/4);
                end
            end
            PopCon(PopCon < 0) = 0;
        end
    end
end