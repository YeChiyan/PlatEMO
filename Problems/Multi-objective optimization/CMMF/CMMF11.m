classdef CMMF11 < PROBLEM
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
                if X(i, 1) > 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ X(i, 1));
                elseif X(i, 1) == 0 || X(i, 1) == -1
                    THETA(i) = 1;
                elseif X(i, 1) < 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ (X(i, 1) + 1));
                end
            end
            
            if D > M
                h = sum((X(:, M+1:D) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            G = zeros(N, M);
            for i = 1 : N
                if X(i, 1) >= 0
                    T(i) = (0.97 - X(i, 1) - X(i, 2))^2 + h(i);
                else
                    T(i) = (1/4 - (X(i, 1) + 1)^2 - 6 * X(i, 2)^2)^2 + h(i);
                end
                G(i, :) = 1 - [1, cumprod(sin(pi/2 * THETA(i)), 2)] .* [cos(pi/2 * THETA(i)), 1];
            end
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            M = obj.M;
            [N, D] = size(X);
            THETA = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) > 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ X(i, 1));
                elseif X(i, 1) == 0 || X(i, 1) == -1
                    THETA(i) = 1;
                elseif X(i, 1) < 0
                    THETA(i) = 2/pi * atan(abs(X(i, 2)) ./ (X(i, 1) + 1));
                end
            end
            
            PopCon = zeros(N, 4);
            for i = 1 : N
                PopCon(i, 1) = -X(i, 2) + eps;
                if X(i, 1) >= 0
                    if X(i, 2) >= X(i, 1)
                        PopCon(i, 2) = -THETA(i) + 2/3;
                    else
                        PopCon(i, 2) = THETA(i) - 1/3;
                    end
                    if X(i, 1)^2 + X(i, 2)^2 <= 0.25
                        PopCon(i, 3) = -X(i, 1)^2 - 9 * X(i, 2)^2 + 1/4;
                        PopCon(i, 4) = X(i, 1)^2 + 4 * X(i, 2)^2 - 1/4;
                    else
                        PopCon(i, 3) = -X(i, 1) - X(i, 2) + 0.96;
                        PopCon(i, 4) = X(i, 1) + X(i, 2) - 0.98;
                    end
                else
                    if X(i, 2) >= 1 + X(i, 1)
                        PopCon(i, 2) = -THETA(i) + 2/3;
                    elseif X(i, 2) < 1 + X(i, 1)
                        PopCon(i, 2) = THETA(i) - 1/3;
                    end
                    if (X(i, 1) + 1)^2 + X(i, 2)^2 <= 0.25
                        PopCon(i, 3) = -(X(i, 1) + 1)^2 - 9 * X(i, 2)^2 + 1/4;
                        PopCon(i, 4) = (X(i, 1) + 1)^2 + 4 * X(i, 2)^2 - 1/4;
                    else
                        PopCon(i, 3) = -(X(i, 1) + 1) - X(i, 2) + 0.96;
                        PopCon(i, 4) = (X(i, 1) + 1) + X(i, 2) - 0.98;
                    end
                end
            end
            PopCon(PopCon < 0) = 0;
        end
    end
end