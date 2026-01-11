classdef CMMF2 < PROBLEM
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
            
            Sx = cumsum(X(:, 1:M).^2, 2, 'reverse');
            THETA = 2/pi * atan(sqrt(Sx(:, 2:end)) ./ abs(X(:, 1:M-1)));
            
            if D > M
                h = sum(100 * ((X(:, M+1:D-1) - OptX).^2 - (X(:, M+2:D) - OptX)).^2 + (X(:, M+1:D-1) - OptX).^2, 2);
            else
                h = zeros(N, 1);
            end
            
            T = zeros(N, 1);
            for i = 1 : N
                if X(i, 1) < 0
                    T(i) = (1/16 - X(i, 1)^2 - 0.25 * X(i, 2)^2)^2 + h(i);
                else
                    T(i) = (0.25 - X(i, 1)^2 - 0.25 * X(i, 2)^2)^2 + h(i);
                end
            end
            
            G = [ones(N, 1), cumprod(sin(pi/2 * THETA), 2)] .* [cos(pi/2 * THETA), ones(N, 1)];
            PopObj = G .* repmat((1+T), 1, M);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj, X)
            M = obj.M;
            [N, D] = size(X);
            
            Sx = cumsum(X(:, 1:M).^2, 2, 'reverse');
            THETA = 2/pi * atan(sqrt(Sx(:, 2:end)) ./ abs(X(:, 1:M-1)));
            
            b = 2.0/4.0;
            PopCon = zeros(N, 3);
            for i = 1 : N
                if X(i, 1) < 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -Sx(i, 1) - 0.04;
                    PopCon(i, 2) = Sx(i, 1) - 0.25;
                    PopCon(i, 3) = THETA(i, 1) - b;
                elseif X(i, 1) < 0 && X(i, 2) < 0
                    PopCon(i, 1) = -Sx(i, 1) - 0.04;
                    PopCon(i, 2) = Sx(i, 1) - 0.25;
                    PopCon(i, 3) = -THETA(i, 1) + b;
                elseif X(i, 1) >= 0 && X(i, 2) >= 0
                    PopCon(i, 1) = -(THETA(i, 1) - b);
                    PopCon(i, 2) = -Sx(i, 1) + 0.36;
                    PopCon(i, 3) = Sx(i, 1) - 1;
                else
                    PopCon(i, 1) = THETA(i, 1) - b;
                    PopCon(i, 2) = -Sx(i, 1) + 0.16;
                    PopCon(i, 3) = Sx(i, 1) - 0.49;
                end
            end
            PopCon(PopCon < 0) = 0;
        end
    end
end