classde f CMMF1 < PROBLEM
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
    end
end