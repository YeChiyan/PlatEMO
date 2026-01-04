classdef MMF11 < PROBLEM
    % <multi> <real> <multimodal>
    % MMF11 test function
    
    properties
        POS; % Pareto optimal set
    end
    methods
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            obj.lower    = [0.1, 0.1];
            obj.upper    = [1.1, 1.1];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            PopObj(:,1) = X(:,1);
            temp1 = (sin(2*pi*X(:,2))).^6;
            temp2 = exp(-2*log10(2).*((X(:,2)-0.1)/0.8).^2);
            g = 2 - temp2.*temp1;
            PopObj(:,2) = g./X(:,1);
        end
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0.1, 1.1, N)';
            R(:,2) = 1 ./ R(:,1); % Approximation assuming g_min = 1
            
            % Generate POS for IGDX
            x1 = linspace(0.1, 1.1, N)';
            % PS exists at y=0.25 and y=0.75 (approx where sin^6 is 1)
            obj.POS = [x1, repmat(0.25, N, 1); x1, repmat(0.75, N, 1)];
        end
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
        function score = CalMetric(obj,metName,Population)
            switch metName
                case 'IGDX'
                    score = feval(metName,Population,obj.POS);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
    end
end