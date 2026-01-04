classdef MMF13 < PROBLEM
    % <multi> <real> <multimodal>
    % MMF13 test function
    
    properties
        POS; % Pareto optimal set
    end
    methods
        function Setting(obj)
            obj.M = 2;
            obj.D = 3;
            obj.lower    = [0.1, 0.1, 0.1];
            obj.upper    = [1.1, 1.1, 1.1];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            PopObj(:,1) = X(:,1);
            g = 2 - exp(-2*log10(2)*((X(:,2)+sqrt(X(:,3))-0.1)/0.8).^2).*(sin(2*pi*(X(:,2)+sqrt(X(:,3))))).^6;
            PopObj(:,2) = g./X(:,1);
        end
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0.1, 1.1, N)';
            R(:,2) = 1 ./ R(:,1);
            
            % Sample POS: y + sqrt(z) = 0.25 or 0.75
            x1 = linspace(0.1, 1.1, 50)';
            y  = linspace(0.1, 0.75^2, 20)'; % y can't exceed 0.75 or something
            pos = [];
            % Solution 1: y + sqrt(z) = 0.25 => sqrt(z) = 0.25 - y >= 0 => y <= 0.25
            y1 = linspace(0.1, 0.25-0.01, 10)';
            z1 = (0.25 - y1).^2;
            % Solution 2: y + sqrt(z) = 0.75 => sqrt(z) = 0.75 - y >= 0 => y <= 0.75
            y2 = linspace(0.1, 0.75-0.01, 20)';
            z2 = (0.75 - y2).^2;
            
            for i = 1:length(x1)
                for j = 1:length(y1)
                    pos = [pos; x1(i), y1(j), z1(j)];
                end
                for j = 1:length(y2)
                    pos = [pos; x1(i), y2(j), z2(j)];
                end
            end
            obj.POS = pos;
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